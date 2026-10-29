// 地形: ぼかし・窪地・安息角・海流・傾斜・地形の生成・川(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 分離型ボックスブラー(累積和を使うので半径に関係なく速い)
void box_blur(std::vector<float>& f,int R,int passes){
    static std::vector<float> tmp;
    tmp.assign(GRID_W*GRID_H,0.f);
    float inv=1.f/(2*R+1);
    for(int p=0;p<passes;p++){
        for(int y=0;y<GRID_H;y++){                    // 東西も端で止める(開いた境界)
            float sum=0.f;
            for(int dx=-R;dx<=R;dx++){
                int nx=std::max(0,std::min(GRID_W-1,dx));
                sum+=f[idx(nx,y)];
            }
            for(int x=0;x<GRID_W;x++){
                tmp[idx(x,y)]=sum*inv;
                int o=std::max(0,std::min(GRID_W-1,x-R));
                int i=std::max(0,std::min(GRID_W-1,x+R+1));
                sum+=f[idx(i,y)]-f[idx(o,y)];
            }
        }
        for(int x=0;x<GRID_W;x++){                    // 南北は端で止める
            float sum=0.f;
            for(int dy=-R;dy<=R;dy++){
                int ny=std::max(0,std::min(GRID_H-1,dy));
                sum+=tmp[idx(x,ny)];
            }
            for(int y=0;y<GRID_H;y++){
                f[idx(x,y)]=sum*inv;
                int oy=std::max(0,std::min(GRID_H-1,y-R));
                int iy=std::max(0,std::min(GRID_H-1,y+R+1));
                sum+=tmp[idx(x,iy)]-tmp[idx(x,oy)];
            }
        }
    }
}

// 窪地を埋める。どの点からも必ず海まで水が流れるようにする(Priority-Flood法)
void fill_depressions(std::vector<float>& z, float sea_lv){
    struct Node { float v; int i; };
    struct Cmp { bool operator()(const Node&a,const Node&b) const { return a.v>b.v; } };
    std::priority_queue<Node,std::vector<Node>,Cmp> pq;
    std::vector<char> done(GRID_W*GRID_H,0);
    for(int i=0;i<GRID_W*GRID_H;i++)
        if(z[i]<sea_lv){ done[i]=1; pq.push({z[i],i}); }
    // 南の端(海)だけを出口にする。北の端は傾動する地塊の頂(分水嶺)で、
    // 区画の外から川は入らず、外へも流れ出ない(条件)。東西の端も同じく、川の出入りは無いとする。北の端を出口にすると、
    // そこは削られずに隆起だけを受ける基準面になり、高さが限りなく伸びる
    for(int x=0;x<GRID_W;x++){
        int i=idx(x,GRID_H-1);
        if(!done[i]){ done[i]=1; pq.push({z[i],i}); }
    }
    while(!pq.empty()){
        Node n=pq.top(); pq.pop();
        int x=n.i%GRID_W, y=n.i/GRID_W;
        for(int d=0;d<4;d++){
            const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
            int nx=x+DX[d], ny=y+DY[d];
            if(nx<0||nx>=GRID_W) continue;
            if(ny<0||ny>=GRID_H) continue;
            int j=idx(nx,ny);
            if(done[j]) continue;
            if(z[j]<n.v+FILL_EPS_M) z[j]=n.v+FILL_EPS_M;
            done[j]=1; pq.push({z[j],j});
        }
    }
}

// 最急降下で下流を決め、上流から流量を積む
void route_flow(const std::vector<float>& zf, const std::vector<float>& rain,
                std::vector<int>& rec, std::vector<float>& acc,
                std::vector<int>& order){
    int N=GRID_W*GRID_H;
    order.resize(N);
    for(int i=0;i<N;i++) order[i]=i;
    std::sort(order.begin(),order.end(),
              [&](int a,int b){ return zf[a]>zf[b]; });   // 高い順
    rec.assign(N,-1);
    acc=rain;
    const int DX[8]={1,-1,0,0,1,1,-1,-1}, DY[8]={0,0,1,-1,1,-1,1,-1};
    for(int i=0;i<N;i++){
        int x=i%GRID_W, y=i/GRID_W;
        int best=-1; float bs=0.f;
        for(int d=0;d<8;d++){
            int nx=x+DX[d], ny=y+DY[d];
            if(nx<0||nx>=GRID_W) continue;
            if(ny<0||ny>=GRID_H) continue;
            int j=idx(nx,ny);
            float len=(d<4)?1.f:1.41421f;
            float s=(zf[i]-zf[j])/len;
            if(s>bs){ bs=s; best=j; }
        }
        rec[i]=best;
    }
    for(int k=0;k<N;k++){                 // 高い順に流量を下流へ渡す
        int i=order[k];
        if(rec[i]>=0) acc[rec[i]]+=acc[i];
    }
}

// 安息角。急な斜面の物質は下の隣へ落ちる(崩れ)。陸の中だけで動かし、質量は保つ
void relax_repose(std::vector<float>& zm,int max_pass){
    const float dmax=REPOSE_SLOPE*CELL_M;
    const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
    for(int pass=0;pass<max_pass;pass++){
        int moved=0;
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(zm[i]<0.f) continue;
            for(int d=0;d<4;d++){
                int nx=x+DX[d], ny=y+DY[d];
                if(nx<0||nx>=GRID_W) continue;
                if(ny<0||ny>=GRID_H) continue;
                int j=idx(nx,ny);
                if(zm[j]<0.f) continue;
                float diff=zm[i]-zm[j];
                if(diff<=dmax) continue;
                float mv=(diff-dmax)*0.5f;
                zm[i]-=mv; zm[j]+=mv; moved++;
            }
        }
        if(!moved) break;
    }
}

void erode_step(std::vector<float>& z, const std::vector<float>& rain,
                const std::vector<float>& uplift, float sea_lv){
    int N=GRID_W*GRID_H;
    for(int i=0;i<N;i++) z[i]+=uplift[i]*EROSION_DT_YR;
    std::vector<float> zf=z;
    fill_depressions(zf,sea_lv);
    static std::vector<int> rec, order;
    static std::vector<float> acc;
    route_flow(zf,rain,rec,acc,order);

    // 陰的スキーム。下流の新しい高さを先に求め、それを使って上流を解く。
    // 大きな時間刻みでも安定するため、収束が速い
    for(int k=N-1;k>=0;k--){          // order は高い順。逆に辿れば下流から
        int i=order[k];
        int r=rec[i];
        if(r<0) continue;             // 河口。海面に固定され、削られない
        if(z[i]<sea_lv) continue;
        if(acc[i]<CHANNEL_AREA) continue;  // 谷頭より上は川が無い
        int x1=i%GRID_W, y1=i/GRID_W, x2=r%GRID_W, y2=r/GRID_W;
        int dx=x1-x2;
        int dy=y1-y2;
        float dist=std::sqrt((float)(dx*dx+dy*dy))*CELL_M;               // m
        float area=acc[i]*CELL_M*CELL_M;                                  // m²
        float c=EROSION_K*std::pow(area,EROSION_M)/dist;                  // /年
        z[i]=(z[i]+EROSION_DT_YR*c*z[r])/(1.f+c*EROSION_DT_YR);
    }

    // 斜面拡散 dz/dt = D ∇²z。土壌の匍行で尾根が丸くなる。
    // 3×3 の平均との差は ∇²z·Δx²/3 なので、1 回に均す割合は 3·D·Δt/Δx²。
    // 安定のため、割合が 0.5 を超えないように刻みを分ける
    const float kd=3.f*HILL_DIFFUSE_M2_YR*EROSION_DT_YR/(CELL_M*CELL_M);
    const int npass=std::max(1,(int)std::ceil(kd/0.5f));
    static std::vector<float> sm;
    for(int p=0;p<npass;p++){
        sm=z;
        box_blur(sm,1,1);
        for(int i=0;i<N;i++) z[i]+=(sm[i]-z[i])*(kd/npass);
    }
    // 匍行より速く隆起する急な斜面は、安息角を超えた所で崩れて下へ移る
    // (しきい値の斜面。Roering et al. 1999、Montgomery & Brandon 2002)
    relax_repose(z,REPOSE_PASS);
}

// 吹送流。区画は大洋の循環を決めるには小さすぎる(海流は区画の外の条件)。
// 区画の中の海面の流れは、卓越風が押す分だけとする。水温の偏差は持たない
void compute_currents(){
    std::fill(sst_anom.begin(),sst_anom.end(),0.f);
    std::fill(cur_u.begin(),cur_u.end(),0.f);
    std::fill(cur_v.begin(),cur_v.end(),0.f);
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
        int i=idx(x,y);
        if(sea[i]) cur_u[i]=(wind[y]>0.f)?1.f:(wind[y]<0.f?-1.f:0.f);   // 向きだけ。速さは BUOY_DRIFT が持つ
    }
}

// 土地の地形を作る。
// 1. 初期の面: 南が低く北が高い緩い傾きに、波長の違う起伏を重ねる
// 2. 区画の外の地殻変動として、北ほど大きく隆起させる(傾動地塊。条件)
// 3. 隆起と川の侵食(ストリームパワー則)と斜面拡散を繰り返し、谷と尾根を刻ませる
// 4. すべて m・年の単位で計算するので、出てきた高さをそのまま使う(写し直さない)
// 5. 安息角より急な斜面を崩し、物質を下へ移す(質量は保つ)
// 東西は開いた境界。区画の外とはつながらない(外の空気の出入りは風の計算で条件として与える)
// 勾配(m/m)。中心差分。崖崩れで標高が変わったら、そのまわりだけ計算し直す
inline void slope_at(int x,int y){
    int xr=std::min(GRID_W-1,x+1), xl=std::max(0,x-1);
    int yu=std::max(0,y-1), yd=std::min(GRID_H-1,y+1);
    float gx=(elev[idx(xr,y)]-elev[idx(xl,y)])/(float)(xr-xl);
    float gy=(elev[idx(x,yd)]-elev[idx(x,yu)])/(float)(yd-yu);
    slope[idx(x,y)]=std::sqrt(gx*gx+gy*gy)*ELEV_M/CELL_M;
}
void slope_around(int ci){
    int x=ci%GRID_W, y=ci/GRID_W;
    for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
        int nx=x+dx, ny=y+dy;
        if(nx<0||nx>=GRID_W||ny<0||ny>=GRID_H) continue;
        slope_at(nx,ny);
    }
}

void init_terrain(){
    const int N=GRID_W*GRID_H;
    // 起伏。オクターブごとに 0〜1 へ揃え、最後に標準偏差で正規化する
    auto octaves=[&](const int* RR,const float* RA,int n){
        std::vector<float> det(N,0.f), layer(N);
        for(int oc=0;oc<n;oc++){
            for(int i=0;i<N;i++) layer[i]=frand(0.f,1.f);
            box_blur(layer,RR[oc],3);
            float a=1e9f,b=-1e9f;
            for(float v:layer){ if(v<a)a=v; if(v>b)b=v; }
            for(int i=0;i<N;i++) det[i]+=((layer[i]-a)/std::max(1e-6f,b-a)-0.5f)*RA[oc];
        }
        double m=0; for(float v:det) m+=v; m/=N;
        double s2=0; for(float v:det){ double d=v-m; s2+=d*d; }
        float sd=(float)std::sqrt(s2/N)+1e-6f;
        for(auto& v:det) v=(float)((v-m)/sd);
        return det;
    };
    const int   RR1[5]={60,30,15,7,3};
    const float RA1[5]={1.f,0.6f,0.35f,0.2f,0.1f};
    std::vector<float> rough=octaves(RR1,RA1,5);
    const int   RR2[2]={40,20};
    const float RA2[2]={1.f,0.5f};
    std::vector<float> vary=octaves(RR2,RA2,2);

    std::vector<float> h(N), up(N), rain(N,1.f);
    for(int y=0;y<GRID_H;y++){
        float s=1.f-(float)y/(GRID_H-1);              // 南端 0、北端 1
        float land=(s-COAST_FRAC)/(1.f-COAST_FRAC);   // 負なら海の下
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            // 初期の面: 沈降する帯の海底の傾き(南の端で SEA_DEPTH_M)をそのまま陸へ延ばした面
            h[i]=SEA_DEPTH_M/COAST_FRAC*(s-COAST_FRAC)+TERRAIN_NOISE_M*rough[i];
            float u=UPLIFT_MM_YR*1e-3f*clamp01(land)*(1.f+UPLIFT_VARY*clampf(vary[i],-1.f,1.f));
            up[i]=std::max(0.f,u);
        }
    }
    for(int it=0;it<EROSION_STEPS;it++) erode_step(h,rain,up,0.f);

    std::vector<float> zm=h;   // すでに m
    relax_repose(zm,200);
    for(int i=0;i<N;i++) elev[i]=SEA_LEVEL+zm[i]/ELEV_M;

    for(int i=0;i<N;i++) sea[i]=(elev[i]<SEA_LEVEL)?1:0;
    // 勾配は m/m で持つ(崖崩れと歩く速さが物理の勾配で判断できるように)
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++) slope_at(x,y);
    elev_smooth=elev;
    box_blur(elev_smooth,4,2);
    // 卓越風は条件。区画の緯度で、惑星の風の帯の式を一度だけ評価し、区画全体で一様にする
    {
        float w=-std::cos(LAT_N*3.14159265f*3.f)*WIND_STRENGTH;
        for(int y=0;y<GRID_H;y++) wind[y]=w;
    }
    // 陸からの距離。沿岸は河川の栄養が届き、沖は届かない
    std::fill(shelf.begin(),shelf.end(),0.f);
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
        int i=idx(x,y);
        if(!sea[i]){ shelf[i]=1.f; continue; }
        float best=0.f;
        for(int dy=-SHELF_RANGE;dy<=SHELF_RANGE;dy++)
        for(int dx=-SHELF_RANGE;dx<=SHELF_RANGE;dx++){
            int nx=x+dx, ny=y+dy;
            if(nx<0||nx>=GRID_W) continue;
            if(ny<0||ny>=GRID_H) continue;
            if(sea[idx(nx,ny)]) continue;
            float d=std::sqrt((float)(dx*dx+dy*dy));
            float v=1.f-d/(SHELF_RANGE+1.f);
            if(v>best) best=v;
        }
        shelf[i]=best;
    }
    // 湧昇: 風が岸から沖へ吹き、表層水を運び去る場所に深層水が上がる
    for(int y=0;y<GRID_H;y++){
        int d=(wind[y]>0)?1:-1;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(!sea[i]) continue;
            int back=x-d;                              // 風上側
            if(back<0||back>=GRID_W) continue;         // 区画の外は外洋とみなす(陸からの離岸風は無い)
            if(!sea[idx(back,y)] && shelf[i]>0.4f) upwell[i]=1;
        }
    }
    // 陸からの近さ。沿岸の海は、陸の気温に引かれる
    std::fill(coast_w.begin(),coast_w.end(),0.f);
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
        int i=idx(x,y);
        if(!sea[i]) continue;
        float best=0.f;
        for(int dy=-COAST_RANGE;dy<=COAST_RANGE;dy++)
        for(int dx=-COAST_RANGE;dx<=COAST_RANGE;dx++){
            int nx=x+dx, ny=y+dy;
            if(nx<0||nx>=GRID_W) continue;
            if(ny<0||ny>=GRID_H) continue;
            if(sea[idx(nx,ny)]) continue;
            float d=std::sqrt((float)(dx*dx+dy*dy));
            float v=1.f-d/(COAST_RANGE+1.f);
            if(v>best) best=v;
        }
        coast_w[i]=best*best;          // 岸に近いほど強く引かれる
    }
}

inline float altitude(int i){ return std::max(0.f, elev[i]-SEA_LEVEL); }

// 各セルから最も低い隣へ水を流し(D8)、上流面積を累積して河川を決める
void compute_rivers(){
    std::vector<int> order(GRID_W*GRID_H);
    for(int i=0;i<GRID_W*GRID_H;i++) order[i]=i;
    std::sort(order.begin(),order.end(),
              [](int a,int b){ return elev[a]>elev[b]; });
    std::fill(flow_acc.begin(),flow_acc.end(),1.f);
    for(int i : order){
        if(sea[i]) continue;
        int x=i%GRID_W, y=i/GRID_W;
        int best=-1; float lowest=elev[i];
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            if(dx==0&&dy==0) continue;
            int nx=x+dx, ny=y+dy;
            if(nx<0||nx>=GRID_W) continue;
            if(ny<0||ny>=GRID_H) continue;
            int ni=idx(nx,ny);
            if(elev[ni]<lowest){ lowest=elev[ni]; best=ni; }
        }
        if(best>=0) flow_acc[best]+=flow_acc[i];   // 下流へ積み上げる
    }
    for(int i=0;i<GRID_W*GRID_H;i++)
        river[i] = (!sea[i] && flow_acc[i]>RIVER_THRESHOLD) ? 1 : 0;
}

