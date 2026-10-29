// 降水・土の水・地下水・気候区分(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 降水: 緯度帯(赤道で湿り、亜熱帯で乾き、中緯度で再び湿る)に地形性降雨を重ねる
// 水蒸気を風で運ぶ。海で蒸発し、斜面を上るか空気が冷えると降る。
// 降水量は緯度ではなく、海からの距離と地形で決まる。
void compute_precip(){
    // 崖崩れで標高が変わるので、使う直前に作り直す
    elev_smooth=elev; box_blur(elev_smooth,4,2);
    for(int y=0;y<GRID_H;y++){
        float u=wind[y];
        int d=(u>0)?1:-1;
        int x0=(d>0)?0:GRID_W-1;
        // 大気の大循環による上昇の強さは条件。区画の緯度で評価する
        float lift=0.5f+0.5f*std::cos(LAT_N*3.14159265f*3.f);
        // 風上の端から入る空気の水蒸気は条件。相対湿度 AIR_IN_RH の空気が、端の気温で入る
        float m=AIR_IN_RH*MOIST_CAP*clamp01((temper[idx(x0,y)]+20.f)/50.f);
        {
            for(int st=0;st<GRID_W;st++){
                int x=x0+d*st;
                int i=idx(x,y);
                if(sea[i]){
                    float sst=std::max(0.f,temper[i]);
                    m += EVAP_COEF*sst*(1.f-m/MOIST_CAP);   // 暖かい海ほど蒸発する
                } else {
                    m += LAND_RECYCLE*soil_water[i];        // 陸からの再循環
                }
                if(m>MOIST_CAP) m=MOIST_CAP;
                // 空気が保持できる量は気温で決まる。冷えれば溢れて降る
                float cap=MOIST_CAP*clamp01((temper[i]+20.f)/50.f);
                // 斜面は広い幅で測る。1セルの凹凸で雨を落とし切らせない
                int up=std::max(0,std::min(GRID_W-1,x-d*ORO_SPAN));
                int span=std::abs(x-up);                  // 風上の端の近くでは測れる幅が縮む
                float rise=(span>0)?std::max(0.f,(elev_smooth[i]-elev_smooth[idx(up,y)])
                                        /(float)span):0.f;
                float rate=RAIN_BASE*lift+OROGRAPHIC*rise;
                float rain=m*rate+std::max(0.f,m-cap)*0.5f;
                if(rain>m) rain=m;
                m-=rain;
                precip[i]=rain;
            }
        }
    }
    // 降水は数十キロの範囲に広がる。行ごとの独立計算による筋も、これで消える
    box_blur(precip,PRECIP_BLUR,2);
}

// 陸の土壌水分の目標(update_water の緩和の釣り合いの値)。乾燥指数も記録する
inline float water_target(int i,float wet){
    // 乾燥限界は気温に比例する。寒い土地は少ない雨でも乾かない
    float thr=ARID_COEF*std::max(0.5f,temper[i]+7.f);
    float A=precip[i]*wet/std::max(1e-9f,thr);
    aridity[i]=A;
    float target=A/(A+ARID_SOFT);          // 飽和させて差をなだらかにする
    if(river[i]) target=std::max(target,RIVER_SUPPLY);
    return target;
}

// 水収支: AW = [(1-RO)*AW + P - ET]+
// ケッペンの乾燥指数にもとづいて土壌水分の目標を決め、そこへ緩やかに近づける
void update_water(int frame){
    float ph=std::sin(2.f*3.14159265f*frame/SEASON_PERIOD);
    for(int y=0;y<GRID_H;y++){
        float wet=1.f+0.45f*LAT_N*ph;             // 雨季と乾季(北半球)
        if(wet<0.f) wet=0.f;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]){
                soil_water[i]=FIELD_CAP;
                // 有機物が沈降して栄養が抜ける。沿岸と湧昇域だけが補われる
                // 帳簿は更新の前後の差で測る。式を分けると FMA の丸めが変わり結果がずれるため、更新の式は元のまま
                float n0=soil_n[i];
                soil_n[i]-=soil_n[i]*N_SINK;
                float n1=soil_n[i];
                soil_n[i]+=N_RUNOFF*shelf[i]*shelf[i];
                g_bud.n_sink+=n1-n0; g_bud.n_runoff+=soil_n[i]-n1;
                if(upwell[i]){ float n2=soil_n[i]; soil_n[i]+=N_UPWELL; g_bud.n_upwell+=soil_n[i]-n2; }
                float pre=soil_n[i];
                if(soil_n[i]<0.f) soil_n[i]=0.f;
                g_bud.n_clamp+=soil_n[i]-pre;
                continue;
            }
            float target=water_target(i,wet);
            float w0=soil_water[i];
            soil_water[i]+=(target-soil_water[i])*WATER_RELAX;
            soil_water[i]=clamp01(soil_water[i]);
            g_bud.w_relax+=soil_water[i]-w0;
        }
    }
}

// 地下水の流れ(ダルシーの式、ブシネスクの近似)。水頭 = 岩盤の高さ + 飽和した厚さ。
// 隣のマスとの間の流れ = K × 上流側の飽和した厚さ × 水頭の差 / Δx × 辺の長さ。海のマスは海面の水頭で一定(帯水層は海につながる)。
// 海の境界は出る向きだけ(陸から海へは流れ出るが、海から陸へは入らない。塩の場が入るまでの扱い)。
// 区画の端の外は、端のマスと一つ内側のマスの水頭の傾きがそのまま続くとする(開いた境界の条件)。
// 地下水面が地表に届いた分は湧き出し、まず根の層の土を満たし、残りは地表を流れ去る
void update_groundwater(int frame){
    if(frame%GW_INTERVAL!=0) return;
    const int N=GRID_W*GRID_H;
    // 辺を通る流量 K·h·(ΔH/Δx)·Δx = K·h·ΔH を、マスの面積 Δx² と比産出率で割って厚さの変化にする
    const float c=GW_K_M*(float)(GW_INTERVAL*DT)/(CELL_M*CELL_M*GW_SY);
    static std::vector<float> head(N), dh(N);
    static std::vector<double> r_sea(GRID_H), r_edge(GRID_H), r_out(GRID_H), r_clamp(GRID_H);
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int i=ya*GRID_W;i<yb*GRID_W;i++)
            head[i]=sea[i]?0.f:(elev[i]-SEA_LEVEL)*ELEV_M-GW_DEPTH_M+gw_h[i]; });
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++){
            double se=0,ed=0;
            for(int x=0;x<GRID_W;x++){
                int i=idx(x,y);
                if(sea[i]){ dh[i]=0.f; continue; }
                float hi=gw_h[i], d=0.f;
                const int dx4[4]={1,-1,0,0}, dy4[4]={0,0,1,-1};
                for(int k=0;k<4;k++){
                    int nx=x+dx4[k], ny=y+dy4[k];
                    if(nx<0||nx>=GRID_W||ny<0||ny>=GRID_H){
                        int ox=x-dx4[k], oy=y-dy4[k];     // 内側の隣(端の外の水頭を外挿する)
                        if(ox<0||ox>=GRID_W||oy<0||oy>=GRID_H) continue;
                        float f=c*hi*(head[i]-head[idx(ox,oy)]);
                        d+=f; ed+=f*GW_W;
                        continue;
                    }
                    int j=idx(nx,ny);
                    float dH=head[j]-head[i];
                    // 海から陸へは入らない(D-2026-10-28-1 の作者の答え。塩の場が入るまで、海水を淡水として作らない)
                    if(sea[j]&&dH>0.f) continue;
                    float up=dH>0.f?gw_h[j]:hi;   // 上流側の飽和した厚さ
                    float f=c*up*dH;
                    d+=f;
                    if(sea[j]) se+=f*GW_W;
                }
                dh[i]=d;
            }
            r_sea[y]=se; r_edge[y]=ed;
        }
    });
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++){
            double out=0,cl=0;
            for(int x=0;x<GRID_W;x++){
                int i=idx(x,y);
                if(sea[i]) continue;
                float h=gw_h[i]+dh[i];
                if(h<0.f){ cl-=h*GW_W; h=0.f; }
                if(h>GW_DEPTH_M){
                    float ex=(h-GW_DEPTH_M)*GW_W;      // 湧き出す水(土の水の単位)
                    float add=std::min(ex,std::max(0.f,FIELD_CAP-soil_water[i]));
                    soil_water[i]+=add;
                    out-=ex-add;
                    h=GW_DEPTH_M;
                }
                gw_h[i]=h;
            }
            r_out[y]=out; r_clamp[y]=cl;
        }
    });
    for(int y=0;y<GRID_H;y++){                 // 足す順を固定する
        g_bud.w_gw_sea+=r_sea[y]; g_bud.w_gw_edge+=r_edge[y];
        g_bud.w_gw_out+=r_out[y]; g_bud.w_gw_clamp+=r_clamp[y];
    }
}

// 年平均気温と降水から気候区分を決める(ケッペンの判定順に従う)
void classify_climate(){
    for(int y=0;y<GRID_H;y++){
        float lat=LAT_N;
        float amp=SEASON_AMP*lat;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]){ climate[i]=0; continue; }
            float Tann=TEMP_EQUATOR+(TEMP_POLE-TEMP_EQUATOR)*lat
                      -LAPSE_RATE*altitude(i);
            temp_annual[i]=Tann;
            float Tw=Tann+amp, Tc=Tann-amp;
            float thr=ARID_COEF*std::max(0.5f,Tann+7.f);
            float A=precip[i]/std::max(1e-9f,thr);
            if(Tw<0.f)       climate[i]=1;   // EF 氷雪
            else if(Tw<10.f) climate[i]=2;   // ET ツンドラ
            else if(A<0.5f)  climate[i]=3;   // BW 砂漠
            else if(A<1.0f)  climate[i]=4;   // BS ステップ
            else if(Tc>=18.f) climate[i]=(A<2.0f)?7:8;  // Aw サバナ / Af 熱帯雨林
            else if(Tc>-3.f)  climate[i]=6;  // C 温帯
            else              climate[i]=5;  // D 亜寒帯
        }
    }
}

