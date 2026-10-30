// 老化・病原体・匂い(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 老化。寿命の後半に入ると、代謝が上がり、繁殖の力が落ちる
inline float senesce(int age,float lifespan,int lmin=LIFE_MIN,int lmax=LIFE_MAX){
    float L=lmin+(lmax-lmin)*lifespan;
    float u=(float)age/L;
    if(u<SENESCE_START) return 0.f;
    float t=(u-SENESCE_START)/(1.f-SENESCE_START);
    return t*t*SENESCE_RATE;
}

inline float saturation(const Plant& p){ return 0.25f+0.75f*(1.f-p.shade_tol); }
inline float opt_temp_p(const Plant& p){ return 26.f-26.f*p.cold_tol; }
inline float kill_temp_p(const Plant& p){
    return 8.f-62.f*p.cold_tol-(p.dormant?DORM_HARDEN:0.f);
}
inline float opt_temp_a(const Animal& a){ return 26.f-26.f*a.cold_tol; }
inline float kill_temp_a(const Animal& a){ return  8.f-62.f*a.cold_tol; }

// 病原体は、その土地で最も多い免疫型を追いかける。
// 型が一致した宿主は病み、稀な型は逃れる(負の頻度依存選択)
inline float immune_match(const float* im,int ci,std::vector<float>* strain){
    float d=0.f;
    for(int k=0;k<N_IMMUNE;k++){ float t=im[k]-strain[k][ci]; d+=t*t; }
    d=std::sqrt(d/N_IMMUNE);
    float z=d/PATH_WIDTH;
    return std::exp(-z*z);          // 近いほど強く感染する
}

void update_pathogens(){
    static std::vector<float> sum_im[N_IMMUNE];
    static std::vector<float> host_n;
    for(int k=0;k<N_IMMUNE;k++) sum_im[k].assign(GRID_W*GRID_H,0.f);
    host_n.assign(GRID_W*GRID_H,0.f);

    // --- 植物 ---
    for(const auto& p:plants){
        if(!p.alive) continue;
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int ci=idx(gx,gy);
        for(int k=0;k<N_IMMUNE;k++) sum_im[k][ci]+=p.im[k];
        host_n[ci]+=1.f;
    }
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        float n=host_n[ci];
        if(n>0.f){
            for(int k=0;k<N_IMMUNE;k++){
                float mean=sum_im[k][ci]/n;
                pstrain[k][ci]+=(mean-pstrain[k][ci])*PATH_ADAPT;   // 多数派へ適応
            }
            pload[ci]+=PATH_GROW*std::min(6.f,n)-PATH_DECAY*pload[ci];
        } else pload[ci]-=PATH_DECAY*3.f*pload[ci];
        if(pload[ci]<0.f) pload[ci]=0.f;
        if(pload[ci]>PATH_LOAD_MAX) pload[ci]=PATH_LOAD_MAX;
    }

    // --- 動物 ---
    for(int k=0;k<N_IMMUNE;k++) sum_im[k].assign(GRID_W*GRID_H,0.f);
    host_n.assign(GRID_W*GRID_H,0.f);
    for(const auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int ci=idx(gx,gy);
        for(int k=0;k<N_IMMUNE;k++) sum_im[k][ci]+=a.im[k];
        host_n[ci]+=1.f;
    }
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        float n=host_n[ci];
        if(n>0.f){
            for(int k=0;k<N_IMMUNE;k++){
                float mean=sum_im[k][ci]/n;
                astrain[k][ci]+=(mean-astrain[k][ci])*PATH_ADAPT;
            }
            aload[ci]+=PATH_GROW*std::min(6.f,n)-PATH_DECAY*aload[ci];
        } else aload[ci]-=PATH_DECAY*3.f*aload[ci];
        if(aload[ci]<0.f) aload[ci]=0.f;
        if(aload[ci]>PATH_LOAD_MAX) aload[ci]=PATH_LOAD_MAX;
    }
}

// 匂いを放ち、拡散させ、風で流す。
// 視覚と違って時間に残り、遮蔽を回り込み、発信者の意図と無関係に漏れる
void update_odor(){
    // --- 放出 ---
    for(const auto& p:plants){
        if(!p.alive||p.seed_wait>0) continue;   // 種子は匂わない
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int i=idx(gx,gy);
        odor[0][i]+=ODOR_PLANT*(p.height+p.tough*1.5f+p.carnivory);
        if(pload[i]>0.5f) odor[0][i]+=ODOR_SICK*0.4f;
    }
    for(const auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int i=idx(gx,gy);
        float hide=(a.depth>0.5f)?0.25f:1.f;      // 潜れば漏れにくい
        odor[0][i]+=ODOR_ANIMAL*(a.body+a.toxin*1.2f)*hide;
        // 発情。繁殖できる状態の個体が、遺伝子 estrus の量だけ出す(出す量と費用の釣り合いは進化に任せる)
        if(estrus_ready(a))
            odor[0][i]+=ODOR_ESTRUS*a.estrus;
        if(aload[i]>0.3f) odor[0][i]+=ODOR_SICK;
    }
    for(const auto& cp:corpses){
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        odor[1][idx(gx,gy)]+=ODOR_ROT*cp.energy;
    }
    for(int i=0;i<GRID_W*GRID_H;i++){
        if(fire[i]>0.05f)            odor[3][i]+=ODOR_SMOKE*fire[i];
        else if(burn_age[i]<BURN_WINDOW) odor[3][i]+=ODOR_SMOKE*0.05f;
    }

    // --- 移流・拡散・散逸 ---
    for(int k=0;k<N_ODOR;k++){
        std::vector<float>& g=odor[k];
        // 1 回の更新(ODOR_INTERVAL tick)で残る割合。指数の形で換算する
        const float keep=(k<3)
            ? (float)std::exp(-(ODOR_K_OH[k]*ODOR_OH_PEAK*g_sun+ODOR_K_O3[k]*ODOR_O3)*60.0*ODOR_INTERVAL*DT)
            : 1.f-ODOR_DECAY;
        // 風で流れる
        parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++){
            float u=wind[y];
            int d=(u>0)?1:-1;
            if(gust_on && gust_row[y]){
                // 突風のある行: 卓越風と突風を重ねた風で、出発点を二次元で戻る(セミラグランジュ法、双線形で補う)。
                // 区画の外から入る空気は、区画の生き物の匂いを持たない(条件)
                const float bx=(std::fabs(u)<1e-4f)?0.f:d*ODOR_SHIFT;
                for(int x=0;x<GRID_W;x++){
                    int i=idx(x,y);
                    float X=x-(bx+gust_u[i]*GUST_ODOR_CELLS), Y=y-gust_v[i]*GUST_ODOR_CELLS;
                    int x0=(int)std::floor(X), y0=(int)std::floor(Y);
                    float fx=X-x0, fy=Y-y0, v=0.f;
                    for(int jy=0;jy<2;jy++) for(int jx=0;jx<2;jx++){
                        int sx=x0+jx, sy=y0+jy;
                        if(sx<0||sx>=GRID_W||sy<0||sy>=GRID_H) continue;
                        v+=g[idx(sx,sy)]*(jx?fx:1.f-fx)*(jy?fy:1.f-fy);
                    }
                    // 地表の層の質量の保存: 広がる流れ(∇·u > 0)では、層から出た空気の分だけ上から空気が降りて入る。
                    // 上の空気は地表の匂いを持たない(条件)ので、濃さは exp(−∇·u·Δt) に薄まる。
                    // 集まる流れ(∇·u < 0)では、空気はその濃さのまま上へ抜けるので、濃さは変わらない(匂いは区画の上へ出る)
                    const int xl=(x>0)?x-1:x, xr=(x<GRID_W-1)?x+1:x, yl=(y>0)?y-1:y, yr=(y<GRID_H-1)?y+1:y;
                    float dv=((gust_u[idx(xr,y)]-gust_u[idx(xl,y)])/(float)(xr-xl)
                             +(gust_v[idx(x,yr)]-gust_v[idx(x,yl)])/(float)(yr-yl))*GUST_ODOR_CELLS;   // 1 回の更新の ∇·u·Δt
                    if(dv>0.f) v*=std::exp(-dv);
                    odor_buf[i]=v;
                }
                continue;
            }
            if(std::fabs(u)<1e-4f){ for(int x=0;x<GRID_W;x++) odor_buf[idx(x,y)]=g[idx(x,y)]; continue; }
            // セミラグランジュ法: 各マスの匂いは、ODOR_SHIFT マス風上にあった空気が運ばれてきたもの。間は線形に補う
            int n=(int)ODOR_SHIFT; float f=ODOR_SHIFT-(float)n;
            for(int x=0;x<GRID_W;x++){
                int a=x-d*n, b=x-d*(n+1);
                // 風上の端から入る空気は、区画の生き物の匂いを持たない(条件)
                float ga=(a>=0&&a<GRID_W)?g[idx(a,y)]:0.f;
                float gb=(b>=0&&b<GRID_W)?g[idx(b,y)]:0.f;
                odor_buf[idx(x,y)]=ga+(gb-ga)*f;
            }
        }
        });
        g.swap(odor_buf);
        // 拡散して散逸する
        parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++) for(int x=0;x<GRID_W;x++){
            float s=0.f; int c=0;
            for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
                int nx=x+dx, ny=y+dy;
                if(nx<0||nx>=GRID_W) continue;
                if(ny<0||ny>=GRID_H) continue;
                s+=g[idx(nx,ny)]; c++;
            }
            int i=idx(x,y);
            float v=g[i]+(s/c-g[i])*ODOR_DIFFUSE;
            v*=keep;
            odor_buf[i]=(v<1e-5f)?0.f:v;
        }
        });
        g.swap(odor_buf);
        double t=0; for(int i=0;i<GRID_W*GRID_H;i++) t+=g[i];
        odor_total[k]=t/(GRID_W*GRID_H);
    }
}

