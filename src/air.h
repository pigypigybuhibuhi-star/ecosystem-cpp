// 気候(気温・湿度)と大気の気体(phase3.cpp から分けた。phase3.cpp が順に #include する)
void update_climate(int frame){
    const float PI=3.14159265f;
    float ph=std::sin(2.f*PI*frame/SEASON_PERIOD);
    // 太陽赤緯。夏の半球へ傾く
    sun_dec = AXIAL_TILT*PI/180.f*ph;
    float sd=std::sin(sun_dec), cd=std::cos(sun_dec);
    // 自転の位相。時間とともに昼の帯が西へ移る
    float rot=2.f*PI*(float)(frame%DAY_PERIOD)/DAY_PERIOD;
    double dsum=0;

    // 緯度と太陽時は区画の中で一様。東西 9 km は経度 0.1° ほどで、太陽時の差は 30 秒に満たない
    float phi=LAT_N*PI/2.f;
    float sp=std::sin(phi), cp=std::cos(phi);
    float base=TEMP_EQUATOR+(TEMP_POLE-TEMP_EQUATOR)*LAT_N;
    float base_sea=SST_EQUATOR+(SST_POLE-SST_EQUATOR)*LAT_N;
    float seas_land=SEASON_AMP*LAT_N*ph;
    float seas_sea =SEASON_AMP*SEA_SEASON*LAT_N*ph;   // 海は熱容量で振れにくい
    float sinh=sp*sd+cp*cd*std::cos(rot);              // 時角は自転の位相だけで決まる
    float sun=(sinh>0.f)?sinh:0.f;                     // 地平線の下なら日射はゼロ
    float dl=clamp01(sinh*5.f);                        // 薄明を挟んだ昼夜
    g_sun=sun;
    g_t_bg=base+seas_land+DIURNAL_LAND*(dl-0.5f);
    for(int y=0;y<GRID_H;y++){
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            solar[i]=sun*SOLAR_GAIN;
            daylit[i]=dl;
            dsum+=dl;
            // 日較差。海は熱容量が大きくほとんど振れない
            float amp=sea[i]?DIURNAL_SEA:DIURNAL_LAND;
            if(sea[i]){
                float t_sea=base_sea+seas_sea+amp*(dl-0.5f)+sst_anom[i];
                // 沿岸の海は、隣接する陸の気温に引かれる
                if(coast_w[i]>0.01f){
                    float t_land=base+seas_land+amp*(dl-0.5f)+sst_anom[i];
                    float pull=COAST_CHILL*coast_w[i];
                    if(t_land<t_sea) t_sea+=(t_land-t_sea)*pull;        // 冷気は強く効く
                    else             t_sea+=(t_land-t_sea)*pull*0.35f;  // 暖気は弱く
                }
                temper[i]=t_sea;
            } else {
                temper[i]=base+seas_land+amp*(dl-0.5f)+sst_anom[i]
                         -LAPSE_RATE*altitude(i);
            }
        }
    }
    // 海氷。凍れば光を遮り、陸の生物には橋になる
    // 氷は一日では張らないし、一日では融けない。厚さとして蓄える
    static std::vector<float> ice_amt(GRID_W*GRID_H,0.f);
    for(int i=0;i<GRID_W*GRID_H;i++){
        if(!sea[i]){ ice[i]=0; ice_amt[i]=0.f; continue; }
        float target = (temper[i]<FREEZE_TEMP) ? 1.f
                     : (temper[i]>THAW_TEMP)   ? 0.f
                     : (THAW_TEMP-temper[i])/(THAW_TEMP-FREEZE_TEMP);
        ice_amt[i]+=(target-ice_amt[i])*ICE_LAG;
        ice[i]=(ice_amt[i]>0.5f)?1:0;
    }
    day_light=(float)(dsum/(GRID_W*GRID_H));
}

inline double grid_sum(const std::vector<float>& g){
    double s=0; for(int i=0;i<GRID_W*GRID_H;i++) s+=g[i]; return s;
}
void update_gases(){
    double s0[2]={grid_sum(co2),grid_sum(o2)};
    // 東西は開いた境界。風上の端からは区画の外の空気(基準の濃度。条件)が入り、
    // 風下の端からは区画の空気が出ていく。その出入りを端の交換として数える
    double edge[2]={0,0};
    for(int y=0;y<GRID_H;y++){
        float u=wind[y];
        if(std::fabs(u)<1e-4f) continue;
        for(int pass=0;pass<2;pass++){
            std::vector<float>& g=(pass==0)?co2:o2;
            const float base=0.f;   // 区画の外の空気は柱の量のまま(ずれ 0)
            // セミラグランジュ法: 各マスの空気は、GAS_SHIFT マス風上にあった空気が運ばれてきたもの。
            // 出発点が区画の外なら、区画の外の空気(基準の濃度)が入る。間は線形に補う
            int d=(u>0)?1:-1;
            int n=(int)GAS_SHIFT; float f=GAS_SHIFT-(float)n;
            double r0=0, r1=0;
            for(int x=0;x<GRID_W;x++) r0+=g[idx(x,y)];
            for(int x=0;x<GRID_W;x++){
                int a=x-d*n, b=x-d*(n+1);
                float ga=(a>=0&&a<GRID_W)?g[idx(a,y)]:base;
                float gb=(b>=0&&b<GRID_W)?g[idx(b,y)]:base;
                gbuf[idx(x,y)]=ga+(gb-ga)*f;
            }
            for(int x=0;x<GRID_W;x++){ g[idx(x,y)]=gbuf[idx(x,y)]; r1+=g[idx(x,y)]; }
            edge[pass]+=r1-r0;   // 区画の中の移しかえは総量を変えない。行の総量の変化は、端から入った分と出た分の差
        }
    }
    for(int pass=0;pass<2;pass++){
        std::vector<float>& g=(pass==0)?co2:o2;
        const float base=0.f;
        float floor_=(pass==0)?-CO2_COL:-O2_COL;   // 全量が 0 を下回らない
        double s1=grid_sum(g), relax=0;
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            float s=0.f; int c=0;
            for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
                int nx=x+dx,ny=y+dy;
                if(nx<0||nx>=GRID_W) continue;
                if(ny<0||ny>=GRID_H) continue;
                s+=g[idx(nx,ny)]; c++;
            }
            float mean=s/c; int i=idx(x,y);
            gbuf[i]=g[i]+(mean-g[i])*GAS_DIFFUSE+(base-g[i])*GAS_RELAX;
            relax+=(base-g[i])*GAS_RELAX;
            if(gbuf[i]<floor_) gbuf[i]=floor_;
        }
        g.swap(gbuf);
        // 移流・拡散は、端の交換を除けば総量を保つ。ずれた分は数値誤差として別に数える
        double trans=(s1-s0[pass]-edge[pass])+(grid_sum(g)-s1-relax);
        if(pass==0){ g_bud.c_relax+=relax; g_bud.c_edge+=edge[pass]; g_bud.c_trans+=trans; }
        else       { g_bud.o_relax+=relax; g_bud.o_edge+=edge[pass]; g_bud.o_trans+=trans; }
    }
}

