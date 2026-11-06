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
    double sp=0; long nl=0;
    for(int i=0;i<GRID_W*GRID_H;i++) if(!sea[i]){ sp+=precip[i]; nl++; }
    precip_land_mean=nl?sp/nl:0.0;
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

// 雨が地面に染み込む。inf は染み込む量(土の水の単位。浸透の上限で切った後)。
// まず根の層を圃場容水量まで満たし、残りは重力で下へ染みて地下水になる。帯水層が満ちていれば染み込めず地表を流れ去る(帳簿の外へ出る)。
// 染みていく水は根の層の水と混ざり、その濃さの塩を地下へ運ぶ(塩の洗い出し)。雨雲の雨と背景の雨が同じ仕組みを通る。
// 返す値は土と地下水に入った水(土の水の単位)
inline float infiltrate(int i,float inf){
    float add=std::min(inf,std::max(0.f,FIELD_CAP-soil_water[i]));
    soil_water[i]+=add;
    float perc=std::max(0.f,std::min((inf-add)/GW_W,GW_DEPTH_M-gw_h[i]));
    gw_h[i]+=perc;
    if(perc>0.f && salt_soil[i]>0.f){
        float pw=perc*GW_W;
        float sl=salt_soil[i]*pw/(soil_water[i]+pw);
        salt_soil[i]-=sl; salt_gw[i]+=sl;
    }
    return add+perc*GW_W;
}

inline float rel_humidity(int i);   // weather_fire.h

// 水収支(D-2026-10-28-2 の作者の答え (a)): 土の水 = 前の値 + 染み込んだ雨 − 土からの蒸発(蒸散は光合成の側で引く)。
// 根の層に収まらない分は地下水へ、帯水層からも溢れれば流れ去る。目標への緩和は置かない
void update_water(int frame){
    float ph=std::sin(2.f*3.14159265f*frame/SEASON_PERIOD);
    for(int y=0;y<GRID_H;y++){
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(!sea[i]) continue;
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
        }
    }
    // 背景の降水: 降水の場の陸の平均を、条件の年降水量 P_BG_ANN_MM に結びつける(場所ごとの比は降水の場が決める)
    const float pscale=precip_land_mean>0.0
        ? (float)(P_BG_ANN_MM/YEAR_MIN*DT/SOIL_W_MM/precip_land_mean) : 0.f;
    const float ks=GW_KS_MM*(float)DT/SOIL_W_MM;           // 浸透の上限(1 tick)
    const float rs=SOLAR_S0*CLEAR_SKY_TAU*CLOUD_SW_RATIO*g_sun;   // 地表の短波 W/m²(晴天 × 雲の比)
    const float lw_cloud=1.35f*CLOUD_SW_RATIO-0.35f;             // FAO-56 の長波の雲の係数
    const bool  day=g_sun>0.f;
    // 可能蒸発 kg/m²/s を 1 tick の土の水の単位へ(1 kg/m² = 1 mm)
    const float e_unit=60.f*(float)DT/SOIL_W_MM;
    static std::vector<double> r_rain(GRID_H), r_evap(GRID_H);
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++){
            float wet=1.f+0.45f*LAT_N*ph;             // 雨季と乾季(北半球)
            if(wet<0.f) wet=0.f;
            double wr=0, we=0;
            for(int x=0;x<GRID_W;x++){
                int i=idx(x,y);
                if(sea[i]) continue;
                float inf=std::min(precip[i]*wet*pscale,ks);
                if(inf>0.f) wr+=infiltrate(i,inf);
                if(!day || soil_water[i]<=0.f) continue;
                // Priestley–Taylor。飽和水蒸気圧はテテンの式、Δ はその傾き、乾湿計定数は標高の気圧から(FAO-56)
                float T=temper[i];
                float es=0.6108f*std::exp(17.27f*T/(T+237.3f));         // kPa
                float del=4098.f*es/((T+237.3f)*(T+237.3f));
                float z=altitude(i)*ELEV_M;
                float pk=101.3f*std::pow((293.f-0.0065f*z)/293.f,5.26f);
                float gam=0.000665f*pk;
                float ea=rel_humidity(i)*es;
                float tk=T+273.16f;
                float rnl=SB_SIGMA*tk*tk*tk*tk*(0.34f-0.14f*std::sqrt(std::max(0.f,ea)))*lw_cloud;
                float rn=(1.f-SURF_ALBEDO)*rs-rnl;
                if(rn<=0.f) continue;                  // 正味放射が負なら蒸発しない(露は置かない)
                float pet=PT_ALPHA*del/(del+gam)*rn/LATENT_VAP*e_unit;
                float beta=std::min(1.f,soil_water[i]/(EVAP_BETA_FRAC*FIELD_CAP));
                float e=std::min(soil_water[i],pet*beta*ground_frac[i]);
                soil_water[i]-=e;
                we-=e;
            }
            r_rain[y]=wr; r_evap[y]=we;
        }
    });
    for(int y=0;y<GRID_H;y++){ g_bud.w_rain+=r_rain[y]; g_bud.w_evap+=r_evap[y]; }   // 足す順を固定する
}

// 地下水の流れ(ダルシーの式、ブシネスクの近似)。水頭 = 岩盤の高さ + 飽和した厚さ。
// 隣のマスとの間の流れ = K × 上流側の飽和した厚さ × 水頭の差 / Δx × 辺の長さ。海のマスは海面の水頭で一定(帯水層は海につながる)。
// 海のマスとの間は両向きに流れる。海から入る水は塩分 SEA_SALT の海水(D-2026-10-28-1 の作者の答え (b))。
// 塩は水と一緒に運ばれる(移流。辺を通る塩 = 水の流量 × 上流側の濃さ)。分散と拡散は入れていない(下の注)。
// 区画の端の外は、端のマスと一つ内側のマスの水頭の傾きがそのまま続くとする(開いた境界の条件)。端の外の水の濃さは端のマスと同じとする。
// 地下水面が地表に届いた分は湧き出し、まず根の層の土を満たし、残りは地表を流れ去る。湧いた水は帯水層の濃さの塩を持って出る
// 注: 帯水層と根の層の間の塩の拡散は、分子拡散 1e-9 m²/s で厚さ 1〜3 m を渡るのに数十年かかる(Freeze & Cherry 1979)ので入れていない
void update_groundwater(int frame){
    if(frame%GW_INTERVAL!=0) return;
    const int N=GRID_W*GRID_H;
    // 辺を通る流量 K·h·(ΔH/Δx)·Δx = K·h·ΔH を、マスの面積 Δx² と比産出率で割って厚さの変化にする
    const float c=GW_K_M*(float)(GW_INTERVAL*DT)/(CELL_M*CELL_M*GW_SY);
    static std::vector<float> head(N), dh(N), conc(N), ds(N);
    static std::vector<double> r_sea(GRID_H), r_edge(GRID_H), r_out(GRID_H), r_clamp(GRID_H);
    static std::vector<double> rs_sea(GRID_H), rs_edge(GRID_H), rs_out(GRID_H), rs_clamp(GRID_H);
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int i=ya*GRID_W;i<yb*GRID_W;i++){
            head[i]=sea[i]?0.f:(elev[i]-SEA_LEVEL)*ELEV_M-GW_DEPTH_M+gw_h[i];
            // 帯水層の水の塩分(g/kg)。水が無ければ 0
            conc[i]=sea[i]?SEA_SALT:(gw_h[i]>0.f?salt_gw[i]/(gw_h[i]*GW_W):0.f);
        } });
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++){
            double se=0,ed=0,sse=0,sed=0;
            for(int x=0;x<GRID_W;x++){
                int i=idx(x,y);
                if(sea[i]){ dh[i]=0.f; ds[i]=0.f; continue; }
                float hi=gw_h[i], d=0.f, s=0.f;
                const int dx4[4]={1,-1,0,0}, dy4[4]={0,0,1,-1};
                for(int k=0;k<4;k++){
                    int nx=x+dx4[k], ny=y+dy4[k];
                    if(nx<0||nx>=GRID_W||ny<0||ny>=GRID_H){
                        int ox=x-dx4[k], oy=y-dy4[k];     // 内側の隣(端の外の水頭を外挿する)
                        if(ox<0||ox>=GRID_W||oy<0||oy>=GRID_H) continue;
                        float f=c*hi*(head[i]-head[idx(ox,oy)]);
                        float fs=f*GW_W*conc[i];
                        d+=f; ed+=f*GW_W;
                        s+=fs; sed+=fs;
                        continue;
                    }
                    int j=idx(nx,ny);
                    float dH=head[j]-head[i];
                    // 上流側の飽和した厚さ。海のマスは帯水層が海につながり、満ちているとする
                    float up=dH>0.f?(sea[j]?GW_DEPTH_M:gw_h[j]):hi;
                    float f=c*up*dH;
                    float fs=f*GW_W*(dH>0.f?conc[j]:conc[i]);   // 上流側の濃さの塩が一緒に動く
                    d+=f; s+=fs;
                    if(sea[j]){ se+=f*GW_W; sse+=fs; }
                }
                dh[i]=d; ds[i]=s;
            }
            r_sea[y]=se; r_edge[y]=ed; rs_sea[y]=sse; rs_edge[y]=sed;
        }
    });
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++){
            double out=0,cl=0,sout=0,scl=0;
            for(int x=0;x<GRID_W;x++){
                int i=idx(x,y);
                if(sea[i]) continue;
                float h=gw_h[i]+dh[i];
                float sg=salt_gw[i]+ds[i];
                if(h<0.f){ cl-=h*GW_W; h=0.f; }
                if(h<=0.f){ scl-=sg; sg=0.f; }          // 水が無くなったマスに残る塩は、水の切り捨てと一緒に消える(0 であるべき)
                else if(sg<0.f){ scl-=sg; sg=0.f; }
                if(h>GW_DEPTH_M){
                    float ex=(h-GW_DEPTH_M)*GW_W;      // 湧き出す水(土の水の単位)
                    float cs=sg/(h*GW_W);              // 湧く水の塩分
                    float add=std::min(ex,std::max(0.f,FIELD_CAP-soil_water[i]));
                    soil_water[i]+=add;
                    out-=ex-add;
                    float s_add=add*cs, s_ex=ex*cs;
                    salt_soil[i]+=s_add;
                    sout-=s_ex-s_add;
                    sg-=s_ex;
                    if(sg<0.f){ scl-=sg; sg=0.f; }
                    h=GW_DEPTH_M;
                }
                gw_h[i]=h; salt_gw[i]=sg;
            }
            r_out[y]=out; r_clamp[y]=cl; rs_out[y]=sout; rs_clamp[y]=scl;
        }
    });
    for(int y=0;y<GRID_H;y++){                 // 足す順を固定する
        g_bud.w_gw_sea+=r_sea[y]; g_bud.w_gw_edge+=r_edge[y];
        g_bud.w_gw_out+=r_out[y]; g_bud.w_gw_clamp+=r_clamp[y];
        g_bud.s_sea+=rs_sea[y]; g_bud.s_edge+=rs_edge[y];
        g_bud.s_out+=rs_out[y]; g_bud.s_clamp+=rs_clamp[y];
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

