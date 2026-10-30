// 天気(突風・対流の雨雲)・燃料の湿り・火(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 枯れた細い燃料の平衡含水率(乾燥質量あたりの比)。Simard (1968) の式。湿度は %、気温は °F で書かれている
inline float emc_simard(float rh,float tc){
    float H=100.f*rh, T=tc*1.8f+32.f, m;
    if(H<10.f)      m=0.03229f+0.281073f*H-0.000578f*H*T;
    else if(H<50.f) m=2.22749f+0.160107f*H-0.01478f*T;
    else            m=21.0606f+0.005565f*H*H-0.00035f*H*T-0.483199f*H;
    return std::max(0.f,m)*0.01f;
}
// 対流の雨雲。日射と水蒸気から生じ、中層の風に流され、約 1 時間で消える。
// 雲の下には雨が降り、周りには冷気が広がる。冷気は雨より広いので、流れてくる雲の前触れになる
struct Storm { float x,y,age,ipeak; };   // 位置(マス。区画の外もありうる)、齢(分)、成熟期の中心の降水強度(mm/h)
std::vector<Storm> storms;
std::mt19937 wx_rng(12344u^0x2545F491u);   // 天気だけの乱数。雲が生じても生き物の乱数列は動かない
std::uniform_real_distribution<float> wx01(0.f,1.f);
std::vector<float> storm_dT(GRID_W*GRID_H,0.f);   // 冷気による気温の下がり(負)。毎 tick 作り直す
bool   storm_dT_on=false;
long   storm_events=0;
// 雷は雲ごとの雨から別の乱数で起こす。雷を足しても雲の乱数列(発生・強さ)は変わらない
std::mt19937 lt_rng(12344u^0x9E3779B9u);
std::vector<int> strikes;             // この tick に区画の陸へ落ちた対地放電のマス。火の更新が着火に使う
long   lightning_events=0;            // 区画の中への対地放電の累計
double storm_rain_mm=0;               // 前の表示からの、陸に降った雨の陸の平均(mm)
// 飽和水蒸気圧の 20 °C との比(マグヌスの式)
inline float es_ratio(float t){ return std::exp(17.62f*t/(243.12f+t)-17.62f*20.f/(243.12f+20.f)); }
// 突風前線の風を雲の一覧から作る。匂いの更新(ODOR_INTERVAL ごと)で使う
void update_gust(){
    if(gust_on){
        std::fill(gust_u.begin(),gust_u.end(),0.f); std::fill(gust_v.begin(),gust_v.end(),0.f);
        std::fill(gust_row.begin(),gust_row.end(),0); gust_on=false;
    }
    const float sc=STORM_COLD_R_M/CELL_M*0.5f, RG=3.f*sc;
    for(const auto& s:storms){
        float I=s.ipeak*std::sin(3.14159265f*s.age/STORM_LIFE);
        if(I<=0.f) continue;
        float dT=STORM_COLD_DT*std::min(2.f,I/STORM_I_MEAN);
        float c=GUST_FROUDE*std::sqrt(9.81f*dT/(273.15f+g_t_bg)*GUST_DEPTH_M);   // 密度流の速さ(m/s)
        int ya=std::max(0,(int)std::floor(s.y-RG)), yb=std::min(GRID_H,(int)std::ceil(s.y+RG)+1);
        int xa=std::max(0,(int)std::floor(s.x-RG)), xb=std::min(GRID_W,(int)std::ceil(s.x+RG)+1);
        if(ya>=yb||xa>=xb) continue;
        gust_on=true;
        parallel_for(ya,yb,[&](int y0,int y1){
            for(int y=y0;y<y1;y++){
                float dy=y+0.5f-s.y; bool any=false;
                for(int x=xa;x<xb;x++){
                    float dx=x+0.5f-s.x, r2=dx*dx+dy*dy;
                    if(r2>=RG*RG || r2<1e-6f) continue;
                    float r=std::sqrt(r2), q=r/sc;
                    float sp=c*q*std::exp(0.5f*(1.f-q*q));
                    int i=idx(x,y);
                    gust_u[i]+=sp*dx/r; gust_v[i]+=sp*dy/r; any=true;
                }
                if(any) gust_row[y]=1;
            }
        });
    }
    if(gust_on) for(int i=0;i<GRID_W*GRID_H;i++){
        if(sea[i]) continue;
        float m=std::sqrt(gust_u[i]*gust_u[i]+gust_v[i]*gust_v[i]);
        if(m>gust_max) gust_max=m;
    }
}
void update_storms(){
    const float d=(wind[GRID_H/2]>0.f)?1.f:-1.f;
    const float R=STORM_R_M/CELL_M, RC=STORM_COLD_R_M/CELL_M;
    // 発生。区画と、一生のうちに区画へ流れ込みうる風上の範囲で起こす
    if(g_replay){                      // 再生中は記録の各件だけを起こす
        ev_replay_due(EV_STORM,[](const EnvEvent& e){
            storms.push_back({(float)e.x,(float)e.y,0.f,e.strength}); },false);
    } else if(g_sun>0.f){
        const float up=STORM_SHIFT*STORM_LIFE/(float)DT;
        float x0=-RC-(d>0.f?up:0.f), x1=GRID_W+RC+(d<0.f?up:0.f);
        float y0=-RC, y1=GRID_H+RC;
        double area=(double)(x1-x0)*(y1-y0)*CELL_M*CELL_M*1e-6;   // km²
        float cc=es_ratio(g_t_bg);
        double e=STORM_INIT*area*g_sun*cc*DT;
        int n=(int)e; double f=e-n;
        if(f>0.0 && wx01(wx_rng)<f) n++;
        for(int k=0;k<n;k++){
            float x=x0+(x1-x0)*wx01(wx_rng), y=y0+(y1-y0)*wx01(wx_rng);
            float ip=-STORM_I_MEAN*cc*std::log(std::max(1e-6f,1.f-wx01(wx_rng)));
            ip=std::min(ip,STORM_I_MAX);
            storms.push_back({x,y,0.f,ip}); storm_events++;
            ev_record(EV_STORM,(int)std::floor(x),(int)std::floor(y),ip);
        }
    }
    // 流れて歳をとる。一生を終えたか、冷気も届かない風下へ抜けたものは消える
    for(auto& s:storms){ s.x+=d*STORM_SHIFT; s.age+=(float)DT; }
    storms.erase(std::remove_if(storms.begin(),storms.end(),[&](const Storm& s){
        return s.age>=STORM_LIFE || (d>0.f? s.x>GRID_W+RC : s.x<-RC); }),storms.end());
    if(storm_dT_on){ std::fill(storm_dT.begin(),storm_dT.end(),0.f); storm_dT_on=false; }
    strikes.clear();
    if(g_replay)                       // 再生中の雷は記録の各件だけ。再生した雷からは着火しない(火の更新が再生の分岐を通る)
        ev_replay_due(EV_LIGHTNING,[](const EnvEvent& e){ strikes.push_back(idx(e.x,e.y)); });
    if(storms.empty()){                // 雲が無くなれば冷気も突風も消える(残った風が源なしに吹き続けないように)
        if(g_frame%ODOR_INTERVAL==0) update_gust();
        return;
    }
    // 雨と冷気。強さは一生の中で sin の形に盛り上がって衰える(積雲期・成熟期・消散期)
    const float sr=R*0.5f, sc=RC*0.5f;          // ガウス形の幅。半径は 2σ
    static std::vector<double> row_w(GRID_H);
    std::fill(row_w.begin(),row_w.end(),0.0);
    for(const auto& s:storms){
        float I=s.ipeak*std::sin(3.14159265f*s.age/STORM_LIFE);
        if(I<=0.f) continue;
        float dTc=-STORM_COLD_DT*std::min(2.f,I/STORM_I_MEAN);
        if(!g_replay){
            // 対地放電。この tick の雨の質量(中心の強さ × ガウス形の面積、半径 2σ の内側)÷ rain yield
            double sm=sr*CELL_M;
            double kg=(double)I/60.0*DT*2.0*3.14159265358979*sm*sm*(1.0-std::exp(-2.0));
            double e=kg/LIGHT_RAIN_YIELD_KG;
            int n=(int)e; double f=e-n;
            if(f>0.0 && wx01(lt_rng)<f) n++;
            for(int k=0;k<n;k++){
                // 落ちる場所は雨の芯と同じ広がり(ボックス・ミュラー)
                float u1=std::max(1e-7f,wx01(lt_rng)), u2=wx01(lt_rng);
                float rr=sr*std::sqrt(-2.f*std::log(u1));
                int gx=(int)std::floor(s.x+rr*std::cos(6.2831853f*u2)), gy=(int)std::floor(s.y+rr*std::sin(6.2831853f*u2));
                if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
                lightning_events++;
                ev_record(EV_LIGHTNING,gx,gy,1.f);     // 強さ 1 = 対地放電 1 回
                if(!sea[idx(gx,gy)]) strikes.push_back(idx(gx,gy));
            }
        }
        float wc=I/60.f*(float)DT/SOIL_W_MM;    // 中心で 1 tick に降る量(土の水の単位)
        int ya=std::max(0,(int)std::floor(s.y-RC)), yb=std::min(GRID_H,(int)std::ceil(s.y+RC)+1);
        int xa=std::max(0,(int)std::floor(s.x-RC)), xb=std::min(GRID_W,(int)std::ceil(s.x+RC)+1);
        if(ya>=yb||xa>=xb) continue;
        storm_dT_on=true;
        parallel_for(ya,yb,[&](int y0,int y1){
            for(int y=y0;y<y1;y++){
                float dy=y+0.5f-s.y; double wsum=0;
                for(int x=xa;x<xb;x++){
                    int i=idx(x,y);
                    if(sea[i]) continue;
                    float dx=x+0.5f-s.x, r2=dx*dx+dy*dy;
                    if(r2>=RC*RC) continue;
                    storm_dT[i]+=dTc*std::exp(-r2/(2.f*sc*sc));
                    if(r2>=R*R) continue;
                    // 浸透の上限を超える強さの分は地表を流れ去る(流れた水は帳簿の外)
                    float inf=std::min(wc*std::exp(-r2/(2.f*sr*sr)),GW_KS_MM*(float)DT/SOIL_W_MM);
                    float add=std::min(inf,std::max(0.f,FIELD_CAP-soil_water[i]));
                    soil_water[i]+=add;
                    // 根の層を満たした残りは下へ染みて地下水になる。帯水層が満ちていれば流れ去る
                    float perc=std::max(0.f,std::min((inf-add)/GW_W,GW_DEPTH_M-gw_h[i]));
                    gw_h[i]+=perc;
                    // 染みていく水は根の層の水と混ざり、その濃さの塩を地下へ運ぶ(塩の洗い出し)
                    if(perc>0.f && salt_soil[i]>0.f){
                        float pw=perc*GW_W;
                        float sl=salt_soil[i]*pw/(soil_water[i]+pw);
                        salt_soil[i]-=sl; salt_gw[i]+=sl;
                    }
                    wsum+=add+perc*GW_W;
                }
                row_w[y]+=wsum;
            }
        });
    }
    double w=0; for(int y=0;y<GRID_H;y++) w+=row_w[y];   // 足す順を固定する
    g_bud.w_storm+=w;
    storm_rain_mm+=w*SOIL_W_MM;
    if(storm_dT_on) for(int i=0;i<GRID_W*GRID_H;i++) temper[i]+=storm_dT[i];
    if(g_frame%ODOR_INTERVAL==0) update_gust();
}

// 燃料の湿りを大気に追従させる。水蒸気圧は日平均の気温と条件の相対湿度で決まり、一日の中では変わらない。
// 昼に気温が上がれば相対湿度が下がり(飽和水蒸気圧はマグヌスの式)、燃料は時定数で乾く
void update_fuel_moisture(){
    parallel_for(0,GRID_H,[&](int ya,int yb){
        for(int y=ya;y<yb;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]) continue;
            float tn=temper[i];
            // 日平均の気温。雨雲の冷気で下がった分は空気の水蒸気圧を変えないので除く(冷えた分だけ相対湿度が上がる)
            float tm=tn-storm_dT[i]-DIURNAL_LAND*(daylit[i]-0.5f);
            float rh=LAND_RH_MEAN*std::exp(17.62f*tm/(243.12f+tm)-17.62f*tn/(243.12f+tn));
            if(rh>1.f) rh=1.f;
            float emc=emc_simard(rh,tn);
            if(fuel_moist[i]<0.f) fuel_moist[i]=emc;          // 起動時と、古いセーブの読み込みの後は平衡の値から
            else fuel_moist[i]+=(emc-fuel_moist[i])*FUEL_M_LAG;
        }
    });
}

void update_fire(int frame){
    update_fuel_moisture();
    // 酸素の分圧で判定する。柱の気体の全量との比で見る
    double o2m=0.0, tot=0.0;
    for(int i=0;i<GRID_W*GRID_H;i++){ o2m+=o2_at(i); tot+=o2_at(i)+co2_at(i)+N2_COL; }
    float frac = (tot>1e-6f) ? o2m/tot : 0.f;    // 0..1
    // 地球の21%を基準に、火の窓(17%〜35%)へ写す
    float o2r = frac/0.21f;
    if(o2r<FIRE_O2_MIN && !g_replay){         // 酸素が足りず火が広がらない(再生の出来事は外から与える)
        std::fill(fire.begin(),fire.end(),0.f);
        return;
    }
    // 酸素が多いほど、湿った燃料まで燃えるようになる。今の大気(21%)で消炎含水率が実測の値になる
    float wet_ok=FIRE_MX*std::min(2.0f,o2r);

    if(g_replay){
        // 再生中: 落雷と燃え広がりを止め、記録の各件だけを置く。
        // 置いた火は隣へ移らず、その場で減衰して消える(連鎖を起こさない)
        for(int i=0;i<GRID_W*GRID_H;i++){
            float f=fire[i]*(1.f-FIRE_DECAY);
            fire_buf[i]=(f<0.01f)?0.f:f;
        }
        ev_replay_due(EV_FIRE,[](const EnvEvent& e){
            int i=idx(e.x,e.y);
            fire_buf[i]=std::max(fire_buf[i],e.strength);
            fire_events++;
            ev_record(EV_FIRE,e.x,e.y,e.strength);
        });
    } else {
    for(int y=0;y<GRID_H;y++){
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]||fuel[i]<FIRE_FUEL_MIN){ fire_buf[i]=0.f; continue; }
            if(fuel_moist[i]>=wet_ok){ fire_buf[i]=fire[i]*(1.f-FIRE_DECAY*2.f); continue; }
            float f=fire[i];
            // 風下へ移る。隣の火勢を集める
            int up=(wind[y]>0)?x-1:x+1;
            float nb=(up>=0&&up<GRID_W)?fire[idx(up,y)]*1.6f:0.f;   // 区画の外から火は入らない
            if(y>0)          nb+=fire[idx(x,y-1)];
            if(y<GRID_H-1)   nb+=fire[idx(x,y+1)];
            if(x<GRID_W-1)   nb+=fire[idx(x+1,y)];
            if(x>0)          nb+=fire[idx(x-1,y)];
            // 湿りによる燃焼の弱まり(Rothermel 1972 の η_M)。消炎含水率で 0
            float r=fuel_moist[i]/wet_ok;
            float dryness=std::max(0.f,1.f-2.59f*r+5.11f*r*r-3.52f*r*r*r);
            f+=nb*FIRE_SPREAD*dryness;
            f*=(1.f-FIRE_DECAY);
            if(f>3.f) f=3.f;
            fire_buf[i]=f;
            // 燃え移りで燃え始めたマスを1件として記録する(落雷による着火は下で記録する)
            if(fire[i]<0.05f && f>=0.05f) ev_record(EV_FIRE,x,y,f);
        }
    }
    // 落雷による着火。雨雲から落ちた対地放電のうち、長い連続電流を伴うものが乾いた燃料に火を付ける
    for(int i:strikes){
        if(fuel[i]<FIRE_FUEL_MIN || fuel_moist[i]>=wet_ok) continue;
        float r=fuel_moist[i]/wet_ok;
        float dryness=std::max(0.f,1.f-2.59f*r+5.11f*r*r-3.52f*r*r*r);
        if(dryness<=0.f) continue;
        if(wx01(lt_rng)<LIGHT_LCC*dryness*o2r){
            float f0=fire_buf[i];
            fire_buf[i]=std::min(3.f,f0+(1.f-FIRE_DECAY));
            fire_events++;
            if(fire[i]<0.05f && f0<0.05f) ev_record(EV_FIRE,i%GRID_W,i/GRID_W,fire_buf[i]);
        }
    }
    }
    fire.swap(fire_buf);

    // 燃焼: 死骸を灰にし、炭素を大気へ戻す
    double area=0;
    for(auto& cp:corpses){
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int i=idx(gx,gy);
        if(fire[i]<0.05f) continue;
        float burn=cp.energy*FIRE_BURN*std::min(1.f,fire[i]);
        cp.energy-=burn;
        { float bn=cp.n*FIRE_BURN*std::min(1.f,fire[i]); cp.n-=bn; g_bud.n_fire-=bn; }   // 燃えた分の窒素は揮発する
        co2_add(i,burn*CO2_PER_ENERGY);
        o2_add(i,-std::min(burn*FIRE_O2_USE,o2_at(i)));
    }
    for(int i=0;i<GRID_W*GRID_H;i++) if(fire[i]>=0.05f) area+=1;
    for(int i=0;i<GRID_W*GRID_H;i++){
        if(fire[i]>=0.05f) burn_age[i]=0;
        else if(burn_age[i]<999999) burn_age[i]++;
    }
    fire_area=area;
}

