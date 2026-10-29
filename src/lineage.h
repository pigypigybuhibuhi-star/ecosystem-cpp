// 系統(観察のための分類)・平均・要約・絞り込み(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ===== 系統 =====
struct Lineage {
    int   id, kind, num;      // kind 0=植物 1=動物 / num は種別ごとの通し番号
    float rep[REP_N];
    int   nrep, birth_tick, death_tick, parent, missing;
    bool  alive;
    int   pop, max_pop;
    float mean[REP_N];        // 表示用の平均形質
    float final_mean[REP_N];  // 絶滅した瞬間の平均。個体が消えても残す
    int   final_pop;          // 絶滅直前の個体数
};
std::vector<Lineage> lineages;
int next_lineage_id = 0;
int next_plant_num = 1, next_animal_num = 1;

inline int plant_vec(const Plant& p, float* v){
    v[0]=p.absorb         *1.0f;
    v[1]=p.max_height/5.f *2.5f;
    v[2]=p.shade_tol      *2.0f;
    v[3]=p.tough          *2.0f;
    v[4]=p.clonal         *1.5f;
    v[5]=p.disperse       *1.0f;
    v[6]=p.cold_tol       *1.5f;
    v[7]=p.aquatic        *3.0f;   // 生息環境は最も重い軸
    v[8]=p.drought_tol    *2.0f;
    return 9;
}
inline int animal_vec(const Animal& a, float* v){
    v[0]=a.speed          *1.0f;
    v[1]=a.turn           *0.5f;
    v[2]=a.reach/5.f      *2.0f;
    v[3]=a.size/4.f       *2.5f;
    v[4]=a.diet           *3.0f;
    v[5]=a.jaw/2.f        *1.5f;
    v[6]=a.cold_tol       *1.5f;
    // 神経網の規模を、結線の強さの平均で代表させる
    { float s=0.f; for(int k=0;k<N_W;k++) s+=std::fabs(a.w[k]);
      v[7]=clamp01(s/(N_W*1.5f))*0.8f; }
    // 神経網の形。深さと幅が、生態的な位置を分ける
    { int nn=0; for(int k=0;k<N_NODE;k++) if(a.ngain[k]>NODE_ON) nn++;
      v[8]=(float)nn/N_NODE *1.5f; }
    v[9]=a.aquatic        *3.0f;
    // 感覚器の構成。何をどこまで見るかは、生態的な位置を強く表す
    float rs=0.f, tn=0.f; int nv=0;
    for(int k=0;k<N_SENSOR;k++)
        if(a.sen[k].range>=SENSOR_MIN){ rs+=a.sen[k].range; tn+=a.sen[k].tune; nv++; }
    v[10]=rs/(SENSOR_MAX*N_SENSOR) *2.0f;
    v[11]=(nv? tn/nv : 0.f)        *2.5f;
    v[12]=a.arboreal               *2.0f;   // 地上性か樹上性かは生態的な位置を分ける
    return 13;
}
inline float vdist(const float* a,const float* b,int n){
    float s=0.f; for(int i=0;i<n;i++){ float d=a[i]-b[i]; s+=d*d; } return std::sqrt(s);
}
inline float depthn(int i){ return std::max(0.f,(SEA_LEVEL-elev[i])/SEA_LEVEL); }
// 水に浸かることへの不適合(植物)。aquatic は水に浸かることへの強さ(通気組織)。塩は別に salt_mismatch で測る
inline float habitat_p(const Plant& p,int ci){
    if(sea[ci]) return (1.f-p.aquatic)*(SHALLOW_EASE+(1.f-SHALLOW_EASE)*depthn(ci));
    return p.aquatic;
}
// 根の層の水の塩分(g/kg)。海は SEA_SALT。土の水が無いところに塩があれば、結晶の塩の上で飽和の濃さとみなす
inline float soil_salinity(int ci){
    if(sea[ci]) return SEA_SALT;
    float w=soil_water[ci], s=salt_soil[ci];
    if(s<=0.f) return 0.f;
    if(w<=1e-6f) return HALITE_SAT;
    return std::min(HALITE_SAT,s/w);
}
// 塩への不適合。塩への強さ salt_tol は、適合溶質を作り浸透圧を合わせられる塩分(海水を 1 とする比)。
// まわりの水がそれより塩辛ければ水を吸えず(生理的な乾燥)、薄ければ溜めた溶質と輸送の仕組みが無駄になる
// (塩生植物の多くは淡水で育ちが落ちる。Flowers & Colmer 2008)。生息環境の不適合と同じく差の二乗で費用にする
inline float salt_mismatch(const Plant& p,int ci){
    return soil_salinity(ci)/SEA_SALT-p.salt_tol;
}
inline float habitat_a(const Animal& a,int ci){
    if(sea[ci]) return (1.f-a.aquatic)*(SHALLOW_EASE+(1.f-SHALLOW_EASE)*depthn(ci));
    return a.aquatic;
}
sf::Color lineage_color(int i){
    float h=std::fmod(0.13f+(float)i*0.61803399f,1.0f)*360.f;
    float s=0.70f,v=1.0f;
    float c=v*s,x=c*(1.f-std::fabs(std::fmod(h/60.f,2.f)-1.f)),m=v-c;
    float r,g,b;
    if(h<60){r=c;g=x;b=0;} else if(h<120){r=x;g=c;b=0;}
    else if(h<180){r=0;g=c;b=x;} else if(h<240){r=0;g=x;b=c;}
    else if(h<300){r=x;g=0;b=c;} else {r=c;g=0;b=x;}
    return sf::Color((int)((r+m)*255),(int)((g+m)*255),(int)((b+m)*255));
}

void update_lineages(int kind,int tick){
    int total = (kind==0)?(int)plants.size():(int)animals.size();
    if(total<=0){
        for(auto& L:lineages)
            if(L.alive && L.kind==kind){ L.alive=false; L.death_tick=tick; }
        return;
    }
    const int D = (kind==0)?9:13;
    std::vector<float> reps;
    float v[14];
    int stride=std::max(1,total/LINEAGE_SAMPLE), seen=0;
    for(int i=0;i<total;i++){
        if((seen++ % stride)!=0) continue;
        if(kind==0) plant_vec(plants[i],v); else animal_vec(animals[i],v);
        int best=-1; float bd=1e9f; int nc=(int)reps.size()/D;
        for(int c=0;c<nc;c++){ float d=vdist(v,&reps[c*D],D); if(d<bd){bd=d;best=c;} }
        if(best<0 || bd>=SPECIES_THRESHOLD)
            for(int k=0;k<D;k++) reps.push_back(v[k]);
    }
    int nc=(int)reps.size()/D;
    int oldN=(int)lineages.size();
    std::vector<int>  cl(nc,-1);
    std::vector<char> claimed(oldN,0);
    for(int c=0;c<nc;c++){
        int best=-1; float bd=1e9f;
        for(int L=0;L<oldN;L++){
            if(!lineages[L].alive || lineages[L].kind!=kind) continue;
            float d=vdist(&reps[c*D],lineages[L].rep,D);
            if(d<bd){ bd=d; best=L; }
        }
        if(best>=0 && bd<SPECIES_THRESHOLD && !claimed[best]){
            cl[c]=best; claimed[best]=1; lineages[best].missing=0;
            for(int k=0;k<D;k++) lineages[best].rep[k]=reps[c*D+k];
        } else {
            Lineage L; L.id=next_lineage_id++; L.kind=kind; L.nrep=D;
            L.num = (kind==0)? next_plant_num++ : next_animal_num++;
            for(int k=0;k<REP_N;k++){ L.rep[k]=0.f; L.mean[k]=0.f; L.final_mean[k]=0.f; }
            L.final_pop=0;
            for(int k=0;k<D;k++) L.rep[k]=reps[c*D+k];
            L.birth_tick=tick; L.death_tick=-1; L.missing=0;
            L.parent=(best>=0)?lineages[best].id:-1;
            L.alive=true; L.pop=0; L.max_pop=0;
            cl[c]=(int)lineages.size(); lineages.push_back(L);
        }
    }
    // 一度見失っただけでは絶滅としない(標本の揺らぎを吸収する)
    for(int L=0;L<oldN;L++)
        if(lineages[L].alive && lineages[L].kind==kind && !claimed[L]){
            lineages[L].missing++;
            if(lineages[L].missing >= LINEAGE_GRACE){
                lineages[L].alive=false; lineages[L].death_tick=tick;
                // 最後に見えていた姿を残す。個体がいなくなっても参照できる
                for(int k=0;k<REP_N;k++) lineages[L].final_mean[k]=lineages[L].mean[k];
                lineages[L].final_pop=lineages[L].pop;
            }
        }

    for(auto& L:lineages) if(L.kind==kind){ L.pop=0; for(int k=0;k<REP_N;k++) L.mean[k]=0.f; }
    for(int i=0;i<total;i++){
        if(kind==0) plant_vec(plants[i],v); else animal_vec(animals[i],v);
        int best=-1; float bd=1e9f;
        for(int c=0;c<nc;c++){ float d=vdist(v,&reps[c*D],D); if(d<bd){bd=d;best=c;} }
        if(best<0) continue;
        Lineage& L = lineages[cl[best]];
        if(kind==0){
            const Plant& p=plants[i];
            plants[i].lin=L.id;
            L.mean[0]+=p.absorb;   L.mean[1]+=p.max_height; L.mean[2]+=p.shade_tol;
            L.mean[3]+=p.tough;    L.mean[4]+=p.clonal;     L.mean[5]+=p.disperse;
            L.mean[6]+=p.cold_tol; L.mean[7]+=p.height; L.mean[8]+=p.aquatic;
            L.mean[9]+=p.drought_tol;
        } else {
            const Animal& a=animals[i];
            animals[i].lin=L.id;
            L.mean[0]+=a.speed; L.mean[1]+=a.turn;  L.mean[2]+=a.reach;
            L.mean[3]+=a.size;  L.mean[4]+=a.diet;  L.mean[5]+=a.jaw;
            L.mean[6]+=a.cold_tol; L.mean[7]+=a.body;
            L.mean[8]+=a.w[0];
            { int nn=0; for(int k=0;k<N_NODE;k++) if(a.ngain[k]>NODE_ON) nn++;
              L.mean[9]+=(float)nn; } L.mean[10]+=a.aquatic;
            { float rs=0.f; int nv=0;
              for(int k=0;k<N_SENSOR;k++)
                  if(a.sen[k].range>=SENSOR_MIN){ rs+=a.sen[k].range; nv++; }
              L.mean[11]+=rs; }
            L.mean[12]+=a.arboreal;
            L.mean[13]+=a.perch;
        }
        L.pop++;
    }
    for(auto& L:lineages){
        if(L.kind!=kind) continue;
        if(L.pop>0){
            for(int k=0;k<REP_N;k++) L.mean[k]/=(float)L.pop;
            // 生きている間、最後の姿を常に控えておく。
            // 絶滅が確定した時点では mean が既に0になっているため
            for(int k=0;k<REP_N;k++) L.final_mean[k]=L.mean[k];
            L.final_pop=L.pop;
        }
        if(L.pop>L.max_pop) L.max_pop=L.pop;
    }
}

// 生きていれば現在の平均、絶滅していれば最後に見えていた姿を返す
inline const float* view_mean(const Lineage& L){
    return L.alive ? L.mean : L.final_mean;
}
// 集団平均。レーダー図に重ねて「どこが尖っているか」を見せる
float g_pavg[13], g_aavg[13];
void compute_averages(){
    for(int k=0;k<13;k++){ g_pavg[k]=0.f; g_aavg[k]=0.f; }
    int np=0;
    for(const auto& p:plants){
        if(!p.alive||p.seed_wait>0) continue;
        g_pavg[0]+=p.max_height/5.f; g_pavg[1]+=p.absorb;    g_pavg[2]+=p.shade_tol;
        g_pavg[3]+=p.tough;          g_pavg[4]+=p.clonal;    g_pavg[5]+=p.disperse;
        g_pavg[6]+=p.cold_tol;       g_pavg[7]+=p.aquatic;   g_pavg[8]+=p.drought_tol;
        g_pavg[9]+=p.climb;          g_pavg[10]+=p.fire_tol;
        g_pavg[11]+=p.store_cap/STORE_MAX; g_pavg[12]+=p.lifespan;
        np++;
    }
    if(np) for(int k=0;k<13;k++) g_pavg[k]/=np;
    int na=0;
    for(const auto& a:animals){
        if(!a.alive) continue;
        float sr=0.f; for(int k=0;k<N_SENSOR;k++) sr+=a.sen[k].range;
        g_aavg[0]+=a.diet;       g_aavg[1]+=a.body/4.f;  g_aavg[2]+=a.reach/5.f;
        g_aavg[3]+=a.speed;      g_aavg[4]+=a.jaw/2.f;   g_aavg[5]+=a.cold_tol;
        g_aavg[6]+=a.aquatic;    g_aavg[7]+=a.detritus;
        g_aavg[8]+=sr/(SENSOR_MAX*N_SENSOR);
        g_aavg[9]+=a.arboreal;   g_aavg[10]+=a.resp;
        g_aavg[11]+=a.brood;     g_aavg[12]+=a.lifespan;
        na++;
    }
    if(na) for(int k=0;k<13;k++) g_aavg[k]/=na;
}
// その個体が何者かを一文で言う。レーダー図が読めない人にも伝わる
std::string plant_summary(const Plant& p){
    std::vector<std::string> t;
    if(p.aquatic>0.5f)      t.push_back("水生");
    if(p.parasite>0.35f)    t.push_back("寄生");
    if(p.carnivory>0.30f)   t.push_back("食虫");
    if(p.climb>0.45f)       t.push_back("つる");
    if(p.max_height>1.5f)   t.push_back("高木");
    else if(p.max_height<0.4f) t.push_back("低木");
    if(p.shade_tol>0.6f)    t.push_back("陰生");
    else if(p.shade_tol<0.2f) t.push_back("陽生");
    if(p.tough>0.5f)        t.push_back("硬い");
    if(p.cold_tol>0.5f)     t.push_back("耐寒");
    if(p.drought_tol>0.5f)  t.push_back("耐乾");
    if(p.n_fix>0.4f)        t.push_back("窒素固定");
    if(p.clonal>0.5f)       t.push_back("地下茎");
    if(t.empty()) t.push_back("ふつうの草");
    std::string s;
    for(size_t i=0;i<t.size()&&i<5;i++){ if(i) s+="・"; s+=t[i]; }
    // 生き方を一言で
    if(p.parasite>0.35f)      s+="  — 他人から奪う型";
    else if(p.climb>0.45f)    s+="  — 他人に登る型";
    else if(p.max_height>1.5f)s+="  — 競争で勝つ型";
    else if(p.brood>0.5f)     s+="  — 数で押す型";
    else if(p.cold_tol>0.5f||p.drought_tol>0.5f) s+="  — 端で耐える型";
    return s;
}
std::string animal_summary(const Animal& a){
    std::vector<std::string> t;
    if(a.aquatic>0.5f)      t.push_back("水生");
    if(a.body>1.5f)         t.push_back("大型");
    else if(a.body<0.5f)    t.push_back("小型");
    if(a.diet>0.6f)         t.push_back("肉食");
    else if(a.diet<0.2f)    t.push_back("草食");
    else                    t.push_back("雑食");
    if(a.detritus>0.4f)     t.push_back("腐肉食");
    if(a.nocturnal>0.5f)    t.push_back("夜行性");
    if(a.burrow>0.4f)       t.push_back("穴掘り");
    if(a.arboreal>0.4f)     t.push_back("樹上性");
    if(a.toxin>0.3f)        t.push_back("有毒");
    if(a.speed>0.8f)        t.push_back("俊足");
    if(t.empty()) t.push_back("ふつうの動物");
    std::string s;
    for(size_t i=0;i<t.size()&&i<5;i++){ if(i) s+="・"; s+=t[i]; }
    if(a.diet>0.6f&&a.speed>0.7f)      s+="  — 追いかける型";
    else if(a.diet>0.6f)               s+="  — 待ち伏せ型";
    else if(a.detritus>0.4f)           s+="  — 掃除する型";
    else if(a.burrow>0.4f)             s+="  — 隠れる型";
    else                               s+="  — 食べて逃げる型";
    return s;
}
// フィルタ
const int FILTER_N = 10;
const char* FILTER_NAME[FILTER_N] = {
    "全て","現存","絶滅","耐寒型","肉食","草食","高木","硬い",
    "水生","陸生"
};
bool filter_on[FILTER_N] = {true,false,false,false,false,false,false,false,false,false};
bool filter_single(const Lineage& L,int f);
// 有効なフィルタ全てを満たすものだけ通す(ANDで合成)
bool lineage_match_multi(const Lineage& L){
    bool any=false;
    for(int f=0;f<FILTER_N;f++){
        if(!filter_on[f]) continue;
        any=true;
        if(!filter_single(L,f)) return false;
    }
    return any ? true : true;
}
bool filter_single(const Lineage& L,int f){
    switch(f){
        case 0: return true;
        case 1: return L.alive;
        case 2: return !L.alive;
        case 3: return view_mean(L)[6] > 0.60f;
        case 4: return L.kind==1 && view_mean(L)[4] > 0.60f;
        case 5: return L.kind==1 && view_mean(L)[4] < 0.30f;
        case 6: return L.kind==0 && view_mean(L)[1] > 1.50f;
        case 7: return L.kind==0 && view_mean(L)[3] > 0.40f;
        case 8: return (L.kind==0? view_mean(L)[8] : view_mean(L)[10]) > 0.60f;
        case 9: return (L.kind==0? view_mean(L)[8] : view_mean(L)[10]) < 0.30f;
    }
    return true;
}

