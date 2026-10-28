// 動物の放流・スリープの防止・セーブの一覧・植物の形の数え上げ(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 遺伝子の適性に合う場所を探して動物を放つ
void release_animals(int per_species,int frame){
    if(last_animals.empty()){ printf("no stored genes\n"); return; }
    int total=0;
    for(const auto& g : last_animals){
        // この遺伝子にとって最も条件の良いマスを、無作為に探して選ぶ
        for(int n=0;n<per_species;n++){
            int bx=-1,by=-1; float bscore=-1e9f;
            for(int t=0;t<400;t++){
                int gx=(int)frand(0.f,(float)GRID_W);
                int gy=(int)frand(0.f,(float)GRID_H);
                int ci=idx(gx,gy);
                bool s=sea[ci];
                // 生息環境の不適合
                float hm = s ? (1.f-g.aquatic)*(SHALLOW_EASE+(1.f-SHALLOW_EASE)*depthn(ci))
                             : g.aquatic;
                if(hm>0.55f) continue;
                // 移動できない場所は除く
                if( s && g.aquatic<0.20f) continue;
                if(!s && g.aquatic>0.85f) continue;
                // 気温
                float kill=8.f-40.f*g.cold_tol;
                if(temper[ci]<kill+3.f) continue;
                float opt=26.f-26.f*g.cold_tol;
                float dT=(temper[ci]-opt)/25.f;
                float tfit=std::exp(-dT*dT);
                // 餌: 届く高さの植物と、食える死骸
                float food=0.f;
                for(int t=pgrid.b(ci);t<pgrid.e(ci);t++)
                    if(plants[pgrid.items[t]].height<=g.reach) food+=1.f;
                float detr=(g.detritus>0.05f)? corpse_density[ci]*0.5f : 0.f;
                float score = tfit*2.f + food + detr - hm*3.f;
                if(score>bscore){ bscore=score; bx=gx; by=gy; }
            }
            if(bx<0) continue;
            Animal a{};   // 脂肪・潜り具合などの状態は 0 から始める
            a.x=bx+frand(0.f,1.f); a.y=by+frand(0.f,1.f);
            a.dir=frand(0.f,6.2832f);
            a.energy=A_DIVIDE_TH*0.7f; a.alive=true; a.age=0; a.eaten=0;
            a.speed=g.speed; a.turn=g.turn; a.reach=g.reach; a.size=g.size;
            a.diet=g.diet;   a.jaw=g.jaw;   a.cold_tol=g.cold_tol;
            a.aquatic=g.aquatic; a.resp=g.resp; a.detritus=g.detritus;
            a.sexual=g.sexual;   a.arboreal=g.arboreal; a.perch=0.f; a.aerobic=1.f;
            a.mut_rate=0.3f;                    a.mother=-1; a.father=-1; a.gen=0; a.lesion=false;
            for(int k=0;k<N_REC;k++) a.rec[k]=0.f;
            a.brood=g.brood;     a.lifespan=g.lifespan;
            a.toxin=g.toxin;     a.tox_resist=g.tox_resist; a.burrow=g.burrow;
            a.nocturnal=g.nocturnal; a.store_fat=g.store_fat; a.repro_alloc=g.repro_alloc;
            a.bite_size=g.bite_size; a.digest=g.digest;
            for(int k=0;k<N_IMMUNE;k++) a.im[k]=g.im[k];
            a.body=g.size*MATURE_FRAC;      // 成熟した状態で放つ
            a.id=next_animal_id++; a.lin=g.lin;
            for(int k=0;k<N_W;k++) a.w[k]=g.w[k];
            for(int k=0;k<N_NODE;k++) a.ngain[k]=g.ngain[k];
            for(int k=0;k<N_SENSOR;k++) a.sen[k]=g.sen[k];
            animals.push_back(a); total++;
            g_bud.c_immig+=(a.energy+a.fat+a.body*CORPSE_RETURN)*CO2_PER_ENERGY;   // 外から放した動物が持ち込む炭素
        }
    }
    printf("=== released %d animals from %d species at t=%d ===\n",
           total,(int)last_animals.size(),frame);
}

// 長時間の放置中に画面や本体が眠らないようにする
#ifdef __APPLE__
static IOPMAssertionID g_sleep_assert = 0;
void prevent_sleep_begin(){
    CFStringRef reason = CFSTR("evosim long run");
    IOReturn r = IOPMAssertionCreateWithName(
        kIOPMAssertionTypeNoDisplaySleep, kIOPMAssertionLevelOn,
        reason, &g_sleep_assert);
    printf(r==kIOReturnSuccess ? "sleep prevented\n" : "sleep prevention failed\n");
}
void prevent_sleep_end(){
    if(g_sleep_assert) IOPMAssertionRelease(g_sleep_assert);
}
#elif defined(_WIN32)
void prevent_sleep_begin(){
    if(SetThreadExecutionState(ES_CONTINUOUS|ES_SYSTEM_REQUIRED|ES_DISPLAY_REQUIRED))
        printf("sleep prevented\n");
    else
        printf("sleep prevention failed\n");
}
void prevent_sleep_end(){
    SetThreadExecutionState(ES_CONTINUOUS);
}
#else
void prevent_sleep_begin(){}
void prevent_sleep_end(){}
#endif

// C++17 では file_time_type から time_t への標準の変換がないため、
// 現在時刻との差分を取って system_clock に移す
std::time_t file_time_to_t(std::filesystem::file_time_type ft){
    auto sys = std::chrono::time_point_cast<std::chrono::system_clock::duration>(
        ft - std::filesystem::file_time_type::clock::now()
           + std::chrono::system_clock::now());
    return std::chrono::system_clock::to_time_t(sys);
}

// 植物の姿。遺伝子の組み合わせから、どの型に属するかを決める
const int N_FORM = 8;
const char* FORM_NAME[N_FORM]={
    "プランクトン","寄生植物","食虫植物","着生植物","つる","竹型","樹木","草本"};
inline int plant_form(const Plant& p){
    if(p.aquatic > 0.50f) return 0;
    if(p.parasite  > 0.35f) return 1;
    if(p.carnivory > 0.30f) return 2;
    if(p.climb>0.45f && p.max_height<0.8f) return 3;   // 登るが幹を持たない
    if(p.climb     > 0.45f) return 4;
    if(p.clonal    > 0.50f) return 5;
    if(p.max_height> 1.50f) return 6;
    return 7;
}
// 型ごとの、世界に現存する最大と最小の背丈
float form_min[N_FORM], form_max[N_FORM];
int   form_cnt[N_FORM];
void scan_forms(){
    for(int k=0;k<N_FORM;k++){ form_min[k]=1e9f; form_max[k]=-1e9f; form_cnt[k]=0; }
    for(const auto& p:plants){
        if(!p.alive) continue;
        int f=plant_form(p);
        if(p.max_height<form_min[f]) form_min[f]=p.max_height;
        if(p.max_height>form_max[f]) form_max[f]=p.max_height;
        form_cnt[f]++;
    }
    for(int k=0;k<N_FORM;k++) if(form_cnt[k]==0){ form_min[k]=0.f; form_max[k]=1.f; }
}

