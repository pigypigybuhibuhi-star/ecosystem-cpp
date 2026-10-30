// セーブ形式 v2(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ===== セーブ形式 v2 =====
// 各レコードの並びを、ファイル冒頭に名前で宣言する。
// 読み込みは名前で対応づけるので、遺伝子を足しても消しても古いセーブが読める。
//   足した遺伝子 → 既定値で埋める
//   消した遺伝子 → 読み飛ばす
//   並べ替え    → 影響なし

template<class T>
struct Fld {
    std::string name;
    std::function<double(const T&)> get;
    std::function<void(T&,double)> set;
    double def;
};
// 普通のメンバ(float / int / bool)
template<class T,class M>
Fld<T> F(const char* n, M T::* p, double d){
    return { n, [p](const T& o){ return (double)(o.*p); },
                [p](T& o,double v){ o.*p=(M)v; }, d };
}
// 配列メンバの一要素
template<class T,class M,size_t N>
Fld<T> FA(const std::string& n, M (T::*p)[N], size_t k, double d){
    return { n, [p,k](const T& o){ return (double)((o.*p)[k]); },
                [p,k](T& o,double v){ (o.*p)[k]=(M)v; }, d };
}

// --- 各レコードの定義。遺伝子を足すときは、ここに一行足すだけ ---
std::vector<Fld<Plant>> plant_schema(){
    std::vector<Fld<Plant>> s={
        F("x",&Plant::x,0), F("y",&Plant::y,0),
        F("energy",&Plant::energy,1), F("height",&Plant::height,0.05),
        F("absorb",&Plant::absorb,0.9), F("max_height",&Plant::max_height,0.5),
        F("shade_tol",&Plant::shade_tol,0.1), F("tough",&Plant::tough,0),
        F("clonal",&Plant::clonal,0), F("disperse",&Plant::disperse,0.2),
        F("cold_tol",&Plant::cold_tol,0.05), F("aquatic",&Plant::aquatic,0),
        F("salt_tol",&Plant::salt_tol,0),
        F("drought_tol",&Plant::drought_tol,0), F("dorm_temp",&Plant::dorm_temp,-15),
        F("dormant",&Plant::dormant,0), F("store_cap",&Plant::store_cap,0),
        F("water_store",&Plant::water_store,0), F("eff_height",&Plant::eff_height,0.05),
        F("fire_tol",&Plant::fire_tol,0), F("climb",&Plant::climb,0),
        F("sexual",&Plant::sexual,0.2), F("brood",&Plant::brood,0.1),
        F("n_fix",&Plant::n_fix,0), F("n_store",&Plant::n_store,0.5),
        F("lifespan",&Plant::lifespan,0.05), F("parasite",&Plant::parasite,0),
        F("carnivory",&Plant::carnivory,0), F("fire_seed",&Plant::fire_seed,0),
        F("repro_alloc",&Plant::repro_alloc,0.5), F("seed_bank",&Plant::seed_bank,0),
        F("defense",&Plant::defense,0), F("seed_wait",&Plant::seed_wait,0),
        F("id",&Plant::id,0), F("genet",&Plant::genet,0),
        F("age",&Plant::age,0), F("lin",&Plant::lin,-1),
    };
    // 後から足した遺伝子。抜けていると読み込み時に未初期化になる
    s.push_back(F("seed_bank",&Plant::seed_bank,0));
    s.push_back(F("buoyancy",&Plant::buoyancy,0.6));
    s.push_back(F("mut_rate",&Plant::mut_rate,0.3));
    s.push_back(F("mother",&Plant::mother,-1));
    s.push_back(F("father",&Plant::father,-1));
    s.push_back(F("gen",&Plant::gen,0));
    for(int k=0;k<N_IMMUNE;k++) s.push_back(FA("im"+std::to_string(k),&Plant::im,k,0.5));
    return s;
}
std::vector<Fld<Animal>> animal_schema(){
    std::vector<Fld<Animal>> s={
        F("x",&Animal::x,0), F("y",&Animal::y,0), F("dir",&Animal::dir,0),
        F("energy",&Animal::energy,1), F("speed",&Animal::speed,0.5),
        F("turn",&Animal::turn,0.2), F("reach",&Animal::reach,0.7),
        F("size",&Animal::size,0.5), F("diet",&Animal::diet,0),
        F("jaw",&Animal::jaw,1), F("cold_tol",&Animal::cold_tol,0.3),
        F("body",&Animal::body,0.1), F("aquatic",&Animal::aquatic,0),
        F("salt_tol",&Animal::salt_tol,0),
        F("resp",&Animal::resp,0.3), F("detritus",&Animal::detritus,0),
        F("sexual",&Animal::sexual,0.2), F("arboreal",&Animal::arboreal,0),
        F("perch",&Animal::perch,0), F("brood",&Animal::brood,0.1),
        F("lifespan",&Animal::lifespan,0.2), F("toxin",&Animal::toxin,0),
        F("tox_resist",&Animal::tox_resist,0), F("burrow",&Animal::burrow,0),
        F("nocturnal",&Animal::nocturnal,0), F("store_fat",&Animal::store_fat,0.2),
        F("repro_alloc",&Animal::repro_alloc,0.5), F("fat",&Animal::fat,0),
        F("depth",&Animal::depth,0),
        F("id",&Animal::id,0), F("age",&Animal::age,0),
        F("eaten",&Animal::eaten,0), F("lin",&Animal::lin,-1),
    };
    for(int k=0;k<N_W;k++) s.push_back(FA("w"+std::to_string(k),&Animal::w,k,0));
    for(int k=0;k<N_NODE;k++)
        s.push_back(FA("ng"+std::to_string(k),&Animal::ngain,k,k<3?1.0:0.05));
    s.push_back(F("mut_rate",&Animal::mut_rate,0.3));
    s.push_back(F("bite_size",&Animal::bite_size,0.5));
    s.push_back(F("digest",&Animal::digest,0.5));
    s.push_back(F("mother",&Animal::mother,-1));
    s.push_back(F("father",&Animal::father,-1));
    s.push_back(F("gen",&Animal::gen,0));
    for(int k=0;k<N_REC;k++) s.push_back(FA("rec"+std::to_string(k),&Animal::rec,k,0));
    s.push_back(F("lesion",&Animal::lesion,0));
    for(int k=0;k<N_SENSOR;k++){
        std::string b="sen"+std::to_string(k)+"_";
        s.push_back({b+"angle",[k](const Animal& a){return (double)a.sen[k].angle;},
                               [k](Animal& a,double v){a.sen[k].angle=(float)v;},0});
        s.push_back({b+"range",[k](const Animal& a){return (double)a.sen[k].range;},
                               [k](Animal& a,double v){a.sen[k].range=(float)v;},0});
        s.push_back({b+"tune", [k](const Animal& a){return (double)a.sen[k].tune;},
                               [k](Animal& a,double v){a.sen[k].tune=(float)v;},0});
    }
    for(int k=0;k<N_IMMUNE;k++) s.push_back(FA("im"+std::to_string(k),&Animal::im,k,0.5));
    return s;
}
std::vector<Fld<Microbe>> microbe_schema(){
    return { F("x",&Microbe::x,0), F("y",&Microbe::y,0),
             F("energy",&Microbe::energy,0.4), F("rate",&Microbe::rate,0.3),
             F("cold_tol",&Microbe::cold_tol,0.2),
             F("id",&Microbe::id,0), F("age",&Microbe::age,0) };
}
std::vector<Fld<Corpse>> corpse_schema(){
    return { F("x",&Corpse::x,0), F("y",&Corpse::y,0),
             F("energy",&Corpse::energy,0), F("origin",&Corpse::origin,1),
             F("n",&Corpse::n,0) };
}
std::vector<Fld<Lineage>> lineage_schema(){
    std::vector<Fld<Lineage>> s={
        F("id",&Lineage::id,0), F("kind",&Lineage::kind,0), F("num",&Lineage::num,0),
        F("nrep",&Lineage::nrep,0), F("birth",&Lineage::birth_tick,0),
        F("death",&Lineage::death_tick,-1), F("parent",&Lineage::parent,-1),
        F("alive",&Lineage::alive,1), F("pop",&Lineage::pop,0),
        F("max_pop",&Lineage::max_pop,0), F("missing",&Lineage::missing,0),
    };
    for(int k=0;k<REP_N;k++) s.push_back(FA("rep"+std::to_string(k),&Lineage::rep,k,0));
    for(int k=0;k<REP_N;k++) s.push_back(FA("mean"+std::to_string(k),&Lineage::mean,k,0));
    for(int k=0;k<REP_N;k++)
        s.push_back(FA("fmean"+std::to_string(k),&Lineage::final_mean,k,0));
    s.push_back(F("final_pop",&Lineage::final_pop,0));
    return s;
}
// マスはインデックスで引く
std::vector<Fld<int>> cell_schema(){
    std::vector<Fld<int>> s;
    auto fv=[&](const std::string& n,std::vector<float>& v,double d){
        s.push_back({n,[&v](const int& i){return (double)v[i];},
                       [&v](int& i,double x){v[i]=(float)x;},d}); };
    auto iv=[&](const std::string& n,std::vector<int>& v,double d){
        s.push_back({n,[&v](const int& i){return (double)v[i];},
                       [&v](int& i,double x){v[i]=(int)x;},d}); };
    fv("elev",elev,0.4); fv("co2_dev",co2,0); fv("o2_dev",o2,0);
    fv("co2_res",co2_res,0); fv("o2_res",o2_res,0);
    fv("soil_n",soil_n,SOIL_N_BASE); fv("soil_water",soil_water,0.4);
    fv("gw_h",gw_h,0); fv("salt_soil",salt_soil,0); fv("salt_gw",salt_gw,0);
    fv("pload",pload,0); fv("aload",aload,0);
    iv("burn_age",burn_age,999999);
    fv("fuel_moist",fuel_moist,-1);
    for(int k=0;k<N_ODOR;k++) fv("odor"+std::to_string(k),odor[k],0);
    for(int k=0;k<N_IMMUNE;k++){
        fv("pstrain"+std::to_string(k),pstrain[k],0.5);
        fv("astrain"+std::to_string(k),astrain[k],0.5);
    }
    return s;
}

// 書き出しと読み込みを同じ関数で記述するための器。
// フィールドの並びが一箇所にしかないので、片方だけ直し忘れることがない
struct SaveW {
    std::ofstream& f;
    template<class T> void operator()(T& v){ f << v << ' '; }
    void operator()(bool& v){ f << (v?1:0) << ' '; }
    void nl(){ f << '\n'; }
};
struct SaveR {
    std::ifstream& f;
    template<class T> void operator()(T& v){ f >> v; }
    void operator()(bool& v){ int t=0; f >> t; v=(t!=0); }
    void nl(){}
};

// --- 書き出しと読み込みの汎用部 ---
template<class T>
void write_schema(std::ostream& f,const char* tag,const std::vector<Fld<T>>& sc){
    f<<"SCHEMA "<<tag<<" "<<sc.size();
    for(auto& d:sc) f<<' '<<d.name;
    f<<'\n';
}
template<class T>
void write_rec(std::ostream& f,const std::vector<Fld<T>>& sc,const T& o){
    for(size_t k=0;k<sc.size();k++){ if(k) f<<' '; f<<sc[k].get(o); }
    f<<'\n';
}
// ファイルの列を、今のフィールドに対応づける
struct ColMap { std::vector<int> col; std::vector<char> seen; };
template<class T>
ColMap build_map(const std::vector<std::string>& names,
                 const std::vector<Fld<T>>& sc,const char* tag){
    ColMap m; m.seen.assign(sc.size(),0);
    std::unordered_map<std::string,int> by;
    for(size_t k=0;k<sc.size();k++) by[sc[k].name]=(int)k;
    std::string extra;
    for(auto& nm:names){
        auto it=by.find(nm);
        int fi=(it==by.end())?-1:it->second;
        m.col.push_back(fi);
        if(fi>=0) m.seen[fi]=1; else extra+=" "+nm;
    }
    std::string miss;
    for(size_t k=0;k<sc.size();k++) if(!m.seen[k]) miss+=" "+sc[k].name;
    if(!miss.empty())  printf("  %s: filled with defaults:%s\n",tag,miss.c_str());
    if(!extra.empty()) printf("  %s: ignored (no longer used):%s\n",tag,extra.c_str());
    return m;
}
template<class T>
void read_rec(std::istream& f,const std::vector<Fld<T>>& sc,const ColMap& m,T& o){
    for(auto& d:sc) d.set(o,d.def);            // まず既定値で埋める
    for(int fi:m.col){ double v; f>>v; if(fi>=0) sc[fi].set(o,v); }
}

// 乱数の内部状態。これが無いと、ロード後の乱数列が元の続きにならず、
// 同じ地点から再開しても別の歴史になる
template<class S> void io_rng(S& s, std::ostream* w, std::istream* r){
    if(w){ *w << rng << '\n'; }
    if(r){ *r >> rng; }
}

// グラフの履歴。世界の状態ではないが、変遷を見せるために残す
template<class S> void io_hist(S& s){
    auto fv=[&](std::vector<float>& v){
        int n=(int)v.size(); s(n);
        v.resize(n);
        for(int i=0;i<n;i++) s(v[i]);
        s.nl();
    };
    auto iv=[&](std::vector<int>& v){
        int n=(int)v.size(); s(n);
        v.resize(n);
        for(int i=0;i<n;i++) s(v[i]);
        s.nl();
    };
    fv(h_maxH); fv(h_reach); fv(h_diet); fv(h_shade); fv(h_tough); fv(h_jaw);
    fv(h_coldP); fv(h_coldA); fv(h_co2); fv(h_o2); fv(h_temp);
    fv(h_soil); fv(h_nut); fv(h_ice); fv(h_path); fv(h_odor); fv(h_sst);
    iv(h_plant); iv(h_animal); iv(h_corpse);
    Series* ss[13]={&s_plant,&s_animal,&s_maxh,&s_reach,&s_temp,
                    &s_co2,&s_o2,&s_soil,&s_nut,
                    &s_shade,&s_diet,&s_coldP,&s_coldA};
    for(int k=0;k<13;k++)
        for(int sc=0;sc<N_SCALE;sc++){
            fv(ss[k]->v[sc]);
            s(ss[k]->acc[sc]); s(ss[k]->cnt[sc]); s.nl();
        }
}

void save_state(const std::string& fn,int frame){
    auto t0=std::chrono::steady_clock::now();
    std::string tmp=fn+".tmp";                 // 書き終えてから差し替える
    {
        std::ofstream f(tmp);
        if(!f){ printf("save failed: %s\n",fn.c_str()); return; }
        f.precision(8);
        auto cs=cell_schema(); auto ps=plant_schema(); auto as=animal_schema();
        auto ks=corpse_schema(); auto ls=lineage_schema();
        auto ms=microbe_schema();
        SaveW sw{f};
        f<<"P3SAVE2\n";
        f<<"HEADER "<<frame<<' '<<next_plant_id<<' '<<next_animal_id<<' '
         <<next_lineage_id<<' '<<next_plant_num<<' '<<next_animal_num<<' '
         <<GRID_W<<' '<<GRID_H<<'\n';
        write_schema(f,"cell",cs);   write_schema(f,"plant",ps);
        write_schema(f,"animal",as); write_schema(f,"microbe",ms);
        write_schema(f,"corpse",ks);
        write_schema(f,"lineage",ls);
        f<<"CELLS "<<GRID_W*GRID_H<<'\n';
        for(int i=0;i<GRID_W*GRID_H;i++) write_rec(f,cs,i);
        f<<"PLANTS "<<plants.size()<<'\n';
        for(auto& p:plants) write_rec(f,ps,p);
        f<<"ANIMALS "<<animals.size()<<'\n';
        for(auto& a:animals) write_rec(f,as,a);
        f<<"MICROBES "<<microbes.size()<<'\n';
        for(auto& m:microbes) write_rec(f,ms,m);
        f<<"CORPSES "<<corpses.size()<<'\n';
        for(auto& c:corpses) write_rec(f,ks,c);
        f<<"LINEAGES "<<lineages.size()<<'\n';
        for(auto& L:lineages) write_rec(f,ls,L);
        f<<"RNG\n"<<rng<<'\n';
        // 天気の状態。雨雲は 1 時間ほどしか続かないが、再開しても同じ続きになるよう乱数ごと残す
        f<<"WEATHER "<<storms.size()<<'\n'<<wx_rng<<'\n';
        for(auto& s:storms){ char b[160]; snprintf(b,sizeof b,"%a %a %a %a\n",s.x,s.y,s.age,s.ipeak); f<<b; }
        f<<"LIGHTNING\n"<<lt_rng<<'\n'<<lightning_events<<'\n';
        // 地震の状態。余震は本震から長く続くので、本震の一覧と乱数を残す
        f<<"QUAKE "<<mainshocks.size()<<'\n'<<qk_rng<<'\n'<<quake_events<<'\n';
        for(auto& s:mainshocks){ char b[160]; snprintf(b,sizeof b,"%ld %a %a\n",s.tick,s.mag,s.x); f<<b; }
        f<<"HISTORY\n";
        io_hist(sw);
        f<<"END\n";
        if(!f){ printf("save failed while writing\n"); return; }
    }
    std::remove(fn.c_str());
    if(std::rename(tmp.c_str(),fn.c_str())!=0){
        printf("save failed: could not replace %s\n",fn.c_str()); return;
    }
    double ms=std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-t0).count();
    printf("saved: %s (plants %d animals %d lineages %d, %.0f ms)\n",fn.c_str(),
           (int)plants.size(),(int)animals.size(),(int)lineages.size(),ms);
}

bool load_state(const std::string& fn,int& frame){
    std::ifstream f(fn);
    if(!f){ printf("load failed (no file): %s\n",fn.c_str()); return false; }
    std::string tag; f>>tag;
    if(tag!="P3SAVE2"){
        printf("load failed: '%s' is the old format and cannot be read\n",fn.c_str());
        return false;
    }
    try {
        SaveR sr{f};
        std::string sec;
        int gw=0,gh=0;
        f>>sec;
        if(sec!="HEADER") throw std::runtime_error("missing HEADER");
        f>>frame>>next_plant_id>>next_animal_id>>next_lineage_id
         >>next_plant_num>>next_animal_num>>gw>>gh;
        if(gw!=GRID_W||gh!=GRID_H) throw std::runtime_error("grid size differs");

        auto cs=cell_schema(); auto ps=plant_schema(); auto as=animal_schema();
        auto ks=corpse_schema(); auto ls=lineage_schema();
        auto ms=microbe_schema();
        ColMap mc,mp,ma,mk,ml,mm;
        auto read_schema=[&](){
            std::string s,kind; int n; f>>s>>kind>>n;
            if(s!="SCHEMA") throw std::runtime_error("expected SCHEMA, got "+s);
            std::vector<std::string> names(n);
            for(auto& x:names) f>>x;
            if(kind=="cell")    mc=build_map(names,cs,"cell");
            else if(kind=="plant")  mp=build_map(names,ps,"plant");
            else if(kind=="animal") ma=build_map(names,as,"animal");
            else if(kind=="microbe")mm=build_map(names,ms,"microbe");
            else if(kind=="corpse") mk=build_map(names,ks,"corpse");
            else if(kind=="lineage")ml=build_map(names,ls,"lineage");
        };
        for(int k=0;k<6;k++) read_schema();

        auto want=[&](const char* name)->long{
            std::string s; long n=-1; f>>s>>n;
            if(!f||s!=name) throw std::runtime_error(std::string("expected ")+name+", got "+s);
            return n;
        };
        long n=want("CELLS");
        for(int i=0;i<n;i++){ int ii=i; read_rec(f,cs,mc,ii); }
        for(int i=0;i<GRID_W*GRID_H;i++) sea[i]=(elev[i]<SEA_LEVEL)?1:0;

        n=want("PLANTS"); plants.clear(); plants.reserve(n);
        for(long i=0;i<n;i++){
            Plant p; read_rec(f,ps,mp,p);
            if(!f) throw std::runtime_error("plant read failed at "+std::to_string(i));
            p.alive=true; plants.push_back(p);
        }
        n=want("ANIMALS"); animals.clear(); animals.reserve(n);
        for(long i=0;i<n;i++){
            Animal a; read_rec(f,as,ma,a);
            if(!f) throw std::runtime_error("animal read failed at "+std::to_string(i));
            a.alive=true; a.aerobic=1.f; animals.push_back(a);
        }
        n=want("MICROBES"); microbes.clear(); microbes.reserve(n);
        for(long i=0;i<n;i++){
            Microbe m; read_rec(f,ms,mm,m);
            m.alive=true; microbes.push_back(m);
        }
        n=want("CORPSES"); corpses.clear(); corpses.reserve(n);
        for(long i=0;i<n;i++){ Corpse c; read_rec(f,ks,mk,c); corpses.push_back(c); }
        n=want("LINEAGES"); lineages.clear(); lineages.reserve(n);
        for(long i=0;i<n;i++){ Lineage L; read_rec(f,ls,ml,L); lineages.push_back(L); }
        f>>sec;
        if(sec=="RNG"){ f>>rng; f>>sec; }
        storms.clear();
        if(sec=="WEATHER"){
            long n; f>>n; f>>wx_rng;
            for(long i=0;i<n;i++){ std::string a,b,c,d; f>>a>>b>>c>>d;
                storms.push_back({std::strtof(a.c_str(),nullptr),std::strtof(b.c_str(),nullptr),
                                  std::strtof(c.c_str(),nullptr),std::strtof(d.c_str(),nullptr)}); }
            f>>sec;
        }
        if(sec=="LIGHTNING"){ f>>lt_rng; f>>lightning_events; f>>sec; }
        mainshocks.clear();
        if(sec=="QUAKE"){
            long n; f>>n; f>>qk_rng; f>>quake_events;
            for(long i=0;i<n;i++){ long t; std::string a,b; f>>t>>a>>b;
                mainshocks.push_back({t,std::strtof(a.c_str(),nullptr),std::strtof(b.c_str(),nullptr)}); }
            f>>sec;
        }
        if(sec=="HISTORY"){
            io_hist(sr);
            f>>sec;
        }
        if(sec!="END") throw std::runtime_error("END marker missing — file truncated?");

        // 地形から導かれるものを作り直す
        elev_smooth=elev; box_blur(elev_smooth,4,2);
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++) slope_at(x,y);
        // 沿岸の近さを作り直す
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
            coast_w[i]=best*best;
        }
        compute_currents();
        compute_rivers();
        update_climate(frame);
        compute_precip();
        classify_climate();
        pgrid.reset();
        for(const auto& p:plants){
            int gx=(int)p.x, gy=(int)p.y;
            if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) pgrid.count(idx(gx,gy));
        }
        pgrid.finalize();
        for(int i=0;i<(int)plants.size();i++){
            int gx=(int)plants[i].x, gy=(int)plants[i].y;
            if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) pgrid.put(idx(gx,gy),i);
        }
        printf("loaded: %s (plants %d animals %d lineages %d)\n",fn.c_str(),
               (int)plants.size(),(int)animals.size(),(int)lineages.size());
        return true;
    } catch(const std::exception& e){
        printf("load failed: %s\n",e.what());
        printf("  -> world cleared.\n");
        plants.clear(); animals.clear(); corpses.clear(); lineages.clear();
        return false;
    }
}

