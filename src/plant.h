// 植物(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 海の浮く植物の始まりの遺伝子。最初に置く株と、風上の端から入ってくる区画の外の株(条件)に使う
void sea_plant_genes(Plant& p){
    p.energy=PLANT_INIT_E; p.height=0.05f; p.alive=true; p.age=0;
    p.absorb=0.92f; p.max_height=0.25f; p.shade_tol=0.10f;
    p.tough=0.0f; p.clonal=0.0f; p.disperse=0.2f; p.cold_tol=0.05f;
    p.aquatic=1.0f; p.salt_tol=1.0f; p.drought_tol=0.0f; p.dorm_temp=-15.0f; p.dormant=false;
    p.store_cap=0.0f; p.water_store=0.0f; p.eff_height=0.05f;
    p.fire_tol=0.0f; p.climb=0.0f; p.sexual=0.2f;
    for(int k=0;k<N_IMMUNE;k++) p.im[k]=frand(0.f,1.f);
    p.brood=0.1f; p.n_fix=0.0f; p.n_store=1.0f; p.lifespan=0.05f;
    p.parasite=0.f; p.carnivory=0.f; p.fire_seed=0.f; p.repro_alloc=0.5f;
    p.seed_bank=0.f; p.def_const=0.5f; p.defense=0.f; p.seed_wait=0; p.buoyancy=0.6f;
    p.mut_rate=0.3f;        p.mother=-1; p.father=-1; p.gen=0;
    p.id=next_plant_id++; p.genet=p.id; p.lin=-1;
}
// 風上の端から入ってくる浮く植物。東西は開いた境界なので、流されて出る流れと同じく、入る流れがある。
// 区画の外の海の浮く植物の密度は条件で、最初に区画の海に置いた株の、海全体の平均の密度と同じとする。
// 1 行の端を 1 分に越える数 = 密度 [株/マス] × 漂流の速さ [マス/分] × 浮く強さ。持ち込む炭素と窒素は c_immig・n_immig に付ける
void plant_inflow(){
    static double dens=-1.0;
    if(dens<0.0){ long ns=0; for(int i=0;i<GRID_W*GRID_H;i++) if(sea[i]) ns++;
                  dens=ns?(double)(PLANT_INIT_N-PLANT_INIT_LAND)/ns:0.0; }
    static std::vector<int> edge; edge.clear();
    for(int y=0;y<GRID_H;y++){
        if(wind[y]==0.f) continue;
        int x0=(wind[y]>0.f)?0:GRID_W-1;   // 風上の端
        int ci=idx(x0,y);
        if(sea[ci]&&!ice[ci]) edge.push_back(ci);
    }
    if(edge.empty()) return;
    const float b0=0.6f;                   // 外の株の浮く強さ(sea_plant_genes と同じ)
    int n=count_by_rate(dens*SURF_DRIFT_CPM*DT0*b0*edge.size());
    for(int k=0;k<n && (int)plants.size()<PLANT_MAX;k++){
        int ci=edge[(int)(dist01(rng)*edge.size())%edge.size()];
        int gx=ci%GRID_W, gy=ci/GRID_W;
        Plant p{};
        // この分の間に端から流れ込んだ距離の中に置く
        float d=dist01(rng)*(float)(SURF_DRIFT_CPM*DT*b0);
        p.x=(gx==0)?d:(float)GRID_W-1e-3f-d; p.y=gy+dist01(rng);
        if(p.x<0.f) p.x=0.f; if(p.x>=GRID_W) p.x=GRID_W-1e-3f;
        sea_plant_genes(p);
        plants.push_back(p);
        g_bud.c_immig+=(p.energy+p.height*GROW_COST)*CO2_PER_ENERGY;
        g_bud.n_immig+=p.n_store+p.height*N_PER_GROWTH;
    }
}
void update_plants(){
    pgrid.reset();
    for(const auto& p:plants){
        if(!p.alive||p.seed_wait>0) continue;   // 種子は光も遮らず、食われもしない
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        pgrid.count(idx(gx,gy));
    }
    pgrid.finalize();
    for(int i=0;i<(int)plants.size();i++){
        if(!plants[i].alive||plants[i].seed_wait>0) continue;
        int gx=(int)plants[i].x, gy=(int)plants[i].y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        pgrid.put(idx(gx,gy),i);
        fuel[idx(gx,gy)]+=plants[i].height*1.2f;
    }
    std::fill(canopy_h.begin(),canopy_h.end(),0.f);
    auto t_a=std::chrono::steady_clock::now();
    // マスごとに独立しているので、マスを分けて並列に処理できる。
    // order はスレッドごとに持つ必要がある(static のままだと競合する)
    parallel_for(0,GRID_W*GRID_H,[&](int cs,int cee){
    std::vector<int> order;            // スレッドごとに持つ。static だと競合する
    double t_gain=0,t_light=0,t_tf=0,t_wf=0,t_co2=0; long t_n=0;
    for(int ci=cs;ci<cee;ci++){
        int b0=pgrid.start[ci], e0=pgrid.start[ci+1];
        if(b0==e0) continue;
        const int* use=&pgrid.items[b0];
        int cnt=e0-b0;
        // このマスで最も高い株。つるはこれを支えにして登る
        float tallest=0.f;
        for(int t=b0;t<e0;t++)
            if(plants[pgrid.items[t]].height>tallest)
                tallest=plants[pgrid.items[t]].height;
        for(int t=b0;t<e0;t++){
            Plant& q=plants[pgrid.items[t]];
            float support=std::max(0.f,tallest-q.height);
            // つるは自立する幹を作らない。支えがなければ地面を這うしかない
            float self=q.height*(1.f-CLIMB_PENALTY*q.climb);
            q.eff_height=self+support*q.climb*CLIMB_GAIN;
            if(q.eff_height>tallest) tallest=q.eff_height;   // つるも林冠を作る
        }
        canopy_h[ci]=tallest;
        if(cnt>1){                               // 1株なら背比べは要らない
            order.assign(use,use+cnt);
            std::sort(order.begin(),order.end(),[](int a,int b){
                return plants[a].eff_height>plants[b].eff_height; });
            use=order.data();
        }
        float remaining=solar[ci];
        if(ice[ci]) remaining *= ICE_LIGHT;     // 氷が光を遮る
        for(int t=0;t<cnt;t++){
            if(remaining<=0.f) break;
            int pi=use[t];
            Plant& p=plants[pi];
            if(p.dormant) continue;
            // 浮く個体は水面近くに留まり、深さに関わらず光を受ける
            float att=1.f;
            if(sea[ci]) att=std::exp(-WATER_ATTEN*depthn(ci)*(1.f-p.buoyancy));
            float got=std::min(remaining*att,saturation(p)*SAT_FLUX);
            remaining -= got/std::max(0.05f,att);
            float dT=(temper[ci]-opt_temp_p(p))/TEMP_WIDTH_P;
            float tf=std::exp(-dT*dT);
            float gain=got*p.absorb*tf;
            // 水: 耐乾性が高いほど乏しい水で足りるが、装置に維持費がかかる
            float wneed = 0.55f*(1.f-0.75f*p.drought_tol);
            // 土壌が乏しければ、体に貯めた水で補う
            float avail_w = soil_water[ci]
                          + p.water_store/std::max(0.05f,wneed)*0.25f;
            float wf = clamp01(avail_w/std::max(0.02f,wneed));
            gain *= wf;
            if(soil_water[ci]<wneed && p.water_store>0.f){
                float draw=std::min(p.water_store,(wneed-soil_water[ci])*STORE_DRAW);
                p.water_store-=draw;
                w_transp_cell[ci]-=draw;                  // 貯水は蒸散で大気へ出る
            }
            { float w0=soil_water[ci];
              soil_water[ci] = std::max(0.f, soil_water[ci]-gain*0.010f);  // 蒸散
              if(!sea[ci]) w_transp_cell[ci]+=soil_water[ci]-w0; }
            // 炭素(リービッヒの最小律)
            float need=gain*CO2_PER_ENERGY;
            float avail=std::min(need,co2_at(ci));
            if(need>1e-9f) gain*=avail/need;
            co2_add(ci,-avail); o2_add(ci,avail*O2_PER_CO2);   // マスごとに一つのスレッドだけが触るので、端数の持ち越しも決定的
            p.energy+=gain;
            t_gain+=gain; t_n++;
            t_light+=got; t_tf+=tf; t_wf+=wf;
            t_co2+=(need>1e-9f? avail/need : 1.f);
        }
    }
    // 各スレッドの集計を、最後に一度だけ足し合わせる
    { static std::mutex dbg_mx;
      std::lock_guard<std::mutex> lk(dbg_mx);
      dbg_gain+=t_gain; dbg_n+=t_n;
      dbg_light+=t_light; dbg_tf+=t_tf; dbg_wf+=t_wf; dbg_co2+=t_co2; }
    });
    { auto t_b=std::chrono::steady_clock::now();
      static double acc=0; static int cnt=0;
      acc+=std::chrono::duration<double,std::milli>(t_b-t_a).count();
      if(++cnt>=60){ printf("  [prof] 光合成ループ %.1f ms\n",acc/cnt); acc=0; cnt=0; } }
    // 寄生: 吸器を同じマスの最も大きな株に刺し、エネルギーを直接奪う。
    // 宿主の tough が防御として働く。吸器が宿主の維管束とつながる面の広さは吸器への投資 parasite に
    // 比例するので、奪う率は parasite × 刺さり具合。閾値は置かず、寄生者も宿主になりうる
    static int para_tick=0;
    if(++para_tick%PARA_INTERVAL==0)
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        int b0=pgrid.start[ci], e0=pgrid.start[ci+1];
        if(e0-b0<2) continue;
        int host=-1; float best=-1.f;
        for(int t=b0;t<e0;t++){
            const Plant& q=plants[pgrid.items[t]];
            if(q.energy>best){ best=q.energy; host=pgrid.items[t]; }
        }
        if(host<0) continue;
        for(int t=b0;t<e0;t++){
            int pi=pgrid.items[t];
            if(pi==host) continue;
            Plant& q=plants[pi];
            if(q.parasite<=0.f) continue;
            float pierce=q.parasite/(q.parasite+plants[host].defense);
            float drain=plants[host].energy*-std::expm1(-PARASITE_DRAIN_K*q.parasite*pierce*(float)(PARA_INTERVAL*DT));
            plants[host].energy-=drain;
            q.energy+=drain*0.85f;
            a_burn(co2.data(),o2.data(),ci,drain*0.15f);   // 同化しなかった 15% は寄生者の呼吸で CO2 に(O2 を 1:1)
            // 窒素も宿主の師管液から奪う。宿主の持つ分と、自分の空きまで
            float nd=std::min({plants[host].n_store,drain*0.02f,std::max(0.f,N_STORE_MAX-q.n_store)});
            plants[host].n_store-=nd; q.n_store+=nd;
        }
    }
    static int genet_tick=0;
    if(++genet_tick%GENET_INTERVAL==0){   // 地下茎の融通は10tickに一度で足りる
        // 各株は地下茎の通りやすさ a = 1−exp(−k·clonal·n) で共通の節につながる。節は蓄えを持たないので、
        // その値は a で重みを付けた平均になり(キルヒホッフの電流則)、Σ a·(節 − E) = 0 で genet の合計が保たれる
        struct GenetSum{ double aE=0, a=0; int n=0; };
        static std::unordered_map<int,GenetSum> genet;
        genet.clear();
        // 閾値は置かない。地下茎への投資 clonal が小さければ通りやすさ a も小さく、融通は連続に小さくなる
        for(const auto& p:plants){
            if(!p.alive||p.clonal<=0.f) continue;
            double a=-std::expm1(-GENET_SHARE_K*p.clonal*(float)(GENET_INTERVAL*DT));
            auto& e=genet[p.genet]; e.aE+=a*p.energy; e.a+=a; e.n++;
        }
        for(auto& p:plants){
            if(!p.alive||p.clonal<=0.f) continue;
            auto it=genet.find(p.genet);
            if(it==genet.end()||it->second.n<2||it->second.a<=0.0) continue;
            double a=-std::expm1(-GENET_SHARE_K*p.clonal*(float)(GENET_INTERVAL*DT));
            double node=it->second.aE/it->second.a;
            p.energy=(float)(p.energy+a*(node-p.energy));
        }
    }
    static int pphase=0;
    pphase=(pphase+1)%PLANT_STRIDE;
    std::vector<Plant> babies;
    int pidx=-1;
    for(auto& p:plants){
        pidx++;
        if(!p.alive) continue;
        // 土中で待つ種子。光合成もせず、食われもせず、ほぼ無費用で耐える
        if(p.seed_wait>0){
            // 待機中の種子は数が多い。8tickに一度だけ、まとめて処理する
            if((pidx+pphase)%SEED_STRIDE!=0) continue;
            p.seed_wait-=SEED_STRIDE;
            // 種子の呼吸。使った炭素は CO2 になり、O2 を 1:1 で使う(体に無いぶんは出せない)
            { int sx=(int)p.x, sy=(int)p.y;
              if(sx>=0&&sx<GRID_W&&sy>=0&&sy<GRID_H)
                  a_burn(co2.data(),o2.data(),idx(sx,sy),std::min(BANK_COST*SEED_STRIDE,std::max(0.f,p.energy))); }
            p.energy-=BANK_COST*SEED_STRIDE;
            // 死んだ種子は朽ちる。体の炭素と窒素は死骸に
            if(dist01(rng)<BANK_PREDATION || p.energy<=0.f){
                p.alive=false;
                corpses.push_back({p.x,p.y,std::max(0.f,p.energy)+p.height*GROW_COST,0,p.n_store+p.height*N_PER_GROWTH});
                continue;
            }
            if(p.seed_wait<=0){
                // 発芽の時だけ環境を見る
                int gx0=(int)p.x, gy0=(int)p.y;
                int ci0=(gx0>=0&&gx0<GRID_W&&gy0>=0&&gy0<GRID_H)?idx(gx0,gy0):0;
                if(soil_water[ci0]<0.25f){ p.seed_wait=SEED_RETRY; continue; }
                p.seed_wait=0; p.height=0.05f;
            }
            continue;
        }
        if((pidx+pphase)%PLANT_STRIDE!=0) continue;
        p.age+=PLANT_STRIDE;   // 体の更新は PLANT_STRIDE tick に一度。その間に過ぎた分だけ歳をとる
        int gx=(int)p.x, gy=(int)p.y;
        int ci=(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H)?idx(gx,gy):0;
        // 休眠の切り替え。覚醒には余裕をもたせ、境界でばたつかせない
        if(p.dormant){ if(temper[ci]>p.dorm_temp+DORM_HYST) p.dormant=false; }
        else         { if(temper[ci]<p.dorm_temp)            p.dormant=true;  }

        if(temper[ci]<kill_temp_p(p)){
            p.alive=false;
            corpses.push_back({p.x,p.y,std::max(0.f,p.energy)+p.height*GROW_COST,0,p.n_store+p.height*N_PER_GROWTH});   // 体の炭素と窒素はすべて死骸に
            continue;
        }
        // 炎。厚い樹皮を持つ個体は耐える
        if(fire[ci]>0.05f){
            // 厚い樹皮でも完全には守れない。炎は幹を焼き、葉を奪う
            float dmg=fire[ci]*(1.f-p.fire_tol*0.65f);
            // 焼けた組織の炭素は燃焼で CO2 になり、死骸の燃焼と同じ割合で O2 を使う
            float h0=p.height;
            float burnt=std::min(dmg*FIRE_DAMAGE,std::max(0.f,p.energy));
            p.energy-=dmg*FIRE_DAMAGE;
            p.height=std::max(0.03f,p.height-dmg*FIRE_SHOOT_LOSS);
            burnt+=(h0-p.height)*GROW_COST;
            g_bud.n_fire-=(h0-p.height)*N_PER_GROWTH;           // 焼けた組織の窒素は揮発する
            if(p.energy<=0.f){
                p.alive=false;
                // 炭として残るぶんは、体の窒素も同じ割合で残す。残りは揮発する
                float nb=p.n_store+p.height*N_PER_GROWTH, kc=0.25f/GROW_COST;
                corpses.push_back({p.x,p.y,p.height*0.25f,0,nb*kc});   // 炭として残るぶん
                g_bud.n_fire-=nb*(1.f-kc);
                burnt+=p.height*(GROW_COST-0.25f);              // 残りの体は燃える
            }
            co2_add(ci,burnt*CO2_PER_ENERGY);
            o2_add(ci,-std::min(burnt*FIRE_O2_USE,o2_at(ci)));
            if(!p.alive) continue;
        }
        // 窒素を取り込む。固定できる個体は大気から、できない個体は土壌から
        if(!p.dormant){
            // 体内が満ちていれば固定も吸収も止まる(窒素を捨てずに、取らない)
            float fix=std::min(N_FIX_RATE*p.n_fix*PLANT_STRIDE,std::max(0.f,N_STORE_MAX-p.n_store));
            p.n_store += fix;
            g_bud.n_fix += fix;
            // 同じマスの株数で分け合う。密集すれば一株あたりの取り分が減る
            int nn2=std::max(1,pgrid.n_at(ci));
            float take=std::min(soil_n[ci],N_UPTAKE*PLANT_STRIDE/(float)nn2);
            take=std::min(take,std::max(0.f,N_STORE_MAX-p.n_store));
            soil_n[ci]-=take; p.n_store+=take;
        }
        if(!p.dormant && daylit[ci]>0.3f && p.height<p.max_height && p.energy>0.5f){
            float g=std::min(GROW_RATE*PLANT_STRIDE,p.max_height-p.height);
            float nn=g*N_PER_GROWTH;
            if(p.n_store<nn) g*=p.n_store/std::max(1e-6f,nn);   // 窒素が律速する
            p.n_store=std::max(0.f,p.n_store-g*N_PER_GROWTH);
            p.height+=g; p.energy-=g*GROW_COST;
            // 建設は呼吸を伴う。出す炭素は体のエネルギーから払う
            float r=std::min(g*GROW_COST*GROW_RESP,std::max(0.f,p.energy));
            p.energy-=r;
            float need=r*CO2_PER_ENERGY;
            o2_add(ci,-std::min(need,o2_at(ci))); co2_add(ci,need);
        }
        float hm=habitat_p(p,ci);
        float hs=salt_mismatch(p,ci);
        float upkeep = RESP_COEF*saturation(p)+HEIGHT_COST*p.height
                 +TOUGH_COST*p.defense+CLONAL_COST*p.clonal+COLD_COST_P*p.cold_tol
                     +HABITAT_COST_P*hm*hm+DROUGHT_COST*p.drought_tol
                     +HABITAT_COST_P*hs*hs;
        // 植物も呼吸する: 維持のために酸素を使い、二酸化炭素を出す。
        // 気体のやりとりは、実際に払った維持費に合わせて下でまとめて行う。ここでは酸欠かどうかだけを見る
        {
            float need = upkeep*PLANT_RESP_CO2*CO2_PER_ENERGY*PLANT_STRIDE
                       * (p.dormant?DORM_UPKEEP:1.f);
            if(need>1e-9f && o2_at(ci)<need*0.5f) upkeep *= 1.3f;   // 酸欠は割高になる
        }
        // 余裕があるとき、組織に水を貯める
        if(!p.dormant && p.water_store<p.store_cap && soil_water[ci]>0.5f){
            float f=std::min(STORE_FILL,p.store_cap-p.water_store);
            f=std::min(f,(soil_water[ci]-0.5f)*STORE_SURPLUS);
            // 貯めた水は土から来る。貯めたぶんだけ土の水が減る(海のマスでは海から入る)
            if(f>0.f){ p.water_store+=f; soil_water[ci]-=f; if(sea[ci]) g_bud.w_sea+=f; }
        }
        // 防御。常に張る分(構成的防御)と、食害の匂いを感じて張る分(誘導防御)の割合は遺伝子 def_const が決める。
        // 常に張れば食われ始めから守れるが費用が続き、匂いで張れば費用は食害のあるときだけだが立ち上がりが遅れる
        float threat=clamp01(odor[2][ci]/1.5f);
        float dtarget=p.tough*(p.def_const+(1.f-p.def_const)*threat);
        p.defense += (dtarget - p.defense)*INDUCE_RATE;
        p.defense -= p.defense*INDUCE_DECAY;
        if(p.defense<0.f) p.defense=0.f;
        if(p.defense>p.tough) p.defense=p.tough;
        // 浮く個体は海流に流される。根を張れないので固定できない
        if(sea[ci] && p.buoyancy>0.4f){
            p.x+=cur_u[ci]*BUOY_DRIFT*p.buoyancy;
            p.y+=cur_v[ci]*BUOY_DRIFT*p.buoyancy;
            // 東西は開いた境界。流されて端を越えた個体は区画の外へ出る。体の炭素と窒素は外へ出たぶんとして帳簿に付ける
            // (貯水は、死んだ株の貯水を戻す処理が区画の外の水として帳簿に付ける)
            if(p.x<0.f || p.x>=GRID_W){
                p.alive=false;
                g_bud.c_immig-=(p.energy+p.height*GROW_COST)*CO2_PER_ENERGY;
                g_bud.n_immig-=p.n_store+p.height*N_PER_GROWTH;
                continue;
            }
            if(p.y<0.f) p.y=0.f;     if(p.y>=GRID_H) p.y=GRID_H-1.f;
        }
        upkeep += BUOY_COST*p.buoyancy
                + STORE_COST*p.store_cap + FIRE_TOL_COST*p.fire_tol
                + CLIMB_COST*p.climb + N_FIX_COST*p.n_fix
                + PARASITE_COST*p.parasite + CARNIVORY_COST*p.carnivory
                + FIRE_SEED_COST*p.fire_seed;
        // 暗所では維持費が落ちる。葉の光化学系を回す必要がないため
        upkeep *= (NIGHT_RESP+(1.f-NIGHT_RESP)*daylit[ci]);
        float sen=senesce(p.age,p.lifespan);
        upkeep *= (1.f+sen);
        upkeep += LIFESPAN_COST*p.lifespan;
        // 自己間引き。混み合うほど、光以外の要因でも消耗する
        {
            int nn=pgrid.n_at(ci);
            if(nn>4) upkeep += (nn-4)*CROWD_COST;
        }
        if(p.dormant) upkeep*=DORM_UPKEEP;   // 葉を落として代謝を絞る
        // 発病。この土地の病原体と型が合ってしまった個体が病む
        if(pload[ci]>0.02f){
            float inf=pload[ci]*immune_match(p.im,ci,pstrain);
            upkeep+=inf*PATH_VIRULENCE;
            infect_p+=inf;
        }
        // 払った維持費の炭素は呼吸で CO2 になる(体に無いぶんは出せない)
        { float c=std::min(upkeep*PLANT_STRIDE,std::max(0.f,p.energy))*PLANT_RESP_CO2*CO2_PER_ENERGY;
          o2_add(ci,-std::min(c,o2_at(ci))); co2_add(ci,c); }
        p.energy-=upkeep*PLANT_STRIDE;
        dbg_cost+=upkeep;      // 更新を飛ばしたぶんをまとめて払う
        if(p.energy<=0.f){
            p.alive=false;
            corpses.push_back({p.x,p.y,p.height*GROW_COST,0,p.n_store+p.height*N_PER_GROWTH});   // 体の炭素と窒素はすべて死骸に
            continue;
        }
        if(!p.dormant && p.energy>=PLANT_DIVIDE_TH
           && (int)plants.size()+(int)babies.size() < PLANT_MAX
           && count_by_rate(1.0)>0){      // 繁殖の判定は段階 3 より前の 1 回の更新に 1 回
            // 有性生殖を選ぶなら、まず相手を見つける必要がある
            const Plant* mate=nullptr;
            if(dist01(rng)<p.sexual){
                int r=(int)SEX_SEARCH;
                for(int dy=-r;dy<=r&&!mate;dy++) for(int dx=-r;dx<=r&&!mate;dx++){
                    int nx=gx+dx, ny=gy+dy;
                    if(nx<0||nx>=GRID_W) continue;
                    if(ny<0||ny>=GRID_H) continue;
                    int ni=idx(nx,ny);
                    for(int t=pgrid.b(ni);t<pgrid.e(ni);t++){
                        const Plant& q=plants[pgrid.items[t]];
                        if(!q.alive||&q==&p||q.lin!=p.lin) continue;
                        mate=&q; break;
                    }
                }
                if(!mate) continue;      // 相手がいなければ、この機会を逃す
            }
            { double v[N_SG_P]; sg_vec_p(p,v);
              for(int k=0;k<N_SG_P;k++){ sg_p_rep[k]+=v[k]; cs_p_rep[k]+=v[k]; }
              sg_p_nr++; cs_p_nr++; }
            int nb=1+(int)(p.brood*(MAX_BROOD-1)+0.5f);
            // 子に渡せるのは投資したぶんまで。多産ほど種子一つあたりの構造に余分がかかり、
            // その余分は作るときの呼吸で CO2 になる
            float share=1.f/(nb*(1.f+(nb-1)*BROOD_TAX));
            float invest=p.energy*(0.20f+0.55f*p.repro_alloc);
            p.energy-=invest;
            { float c=invest*(1.f-nb*share)*CO2_PER_ENERGY;
              o2_add(ci,-std::min(c,o2_at(ci))); co2_add(ci,c); }
            for(int bi=0;bi<nb;bi++){
            Plant c=p;
            c.energy=invest*share;
            // 組み換え: 各形質を、どちらかの親から受け取る
            if(mate){
                #define PICK(F) c.F = (dist01(rng)<0.5f)? p.F : mate->F
                PICK(absorb); PICK(max_height); PICK(shade_tol); PICK(tough);
                PICK(clonal); PICK(disperse);   PICK(cold_tol);  PICK(aquatic); PICK(salt_tol);
                PICK(drought_tol); PICK(dorm_temp); PICK(store_cap);
                PICK(fire_tol); PICK(climb); PICK(sexual); PICK(brood);
                PICK(n_fix); PICK(lifespan); PICK(parasite);
                PICK(carnivory); PICK(fire_seed); PICK(repro_alloc); PICK(seed_bank);
                PICK(buoyancy); PICK(def_const);
                #undef PICK
                // 免疫は座位ごとに独立に受け継ぐ。組み合わせで新しい型が生まれる
                for(int k=0;k<N_IMMUNE;k++)
                    c.im[k]=(dist01(rng)<0.5f)? p.im[k] : mate->im[k];
            }
            bool stay=(dist01(rng)<p.clonal);
            float d=stay?1.5f:(1.5f+SEED_FAR*p.disperse);
            float ang=frand(0.f,6.2832f);
            c.x=p.x+std::cos(ang)*d; c.y=p.y+std::sin(ang)*d;
            if(c.x<0.f) c.x=0.f;     if(c.x>=GRID_W) c.x=GRID_W-1.f;
            if(c.y<0.f) c.y=0.f;     if(c.y>=GRID_H) c.y=GRID_H-1.f;
            int cgx=(int)c.x, cgy=(int)c.y;
            if(cgx<0||cgx>=GRID_W||cgy<0||cgy>=GRID_H) continue;   // この子は失われる
            // 火後発芽: 種子ごとに確率 fire_seed で、焼け跡でしか芽吹かない種子になる
            // (球果の閉じ方の度合いは連続で、部分的な閉鎖が多い。Lamont et al. 1991)。
            // 焼け跡にしか芽吹かない代わり、競争相手がいない。焼けて間もないマスでは乱数を引かない
            if(burn_age[idx(cgx,cgy)]>BURN_WINDOW && p.fire_seed>0.f && dist01(rng)<p.fire_seed){
                corpses.push_back({c.x,c.y,c.energy,0});   // 芽吹かない種子は朽ちる
                continue;
            }
            c.height=0.05f; c.age=0; c.lin=p.lin;
            c.energy-=c.height*GROW_COST;   // 芽生えの体は種子のエネルギーで作る
            if(c.energy<=0.f){               // 小さすぎる種子は芽生えを作れずに朽ちる
                corpses.push_back({c.x,c.y,c.energy+c.height*GROW_COST,0});
                continue;
            }
            // 芽生えの体の窒素は親の貯えから渡す。足りなければ芽生えを作れずに朽ちる
            { float cn=c.height*N_PER_GROWTH;
              if(p.n_store<cn){ corpses.push_back({c.x,c.y,c.energy+c.height*GROW_COST,0}); continue; }
              p.n_store-=cn; }
            c.defense=0.f;
            // 土中で待つか、すぐ芽吹くか
            c.seed_wait = (dist01(rng)<p.seed_bank)
                        ? (int)(frand((float)BANK_MIN,(float)BANK_MAX)) : 0;
            c.id=next_plant_id++;
            c.genet=stay?p.genet:c.id;
            // 変異は、組み換えで両親から選んだ値に掛ける
            c.absorb    =mutate(c.absorb,    0.05f,1.0f);
            c.max_height=mutate(c.max_height,0.05f,5.0f);
            c.shade_tol =mutate(c.shade_tol, 0.0f,1.0f);
            c.tough     =mutate(c.tough,     0.0f,1.0f);
            c.clonal    =mutate(c.clonal,    0.0f,1.0f);
            c.disperse  =mutate(c.disperse,  0.0f,1.0f);
            c.cold_tol  =mutate(c.cold_tol,  0.0f,1.0f);
            c.aquatic     =mutate(c.aquatic,     0.0f,1.0f);
            c.salt_tol    =mutate(c.salt_tol,    0.0f,1.0f);
            c.def_const   =mutate(c.def_const,   0.0f,1.0f);
            c.drought_tol =mutate(c.drought_tol, 0.0f,1.0f);
            c.dorm_temp   =mutate(c.dorm_temp,  -20.0f,20.0f);
            c.store_cap   =mutate(c.store_cap,   0.0f,STORE_MAX);
            c.fire_tol    =mutate(c.fire_tol,    0.0f,1.0f);
            c.climb       =mutate(c.climb,       0.0f,1.0f);
            c.sexual      =mutate(c.sexual,      0.0f,1.0f);
            for(int k=0;k<N_IMMUNE;k++) c.im[k]=mutate(c.im[k],0.0f,1.0f);
            c.brood   =mutate(c.brood,   0.0f,1.0f);
            c.n_fix   =mutate(c.n_fix,   0.0f,1.0f);
            c.lifespan   =mutate(c.lifespan,   0.0f,1.0f);
            c.parasite   =mutate(c.parasite,   0.0f,1.0f);
            c.carnivory  =mutate(c.carnivory,  0.0f,1.0f);
            c.fire_seed  =mutate(c.fire_seed,  0.0f,1.0f);
            c.repro_alloc=mutate(c.repro_alloc,0.0f,1.0f);
            c.seed_bank  =mutate(c.seed_bank,  0.0f,1.0f);
            c.buoyancy   =mutate(c.buoyancy,   0.0f,1.0f);
            c.mut_rate   =mutate(p.mut_rate,   0.0f,1.0f);
            c.mother=p.id; c.father=mate?mate->id:-1; c.gen=p.gen+1;
            c.n_store =0.f;
            c.dormant     =p.dormant;
            c.water_store =0.f;
            babies.push_back(c);
            }
        }
    }
    for(auto& b:babies) plants.push_back(b);
    plant_inflow();
    // 死んだ株の貯水は、そのマスの土へ戻る(海のマスでは海へ)
    for(const auto& p:plants){
        if(p.alive||p.water_store<=0.f) continue;
        int gx=(int)p.x, gy=(int)p.y;
        int ci=(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H)?idx(gx,gy):-1;
        if(ci>=0&&!sea[ci]) soil_water[ci]+=p.water_store;
        else g_bud.w_sea-=p.water_store;
    }
    plants.erase(std::remove_if(plants.begin(),plants.end(),
        [](const Plant& p){ return !p.alive; }),plants.end());
}



