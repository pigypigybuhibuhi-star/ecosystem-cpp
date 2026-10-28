// 微生物(phase3.cpp から分けた。phase3.cpp が順に #include する)
void update_microbes(){
    static int mphase=0;
    mphase=(mphase+1)%MICROBE_STRIDE;
    // update_animals の最後で小さな死骸を消しているので、そこで作った cgrid の添字はずれている。使う直前に作り直す
    rebuild_cgrid();
    std::vector<Microbe> babies;
    // 添字で飛ばす。全件を舐めるのをやめ、該当するものだけに触れる
    int n_mic=(int)microbes.size();
    for(int midx=mphase;midx<n_mic;midx+=MICROBE_STRIDE){
        Microbe& m=microbes[midx];
        if(!m.alive) continue;
        m.age+=MICROBE_STRIDE;
        int gx=(int)m.x, gy=(int)m.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H){ m.alive=false; continue; }
        int ci=idx(gx,gy);
        // 気温。低温では活動が落ちる
        float kill=-8.f-30.f*m.cold_tol;   // 胞子は凍結にも耐える
        // 死んだ細胞の炭素は、他の微生物にすぐ分解されて大気へ戻る
        if(temper[ci]<kill){ a_burn(co2.data(),o2.data(),ci,m.energy); m.alive=false; d_micro++; continue; }
        float dT=(temper[ci]-(24.f-20.f*m.cold_tol))/22.f;
        float act=std::exp(-dT*dT);

        // 死骸を分解する。速いほど多く取れるが、酵素の維持費も上がる
        int cb=cgrid.b(ci), ce=cgrid.e(ci);
        if(cb<ce && frand(0.f,1.f)<MICROBE_TURNOVER){ a_burn(co2.data(),o2.data(),ci,std::max(0.f,m.energy)); m.alive=false; d_micro++; continue; }
        // 分解の速さは基質の量に比例する(一次反応)。マスにある死骸のすべてに同じ割合で取りつく
        // 1 回の更新(STRIDE tick)で分解される割合は、指数の形で換算する
        float dfrac=(cb<ce)? -std::expm1(-MICROBE_DECOMP_K*m.rate*act*(float)(MICROBE_STRIDE*DT)) : 0.f;
        // 取り込みの上限。マスの死骸から取れる量が上限を超えるときは、すべての死骸から取る割合を同じ比で下げる
        if(cb<ce){
            float avail=0.f; for(int t=cb;t<ce;t++) avail+=corpses[cgrid.items[t]].energy;
            float vcap=MICROBE_VMAX*act*std::max(0.f,m.energy)*(float)(MICROBE_STRIDE*DT);
            g_mic_upd++;
            if(avail*dfrac>vcap){ dfrac=vcap/avail; g_mic_cap++; }
        }
        for(int t=cb;t<ce;t++){
            int pk=cgrid.items[t];
            float dec=corpses[pk].energy*dfrac;
            // 分解した割合だけ、死骸の窒素も無機化して土へ出る
            float dn=corpses[pk].n*std::min(1.f,dec/std::max(1e-9f,corpses[pk].energy));
            corpses[pk].n-=dn;
            corpses[pk].energy-=dec;
            m.energy+=dec*MICROBE_EFF;
            // 残りは大気と土へ。体にしなかった分は呼吸で酸化されるので、炭素は全量が CO2 になる
            // 好気的な呼吸なので O2 を 1:1 で使う
            a_burn(co2.data(),o2.data(),ci,dec*(1.f-MICROBE_EFF));
            soil_n[ci]+= dn;
        }
        { float spend=(MICROBE_COST*(0.5f+m.rate)*std::max(0.f,m.energy)+m.cold_tol*MICROBE_COLD_COST*std::max(0.f,m.energy)/MICROBE_REF)*MICROBE_STRIDE;
          a_burn(co2.data(),o2.data(),ci,std::min(spend,std::max(0.f,m.energy)));   // 維持の呼吸(O2 を 1:1)
          m.energy-=spend; }
        if(m.energy<=0.f){ m.alive=false; d_micro++; continue; }
        if(m.energy>=MICROBE_DIV && (int)microbes.size()<800000){
            m.energy*=0.5f;
            Microbe c=m;
            c.energy=m.energy;
            float a=frand(0.f,6.2832f), d=-MICROBE_DISP*std::log(1.f-frand(0.f,1.f));   // 指数分布。1-u は (0,1]
            c.x=m.x+std::cos(a)*d; c.y=m.y+std::sin(a)*d;
            if(c.x<0.f) c.x=0.f;     if(c.x>=GRID_W) c.x=GRID_W-1.f;
            if(c.y<0.f) c.y=0.f;     if(c.y>=GRID_H) c.y=GRID_H-1.f;
            c.rate    =mutate(m.rate,    0.05f,1.0f);
            c.cold_tol=mutate(m.cold_tol,0.0f,1.0f);
            c.id=next_microbe_id++; c.age=0;
            babies.push_back(c);
        }
    }
    for(auto& b:babies) microbes.push_back(b);
    microbes.erase(std::remove_if(microbes.begin(),microbes.end(),
        [](const Microbe& m){ return !m.alive; }),microbes.end());
    // 胞子の定着。STRIDE tick に 1 回、微生物のいない死骸のマスを確かめる。
    // 胞子の型は、空気に混ざった土地の微生物の胞子から一つ(生きている個体を一つ選び、分裂と同じく変異させる)。
    // 群落の体の炭素はそのマスの死骸を食って作ったもので、分解と同じく、食った炭素の EFF が体に、残りが CO2 に、窒素は土に出る
    if(mphase==0 && !microbes.empty() && !corpses.empty()){
        static std::vector<unsigned char> occ;
        occ.assign(GRID_W*GRID_H,0);
        for(const auto& m:microbes){ int gx=(int)m.x, gy=(int)m.y; if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) occ[idx(gx,gy)]=1; }
        int n_src=(int)microbes.size();
        for(int i=0;i<(int)corpses.size();i++){
            int gx=(int)corpses[i].x, gy=(int)corpses[i].y;
            if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
            int ci=idx(gx,gy);
            int cb=cgrid.b(ci), ce=cgrid.e(ci);
            if(cb>=ce || cgrid.items[cb]!=i || occ[ci]) continue;   // マスの最初の死骸でだけ、一度だけ確かめる
            if((int)microbes.size()>=800000) break;
            const Microbe& src=microbes[(size_t)frand(0.f,(float)n_src-0.001f)];
            float dT=(temper[ci]-(24.f-20.f*src.cold_tol))/22.f;
            if(temper[ci]<-8.f-30.f*src.cold_tol) continue;
            float act=std::exp(-dT*dT);
            float pr=-std::expm1(act*std::log1p(-MICROBE_COLONIZE));   // 活動が低いと育つのが遅れ、率が act 倍になる
            if(frand(0.f,1.f)>=pr) continue;
            float avail=0.f; for(int t=cb;t<ce;t++) avail+=corpses[cgrid.items[t]].energy;
            if(avail<=0.f) continue;
            float g=std::min(MICROBE_COLONY_E,avail*MICROBE_EFF);
            float frac=g/MICROBE_EFF/avail;
            for(int t=cb;t<ce;t++){
                int pk=cgrid.items[t];
                float dec=corpses[pk].energy*frac;
                float dn=corpses[pk].n*std::min(1.f,frac);
                corpses[pk].n-=dn; corpses[pk].energy-=dec;
                a_burn(co2.data(),o2.data(),ci,dec*(1.f-MICROBE_EFF));
                soil_n[ci]+=dn;
            }
            Microbe c=src;
            c.x=gx+frand(0.f,0.999f); c.y=gy+frand(0.f,0.999f);
            c.energy=g; c.alive=true;
            c.rate    =mutate(src.rate,    0.05f,1.0f);
            c.cold_tol=mutate(src.cold_tol,0.0f,1.0f);
            c.id=next_microbe_id++; c.age=0;
            microbes.push_back(c);
            occ[ci]=1; g_mic_col++;
        }
    }
    // 分解者が絶えた場合、死骸の山から再び湧く
    if(microbes.size()<50 && corpses.size()>1000){
        for(int t=0;t<200;t++){
            const Corpse& cp=corpses[(size_t)frand(0.f,(float)corpses.size()-0.001f)];
            Microbe m; m.x=cp.x; m.y=cp.y;
            m.energy=0.4f; m.alive=true; m.age=0;
            m.rate=frand(0.1f,0.8f); m.cold_tol=frand(0.f,0.6f);
            m.id=next_microbe_id++;
            microbes.push_back(m);
            g_bud.c_immig+=m.energy*CO2_PER_ENERGY;   // 外から飛んできた胞子が持ち込む炭素
        }
    }
}

// 山火事。落雷で着火し、乾いた燃料の上を風下へ広がる。
// 埋もれた炭素を大気に戻し、焼け跡という新しい場所を開く。

