// 崖崩れ(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 崖崩れ。急斜面の表土が滑り、その場の生物を巻き込んで斜面の下へ運ぶ
// 崩土は最も急な下りを伝って流れ、出発点からの見通し角が H/L を下回るところで止まる(少なくとも一つ下のマスまでは動く)
void slide_at(int ci){
    slide_events++;
    ev_record(EV_SLIDE,ci%GRID_W,ci/GRID_W,slope[ci]);
    // 到達点を決める
    float z0=elev[ci];
    int cur=ci; float run=0.f;
    for(int st=0;st<SLIDE_RUN_MAX;st++){
        int x=cur%GRID_W, y=cur/GRID_W, lo=-1; float le=elev[cur], ld=0.f;
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            if(!dx&&!dy) continue;
            int nx=x+dx, ny=y+dy;
            if(nx<0||nx>=GRID_W) continue;
            if(ny<0||ny>=GRID_H) continue;
            int j=idx(nx,ny); if(elev[j]<le){ le=elev[j]; lo=j; ld=(dx&&dy)?1.41421356f:1.f; } }
        if(lo<0) break;                               // 窪地。ここに溜まる
        float r2=run+ld*CELL_M;
        if(st>0 && (z0-le)*ELEV_M<SLIDE_HL*r2) break; // 見通し角が H/L を下回る
        cur=lo; run=r2;
        if(sea[cur]) break;                           // 海に入ったらそこで沈む
    }
    float dep_x=(cur%GRID_W)+0.5f, dep_y=(cur/GRID_W)+0.5f;
    for(int t=pgrid.b(ci);t<pgrid.e(ci);t++){
        Plant& p=plants[pgrid.items[t]];
        if(!p.alive) continue;
        p.alive=false;
        // 体の炭素(建設に払った GROW_COST ぶん)と残りのエネルギーは死骸として崩土と一緒に運ばれる
        corpses.push_back({dep_x,dep_y,std::max(0.f,p.energy)+p.height*GROW_COST,0,p.n_store+p.height*N_PER_GROWTH});
    }
    // 巻き込まれた動物は埋まって死ぬ(潜っていても、木の上にいても、立っている表土ごと流される)
    for(auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H||idx(gx,gy)!=ci) continue;
        a.alive=false; d_slide++;
        life_note(a.id,g_frame,"崖崩れ",(float)a.age);
        corpses.push_back({dep_x,dep_y,a.body*CORPSE_RETURN+std::max(0.f,a.energy)+a.fat,1});
    }
    if(cur==ci) return;                               // 動けない(下りが無い)なら地形も表土も変わらない
    // 表土の窒素と土砂を到達点へ運ぶ。崩れた跡は SLIDE_DEPTH_M だけ低くなり、傾斜が緩む
    float moved=soil_n[ci]*0.65f;
    soil_n[ci]-=moved; soil_n[cur]+=moved;
    float dz=SLIDE_DEPTH_M/ELEV_M;
    elev[ci]-=dz;
    if(sea[cur]) elev[cur]=std::min(elev[cur]+dz,SEA_LEVEL-1e-6f);  // 海底に沈む。海面より上に出る分は沖へ散る(土砂は帳簿の量ではない)
    else         elev[cur]+=dz;
    if(!sea[ci] && elev[ci]<SEA_LEVEL) elev[ci]=SEA_LEVEL;           // 陸のまま(海の判定は起動時の地形で決まる)
    slope_around(ci); slope_around(cur);
}
// 無限斜面の安全率 FS = c'/(γ z sinα cosα) + (1 − m γw/γ) tanφ'/tanα。m は崩れる厚さのうち地下水で飽和した割合
// (帯水層は地表から GW_DEPTH_M の厚さで、飽和した厚さ gw_h が GW_DEPTH_M − SLIDE_DEPTH_M を超えた分が崩れる層に入る)
inline double slope_fs(int ci){
    double s=slope[ci];
    double h=std::sqrt(1.0+s*s), sa=s/h, ca=1.0/h;
    double m=std::min(1.0,std::max(0.0,(double)(gw_h[ci]-(GW_DEPTH_M-SLIDE_DEPTH_M))/SLIDE_DEPTH_M));
    return SLOPE_C_KPA/(SLOPE_GAMMA*SLIDE_DEPTH_M)/(sa*ca)+(1.0-m*WATER_GAMMA/SLOPE_GAMMA)*REPOSE_SLOPE/s;
}
void update_slides(){
    // 崖崩れ。地下水が崩れる層に上がり、安全率が 1 を下回った斜面が滑り、その場の生物を押し流す
    if(g_replay){                      // 再生中は記録の各件だけを起こす
        ev_replay_due(EV_SLIDE,[](const EnvEvent& e){ slide_at(idx(e.x,e.y)); });
        return;
    }
    const float wet=GW_DEPTH_M-SLIDE_DEPTH_M;
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        if(sea[ci]||gw_h[ci]<=wet||slope[ci]<=0.f) continue;   // 乾いた斜面は FS > 1(粘着力と摩擦だけで立つ)
        if(slope_fs(ci)>=1.0) continue;
        slide_at(ci);
    }
}

// 索引は、生物が生まれたり死んだりするたびに古くなる。
// 使う直前に作り直す。添字のずれが範囲外アクセスになるのを防ぐ
void rebuild_cgrid(){
    cgrid.reset();
    for(const auto& c:corpses){
        int gx=(int)c.x, gy=(int)c.y;
        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) cgrid.count(idx(gx,gy));
    }
    cgrid.finalize();
    for(int i=0;i<(int)corpses.size();i++){
        int gx=(int)corpses[i].x, gy=(int)corpses[i].y;
        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) cgrid.put(idx(gx,gy),i);
    }
}
