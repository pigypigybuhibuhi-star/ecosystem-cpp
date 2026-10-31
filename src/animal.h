// 動物(phase3.cpp から分けた。phase3.cpp が順に #include する)
void update_animals(){
    agrid.reset();
    for(const auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        agrid.count(idx(gx,gy));
    }
    agrid.finalize();
    for(int i=0;i<(int)animals.size();i++){
        if(!animals[i].alive) continue;
        int gx=(int)animals[i].x, gy=(int)animals[i].y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        agrid.put(idx(gx,gy),i);
    }
    cgrid.reset();
    for(const auto& cp:corpses){
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        cgrid.count(idx(gx,gy));
    }
    cgrid.finalize();
    for(int i=0;i<(int)corpses.size();i++){
        int gx=(int)corpses[i].x, gy=(int)corpses[i].y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        cgrid.put(idx(gx,gy),i);
    }
    std::fill(fuel.begin(),fuel.end(),0.f);
    for(int ci=0;ci<GRID_W*GRID_H;ci++)
        corpse_density[ci]=(float)cgrid.n_at(ci);
    for(const auto& cp:corpses){
        if(cp.origin!=0) continue;            // 燃えるのは枯れた植物
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        fuel[idx(gx,gy)]+=cp.energy;
    }
    // 食べたのに同化しなかった分(糞)と、齧るときに落とした葉は、そのマスの死骸(分解物)になる。
    // 食べ物の炭素 = 同化 + 糞・落とした葉(質量保存)。マスごとに貯めて、更新の終わりに死骸へ移す
    static std::vector<float> dung(GRID_W*GRID_H,0.f);
    static std::vector<int>   dung_cells;
    // 動物は窒素の在庫を持たない。食べた窒素は同じマスの糞(分解物)にすべて出す
    static std::vector<float> dung_n(GRID_W*GRID_H,0.f);
    auto egest=[&](int c,float amt,float nn=0.f){
        if(amt<=0.f&&nn<=0.f) return;
        if(dung[c]==0.f&&dung_n[c]==0.f) dung_cells.push_back(c);
        if(amt>0.f) dung[c]+=amt;
        if(nn>0.f)  dung_n[c]+=nn;
    };
    std::vector<Animal> babies;
    int n0=(int)animals.size();
    static int aphase=0;
    aphase=(aphase+1)%ANIMAL_STRIDE;
    // 2tickに分けて処理する。個体数が多いので走査自体を半分にする
    for(int i=aphase;i<n0;i+=ANIMAL_STRIDE){
        Animal& a=animals[i];
        if(!a.alive) continue;
        a.age+=ANIMAL_STRIDE;
        // それぞれの感覚器が、自分の向いた方角に、自分の射程だけ探る
        float in[N_IN];
        float sensor_span=0.f;
        for(int k=0;k<N_SENSOR;k++){
            in[k]=0.f;
            const Animal::Sensor& S=a.sen[k];
            if(S.range<SENSOR_MIN) continue;      // 器官として成立していない
            sensor_span+=S.range;
            float ang=a.dir+S.angle;
            float cs=std::cos(ang), sn=std::sin(ang);
            for(int s=1;s<=RAY_STEPS;s++){
                float d=S.range*s/(float)RAY_STEPS;
                int rx=(int)(a.x+cs*d), ry=(int)(a.y+sn*d);
                if(rx<0||rx>=GRID_W) continue;
                if(ry<0||ry>=GRID_H) continue;
                int ri=idx(rx,ry);
                int band=(int)(S.tune*N_TUNE); if(band>N_TUNE-1) band=N_TUNE-1;
                float sig=0.f;
                switch(band){
                case 0: if(!pgrid.empty_at(ri)) sig=1.f; break;      // 植物
                case 1:                                              // 自分より小さい動物
                    { int e=std::min(agrid.e(ri),agrid.b(ri)+4);     // 先頭数体だけ見る
                      for(int t=agrid.b(ri);t<e;t++){
                        int j=agrid.items[t];
                        if(j!=i && animals[j].body < a.body*0.9f){ sig=1.f; break; }
                      } } break;
                case 2:                                              // 自分より大きい動物
                    { int e=std::min(agrid.e(ri),agrid.b(ri)+4);
                      for(int t=agrid.b(ri);t<e;t++){
                        int j=agrid.items[t];
                        if(j!=i && animals[j].body > a.body*1.1f){ sig=1.f; break; }
                      } } break;
                case 3: if(!cgrid.empty_at(ri)) sig=1.f; break;      // 死骸
                case 4: sig=std::min(1.f,fire[ri]); break;           // 炎
                case 5: sig=clamp01(canopy_h[ri]/3.f); break;        // 林冠の高さ
                case 6: sig=clamp01(aload[ri]/PATH_LOAD_MAX); break; // 病原体
                case 11:                                             // 似た体。群れの前提
                    // 相手の体の表面の化学的な手がかり(免疫の型 im)が自分とどれだけ似ているかを、
                    // 自分の型と比べて感じ取る(表現型の照合。Lacy & Sherman 1983)。系統 lin は使わない
                    { float c=0.f; int e=std::min(agrid.e(ri),agrid.b(ri)+6);
                      for(int t=agrid.b(ri);t<e;t++){
                          int j=agrid.items[t];
                          if(j==i) continue;
                          float dm=0.f;
                          for(int m=0;m<N_IMMUNE;m++) dm+=std::fabs(animals[j].im[m]-a.im[m]);
                          c+=1.f-dm/N_IMMUNE;
                      }
                      sig=std::min(1.f,c*0.45f); }
                    break;
                case 12: sig=sea[ri]?1.f:0.f; break;                 // 水際
                // 嗅覚。濃さそのものを読む。遮蔽を越え、時間に残った跡も追える
                default: sig=clamp01(odor[band-7][ri]/ODOR_SCALE); break;
                }
                if(sig>0.02f){ in[k]=sig*(1.f-d/S.range); break; }
            }
        }
        int sci=idx(std::max(0,std::min(GRID_W-1,(int)a.x)),
                    std::max(0,std::min(GRID_H-1,(int)a.y)));
        in[N_SENSOR+0]=std::min(1.f,a.energy/A_DIVIDE_TH);
        in[N_SENSOR+1]=clamp01((temper[sci]+30.f)/70.f);
        in[N_SENSOR+2]=soil_water[sci];
        in[N_SENSOR+3]=std::min(1.f,a.body/3.f);   // 自分の体格。戦うか逃げるかを分ける
        in[N_SENSOR+4]=1.f;                        // バイアス
        // 再帰入力。tanh を通しているので本来は有界だが、念のため縛る
        for(int k=0;k<N_REC;k++){
            float r=a.lesion?0.f:a.rec[k];   // 陰性対照は再帰の入力を切る(試験器の lesion と同じ)
            if(!(r>-1.5f&&r<1.5f)) r=0.f;   // NaN もここで潰れる
            in[N_SENSOR+5+k]=r;
        }
        // 神経網。計算は brain_forward(試験器と同じもの)
        float out[N_OUT];
        int   n_node=0, n_syn=0;
        brain_forward(a,in,out,n_node,n_syn);
        float o0=out[0], o1=out[1];
        float want_eat   =out[2];   // 採食するか
        float want_attack=out[3];   // 襲うか
        for(int k=0;k<N_REC;k++) a.rec[k]=out[4+k];   // 次のtickへ持ち越す
        if(g_track_f||g_trig_f) field_record(a,g_frame,in,out);
        // 神経組織は高くつく。ニューロンと結線の数に比例する
        float brain_use=n_node*(NODE_COST/SYN_COST)+n_syn;
        float turn_mult=(o0+1.f)*0.5f, thrust_mult=(o1+1.f)*0.5f;

        int hci=idx(std::max(0,std::min(GRID_W-1,(int)a.x)),
                    std::max(0,std::min(GRID_H-1,(int)a.y)));
        if(a.body<a.size && a.energy>1.0f){
            float g=std::min(BODY_GROW_RATE,a.size-a.body);
            a.body+=g; a.energy-=g*BODY_GROW_COST;
            // 組織に残るのは CORPSE_RETURN ぶん。残りは合成の呼吸で CO2 になる
            a_burn(co2.data(),o2.data(),hci,g*(BODY_GROW_COST-CORPSE_RETURN));
        }

        // 枝が支えられる高さ。重い体は高くまで登れない
        float bearable = canopy_h[hci]*BRANCH_SUPPORT/std::max(0.15f,a.body);
        float reachable = std::min(canopy_h[hci],bearable);

        // 登り降り。登攀が巧いほど安く登れる
        float want = reachable*a.arboreal;
        if(want>a.perch+0.02f){
            float up=std::min(want-a.perch,CLIMB_RATE0+CLIMB_RATE1*a.arboreal);
            a.perch+=up;
            float e=up*CLIMB_UP_COST/(0.25f+a.arboreal);
            a_burn(co2.data(),o2.data(),hci,std::min(e,a.energy));
            a.energy-=e;
        } else if(a.perch>reachable){
            // 支えが失われた。滑り落ちる
            float drop=a.perch-reachable;
            a.perch=reachable;
            if(drop>0.5f){
                // 傷を治すのに使うエネルギー。これも呼吸で CO2 になる
                float e=drop*FALL_DAMAGE*(1.f-a.arboreal*0.7f);
                a_burn(co2.data(),o2.data(),hci,std::min(e,a.energy));
                a.energy-=e;
            }
        } else if(a.perch>want){
            a.perch=std::max(want,a.perch-DESCEND_RATE);    // 降りるのは安い
        }
        if(a.perch<0.f) a.perch=0.f;

        float terr;
        if(a.perch>0.25f){
            // 樹上。地面の起伏に縛られない。登攀が巧いほど速く渡れる
            terr = 0.35f+0.75f*a.arboreal;
        } else if(sea[hci]) terr = (0.25f + WATER_SPEED*a.aquatic);   // 泳ぐ
        else         terr = (1.f - 0.75f*a.aquatic)                    // 陸を這う
                          / (1.f + SLOPE_COEF*slope[hci]);             // 急坂は遅い
        if(terr<0.05f) terr=0.05f;

        // 酸素の分圧。柱の気体の全量との比で見る
        float o2f = o2_at(hci)/std::max(1e-6f,o2_at(hci)+co2_at(hci)+N2_COL);
        if(sea[hci]) o2f *= DISSOLVED_O2;              // 水に溶ける酸素は僅か
        o2f *= std::exp(-altitude(hci)*ALT_O2_DECAY);  // 高所は空気が薄い
        float o2rel = o2f/0.21f;                       // 地球の21%を1.0とする
        // 呼吸効率が高いほど、薄い酸素でも動ける
        float thr_half = 0.55f*(1.f-0.80f*a.resp);
        float o2fit = std::min(1.f, o2rel*(1.f+thr_half)/(o2rel+thr_half));
        a.aerobic = o2fit;

        a.dir+=frand(-a.turn*turn_mult,a.turn*turn_mult)*TURN_SCALE;   // 1 回の更新ぶんの向きのゆらぎ
        float sp=a.speed*A_MAX_SPEED*thrust_mult*terr*o2fit;
        float nx=a.x+std::cos(a.dir)*sp, ny=a.y+std::sin(a.dir)*sp;
        if(nx<0.f){ nx=0.f; a.dir=3.14159265f-a.dir; }       // 東西の端も南北と同じく区画の縁
        if(nx>=GRID_W){ nx=GRID_W-1.f; a.dir=3.14159265f-a.dir; }
        if(ny<0.f){ ny=0.f; a.dir=-a.dir; }
        if(ny>=GRID_H){ ny=GRID_H-1.f; a.dir=-a.dir; }
        int ngx=(int)nx, ngy=(int)ny;
        bool blocked=false;
        if(ngx>=0&&ngx<GRID_W&&ngy>=0&&ngy<GRID_H){
            int ni=idx(ngx,ngy);
            // 水と陸の境で止める閾値は置かない。泳ぐ・這う速さ(terr)と浸かることへの不適合の費用(habitat_a)が
            // aquatic に連続に従うので、入るかどうかは神経網の向きと、その場の速さと費用の結果で決まる
            // 樹上を渡る: 隣の林冠が近い高さでないと渡れず、降りるしかない
            if(a.perch>0.25f && canopy_h[ni] < a.perch-PERCH_GAP){
                a.perch=std::max(0.f,canopy_h[ni]);         // 木が途切れた。降りる
                float e=0.04f*(1.f-a.arboreal*0.6f);
                a_burn(co2.data(),o2.data(),ni,std::min(e,a.energy));
                a.energy-=e;
            }
        }
        if(blocked) a.dir+=3.14159f;
        else { a.x=nx; a.y=ny; }
        int gx=(int)a.x, gy=(int)a.y;
        int ci=(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H)?idx(gx,gy):0;
        // 潜り具合。地中は気温が和らぐ
        a.depth += (a.burrow - a.depth)*BURROW_RATE;
        float felt = temper[ci];
        if(a.depth>0.05f){
            float year_mean=TEMP_EQUATOR
                +(TEMP_POLE-TEMP_EQUATOR)*LAT_N - LAPSE_RATE*altitude(ci);
            felt = temper[ci]+(year_mean-temper[ci])*BURROW_BUFFER*a.depth;
        }
        if(felt<kill_temp_a(a)){
            a.alive=false; d_cold++;
            corpses.push_back({a.x,a.y,a.body*CORPSE_RETURN+std::max(0.f,a.energy)+a.fat,1});
            continue;
        }
        float dT=(felt-opt_temp_a(a))/25.f;
        float act=std::exp(-dT*dT)*o2fit;   // 酸素が薄ければ動きも鈍る
        // 昼夜。夜行性なら夜に本領を発揮する
        float night=1.f-daylit[ci];      // その土地が今、昼か夜か
        act *= (0.35f + 0.65f*(a.nocturnal*night + (1.f-a.nocturnal)*daylit[ci]));
        // 潜っていれば地表のものを食えない
        if(a.depth>0.5f) act*=0.15f;

        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H){
            // 三つの消化能力。合計に上限があるため、全部は伸ばせない
            float cap_sum = a.diet + a.detritus*0.5f;   // 分解は草食と干渉しにくい
            if(cap_sum > 1.f) cap_sum = 1.f;
            float plant_dig = 1.f - cap_sum*cap_sum;
            float detr_dig  = a.detritus;
            float meat_dig  = a.diet;
            if(plant_dig>0.01f && want_eat>ACT_THRESH){
                // 摂食量は体格に比例する。大きい体は多く処理できる
                int max_bites0=1+(int)(a.body*2.0f);   // 段階 3 より前の 1 回の更新あたり
                if(max_bites0>6) max_bites0=6;
                int max_bites=count_by_rate(max_bites0);
                int bites=0;
                for(int t=pgrid.b(ci);t<pgrid.e(ci)&&bites<max_bites;t++){
                    int pi=pgrid.items[t];
                    Plant& p=plants[pi];
                    if(!p.alive||p.eff_height>a.perch+a.reach) continue;
                    float chew=a.jaw/(a.jaw+p.defense*2.0f);
                    // 一口が大きいほど速いが、噛み砕きが雑になって効率が落ちる
                    float bs=0.4f+1.2f*a.bite_size;
                    float bite=p.energy*A_BITE*bs*chew*act;
                    p.energy-=bite;
                    // 齧られれば葉を失う。光合成能力が落ち、回復に時間がかかる。
                    // これがあって初めて、その場が枯れ、探す必要が生まれる
                    float h0=p.height;
                    p.height=std::max(0.03f,p.height-bite*GRAZE_DEFOLIATE*bs);
                    float eff=(1.f-A_EAT_LOSS)*(0.6f+0.8f*a.digest)/std::max(0.5f,bs);
                    float gain=std::min(1.f,eff*plant_dig);   // 同化は食べた量を超えない
                    a_burn(co2.data(),o2.data(),ci,a_repay(a.energy,bite*gain));   // 登り降りの借りを返す(O2 を 1:1)
                    a.energy+=bite*gain;
                    // 同化しなかった分は糞に、落とした葉の組織(建設に GROW_COST を払った炭素)は落葉に
                    egest(ci,bite*(1.f-gain)+(h0-p.height)*GROW_COST,(h0-p.height)*N_PER_GROWTH);
                    odor[2][ci]+=ODOR_WOUND*bite;   // 傷ついた植物が匂いを放つ
                    a.eaten++; bites++;
                    if(p.energy<=0.f){ p.alive=false; egest(ci,p.height*GROW_COST,p.n_store+p.height*N_PER_GROWTH); }   // 残った茎も分解物に
                }
            }
            // 死骸を漁る: 植物由来は分解能力で、動物由来は肉食能力で処理する
            // 死骸を漁る。分解者は栄養価の低い落葉を大量に処理する
            {
                int cb=cgrid.b(ci), ce=cgrid.e(ci);
                if(cb<ce && a.perch<0.5f && want_eat>ACT_THRESH){
                    int nd=ce-cb;
                    int tries=count_by_rate(1+(int)(3.f*detr_dig));   // 能力が高いほど多く処理する
                    for(int t=0;t<tries;t++){
                        int pk=cgrid.items[cb+(int)(frand(0.f,(float)nd-0.001f))];
                        bool is_plant=(corpses[pk].origin==0);
                        float cap =is_plant?detr_dig:meat_dig;
                        float loss=is_plant?DETRITUS_LOSS:MEAT_LOSS;
                        if(cap<=0.005f) continue;
                        float bite=corpses[pk].energy*SCAV_EFF*cap;
                        float bn=corpses[pk].n*std::min(1.f,SCAV_EFF*cap);
                        corpses[pk].energy-=bite; corpses[pk].n-=bn;
                        a_burn(co2.data(),o2.data(),ci,a_repay(a.energy,bite*(1.f-loss)));
                        a.energy+=bite*(1.f-loss); a.eaten++;
                        egest(ci,bite*loss,bn);   // 同化しなかった分は糞。窒素は糞に
                    }
                }
            }
            // 生きた動物を襲う
            if(meat_dig>0.02f && want_attack>ACT_THRESH){
                int bx=(int)a.x, by=(int)a.y, pk=-1;
                for(int dy=-1;dy<=1&&pk<0;dy++) for(int dx=-1;dx<=1&&pk<0;dx++){
                    int nx=bx+dx, ny=by+dy;
                    if(nx<0||nx>=GRID_W) continue;
                    if(ny<0||ny>=GRID_H) continue;
                    int ni=idx(nx,ny);
                    for(int t=agrid.b(ni);t<agrid.e(ni);t++){
                        int j=agrid.items[t];
                        if(j==i||!animals[j].alive) continue;
                        if(a.body<=animals[j].body*PRED_RATIO) continue;
                        if(std::fabs(a.perch-animals[j].perch)>=0.6f) continue;
                        if(animals[j].depth>0.5f) continue;      // 穴に逃げ込まれた
                        pk=j; break;
                    }
                }
                if(pk>=0){
                    Animal& v=animals[pk];
                    // わずかな肉食性でも噛めるようにする。線形だと谷を渡れない
                    float bite=v.energy*PRED_EFF*std::sqrt(meat_dig)*act;
                    v.energy-=bite;
                    // 毒。解毒できなければ、食っても損をする
                    float tox=v.toxin*(1.f-a.tox_resist);
                    float gain=bite*(1.f-MEAT_LOSS)*(1.f-tox);
                    a_burn(co2.data(),o2.data(),ci,a_repay(a.energy,gain));
                    a.energy+=gain;
                    egest(ci,bite-gain);   // 同化しなかった分(毒で吸収できない分を含む)は糞
                    life_note(a.id,g_frame,"捕食",bite);
                    { float e=v.toxin*TOXIN_HARM*(1.f-a.tox_resist);   // 解毒に使うエネルギー
                      a_burn(co2.data(),o2.data(),ci,std::min(e,a.energy));
                      a.energy-=e; }
                    a.eaten++;
                    if(v.energy<=0.f){
                        v.alive=false; d_preyed++;
                        corpses.push_back({v.x,v.y,v.body*CORPSE_RETURN+std::max(0.f,v.energy)+v.fat,1});
                    }
                }
            }
        }
        float hma=habitat_a(a,ci);
        float hsa=salt_mismatch_a(a,ci);
        // 食虫植物に捕らえられる。土壌の窒素が乏しい場所ほど罠が育つ
        if(a.depth<0.2f){
            for(int t=pgrid.b(ci);t<pgrid.e(ci);t++){
                Plant& tp=plants[pgrid.items[t]];
                if(!tp.alive||tp.carnivory<0.12f) continue;
                float catch_p=tp.carnivory*CARNIVORY_CATCH/(1.f+a.body*1.5f);
                if(dist01(rng)>=catch_p) continue;
                // 体のエネルギーの 55% と脂肪の 40% を奪われる(脂肪は脂肪からだけ引く)
                float te=std::max(0.f,a.energy)*0.55f, tf=a.fat*0.4f;
                float taken=te+tf;
                a.energy-=te; a.fat-=tf;
                tp.energy  += taken*0.70f;                 // 直接エネルギーになる
                a_burn(co2.data(),o2.data(),ci,taken*0.30f);   // 残りは消化の呼吸で CO2
                { float n0=tp.n_store;
                  tp.n_store  = std::min(N_STORE_MAX,tp.n_store+taken*CARNIVORY_N);
                  g_bud.n_carn+=tp.n_store-n0; }
                a.eaten--;
                break;
            }
            if(a.energy<=0.f){
                a.alive=false; d_preyed++;
                corpses.push_back({a.x,a.y,a.body*CORPSE_RETURN+a.fat,1});
                continue;
            }
        }
        float upkeep=a.reach*A_REACH_COST
                   +a.body*A_SIZE_COST+a.jaw*JAW_COST+COLD_COST_A*a.cold_tol
                   +HABITAT_COST_A*hma*hma
                   +HABITAT_COST_A*hsa*hsa
                   +a.detritus*DETRITUS_COST+sensor_span*SENSOR_COST
                   +brain_use*SYN_COST
                   +a.sexual*SEX_COST+a.arboreal*ARBOREAL_COST
                   +(1.f-a.mut_rate)*MUT_RATE_COST   // 精密な修復は高くつく
                   +a.digest*DIGEST_COST                   // 長い腸は重い
                   +a.toxin*TOXIN_COST+a.tox_resist*TOXRES_COST
                   +a.burrow*BURROW_COST+a.nocturnal*NOCTURNAL_COST
                   +a.fat*FAT_COST;
        // 基礎代謝と形質の維持費は組織の炭素に比例する。歩く費用は別に足す
        float maint=A_SMR*a.body*CORPSE_RETURN*(1.f+upkeep/A_TRAIT_REF);   // 体の維持費(1 回の更新)
        float metab=maint+sp*A_MOVE_COST;
        float sen=senesce(a.age,a.lifespan,LIFE_MIN_A,LIFE_MAX_A);
        metab *= (1.f+sen);
        metab += LIFESPAN_COST*a.lifespan;
        if(estrus_ready(a)) metab += ESTRUS_COST*a.estrus*a.body;   // 性フェロモンの合成
        // 発病
        if(aload[ci]>0.02f){
            float inf=aload[ci]*immune_match(a.im,ci,astrain);
            metab+=inf*PATH_VIRULENCE;
            infect_a+=inf;
        }
        // 呼吸: 酸素の分圧で効率が決まる。水中と高所は薄い
        float o2_need = metab*O2_PER_METAB*(1.f-RESP_EFFICIENCY*a.resp);
        float avail   = std::min(o2_need,o2_at(ci));
        float eff=(o2_need>1e-9f)?avail/o2_need:1.f;
        if(eff<0.15f) eff=0.15f;
        // 使ったエネルギーは糖の酸化(または発酵)でまかなう。体にある分だけが燃料になる。
        // 好気呼吸は CH2O + O2 → CO2 + H2O(呼吸商 1)なので、CO2 に出る炭素と同じ量の O2 を使う。
        // マスの O2 で足りない分は発酵で、炭素は乳酸などの有機物のまま体の外(糞)に出る
        { float spent=metab/eff + a.resp*RESP_COST;
          float fc=std::min(spent,std::max(0.f,a.energy))*CO2_PER_ENERGY;
          float ox=std::min(fc,std::max(0.f,o2_at(ci)));
          o2_add(ci,-ox); co2_add(ci,ox);
          egest(ci,(fc-ox)/CO2_PER_ENERGY);
          a.energy -= spent; }
        // 余剰を脂肪に回す。飢えれば取り崩す。どちらも繁殖の閾値ではなく体の維持費の日数で測る(D-2026-10-24-3)。
        // 糖の貯え GLYCO_DAYS 日分を超えた余りを、脂肪が store_fat × FAT_DAYS_MAX 日分になるまで回す
        float day_cost=maint*UPD_PER_DAY;
        float glyco=(float)GLYCO_DAYS*day_cost;
        float fat_goal=a.store_fat*(float)FAT_DAYS_MAX*day_cost;
        if(a.energy>glyco && a.fat<fat_goal){
            float put=std::min((a.energy-glyco)*FAT_RATE,
                               fat_goal-a.fat);
            a.energy-=put; a.fat+=put*0.85f;
            a_burn(co2.data(),o2.data(),ci,put*0.15f);   // 脂肪の合成の呼吸
        } else if(a.energy<glyco && a.fat>0.f){
            float take=std::min(a.fat,FAT_DRAW);
            // 代謝の借りを脂肪で返す。借りの分は O2 をまだ使っていないので、O2 を 1:1 で使って CO2 に出す
            a_burn(co2.data(),o2.data(),ci,a_repay(a.energy,take));
            a.fat-=take; a.energy+=take;
        }
        if(a.energy<=0.f){
            a.alive=false; d_starve++;
            life_note(a.id,g_frame,"餓死",(float)a.age);
            corpses.push_back({a.x,a.y,a.body*CORPSE_RETURN+std::max(0.f,a.energy)+a.fat,1});
            continue;
        }
        if((int)animals.size()+(int)babies.size() < ANIMAL_MAX
           && a.energy>=A_DIVIDE_TH*(1.4f-0.6f*a.repro_alloc)
           && a.body>=a.size*MATURE_FRAC
           && count_by_rate(1.0)>0){      // 繁殖の判定は段階 3 より前の 1 回の更新に 1 回
            const Animal* mate=nullptr;
            if(dist01(rng)<a.sexual){
                int r=(int)SEX_SEARCH;
                int bx=(int)a.x, by=(int)a.y;
                for(int dy=-r;dy<=r&&!mate;dy++) for(int dx=-r;dx<=r&&!mate;dx++){
                    int nx=bx+dx, ny=by+dy;
                    if(nx<0||nx>=GRID_W) continue;
                    if(ny<0||ny>=GRID_H) continue;
                    int ni=idx(nx,ny);
                    for(int t=agrid.b(ni);t<agrid.e(ni);t++){
                        const Animal& q=animals[agrid.items[t]];
                        if(!q.alive||&q==&a||q.lin!=a.lin) continue;
                        if(q.body<q.size*MATURE_FRAC) continue;
                        mate=&q; break;
                    }
                }
                if(!mate) continue;
            }
            // 何匹産むか。多く産めば一匹あたりの持ち分が減る
            { double v[N_SG_A]; sg_vec_a(a,v);
              for(int k=0;k<N_SG_A;k++){ sg_a_rep[k]+=v[k]; cs_a_rep[k]+=v[k]; }
              sg_a_nr++; cs_a_nr++; }
            int nb=1+(int)(a.brood*(MAX_BROOD-1)+0.5f);
            // 子に渡せるのは投資したぶんまで。多産の余分は作るときの呼吸で CO2 になる
            float share=1.f/(nb*(1.f+(nb-1)*BROOD_TAX));
            float invest=a.energy*(0.20f+0.55f*a.repro_alloc);
            a.energy-=invest;
            { int bx=(int)a.x, by=(int)a.y;
              int bi0=(bx>=0&&bx<GRID_W&&by>=0&&by<GRID_H)?idx(bx,by):0;
              float c=invest*(1.f-nb*share)*CO2_PER_ENERGY;
              o2_add(bi0,-std::min(c,o2_at(bi0))); co2_add(bi0,c); }
            for(int bi=0;bi<nb;bi++){
            // 子の体 0.10 を作る費用は子の取り分から払う。足りない卵は育たずに死骸になる
            float cE=invest*share, body0=0.10f*BODY_GROW_COST;
            if(cE<=body0){ corpses.push_back({a.x,a.y,cE,1}); continue; }
            { int bx=(int)a.x, by=(int)a.y;
              int bi0=(bx>=0&&bx<GRID_W&&by>=0&&by<GRID_H)?idx(bx,by):0;
              a_burn(co2.data(),o2.data(),bi0,0.10f*(BODY_GROW_COST-CORPSE_RETURN)); }
            Animal c=a;
            c.energy=cE-body0;
            if(mate){
                #define PICK(F) c.F = (dist01(rng)<0.5f)? a.F : mate->F
                PICK(speed); PICK(turn); PICK(reach); PICK(size); PICK(diet);
                PICK(jaw); PICK(cold_tol); PICK(aquatic); PICK(salt_tol); PICK(resp);
                PICK(detritus); PICK(sexual); PICK(arboreal); PICK(brood);
                PICK(lifespan); PICK(toxin); PICK(tox_resist);
                PICK(burrow); PICK(nocturnal); PICK(store_fat); PICK(repro_alloc);
                PICK(estrus);
                #undef PICK
                for(int k=0;k<N_W;k++)
                    c.w[k]=(dist01(rng)<0.5f)? a.w[k] : mate->w[k];
                for(int k=0;k<N_NODE;k++)
                    c.ngain[k]=(dist01(rng)<0.5f)? a.ngain[k] : mate->ngain[k];
                for(int k=0;k<N_SENSOR;k++)
                    c.sen[k]=(dist01(rng)<0.5f)? a.sen[k] : mate->sen[k];
                for(int k=0;k<N_IMMUNE;k++)
                    c.im[k]=(dist01(rng)<0.5f)? a.im[k] : mate->im[k];
            }
            c.x=a.x+frand(-2.f,2.f); c.y=a.y+frand(-2.f,2.f);
            if(c.x<0.f) c.x=0.f;     if(c.x>=GRID_W) c.x=GRID_W-1.f;
            if(c.y<0.f) c.y=0.f;     if(c.y>=GRID_H) c.y=GRID_H-1.f;
            c.dir=frand(0.f,6.2832f); c.age=0; c.eaten=0;
            c.body=0.10f; c.lin=a.lin;
            c.id=next_animal_id++;
            // 変異は、組み換えで両親から選んだ値に掛ける
            c.speed   =mutate(c.speed,   0.0f,1.0f);
            c.turn    =mutate(c.turn,    0.0f,1.0f);
            c.reach   =mutate(c.reach,   0.05f,5.0f);
            c.size    =mutate(c.size,    0.10f,4.0f);
            c.diet    =mutate(c.diet,    0.0f,1.0f);
            c.jaw     =mutate(c.jaw,     0.0f,2.0f);
            c.cold_tol=mutate(c.cold_tol,0.0f,1.0f);
            c.aquatic =mutate(c.aquatic, 0.0f,1.0f);
            c.salt_tol=mutate(c.salt_tol,0.0f,1.0f);
            c.resp     =mutate(c.resp,     0.0f,1.0f);
            c.detritus =mutate(c.detritus, 0.0f,1.0f);
            c.sexual   =mutate(c.sexual,   0.0f,1.0f);
            c.arboreal =mutate(c.arboreal, 0.0f,1.0f);
            c.perch    =0.f;                       // 子は地面から始める
            c.aerobic  =1.f;
            for(int k=0;k<N_IMMUNE;k++) c.im[k]=mutate(c.im[k],0.0f,1.0f);
            // 変異率そのものが遺伝する。環境が荒れた後に上がりうる
            float mr=0.4f+a.mut_rate*MUT_RATE_MAX;
            c.mut_rate =mutate_r(a.mut_rate, 0.f,1.f,mr);
            c.bite_size=mutate_r(a.bite_size,0.f,1.f,mr);
            c.digest   =mutate_r(a.digest,   0.f,1.f,mr);
            for(int k=0;k<N_REC;k++) c.rec[k]=0.f;   // 記憶は引き継がない
            c.lesion=a.lesion;                         // 対照の印は変異させず、母からそのまま写す
            c.mother=a.id; c.father=mate?mate->id:-1; c.gen=a.gen+1;
            life_note(a.id,g_frame,mate?"有性繁殖":"無性繁殖",(float)nb);
            c.brood   =mutate(c.brood,   0.0f,1.0f);
            c.lifespan   =mutate(c.lifespan,   0.0f,1.0f);
            c.toxin      =mutate(c.toxin,      0.0f,1.0f);
            c.tox_resist =mutate(c.tox_resist, 0.0f,1.0f);
            c.burrow     =mutate(c.burrow,     0.0f,1.0f);
            c.nocturnal  =mutate(c.nocturnal,  0.0f,1.0f);
            c.store_fat  =mutate(c.store_fat,  0.0f,1.0f);
            c.repro_alloc=mutate(c.repro_alloc,0.0f,1.0f);
            c.estrus     =mutate(c.estrus,     0.0f,1.0f);
            c.fat=0.f; c.depth=0.f;
            for(int k=0;k<N_W;k++) c.w[k]=mutate(c.w[k],-3.0f,3.0f);
            // ニューロンの有効・無効。稀に大きく飛んで、層の追加や削除が起きる
            for(int k=0;k<N_NODE;k++)
                c.ngain[k]=(dist01(rng)<0.015f) ? frand(0.f,1.6f)
                                                : mutate(c.ngain[k],0.f,1.6f);
            for(int k=0;k<N_SENSOR;k++){
                c.sen[k].angle=mutate(c.sen[k].angle,-3.1416f,3.1416f);
                c.sen[k].range=mutate(c.sen[k].range, 0.0f,SENSOR_MAX);
                // 受容体の特異性は、稀に大きく変わる
                c.sen[k].tune = (dist01(rng)<0.02f) ? frand(0.f,1.f)
                                                    : mutate(c.sen[k].tune,0.0f,1.0f);
            }
            babies.push_back(c);
            }
        }
    }
    for(auto& b:babies) animals.push_back(b);
    if(g_track_f||g_trig_f) for(const auto& a:animals) if(!a.alive) field_death(a,g_frame);
    animals.erase(std::remove_if(animals.begin(),animals.end(),
        [](const Animal& a){ return !a.alive; }),animals.end());
    // 糞と落とした葉を、そのマスの植物由来の死骸に足す。無ければマスの中心に新しく置く
    for(int c:dung_cells){
        int pk=-1;
        for(int t=cgrid.b(c);t<cgrid.e(c);t++)
            if(corpses[cgrid.items[t]].origin==0){ pk=cgrid.items[t]; break; }
        if(pk>=0){ corpses[pk].energy+=dung[c]; corpses[pk].n+=dung_n[c]; }
        else corpses.push_back({(c%GRID_W)+0.5f,(c/GRID_W)+0.5f,dung[c],0,dung_n[c]});
        dung[c]=0.f; dung_n[c]=0.f;
    }
    dung_cells.clear();
    // 分解: 死骸が朽ちると、炭素は大気へ戻る
    for(auto& cp:corpses){
        float rot = cp.energy*ROT_RATE;
        cp.energy -= rot;
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H){
            int ri=idx(gx,gy);
            a_burn(co2.data(),o2.data(),ri,rot);   // 風化も酸化。O2 を 1:1 で使う
            // 朽ちた割合だけ、死骸の窒素が無機化する。海では一部だけが有光層に残り、残りは沈む
            float rn=cp.n*ROT_RATE;
            cp.n-=rn;
            soil_n[ri]+= rn*(sea[ri]?0.45f:1.f);
            if(sea[ri]) g_bud.n_sink-=rn*0.55f;
        }
    }
    // 小さくなった死骸は分解し尽くされたとみなし、残りの炭素を CO2 に出してから消す
    for(const auto& cp:corpses){
        if(cp.energy>=0.02f) continue;
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H){
            int ri=idx(gx,gy);
            soil_n[ri]+=cp.n;                                         // 残りの窒素も無機化する(炭素の無い糞の窒素も)
            if(cp.energy<=0.f) continue;
            co2_add(ri,cp.energy*CO2_PER_ENERGY);
            o2_add(ri,-std::min(cp.energy*CO2_PER_ENERGY,o2_at(ri)));   // 好気的な分解は O2 を 1:1 で使う
        }
    }
    corpses.erase(std::remove_if(corpses.begin(),corpses.end(),
        [](const Corpse& c){ return c.energy<0.02f; }),corpses.end());
}

