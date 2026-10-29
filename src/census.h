// 研究の測定の道具 5: 全体の記録(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ===== 研究の測定の道具 5: 全体の記録 =====
// --census N: N tick ごとに census.csv へ、個体数・遺伝子の分布・選択勾配・環境の出来事を書く。
// 条件どうしを比べるため、どの実行でも同じ列・同じ間隔で書く(1 行 1 値: t,group,name,stat,value)。
// 書き出すだけで、乱数を使わず、世界の状態にも触れない。だから記録の有無で結果は変わらない
int   g_census_every=0;
FILE* g_census_f=nullptr;
static void census_open(){
    if(g_census_every<=0) return;
    g_census_f=fopen("census.csv","w");
    if(g_census_f) fprintf(g_census_f,"t,group,name,stat,value\n");
}
static void census_close(){ if(g_census_f){ fclose(g_census_f); g_census_f=nullptr; } }
// 一つの遺伝子の分布: 平均・標準偏差・10/50/90 百分位
static void census_dist(int t,const char* grp,const char* name,std::vector<double>& v){
    FILE* f=g_census_f; size_t n=v.size();
    double m=0,var=0;
    for(double x:v) m+=x;
    if(n) m/=n;
    for(double x:v) var+=(x-m)*(x-m);
    double sd=n>1?std::sqrt(var/n):0.0;
    fprintf(f,"%d,%s,%s,mean,%.6g\n%d,%s,%s,sd,%.6g\n",t,grp,name,m,t,grp,name,sd);
    const int Q[3]={10,50,90};
    for(int q:Q){
        double x=0;
        if(n){ size_t k=std::min(n-1,(size_t)(q*(n-1)/100.0+0.5));
               std::nth_element(v.begin(),v.begin()+k,v.end()); x=v[k]; }
        fprintf(f,"%d,%s,%s,p%d,%.6g\n",t,grp,name,q,x);
    }
}
static void census_write(int t){
    FILE* f=g_census_f;
    if(!f) return;
    // 個体数
    int np=0,ns=0,na=0,nm=0,lp=0,la=0;
    for(const auto& p:plants) if(p.alive){ if(p.seed_wait>0) ns++; else np++; }
    for(const auto& a:animals) if(a.alive) na++;
    for(const auto& m:microbes) if(m.alive) nm++;
    for(const auto& L:lineages) if(L.alive){ if(L.kind==0) lp++; else la++; }
    fprintf(f,"%d,pop,plant,n,%d\n%d,pop,seed,n,%d\n%d,pop,animal,n,%d\n%d,pop,microbe,n,%d\n%d,pop,corpse,n,%d\n",
            t,np,t,ns,t,na,t,nm,t,(int)corpses.size());
    fprintf(f,"%d,pop,lineage_plant,n,%d\n%d,pop,lineage_animal,n,%d\n",t,lp,t,la);
    fprintf(f,"%d,pop,repro_plant,n,%ld\n%d,pop,repro_animal,n,%ld\n",t,cs_p_nr,t,cs_a_nr);
    fprintf(f,"%d,pop,death_starve,cum,%ld\n%d,pop,death_cold,cum,%ld\n%d,pop,death_preyed,cum,%ld\n%d,pop,death_microbe,cum,%ld\n%d,pop,death_slide,cum,%ld\n",
            t,d_starve,t,d_cold,t,d_preyed,t,d_micro,t,d_slide);
    // 遺伝子の分布と選択勾配(繁殖した親の平均 − 集団の平均。sd で割った値も出す)
    std::vector<double> col;
    double v[std::max(N_SG_A,N_SG_P)];
    for(int k=0;k<N_SG_A;k++){
        col.clear();
        for(const auto& a:animals) if(a.alive){ sg_vec_a(a,v); col.push_back(v[k]); }
        double m=0; for(double x:col) m+=x;
        if(!col.empty()) m/=col.size();
        double var=0; for(double x:col) var+=(x-m)*(x-m);
        double sd=col.size()>1?std::sqrt(var/col.size()):0.0;
        census_dist(t,"geneA",SG_A_NAME[k],col);
        if(cs_a_nr>0 && !col.empty()){
            double d=cs_a_rep[k]/cs_a_nr-m;
            fprintf(f,"%d,selA,%s,diff,%.6g\n%d,selA,%s,std,%.6g\n",t,SG_A_NAME[k],d,t,SG_A_NAME[k],sd>0?d/sd:0.0);
        }
    }
    for(int k=0;k<N_SG_P;k++){
        col.clear();
        for(const auto& p:plants) if(p.alive&&p.seed_wait<=0){ sg_vec_p(p,v); col.push_back(v[k]); }
        double m=0; for(double x:col) m+=x;
        if(!col.empty()) m/=col.size();
        double var=0; for(double x:col) var+=(x-m)*(x-m);
        double sd=col.size()>1?std::sqrt(var/col.size()):0.0;
        census_dist(t,"geneP",SG_P_NAME[k],col);
        if(cs_p_nr>0 && !col.empty()){
            double d=cs_p_rep[k]/cs_p_nr-m;
            fprintf(f,"%d,selP,%s,diff,%.6g\n%d,selP,%s,std,%.6g\n",t,SG_P_NAME[k],d,t,SG_P_NAME[k],sd>0?d/sd:0.0);
        }
    }
    col.clear(); for(const auto& m:microbes) if(m.alive) col.push_back(m.rate);
    census_dist(t,"geneM","rate",col);
    col.clear(); for(const auto& m:microbes) if(m.alive) col.push_back(m.cold_tol);
    census_dist(t,"geneM","cold",col);
    // 環境の出来事(--record と同じ種類。前の記録からの件数と強さの和)と、外から与える条件
    for(int k=0;k<EV_N;k++)
        fprintf(f,"%d,event,%s,n,%ld\n%d,event,%s,strength_sum,%.6g\n",t,EV_NAME[k],g_cs_ev_n[k],t,EV_NAME[k],g_cs_ev_s[k]);
    double tl=0; int nl=0;
    for(int i=0;i<GRID_W*GRID_H;i++) if(!sea[i]){ tl+=temper[i]; nl++; }
    fprintf(f,"%d,env,land_temp,mean,%.6g\n%d,env,day_light,mean,%.6g\n%d,env,sun_dec_deg,value,%.6g\n",
            t,nl?tl/nl:0.0,t,(double)day_light,t,sun_dec*180.0/3.14159265);
    fflush(f);
    for(int k=0;k<N_SG_A;k++) cs_a_rep[k]=0;
    for(int k=0;k<N_SG_P;k++) cs_p_rep[k]=0;
    cs_a_nr=cs_p_nr=0;
    for(int k=0;k<EV_N;k++){ g_cs_ev_n[k]=0; g_cs_ev_s[k]=0; }
}

