// 研究の測定の道具: 神経網の書き出しと試験器(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ===== 研究の測定の道具: 神経網の書き出しと試験器 =====
// 書き出しは乱数を使わず、世界の状態にも触れず、標準出力にも何も書かない。
// だからシミュレーションの結果とログは、書き出しの有無で変わらない
int g_dump_every=0;     // --dump-brains 間隔(tick)。0 なら書き出さない
int g_dump_max=1000;    // --dump-max 一回に書き出す動物の上限(生きている順に等間隔で選ぶ)
std::string g_series_path;   // --assay-series ファイル。書き出すたびに、同じ動物を試験器にかけて 1 行ずつ足す
static void assay_series(int t,const std::vector<int>& pick);

void dump_brains(int t){
    std::vector<int> alive_idx;
    for(int i=0;i<(int)animals.size();i++) if(animals[i].alive) alive_idx.push_back(i);
    if(alive_idx.empty()) return;
    int stride=std::max(1,(int)((alive_idx.size()+g_dump_max-1)/g_dump_max));
    int n=0; for(size_t q=0;q<alive_idx.size();q+=stride) n++;
    char path[64]; snprintf(path,sizeof(path),"brains_t%08d.txt",t);
    FILE* f=fopen(path,"w");
    if(!f) return;
    fprintf(f,"# brains t=%d n=%d N_SENSOR=%d N_TUNE=%d N_IN=%d N_OUT=%d N_REC=%d MAX_LAYER=%d MAX_NODE=%d N_W=%d\n",
            t,n,N_SENSOR,N_TUNE,N_IN,N_OUT,N_REC,MAX_LAYER,MAX_NODE,N_W);
    fprintf(f,"# id gen lin body size | sensor(angle range tune) x%d | ngain x%d | w x%d  (16進の浮動小数で正確に保存)\n",
            N_SENSOR,N_NODE,N_W);
    std::vector<int> pick;
    for(size_t q=0;q<alive_idx.size();q+=stride){
        const Animal& a=animals[alive_idx[q]];
        pick.push_back(alive_idx[q]);
        fprintf(f,"%d %d %d %a %a",a.id,a.gen,a.lin,a.body,a.size);
        for(int k=0;k<N_SENSOR;k++) fprintf(f," %a %a %a",a.sen[k].angle,a.sen[k].range,a.sen[k].tune);
        for(int k=0;k<N_NODE;k++)   fprintf(f," %a",a.ngain[k]);
        for(int k=0;k<N_W;k++)      fprintf(f," %a",a.w[k]);
        fprintf(f,"\n");
    }
    fclose(f);
    if(!g_series_path.empty()) assay_series(t,pick);
}

// 試験器(--assay ファイル)。書き出した神経網に標準の刺激を与え、
// Anderson & Adolphs (2014) の四つの性質(価・持続・強度・一般化)を測る。
// 世界は作らず、神経網だけを brain_forward で動かす。1 段 = 動物の 1 回の更新(ANIMAL_STRIDE tick)
struct AssayStim { const char* name; int band; int kind; };   // kind: +1 脅威, -1 餌, 0 その他
static const AssayStim ASSAY_STIM[]={
    {"predator",2,+1}, {"fire",4,+1}, {"smoke",10,+1},
    {"plant",0,-1},    {"corpse",3,-1},
    {"prey",1,0},      {"kin",11,0},
};
const int   N_ASSAY_STIM=(int)(sizeof(ASSAY_STIM)/sizeof(ASSAY_STIM[0]));
const int   ASSAY_PRE=60, ASSAY_ON=20, ASSAY_OFF=120;
const float ASSAY_EPS=0.02f;      // これより小さい再帰回路の変化は「反応なし」
const int   ASSAY_PERSIST_MIN=3;  // 持続とみなす段数(刺激が去った後も変化が 10% より大きく残る段)
const float ASSAY_FANOUT_MIN=0.05f;

struct AssayRun { float r0[N_REC], ron[N_REC], o0[4], post_out[4]; int persist; float dev; };

static inline int assay_band(float tune){ int b=(int)(tune*N_TUNE); if(b>N_TUNE-1) b=N_TUNE-1; if(b<0) b=0; return b; }

// 刺激は感覚器の向きによらず、合う感覚器すべてに届くとみなす(単純化)。
// 世界では最初のレイの段で見つけたとき入力は 信号×(1−1/3) なので、強さ s を 2/3·s として与える
static void assay_inputs(const Animal& a,int band,float s,const float* rec,bool lesion,float* in){
    for(int k=0;k<N_SENSOR;k++){
        in[k]=0.f;
        if(a.sen[k].range<SENSOR_MIN) continue;
        if(band>=0 && assay_band(a.sen[k].tune)==band) in[k]=s*(2.f/3.f);
    }
    in[N_SENSOR+0]=0.5f;                       // エネルギー: 繁殖の閾値の半分
    in[N_SENSOR+1]=clamp01((15.f+30.f)/70.f);  // 気温 15 °C
    in[N_SENSOR+2]=0.4f;                       // 土の水
    in[N_SENSOR+3]=std::min(1.f,a.body/3.f);
    in[N_SENSOR+4]=1.f;                        // バイアス
    for(int k=0;k<N_REC;k++){                  // 陰性対照(lesion)は再帰の入力を 0 に切る
        float r=lesion?0.f:rec[k];
        if(!(r>-1.5f&&r<1.5f)) r=0.f;
        in[N_SENSOR+5+k]=r;
    }
}

static AssayRun assay_run(const Animal& a,int band,float s,bool lesion){
    AssayRun R{};
    float rec[N_REC]={0}, in[N_IN], out[N_OUT]; int nn=0,ns=0;
    auto step=[&](int b,float sv){
        assay_inputs(a,b,sv,rec,lesion,in);
        brain_forward(a,in,out,nn,ns);
        for(int k=0;k<N_REC;k++) rec[k]=out[4+k];
    };
    for(int t=0;t<ASSAY_PRE;t++) step(-1,0.f);          // 刺激なしで落ち着かせる
    for(int k=0;k<N_REC;k++) R.r0[k]=rec[k];
    for(int o=0;o<4;o++) R.o0[o]=out[o];
    for(int t=0;t<ASSAY_ON;t++) step(band,s);           // 刺激を与える
    float d2=0.f;
    for(int k=0;k<N_REC;k++){ R.ron[k]=rec[k]; d2+=(rec[k]-R.r0[k])*(rec[k]-R.r0[k]); }
    R.dev=std::sqrt(d2);
    bool back=false;
    for(int t=0;t<ASSAY_OFF;t++){                       // 刺激を止める
        step(-1,0.f);
        float e2=0.f;
        for(int k=0;k<N_REC;k++) e2+=(rec[k]-R.r0[k])*(rec[k]-R.r0[k]);
        if(!back){ if(R.dev>ASSAY_EPS && std::sqrt(e2)>0.1f*R.dev) R.persist++; else back=true; }
        // 刺激が去った直後の 1 段を除き、続く 10 段の行動の出力のずれ。
        // 感覚の入力はもう 0 なので、ずれは再帰回路(内部の状態)を通してしか起きない
        if(t>=1 && t<11) for(int o=0;o<4;o++) R.post_out[o]+=std::fabs(out[o]-R.o0[o])/10.f;
    }
    return R;
}

static float assay_cos(const float* a,const float* b){
    float ab=0,aa=0,bb=0;
    for(int k=0;k<N_REC;k++){ ab+=a[k]*b[k]; aa+=a[k]*a[k]; bb+=b[k]*b[k]; }
    if(aa<1e-12f||bb<1e-12f) return 0.f;
    return ab/std::sqrt(aa*bb);
}

struct AssayResult {
    bool settled; int sens_mask, n_threat_resp;
    float vcos, scorr, fanin; int persist, fanout;
    bool v,p,s,g;
    float dev[N_ASSAY_STIM]; int pers[N_ASSAY_STIM];
};

static AssayResult assay_brain(const Animal& a,bool lesion){
    AssayResult X{};
    for(int k=0;k<N_SENSOR;k++) if(a.sen[k].range>=SENSOR_MIN) X.sens_mask|=1<<assay_band(a.sen[k].tune);
    AssayRun sham=assay_run(a,-1,0.f,lesion);            // 刺激なしでも揺れ続ける神経網は除く
    X.settled = sham.dev<=ASSAY_EPS;
    AssayRun R[N_ASSAY_STIM];
    float vthr[N_REC]={0}, vfood[N_REC]={0}; int nthr=0,nfood=0;
    int strongest=-1; float sdev=0.f;
    for(int i=0;i<N_ASSAY_STIM;i++){
        R[i]=assay_run(a,ASSAY_STIM[i].band,1.f,lesion);
        X.dev[i]=R[i].dev; X.pers[i]=R[i].persist;
        if(R[i].dev<=ASSAY_EPS) continue;
        if(ASSAY_STIM[i].kind>0){
            for(int k=0;k<N_REC;k++) vthr[k]+=R[i].ron[k]-R[i].r0[k];
            nthr++;
            if(R[i].dev>sdev){ sdev=R[i].dev; strongest=i; }
            if(R[i].persist>X.persist) X.persist=R[i].persist;
        } else if(ASSAY_STIM[i].kind<0){
            for(int k=0;k<N_REC;k++) vfood[k]+=R[i].ron[k]-R[i].r0[k];
            nfood++;
        }
    }
    X.n_threat_resp=nthr;
    // 価: 脅威と餌で、再帰回路の値が逆向きに動くか
    X.vcos = (nthr&&nfood)? assay_cos(vthr,vfood) : 0.f;
    X.v = nthr>0 && nfood>0 && X.vcos<-0.5f;
    // 持続: 刺激が去った後も変化が残るか
    X.p = X.persist>=ASSAY_PERSIST_MIN;
    // 強度: 最も強く反応した脅威で、刺激を強めるほど変化が大きくなるか
    X.scorr=0.f;
    if(strongest>=0){
        const float lv[4]={0.25f,0.5f,0.75f,1.f};
        float dv[4], md=0,ms=0;
        for(int q=0;q<4;q++){ dv[q]=assay_run(a,ASSAY_STIM[strongest].band,lv[q],lesion).dev; md+=dv[q]/4; ms+=lv[q]/4; }
        float cxy=0,cxx=0,cyy=0;
        for(int q=0;q<4;q++){ cxy+=(lv[q]-ms)*(dv[q]-md); cxx+=(lv[q]-ms)*(lv[q]-ms); cyy+=(dv[q]-md)*(dv[q]-md); }
        X.scorr = (cyy>1e-12f)? cxy/std::sqrt(cxx*cyy) : 0.f;
    }
    X.s = strongest>=0 && X.scorr>0.9f;
    // 一般化: 違う脅威で同じ向きの変化が起きるか(入り口)、その変化が複数の行動を変えるか(出口)
    X.fanin=0.f; int np=0;
    for(int i=0;i<N_ASSAY_STIM;i++) for(int j=i+1;j<N_ASSAY_STIM;j++){
        if(ASSAY_STIM[i].kind<=0||ASSAY_STIM[j].kind<=0) continue;
        if(R[i].dev<=ASSAY_EPS||R[j].dev<=ASSAY_EPS) continue;
        float di[N_REC],dj[N_REC];
        for(int k=0;k<N_REC;k++){ di[k]=R[i].ron[k]-R[i].r0[k]; dj[k]=R[j].ron[k]-R[j].r0[k]; }
        X.fanin+=assay_cos(di,dj); np++;
    }
    if(np) X.fanin/=np;
    X.fanout=0;
    if(strongest>=0) for(int o=0;o<4;o++) if(R[strongest].post_out[o]>ASSAY_FANOUT_MIN) X.fanout++;
    X.g = np>0 && X.fanin>0.5f && X.fanout>=2;
    return X;
}

// 陽性対照: 恐れに似た回路を手で組んだ神経網。
// 捕食者・火・煙で再帰回路 0 が上がり、餌で下がり、自分への結合 0.9 でゆっくり戻り、
// その値が旋回・前進を上げ、採食を下げる。四つの性質がすべて検出されるはず
static Animal assay_positive_control(){
    Animal a{};
    a.body=0.5f; a.size=0.5f; a.id=-1; a.gen=0; a.lin=-1;
    const int bands[N_SENSOR]={2,4,10,0,3};
    for(int k=0;k<N_SENSOR;k++){ a.sen[k].angle=0.f; a.sen[k].range=10.f; a.sen[k].tune=(bands[k]+0.5f)/N_TUNE; }
    for(int k=0;k<N_NODE;k++) a.ngain[k]=0.f;    // 隠れ層なし。直結だけで組む
    for(int k=0;k<N_W;k++) a.w[k]=0.f;
    float* wd=&a.w[W_L0+W_LL+W_OUT];
    const int R0=N_SENSOR+5;                      // 再帰回路 0 の入力の位置
    wd[4*N_IN+0]= 3.0f; wd[4*N_IN+1]= 3.0f; wd[4*N_IN+2]= 3.0f;   // 脅威 → 上がる
    wd[4*N_IN+3]=-3.0f; wd[4*N_IN+4]=-3.0f;                        // 餌 → 下がる
    wd[4*N_IN+R0]=0.9f;                                            // ゆっくり戻る
    wd[0*N_IN+R0]= 1.0f;    // 旋回
    wd[1*N_IN+R0]= 2.0f;    // 前進(逃げる)
    wd[2*N_IN+R0]=-2.0f;    // 採食をやめる
    return a;
}

static void assay_print(const char* label,const AssayResult& X,int id){
    printf("  %-10s id %-7d settled %d  価 %s(cos %+.2f)  持続 %s(%d 段)  強度 %s(相関 %+.2f)  一般化 %s(入口 %+.2f 出口 %d)\n",
           label,id,X.settled?1:0,
           X.v?"○":"×",X.vcos, X.p?"○":"×",X.persist, X.s?"○":"×",X.scorr, X.g?"○":"×",X.fanin,X.fanout);
}

// 書き出しファイル(--dump-brains)を読む
static bool assay_load_dump(const std::string& path,std::vector<Animal>& out,std::string& info){
    std::ifstream fin(path);
    if(!fin){ info="開けない"; return false; }
    std::string line;
    while(std::getline(fin,line)){
        if(line.empty()) continue;
        if(line[0]=='#'){
            int t,n,ns,nt,ni,no,nr,ml,mn,nw;
            if(sscanf(line.c_str(),"# brains t=%d n=%d N_SENSOR=%d N_TUNE=%d N_IN=%d N_OUT=%d N_REC=%d MAX_LAYER=%d MAX_NODE=%d N_W=%d",
                      &t,&n,&ns,&nt,&ni,&no,&nr,&ml,&mn,&nw)==10){
                if(ns!=N_SENSOR||nt!=N_TUNE||ni!=N_IN||no!=N_OUT||nr!=N_REC||ml!=MAX_LAYER||mn!=MAX_NODE||nw!=N_W){
                    info="神経網の形がこのプログラムと違う"; return false;
                }
                info="t="+std::to_string(t);
            }
            continue;
        }
        Animal a{};
        const char* p=line.c_str(); char* e;
        a.id=(int)strtol(p,&e,10); p=e; a.gen=(int)strtol(p,&e,10); p=e; a.lin=(int)strtol(p,&e,10); p=e;
        a.body=(float)strtod(p,&e); p=e; a.size=(float)strtod(p,&e); p=e;
        for(int k=0;k<N_SENSOR;k++){
            a.sen[k].angle=(float)strtod(p,&e); p=e;
            a.sen[k].range=(float)strtod(p,&e); p=e;
            a.sen[k].tune =(float)strtod(p,&e); p=e;
        }
        for(int k=0;k<N_NODE;k++){ a.ngain[k]=(float)strtod(p,&e); p=e; }
        for(int k=0;k<N_W;k++){ a.w[k]=(float)strtod(p,&e); p=e; }
        out.push_back(a);
    }
    return true;
}

// セーブ(P3SAVE2)から動物の神経網だけを読む。列は名前で拾うので、古い版のセーブでも
// 神経網の列(sen*_angle/range/tune, ng*, w*)と body・size・id・gen・lin があれば読める。
// max 匹を超える分は、並び順に等間隔で選ぶ
static bool assay_load_save(const std::string& path,int max_n,std::vector<Animal>& out,std::string& info){
    std::ifstream fin(path);
    if(!fin){ info="開けない"; return false; }
    std::string line;
    if(!std::getline(fin,line) || line.rfind("P3SAVE2",0)!=0){ info="P3SAVE2 ではない"; return false; }
    std::vector<std::string> names; long frame=-1;
    while(std::getline(fin,line)){
        if(line.rfind("HEADER ",0)==0){ frame=strtol(line.c_str()+7,nullptr,10); continue; }
        if(line.rfind("SCHEMA animal ",0)==0){
            const char* p=line.c_str()+14; char* e; long n=strtol(p,&e,10); p=e;
            names.clear();
            for(long k=0;k<n;k++){
                while(*p==' ') p++;
                const char* q=p; while(*q && *q!=' ') q++;
                names.emplace_back(p,q-p); p=q;
            }
            continue;
        }
        if(line.rfind("ANIMALS ",0)==0) break;
    }
    if(names.empty()){ info="動物の SCHEMA が無い"; return false; }
    long n_anim=strtol(line.c_str()+8,nullptr,10);
    // 列の対応: 0..N_W-1 は w、N_W..N_W+N_NODE-1 は ng、その後に sen と基本の列
    std::vector<int> col(names.size(),-1);
    const int C_SEN=N_W+N_NODE, C_BODY=C_SEN+3*N_SENSOR, C_SIZE=C_BODY+1, C_ID=C_BODY+2, C_GEN=C_BODY+3, C_LIN=C_BODY+4, C_N=C_BODY+5;
    std::vector<char> got(C_N,0);
    for(size_t c=0;c<names.size();c++){
        const std::string& nm=names[c];
        int k;
        if(sscanf(nm.c_str(),"w%d",&k)==1 && nm=="w"+std::to_string(k) && k>=0&&k<N_W) col[c]=k;
        else if(sscanf(nm.c_str(),"ng%d",&k)==1 && nm=="ng"+std::to_string(k) && k>=0&&k<N_NODE) col[c]=N_W+k;
        else if(nm.rfind("sen",0)==0 && sscanf(nm.c_str(),"sen%d_",&k)==1 && k>=0&&k<N_SENSOR){
            std::string tail=nm.substr(nm.find('_')+1);
            if(tail=="angle") col[c]=C_SEN+3*k; else if(tail=="range") col[c]=C_SEN+3*k+1; else if(tail=="tune") col[c]=C_SEN+3*k+2;
        }
        else if(nm=="body") col[c]=C_BODY; else if(nm=="size") col[c]=C_SIZE;
        else if(nm=="id") col[c]=C_ID; else if(nm=="gen") col[c]=C_GEN; else if(nm=="lin") col[c]=C_LIN;
        if(col[c]>=0) got[col[c]]=1;
    }
    int missing=0; for(int k=0;k<C_BODY;k++) if(!got[k]) missing++;
    if(missing){ info="神経網の列が "+std::to_string(missing)+" 個足りない(神経網の形が違う)"; return false; }
    long stride=std::max(1L,(n_anim+max_n-1)/max_n);
    std::vector<double> v(C_N);
    for(long i=0;i<n_anim;i++){
        if(!std::getline(fin,line)){ info="動物の行が途中で終わった"; return false; }
        if(i%stride) continue;
        std::fill(v.begin(),v.end(),0.0);
        const char* p=line.c_str(); char* e;
        for(size_t c=0;c<names.size();c++){
            double x=strtod(p,&e); if(e==p) break; p=e;
            if(col[c]>=0) v[col[c]]=x;
        }
        Animal a{};
        for(int k=0;k<N_W;k++) a.w[k]=(float)v[k];
        for(int k=0;k<N_NODE;k++) a.ngain[k]=(float)v[N_W+k];
        for(int k=0;k<N_SENSOR;k++){ a.sen[k].angle=(float)v[C_SEN+3*k]; a.sen[k].range=(float)v[C_SEN+3*k+1]; a.sen[k].tune=(float)v[C_SEN+3*k+2]; }
        a.body=(float)v[C_BODY]; a.size=(float)v[C_SIZE];
        a.id=(int)v[C_ID]; a.gen=(int)v[C_GEN]; a.lin=(int)v[C_LIN];
        out.push_back(a);
    }
    info="t="+std::to_string(frame)+" 動物 "+std::to_string(n_anim)+" 匹から "+std::to_string(out.size())+" 匹";
    return true;
}

static int run_assay(const std::string& path,const std::string& save_path,int max_n,const std::string& out_path){
    printf("assay: 1 段 = %d tick = %.1f 分。前 %d 段・刺激 %d 段・後 %d 段\n",
           ANIMAL_STRIDE,(double)(ANIMAL_STRIDE*DT),ASSAY_PRE,ASSAY_ON,ASSAY_OFF);
    // 陽性対照と、その陰性対照(再帰を切ったもの)
    Animal pc=assay_positive_control();
    printf("対照:\n");
    assay_print("陽性",assay_brain(pc,false),-1);
    assay_print("陽性-切断",assay_brain(pc,true),-1);
    if(path=="control") return 0;

    std::vector<Animal> pop; std::string info; bool ok;
    if(!save_path.empty()) ok=assay_load_save(save_path,max_n,pop,info);
    else                   ok=assay_load_dump(path,pop,info);
    if(!ok){ fprintf(stderr,"assay: %s: %s\n",(save_path.empty()?path:save_path).c_str(),info.c_str()); return 2; }
    printf("読み込み: %s\n",info.c_str());

    FILE* fo=fopen(out_path.empty()? "assay.csv" : out_path.c_str(),"w");
    if(!fo){ fprintf(stderr,"assay: cannot open output\n"); return 2; }
    fprintf(fo,"id,gen,lin,mode,settled,sens_mask,n_threat_resp,valence_cos,persist_steps,scal_corr,fanin_cos,fanout,v,p,s,g,all");
    for(int i=0;i<N_ASSAY_STIM;i++) fprintf(fo,",dev_%s",ASSAY_STIM[i].name);
    fprintf(fo,"\n");
    int cnt[2][5]={{0}}; int nset[2]={0,0};
    for(const Animal& a:pop){
        for(int m=0;m<2;m++){
            AssayResult X=assay_brain(a,m==1);
            bool all=X.v&&X.p&&X.s&&X.g;
            fprintf(fo,"%d,%d,%d,%s,%d,%d,%d,%.4f,%d,%.4f,%.4f,%d,%d,%d,%d,%d,%d",
                    a.id,a.gen,a.lin,m?"lesion":"intact",X.settled?1:0,X.sens_mask,X.n_threat_resp,
                    X.vcos,X.persist,X.scorr,X.fanin,X.fanout,X.v,X.p,X.s,X.g,all);
            for(int i=0;i<N_ASSAY_STIM;i++) fprintf(fo,",%.4f",X.dev[i]);
            fprintf(fo,"\n");
            if(!X.settled) continue;
            nset[m]++;
            cnt[m][0]+=X.v; cnt[m][1]+=X.p; cnt[m][2]+=X.s; cnt[m][3]+=X.g; cnt[m][4]+=all;
        }
    }
    fclose(fo);
    printf("神経網 %d 個(落ち着いたもの: そのまま %d / 再帰を切ったもの %d)\n",(int)pop.size(),nset[0],nset[1]);
    const char* nm[5]={"価","持続","強度","一般化","四つすべて"};
    for(int q=0;q<5;q++)
        printf("  %-10s そのまま %5d (%5.1f%%)   再帰を切ったもの %5d (%5.1f%%)\n",nm[q],
               cnt[0][q],nset[0]?100.0*cnt[0][q]/nset[0]:0.0, cnt[1][q],nset[1]?100.0*cnt[1][q]/nset[1]:0.0);
    printf("結果: %s\n",out_path.empty()? "assay.csv" : out_path.c_str());
    return 0;
}

// 測定の道具 6: 進化の時間変化(--assay-series ファイル、--dump-brains と一緒に使う)。
// 書き出した動物を、そのまま・再帰を切ったものの両方で試験器にかけ、四つの性質が検出された
// 割合を 1 回 1 行ずつ足す。試験器は乱数を使わず、世界の状態を読むだけで、標準出力にも書かない。
// だからシミュレーションの結果とログは変わらない(試験器の計算の時間だけ、その tick が遅くなる)
static void assay_series(int t,const std::vector<int>& pick){
    static FILE* fs=nullptr;
    if(!fs){
        fs=fopen(g_series_path.c_str(),"w");
        if(!fs) return;
        fprintf(fs,"t,mode,n,settled,v,p,s,g,all,v_pct,p_pct,s_pct,g_pct,all_pct,persist_mean\n");
    }
    for(int m=0;m<2;m++){
        int cnt[5]={0}, nset=0; double pers=0;
        for(int i:pick){
            AssayResult X=assay_brain(animals[i],m==1);
            if(!X.settled) continue;
            nset++;
            cnt[0]+=X.v; cnt[1]+=X.p; cnt[2]+=X.s; cnt[3]+=X.g; cnt[4]+=(X.v&&X.p&&X.s&&X.g);
            pers+=X.persist;
        }
        fprintf(fs,"%d,%s,%d,%d",t,m?"lesion":"intact",(int)pick.size(),nset);
        for(int q=0;q<5;q++) fprintf(fs,",%d",cnt[q]);
        for(int q=0;q<5;q++) fprintf(fs,",%.2f",nset?100.0*cnt[q]/nset:0.0);
        fprintf(fs,",%.3f\n",nset?pers/nset:0.0);
    }
    fflush(fs);
}

