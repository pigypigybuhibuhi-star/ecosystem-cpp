// 研究の測定の道具 4: 野外の記録(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ===== 研究の測定の道具 4: 野外の記録 =====
// 書き出すだけで、乱数を使わず、世界の状態にも触れない。だから記録の有無で結果は変わらない。
// (a) --track N [--track-from T]: t>=T に生まれた(放された)最初の N 匹を、生まれてから死ぬまで毎回の更新で track.csv に書く
// (b) --trigger-log W [--trigger-max M]: 脅威(大きい動物・炎・煙)か餌(植物・死骸)の感覚器の入力が 0 から立ち上がったとき、
//     その前 W 回と後 W 回の更新を trigger.csv に書く。一つの個体の記録中に次のきっかけが来ても、記録は延ばさない
int  g_track_n=0, g_track_from=0, g_track_taken=0;
int  g_pc_n=0, g_lesion_n=0;                   // --positive-control N --no-recurrent N: 最初に放す動物のうち対照にする数
std::unordered_map<int,char> g_track_ids;      // 追っている id
FILE* g_track_f=nullptr;
int  g_trig_w=0, g_trig_max=2000, g_trig_ep=0;
FILE* g_trig_f=nullptr;
const int TRIG_K=4+N_SENSOR+N_OUT;            // 1 回の更新の記録: t x y 体力 | 入力(感覚器) | 出力
struct TrigBuf { std::vector<float> ring; int n=0, post=-1, ep=-1, band=-1; unsigned prev=0; };
std::unordered_map<int,TrigBuf> g_trig;
static inline bool trig_band(int b){ return b==2||b==4||b==10||b==0||b==3; }
static inline int sensor_band(float tune){ int b=(int)(tune*N_TUNE); if(b>N_TUNE-1) b=N_TUNE-1; if(b<0) b=0; return b; }
static void field_frame(const Animal& a,int t,const float* in,const float* out,float* fr){
    fr[0]=(float)t; fr[1]=a.x; fr[2]=a.y; fr[3]=a.energy;
    for(int k=0;k<N_SENSOR;k++) fr[4+k]=in[k];
    for(int o=0;o<N_OUT;o++) fr[4+N_SENSOR+o]=out[o];
}
static void field_row(FILE* f,const float* fr){
    fprintf(f,"%d,%.3f,%.3f,%.5f",(int)fr[0],fr[1],fr[2],fr[3]);
    for(int k=4;k<TRIG_K;k++) fprintf(f,",%.4f",fr[k]);
    fprintf(f,"\n");
}
static void field_record(const Animal& a,int t,const float* in,const float* out){
    float fr[TRIG_K];
    field_frame(a,t,in,out,fr);
    if(g_track_f){
        bool on=g_track_ids.count(a.id)>0;
        if(!on && g_track_taken<g_track_n && t>=g_track_from && a.age<=ANIMAL_STRIDE){
            g_track_ids[a.id]=1; g_track_taken++; on=true;
        }
        if(on){
            fprintf(g_track_f,"upd,%d,%d,%d,%d,%.4f,%d,",a.id,a.gen,a.lin,a.mother,a.body,a.lesion?1:0);
            for(int k=0;k<N_SENSOR;k++) fprintf(g_track_f,"%d,",a.sen[k].range<SENSOR_MIN?-1:sensor_band(a.sen[k].tune));
            field_row(g_track_f,fr);
        }
    }
    if(g_trig_f){
        TrigBuf& B=g_trig[a.id];
        if(B.ring.empty()) B.ring.assign((size_t)g_trig_w*TRIG_K,0.f);
        std::copy(fr,fr+TRIG_K,&B.ring[(size_t)(B.n%g_trig_w)*TRIG_K]);
        B.n++;
        unsigned cur=0; int first=-1;
        for(int k=0;k<N_SENSOR;k++){
            if(a.sen[k].range<SENSOR_MIN || in[k]<=0.f) continue;
            int b=sensor_band(a.sen[k].tune);
            if(!trig_band(b)) continue;
            cur|=1u<<b;
            if(!(B.prev&(1u<<b)) && first<0) first=b;
        }
        B.prev=cur;
        if(B.post>0){
            fprintf(g_trig_f,"%d,%d,%d,%d,",B.ep,B.band,a.id,g_trig_w+1-B.post);
            field_row(g_trig_f,fr);
            B.post--;
        } else if(first>=0 && g_trig_ep<g_trig_max){
            B.ep=g_trig_ep++; B.band=first; B.post=g_trig_w;
            int m=std::min(B.n,g_trig_w);                  // 前の記録(今の更新を含む)を古い順に
            for(int q=m-1;q>=0;q--){
                int slot=(B.n-1-q)%g_trig_w;
                fprintf(g_trig_f,"%d,%d,%d,%d,",B.ep,B.band,a.id,-q);
                field_row(g_trig_f,&B.ring[(size_t)slot*TRIG_K]);
            }
        }
    }
}
// 死んだ個体を記録から外す(消す前に呼ぶ)
static void field_death(const Animal& a,int t){
    if(g_track_f && g_track_ids.erase(a.id))
        fprintf(g_track_f,"death,%d,%d,%d,%d,%.4f,%d,,,,,,%d,%.3f,%.3f,%.5f\n",
                a.id,a.gen,a.lin,a.mother,a.body,a.lesion?1:0,t,a.x,a.y,a.energy);
    if(g_trig_f) g_trig.erase(a.id);
}
static void field_open(){
    if(g_track_n>0){
        g_track_f=fopen("track.csv","w");
        if(g_track_f){
            fprintf(g_track_f,"event,id,gen,lin,mother,body,lesion");
            for(int k=0;k<N_SENSOR;k++) fprintf(g_track_f,",band%d",k);
            fprintf(g_track_f,",t,x,y,energy");
            for(int k=0;k<N_SENSOR;k++) fprintf(g_track_f,",in%d",k);
            fprintf(g_track_f,",turn,thrust,eat,attack");
            for(int k=0;k<N_REC;k++) fprintf(g_track_f,",rec%d",k);
            fprintf(g_track_f,"\n");
        }
    }
    if(g_trig_w>0){
        g_trig_f=fopen("trigger.csv","w");
        if(g_trig_f){
            fprintf(g_trig_f,"episode,band,id,step,t,x,y,energy");
            for(int k=0;k<N_SENSOR;k++) fprintf(g_trig_f,",in%d",k);
            fprintf(g_trig_f,",turn,thrust,eat,attack");
            for(int k=0;k<N_REC;k++) fprintf(g_trig_f,",rec%d",k);
            fprintf(g_trig_f,"\n");
        }
    }
}
static void field_close(){
    if(g_track_f){ fclose(g_track_f); g_track_f=nullptr; }
    if(g_trig_f){ fclose(g_trig_f); g_trig_f=nullptr; }
}

