// 共通の状態: 帳簿・乱数・変異・植物と動物と微生物と死骸の構造体・環境の出来事の記録と再生・測定の時間(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ==========================================================

float SCREEN_W = 1280.f, SCREEN_H = 720.f;

std::vector<float> elev(GRID_W*GRID_H, 0.f);
std::vector<char>  sea (GRID_W*GRID_H, 0);
std::vector<float> solar(GRID_W*GRID_H, 0.f);
std::vector<float> temper(GRID_W*GRID_H, 0.f);
std::vector<float> co2(GRID_W*GRID_H, CO2_BASE);
std::vector<float> o2 (GRID_W*GRID_H, O2_BASE);
// 気体の場は float で、1 マスに約 8 ある。その刻み(約 1e-6)より小さい出入り(死骸1つの1分ぶんの分解など)を
// そのまま足すと丸めで消え、質量が保存されない。足しきれなかった端数をマスごとに持ち越し、次に足す
std::vector<float> co2_res(GRID_W*GRID_H, 0.f);
std::vector<float> o2_res (GRID_W*GRID_H, 0.f);
inline void gas_add(std::vector<float>& g,std::vector<float>& r,int c,float d){
    float w=r[c]+d, old=g[c], nv=old+w;
    r[c]=w-(nv-old); g[c]=nv;
}
inline void co2_add(int c,float d){ gas_add(co2,co2_res,c,d); }
inline void o2_add (int c,float d){ gas_add(o2, o2_res, c,d); }
std::vector<float> gbuf(GRID_W*GRID_H, 0.f);
std::vector<float> wind(GRID_H, 0.f);
std::vector<float> slope(GRID_W*GRID_H, 0.f);
std::vector<float> precip(GRID_W*GRID_H, 0.f);
std::vector<float> soil_water(GRID_W*GRID_H, 0.f);
std::vector<float> soil_n(GRID_W*GRID_H, SOIL_N_BASE);
std::vector<float> gw_h(GRID_W*GRID_H, 0.f);   // 陸の帯水層の飽和した厚さ(m、0〜GW_DEPTH_M)。海のマスは使わない
// 塩の在庫(陸のマスだけ。海は塩分 SEA_SALT で一定の外の溜め)。単位は「g/kg × 土の水の単位」で、水の量で割ると濃さ g/kg になる
std::vector<float> salt_soil(GRID_W*GRID_H, 0.f);   // 根の層の土の水に溶けた塩(と、飽和を超えて結晶になった塩)
std::vector<float> salt_gw(GRID_W*GRID_H, 0.f);     // 帯水層の水に溶けた塩
// 収支の帳簿。系の外とのやりとり(明示した出入り)を種類ごとに積算する。
// 読むだけで物理は変えない。BUDGET_EVERY ごとに在庫の変化と並べて出力し、0 に戻す
const int BUDGET_EVERY = 600;
struct Budget {
    double c_relax=0, c_trans=0;     // 大気CO2: 基準値への緩和(海・岩石との交換とみなす) / 輸送の数値誤差と下限での切り捨て
    double o_relax=0, o_trans=0;     // 大気O2: 同上
    double c_edge=0, o_edge=0;       // 東西の端での区画の外の空気との交換(風が運ぶ出入り)
    double n_sink=0, n_runoff=0, n_upwell=0, n_clamp=0;   // 海の窒素: 沈降・河川流入・湧昇・上下限の切り捨て
    double n_fix=0;                  // 大気の N2 からの窒素固定(大気の N2 は帳簿の外)
    double n_fire=0;                 // 燃えた組織の窒素は揮発して大気へ出る(N2・NOx。大気は帳簿の外)
    double n_carn=0;                 // 食虫植物が獲物から得る窒素(動物は窒素の在庫を持たないので、帳簿の外から入る)
    double w_relax=0;                // 陸の土壌水: 降水・河川・流出による目標への緩和(上下限込み)
    double w_transp=0;               // 蒸散: 陸の土の水と植物の貯水から大気へ出る水(大気の水蒸気は帳簿の外)
    double w_sea=0;                  // 海のマスの植物と海との水のやりとり(海の水は帳簿の外)
    double w_storm=0;                // 雨雲から陸の土と地下水に染み込んだ雨(雲の水は大気の水蒸気から来る。大気は帳簿の外)
    double w_gw_sea=0;               // 地下水と海とのやりとり(陸へ入る向きが正。海の水は帳簿の外)
    double w_gw_edge=0;              // 区画の端を通る地下水の出入り(端の外は同じ傾きが続くとする条件)
    double w_gw_out=0;               // 地下水面が地表に出て湧き、土に収まらず地表を流れ去る分(負)
    double w_gw_clamp=0;             // 飽和した厚さが負になるのを切り捨てた分(陽解法の破れ。0 であるべき)
    double s_sea=0, s_edge=0;        // 塩: 地下水と一緒に海・区画の端と出入りする塩(陸へ入る向きが正。海の塩は帳簿の外)
    double s_out=0, s_clamp=0;       // 塩: 湧いて地表を流れ去る水に溶けた塩(負) / 水の切り捨てと一緒に消した塩(0 であるべき)
    double c_immig=0, n_immig=0;     // 土地の外との生き物の出入り(動物の放流・微生物の胞子・浮遊植物が持ち込む分、流されて端から出る浮遊植物の分は負)
} g_bud;
// 光合成ループは並列なので、蒸散の帳簿はマスごとに積み、帳簿を出すときに順に足す(足す順を固定して決定性を保つ)
std::vector<double> w_transp_cell(GRID_W*GRID_H, 0.0);
std::vector<float> shelf(GRID_W*GRID_H, 0.f);   // 陸への近さ(1=岸 0=外洋)
std::vector<char>  upwell(GRID_W*GRID_H, 0);
std::vector<char>  ice(GRID_W*GRID_H, 0);
std::vector<float> coast_w(GRID_W*GRID_H, 0.f);   // 陸からの近さ(1=岸 0=外洋)
std::vector<int>   burn_age(GRID_W*GRID_H, 999999);   // 最後に焼けてからの経過
float day_light = 1.f;
// 表示の状態と、地図ごとのグラフの履歴
bool show_night = true;      // 昼夜の陰影を出すか(Nキー)
int  g_last_view = -1;       // 背景を描き直させるための印
std::vector<float> h_soil, h_nut, h_ice, h_path, h_odor, h_sst;
std::vector<float> daylit(GRID_W*GRID_H, 1.f);   // マスごとの昼夜(0=夜 1=昼)
float sun_dec = 0.f;                   // 今の太陽赤緯(ラジアン)
float g_sun   = 0.f;                   // 今の太陽高度の正弦(地平線の下は 0)。区画の中で一様
float g_t_bg  = 0.f;                   // 海抜 0 の陸の今の気温。区画の外の背景の気温(条件)としても使う
std::vector<float> odor[N_ODOR], odor_buf;
// 雨雲の突風前線による地上 10 m の風(m/s、東向きと南向きの成分)。卓越風に重ねる。匂いの更新の tick に雲から作り直す
std::vector<float> gust_u(GRID_W*GRID_H,0.f), gust_v(GRID_W*GRID_H,0.f);
std::vector<char>  gust_row(GRID_H,0);   // 突風のある行。無い行は卓越風だけで運ぶ
bool   gust_on=false;
float  gust_max=0.f;                     // 前の表示からの、陸の突風の最大(m/s)
double odor_total[N_ODOR]={0,0,0,0};
std::vector<float> flow_acc(GRID_W*GRID_H, 0.f);
std::vector<char>  river(GRID_W*GRID_H, 0);
std::vector<float> temp_annual(GRID_W*GRID_H, 0.f);
std::vector<float> aridity(GRID_W*GRID_H, 0.f);
std::vector<unsigned char> climate(GRID_W*GRID_H, 0);
std::vector<float> elev_smooth(GRID_W*GRID_H, 0.f);
std::vector<float> sst_anom(GRID_W*GRID_H, 0.f);  // 海流による水温偏差
std::vector<float> cur_u(GRID_W*GRID_H, 0.f);     // 海流の東西成分
std::vector<float> cur_v(GRID_W*GRID_H, 0.f);     // 海流の南北成分
// 0=海 1=EF氷雪 2=ETツンドラ 3=BW砂漠 4=BSステップ 5=D亜寒帯 6=C温帯 7=Awサバナ 8=Af熱帯雨林
inline int idx(int x,int y){ return y*GRID_W + x; }

std::mt19937 rng(12344);
// 夜間の自動実行用。起動引数で上書きする
long   g_stop_at  = -1;     // この tick に達したら終了する(-1 で無制限)
bool   g_headless = false;  // 画面を描かずに回す。描画は状態を変えないので結果は同じ
double g_fps_time = 0.0;    // t>=50000 の経過時間の合計
long   g_fps_n    = 0;      // t>=50000 のフレーム数
std::uniform_real_distribution<float> dist01(0.f,1.f);
inline float frand(float a,float b){ return a + (b-a)*dist01(rng); }
// 1 回の更新あたり、平均 n0·STEP_RATIO 回。整数部はそのまま、端数は確率で 1 回足す。
// 端数が 0 のときは乱数を引かない(段階 3.5 では DT = DT0 なので常に引かない)
inline int count_by_rate(double n0){
    double e=n0*STEP_RATIO;
    int n=(int)e; double f=e-n;
    if(f>0.0 && dist01(rng)<f) n++;
    return n;
}
inline float clampf(float v,float lo,float hi){ return v<lo?lo:(v>hi?hi:v); }
inline float clamp01(float v){ return v<0.f?0.f:(v>1.f?1.f:v); }
inline float mutate(float v,float lo,float hi){
    float d=frand(-MUT_STRENGTH,MUT_STRENGTH)*(hi-lo);
    v+=d; if(v<lo)v=lo; if(v>hi)v=hi; return v;
}
// 変異率を個体ごとに変えられる版。修復機構の精度が遺伝子になる
inline float mutate_r(float v,float lo,float hi,float rate){
    float d=frand(-MUT_STRENGTH,MUT_STRENGTH)*(hi-lo)*rate;
    v+=d; if(v<lo)v=lo; if(v>hi)v=hi; return v;
}

struct Plant {
    float x,y,energy,height;
    bool  alive;
    float absorb,max_height,shade_tol,tough,clonal,disperse,cold_tol,aquatic,drought_tol;
    float salt_tol;         // 遺伝子: 塩への強さ(適合溶質で浸透圧を合わせられる塩分。海水を 1 とする比)。aquatic は水に浸かることへの強さ
    float dorm_temp;        // 遺伝子: この気温を下回ると休眠に入る
    float store_cap;        // 遺伝子: 貯水できる量(多肉化)
    float fire_tol;         // 遺伝子: 厚い樹皮。炎に耐える
    float climb;            // 遺伝子: 他の植物に登る。幹を作らず高さを得る
    float sexual;           // 遺伝子: 有性生殖に頼る度合い
    float brood;            // 遺伝子: 一度に作る種子の数
    float n_fix;            // 遺伝子: 共生菌による窒素固定
    float n_store;          // 状態: 体内の窒素
    float lifespan;         // 遺伝子: 寿命の長さ
    float parasite;         // 遺伝子: 吸器で隣の植物から奪う
    float carnivory;        // 遺伝子: 虫を捕らえる。窒素とエネルギーを得る
    float fire_seed;        // 遺伝子: 焼け跡でだけ芽吹く休眠種子
    float repro_alloc;      // 遺伝子: 余剰を繁殖に回す割合
    float seed_bank;        // 遺伝子: 種子を土中で待たせる傾向
    float buoyancy;
    float mut_rate;
    int   mother, father;   // 本当の親。無性なら father=-1
    int   gen;              // 世代数
    float im[N_IMMUNE];
    float water_store;      // 状態: 実際に貯めている水
    float eff_height;       // 状態: つるで得た分を含む実効的な高さ
    float defense;          // 状態: 今張っている防御(0〜tough)
    int   seed_wait;        // 状態: 土中で発芽を待つ残り時間
    bool  dormant;          // 状態: 休眠しているか
    int   id,genet,age,lin;
};
std::vector<Plant> plants;
// マスごとの索引。入れ子のvectorをやめ、連続した配列で持つ
struct Grid {
    std::vector<int> start, items, cur;
    Grid(): start(GRID_W*GRID_H+1,0), cur(GRID_W*GRID_H,0) {}
    inline void reset(){ std::fill(start.begin(),start.end(),0); }
    inline void count(int ci){ start[ci+1]++; }
    void finalize(){
        for(int i=0;i<GRID_W*GRID_H;i++) start[i+1]+=start[i];
        items.resize(start[GRID_W*GRID_H]);
        for(int i=0;i<GRID_W*GRID_H;i++) cur[i]=start[i];
    }
    inline void put(int ci,int v){ items[cur[ci]++]=v; }
    inline int  b(int ci) const { return start[ci]; }
    inline int  e(int ci) const { return start[ci+1]; }
    inline bool empty_at(int ci) const { return start[ci]==start[ci+1]; }
    inline int  n_at(int ci) const { return start[ci+1]-start[ci]; }
};
Grid pgrid, agrid, cgrid;
std::vector<std::vector<int>> plant_grid(GRID_W*GRID_H);
int next_plant_id = 0;

struct Animal {
    float x,y,dir,energy;
    bool  alive;
    float speed,turn,reach,size,diet,jaw,cold_tol,body,aquatic,resp,detritus,sexual;
    float salt_tol;         // 遺伝子: 塩への強さ(浸透圧を調節できる塩分。海水を 1 とする比)。aquatic は水中での呼吸と泳ぎ
    float arboreal;   // 遺伝子: 登攀に適した体
    float brood;      // 遺伝子: 一度に産む数。多産少投資か、少産多投資か
    float lifespan;   // 遺伝子: 寿命の長さ
    float toxin;      // 遺伝子: 毒を持つ
    float tox_resist; // 遺伝子: 毒に耐える
    float burrow;     // 遺伝子: 穴を掘る
    float nocturnal;  // 遺伝子: 夜に活動する
    float store_fat;  // 遺伝子: 脂肪の目標(体の維持費の store_fat × FAT_DAYS_MAX 日分)
    float repro_alloc;// 遺伝子: 余剰を繁殖に回す割合
    float fat;        // 状態: 蓄えた脂肪
    float depth;      // 状態: 今の潜り具合(0=地表)
    float perch;      // 状態: 今いる高さ(0=地上)
    float im[N_IMMUNE];
    float aerobic;    // 状態: 今得られている酸素の充足度(表示用)
    // 感覚器。角度・射程・何に反応するかが、すべて遺伝子
    struct Sensor { float angle, range, tune; };  // tune: 0=植物 .5=動物 1=死骸
    Sensor sen[N_SENSOR];
    float w[N_W];
    float ngain[N_NODE];
    float mut_rate;   // 遺伝子: 自分の変異率。環境が荒れると上がりうる
    float bite_size;  // 遺伝子: 一口の大きさ。大きいと速いが、咀嚼が雑になる
    float digest;     // 遺伝子: 消化効率。腸が長いほど高いが、体が重くなる
    float rec[N_REC];
    // 本当の系譜。形質による分類(lin)とは別に、実際の親子関係を持つ
    int   mother, father;   // 親の id。無性なら father=-1
    int   gen;              // 世代数
    int   id,age,eaten,lin;
    bool  lesion;           // 対照(--no-recurrent): 再帰の入力を 0 に切る。遺伝子ではなく、子にそのまま写る
};
std::vector<Animal> animals;
// 分解者。植物でも動物でもない第三の界。
// 動かず、死骸だけを食い、炭素を大気へ、窒素を土へ返す
struct Microbe {
    float x,y,energy;
    float rate;      // 遺伝子: 分解の速さ。速いと独占できないが、取りこぼさない
    float cold_tol;  // 遺伝子: 低温でも働けるか
    bool  alive;
    int   id,age;
};
std::vector<Microbe> microbes;
int next_microbe_id=0;
long d_micro=0;
long g_mic_upd=0, g_mic_cap=0;   // 基質のあるマスでの微生物の更新の数と、そのうち取り込みが上限で決まった数
long g_mic_col=0;                // 胞子が死骸のマスに定着した数(表示だけ)
int  g_frame=0;   // 現在の tick。関数の引数に足さずに参照できるようにする

// 選択した個体の一生を記録する。進化の観察とデバッグの両方に効く
struct LifeEvent { int tick; char kind[24]; float value; };
std::vector<LifeEvent> life_log;
int  life_target=-1;          // 追跡中の動物の id
void life_note(int id,int tick,const char* kind,float v=0.f){
    if(id!=life_target) return;
    if(life_log.size()>200) return;
    LifeEvent e; e.tick=tick; e.value=v;
    snprintf(e.kind,24,"%s",kind);
    life_log.push_back(e);
}
// 選択勾配を測るための遺伝子ベクトル。
// 繁殖した親と集団全体で平均を比べれば、どちらへ押されているかが分かる
inline void sg_vec_a(const Animal& a,double* v){
    v[0]=a.diet;  v[1]=a.detritus; v[2]=a.reach;    v[3]=a.size;
    v[4]=a.speed; v[5]=a.jaw;      v[6]=a.cold_tol; v[7]=a.aquatic;
    v[8]=a.resp;  v[9]=a.arboreal; v[10]=a.sexual;  v[11]=a.brood;
    v[12]=a.salt_tol;
}
inline void sg_vec_p(const Plant& p,double* v){
    v[0]=p.max_height; v[1]=p.shade_tol; v[2]=p.tough;      v[3]=p.clonal;
    v[4]=p.disperse;   v[5]=p.cold_tol;  v[6]=p.aquatic;    v[7]=p.drought_tol;
    v[8]=p.climb;      v[9]=p.parasite;  v[10]=p.carnivory; v[11]=p.brood;
    v[12]=p.salt_tol;
}
std::vector<std::vector<int>> animal_grid(GRID_W*GRID_H);
int next_animal_id = 0;

struct Corpse { float x,y,energy; int origin; float n=0.f; };   // 0=植物由来 1=動物由来。n は体に含まれる窒素
std::vector<Corpse> corpses;
std::vector<std::vector<int>> corpse_grid(GRID_W*GRID_H);
std::vector<float> corpse_density(GRID_W*GRID_H,0.f);
std::vector<float> fuel(GRID_W*GRID_H,0.f);
std::vector<float> fuel_moist(GRID_W*GRID_H,-1.f);  // 燃料の含水率。負は未設定(最初の更新で平衡の値を置く)
std::vector<float> fire(GRID_W*GRID_H,0.f);
std::vector<float> fire_buf(GRID_W*GRID_H,0.f);
std::vector<float> canopy_h(GRID_W*GRID_H,0.f);   // そのマスで最も高い植物の高さ
// 病原体: 型(免疫座位ごと)と密度。植物と動物で別の病がつく
std::vector<float> pstrain[N_IMMUNE], astrain[N_IMMUNE];
std::vector<float> pload(GRID_W*GRID_H,0.f), aload(GRID_W*GRID_H,0.f);
double infect_p=0, infect_a=0;
long fire_events=0; double fire_area=0;
long slide_events=0;

// 環境の出来事の記録と再生。生き物と無関係に起きる離散的な出来事を
// 1件ずつ tick・位置・種類・強さで書き出し、時刻と場所を混ぜて再生できるようにする。
// 書式: 1行1件 "tick 種類 x y 強さ"。'#' で始まる行は注記
enum EvKind { EV_FIRE=0, EV_SLIDE=1, EV_STORM=2, EV_LIGHTNING=3, EV_QUAKE=4, EV_N };
const char* EV_NAME[EV_N]={"fire","slide","storm","lightning","quake"};
struct EnvEvent { long tick; int kind; int x,y; float strength; };
FILE* g_rec=nullptr;                       // --record の出力先
bool  g_replay=false;                      // --replay 中は該当種類を物理から起こさない
std::vector<EnvEvent> g_replay_ev;         // 混ぜたあと tick 順に並べた出来事
size_t g_replay_pos=0;
long   g_replay_done[EV_N]={0,0,0,0,0};

long   g_cs_ev_n[EV_N]={0,0,0,0,0};        // 全体の記録(--census)用: 前の記録からの件数と強さの和
double g_cs_ev_s[EV_N]={0,0,0,0,0};
void ev_record(int kind,int x,int y,float s){
    g_cs_ev_n[kind]++; g_cs_ev_s[kind]+=s;
    if(!g_rec) return;
    fprintf(g_rec,"%d %s %d %d %.4f\n",g_frame,EV_NAME[kind],x,y,s);
}
// 時刻と場所を、それぞれ同じ種類の出来事の中で入れ替える。
// 種類と強さは出来事に付いたまま残す。種類をまたいで混ぜると、
// 海に崖崩れが起きるような、記録に無い組み合わせになるため
bool ev_load_replay(const char* path,unsigned shuffle_seed){
    std::ifstream f(path);
    if(!f){ fprintf(stderr,"replay: cannot open %s\n",path); return false; }
    std::vector<EnvEvent> ev[EV_N];
    std::string line;
    while(std::getline(f,line)){
        if(line.empty()||line[0]=='#') continue;
        char kn[32]; EnvEvent e;
        if(sscanf(line.c_str(),"%ld %31s %d %d %f",&e.tick,kn,&e.x,&e.y,&e.strength)!=5) continue;
        e.kind=-1;
        for(int k=0;k<EV_N;k++) if(!strcmp(kn,EV_NAME[k])) e.kind=k;
        if(e.kind<0) continue;
        ev[e.kind].push_back(e);
    }
    std::mt19937 srng(shuffle_seed);       // 本体の乱数とは別にする
    g_replay_ev.clear();
    for(int k=0;k<EV_N;k++){
        auto& v=ev[k];
        std::vector<long> ticks; std::vector<std::pair<int,int>> pos;
        for(auto& e:v){ ticks.push_back(e.tick); pos.push_back({e.x,e.y}); }
        std::shuffle(ticks.begin(),ticks.end(),srng);
        std::shuffle(pos.begin(),pos.end(),srng);
        for(size_t i=0;i<v.size();i++){
            v[i].tick=ticks[i]; v[i].x=pos[i].first; v[i].y=pos[i].second;
            g_replay_ev.push_back(v[i]);
        }
        printf("replay: %zu %s events\n",v.size(),EV_NAME[k]);
    }
    std::stable_sort(g_replay_ev.begin(),g_replay_ev.end(),
        [](const EnvEvent& a,const EnvEvent& b){ return a.tick<b.tick; });
    g_replay_pos=0;
    return true;
}
// この tick に起こす再生の出来事を、種類ごとに呼び出し側へ渡す。
// 雨雲は区画の外で生じて流れ込むので、inside=false で区画の外の位置も渡す
template<class F> void ev_replay_due(int kind,F apply,bool inside=true){
    for(size_t i=g_replay_pos;i<g_replay_ev.size() && g_replay_ev[i].tick<=g_frame;i++){
        const EnvEvent& e=g_replay_ev[i];
        if(e.kind!=kind || e.tick<g_frame) continue;
        if(inside && (e.x<0||e.x>=GRID_W||e.y<0||e.y>=GRID_H)) continue;
        apply(e); g_replay_done[kind]++;
    }
}
void ev_replay_advance(){
    while(g_replay_pos<g_replay_ev.size() && g_replay_ev[g_replay_pos].tick<=g_frame) g_replay_pos++;
}

// 選択勾配。繁殖した個体と集団全体で遺伝子の平均を比べる。
// 差がそのまま「その形質が今どちらへ押されているか」になる。
// 収支を手計算せずに、効いているかどうかを測定で言えるようにする
const int N_SG_A = 13;
const char* SG_A_NAME[N_SG_A]={
    "diet","detr","reach","size","speed","jaw",
    "cold","aqua","resp","arbor","sexual","brood","salt"};
const int N_SG_P = 13;
const char* SG_P_NAME[N_SG_P]={
    "maxH","shade","tough","clonal","disp","cold",
    "aqua","dry","climb","paras","carni","brood","salt"};
double sg_a_pop[N_SG_A]={0}, sg_a_rep[N_SG_A]={0};
double sg_p_pop[N_SG_P]={0}, sg_p_rep[N_SG_P]={0};
long   sg_a_np=0, sg_a_nr=0, sg_p_np=0, sg_p_nr=0;
// 全体の記録(--census)用の、繁殖した親の遺伝子の和。ログの表示とは別の間隔で区切るので別に持つ
double cs_a_rep[N_SG_A]={0}, cs_p_rep[N_SG_P]={0};
long   cs_a_nr=0, cs_p_nr=0;
std::vector<float> h_maxH,h_reach,h_diet,h_shade,h_tough,h_jaw,h_coldP,h_coldA,h_co2,h_o2,h_temp;
std::vector<int>   h_plant,h_animal,h_corpse;

// 時間の尺を4段持つ。粗い側は平均を取って別に貯める
const int   N_SCALE = 4;
const int   SCALE_EVERY[N_SCALE] = {1,10,100,1000};   // 何回に一度記録するか
const char* SCALE_NAME[N_SCALE]  = {"1千","1万","10万","全期間"};
const int   SCALE_KEEP = 600;                          // 各段が覚える点の数
struct Series {
    std::vector<float> v[N_SCALE];
    double acc[N_SCALE]={0,0,0,0};
    int    cnt[N_SCALE]={0,0,0,0};
    void push(float x){
        for(int s=0;s<N_SCALE;s++){
            acc[s]+=x; cnt[s]++;
            if(cnt[s]>=SCALE_EVERY[s]){
                v[s].push_back((float)(acc[s]/cnt[s]));
                acc[s]=0; cnt[s]=0;
                if((int)v[s].size()>SCALE_KEEP) v[s].erase(v[s].begin());
            }
        }
    }
};
// 拡大表示用の系列
Series s_corpse, s_defense, s_tough, s_jaw, s_micro,
       s_plant, s_animal, s_maxh, s_reach, s_temp, s_co2, s_o2, s_soil, s_nut,
       s_shade, s_diet, s_coldP, s_coldA;
int peak_plant=1, peak_animal=1, peak_corpse=1;

// 各処理の所要時間(移動平均, ミリ秒)
const int PROF_N = 9;
const char* PROF_NAME[PROF_N] = {
    "climate","precip","water","animals","plants","gases","lineage","bg","draw"
};
double prof[PROF_N]={0};
std::chrono::steady_clock::time_point prof_t0;
inline void prof_begin(){ prof_t0=std::chrono::steady_clock::now(); }
inline void prof_end(int k){
    auto t1=std::chrono::steady_clock::now();
    double ms=std::chrono::duration<double,std::milli>(t1-prof_t0).count();
    prof[k]=prof[k]*0.92+ms*0.08;     // 移動平均で安定させる
}
long d_starve=0, d_cold=0, d_preyed=0, d_slide=0;
// 植物の収支を測る。推測ではなく実測で判断するため
double dbg_gain=0, dbg_cost=0; long dbg_n=0;
double dbg_light=0, dbg_tf=0, dbg_wf=0, dbg_co2=0;

// 最後に生きていた動物の遺伝子を保存しておく(再投入用)
struct AnimalGene {
    float speed,turn,reach,size,diet,jaw,cold_tol,aquatic,resp,detritus,sexual,arboreal;
    float salt_tol;
    float brood,lifespan;
    float toxin,tox_resist,burrow,nocturnal,store_fat,repro_alloc,bite_size,digest;
    float im[N_IMMUNE];
    Animal::Sensor sen[N_SENSOR];
    float w[N_W];
    float ngain[N_NODE];
    int   lin;
};
std::vector<AnimalGene> last_animals;

void remember_animals(){
    if(animals.empty()) return;
    // 種ごとに1体ずつ、個体数の多い順に最大5種
    std::unordered_map<int,int> best;      // lin -> index
    for(int i=0;i<(int)animals.size();i++){
        int L=animals[i].lin;
        auto it=best.find(L);
        if(it==best.end() || animals[i].body>animals[it->second].body) best[L]=i;
    }
    std::vector<int> pick;
    for(auto& kv:best) pick.push_back(kv.second);
    if(pick.size()>5) pick.resize(5);
    last_animals.clear();
    for(int i:pick){
        const Animal& a=animals[i];
        AnimalGene g;
        g.speed=a.speed; g.turn=a.turn; g.reach=a.reach; g.size=a.size;
        g.diet=a.diet;   g.jaw=a.jaw;   g.cold_tol=a.cold_tol;
        g.aquatic=a.aquatic; g.resp=a.resp; g.detritus=a.detritus;
        g.salt_tol=a.salt_tol;
        g.sexual=a.sexual;   g.arboreal=a.arboreal; g.lin=a.lin;
        g.brood=a.brood;     g.lifespan=a.lifespan;
        g.toxin=a.toxin;     g.tox_resist=a.tox_resist; g.burrow=a.burrow;
        g.nocturnal=a.nocturnal; g.store_fat=a.store_fat; g.repro_alloc=a.repro_alloc;
        g.bite_size=a.bite_size; g.digest=a.digest;
        for(int k=0;k<N_IMMUNE;k++) g.im[k]=a.im[k];
        for(int k=0;k<N_W;k++) g.w[k]=a.w[k];
        for(int k=0;k<N_NODE;k++) g.ngain[k]=a.ngain[k];
        for(int k=0;k<N_SENSOR;k++) g.sen[k]=a.sen[k];
        last_animals.push_back(g);
    }
}

