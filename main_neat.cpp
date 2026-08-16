#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <filesystem>
#include <functional>
#ifdef _OPENMP
#include <omp.h>
#endif
// 計測用
double t_prey_nn=0,t_pred_nn=0,t_move=0,t_eat=0,t_field=0,t_draw=0,t_grid=0;
int prof_frames=0;
auto now=[](){ return std::chrono::high_resolution_clock::now(); };
auto ms=[](auto a,auto b){ return std::chrono::duration<double,std::milli>(b-a).count(); };
// 画面と世界を分離
float SCREEN_W=1920.f, SCREEN_H=1020.f;
const float WIDTH=6000.f, HEIGHT=6000.f;   // 世界
const int NUM_PREY=200, NUM_PREDATOR=50, NUM_PLANT=3000;
const float SPEED=2.0f;
const int CELL_SIZE=200;
const int GRID_COLS=(int)WIDTH/CELL_SIZE, GRID_ROWS=(int)HEIGHT/CELL_SIZE;
const float PREY_VISION=150.f, PREDATOR_VISION=300.f;
const int NUM_RAYS=12;
const float TWO_PI=6.28318530718f, PI=3.14159265f;
const float PREDATOR_FOV=60.f*PI/180.f;
const float EAT_DISTANCE=12.f, EAT_GAIN=30.f;
const float METABOLISM=0.05f, INIT_ENERGY=100.f, MAX_ENERGY=120.f;
const float REPRO_ENERGY=100.f;
const int PREY_REPRO_TIME=120, PREDATOR_REPRO_TIME=300;
const float MUT_RATE=0.1f;
const int HALL_SIZE=20;
const int PLANT_REPRO_TIME=300, PLANT_CHILDREN=6, PLANT_MAX=50000;
const int PLANT_MIN=150;   // これ以下になったらランダム補充(絶滅防止)
const float PLANT_SPREAD=100.f;
const float MEAT_ENERGY_STARVE=12.f;
const float MEAT_EAT_RATE=1.0f, MEAT_ROT_RATE=0.05f, MEAT_RADIUS_SCALE=0.2f;
const int NUT_CELL=20;
const int NUT_COLS=(int)WIDTH/NUT_CELL, NUT_ROWS=(int)HEIGHT/NUT_CELL;
const float NUT_DIFFUSION=0.1f, NUT_DECAY=0.999f;
const float NUT_FULL=5.0f, NUT_CONSUME=20.f, PLANT_BASE_SURVIVAL=0.10f;
const float BASE_SIZE=6.f, SIZE_VARIATION=0.5f;
const float MAXENERGY_COEF=20.f, METABOLISM_COEF=0.003f, MOVE_COEF=0.003f;
const float MEAT_SIZE_COEF=5.f, REPRO_COST_COEF=8.f, MUT_RATE_GENE=0.05f;
const float BRAIN_COST_COEF=0.0001f;  // 脳(接続)1本あたりの代謝コスト
const float FEAR_RISE=0.02f, FEAR_DECAY=0.01f, AFFINITY_RISE=0.01f, AFFINITY_DECAY=0.01f;
const int PH_CELL=20;
const int PH_COLS=(int)WIDTH/PH_CELL, PH_ROWS=(int)HEIGHT/PH_CELL;
const float PH_DIFFUSION=0.1f, PH_DECAY=0.99f;
const float STRESS_RISE=0.02f, STRESS_FALL=0.01f, STRESS_MIX=0.005f, STRESS_DECAY=0.005f;
const float PHERO_THRESHOLD=0.05f, STRESS_REPRO_PENALTY=0.5f;
const float STRESS_PREDATOR_SIGHT=0.02f, MUT_STRESS_FACTOR=3.0f;
const int HIST_MAX=6000;

const int PREY_IN=25, PRED_IN=13, N_H1=32, N_H2=16, N_OUT=2;

std::mt19937 rng(12345);
std::uniform_real_distribution<float> dist01(0.f,1.f);
std::normal_distribution<float> gauss(0.f,1.f);
float frand(float lo,float hi){return lo+dist01(rng)*(hi-lo);}
float clamp01(float v){return v<0?0:(v>1?1:v);}
int next_id=0;
int g_innov_counter=1000;
int g_next_node_id=100;

// カメラ
float cam_x=WIDTH/2, cam_y=HEIGHT/2, zoom=SCREEN_H/HEIGHT;
float w2s_x(float wx){return (wx-cam_x)*zoom+SCREEN_W/2;}
float w2s_y(float wy){return (wy-cam_y)*zoom+SCREEN_H/2;}
float s2w_x(float sx){return (sx-SCREEN_W/2)/zoom+cam_x;}
float s2w_y(float sy){return (sy-SCREEN_H/2)/zoom+cam_y;}
bool on_screen(float wx,float wy,float margin=30.f){
    float sx=w2s_x(wx),sy=w2s_y(wy);
    return sx>=-margin&&sx<=SCREEN_W+margin&&sy>=-margin&&sy<=SCREEN_H+margin;
}

// ===== NEAT: センサー・アクチュエータ・脳の3本立て =====
struct Sensor {
    int target;    // 0=植物, 1=敵, 2=同種, 3=肉, 4=内部(energy等)
    float angle;   // facing相対の角度(内部センサーは0)
    float range;   // 視線の長さ(内部センサーは0)
};
struct Actuator {
    int type;      // 0=vx, 1=vy  (将来: 2=フェロモン放出 等)
};
struct NodeGene {
    int id;        // ノード固有番号
    int type;      // 0=入力, 1=出力, 2=中間(hidden)
};
struct ConnGene {
    int in_node;   // 始点ノードid
    int out_node;  // 終点ノードid
    float weight;
    bool enabled;
    int innov;     // イノベーション番号(段階5の交叉で使う。今は通し番号)
};
struct Genome {
    std::vector<Sensor> sensors;      // 感覚器官
    std::vector<Actuator> actuators;  // 出力器官
    std::vector<NodeGene> nodes;      // 脳のノード
    std::vector<ConnGene> conns;      // 脳の接続
    std::vector<int> eval_order;                              // 計算順序(非入力ノード)
    std::vector<std::vector<std::pair<int,float>>> incoming;  // 各ノードへの入力接続
    int cached_n_in=0;
    int cached_active_conns=0;   // 有効接続数(脳コスト計算用)
    // 計算順序と入力接続リストを前計算(構造が確定した時に1回呼ぶ)
    void finalize(){
        int n_in=(int)sensors.size();
        cached_n_in=n_in;
        cached_active_conns=0;
        for(const auto& c:conns) if(c.enabled) cached_active_conns++;
        int max_id=0;
        for(const auto& n:nodes) if(n.id>max_id) max_id=n.id;
        for(const auto& c:conns){ if(c.in_node>max_id)max_id=c.in_node; if(c.out_node>max_id)max_id=c.out_node; }  // ← 接続の参照先も見る
        incoming.assign(max_id+1, {});
        for(const auto& c:conns) if(c.enabled) incoming[c.out_node].push_back({c.in_node, c.weight});
        std::vector<int> deg(max_id+1,0);
        for(const auto& c:conns) if(c.enabled) deg[c.out_node]++;
        std::vector<int> q; for(int i=0;i<n_in;i++) q.push_back(i);
        size_t qi=0; eval_order.clear();
        while(qi<q.size()){
            int nid=q[qi++];
            for(const auto& c:conns){ if(!c.enabled)continue; if(c.in_node==nid){ deg[c.out_node]--; if(deg[c.out_node]==0){ q.push_back(c.out_node); eval_order.push_back(c.out_node);} } }
        }
    }
    // 順伝播: 計算順序対応版
    void forward(const std::vector<float>& sv, std::vector<float>& output) const {
        if(eval_order.empty() && incoming.empty()){
            // finalizeされてない → 出力ゼロで返す(クラッシュ回避)
            output.assign((int)actuators.size(), 0.f);
            return;
        }
        int n_in=cached_n_in;
        int n_out=(int)actuators.size();
        int max_id=(int)incoming.size()-1;
        static thread_local std::vector<float> node_val;
        node_val.assign(max_id+1, 0.f);
        for(int i=0;i<n_in;i++) node_val[i]=sv[i];
        for(int nid : eval_order){
            float s=0.f;
            for(const auto& pr : incoming[nid]) s += node_val[pr.first]*pr.second;
            node_val[nid]=std::tanh(s);
        }
        output.assign(n_out, 0.f);
        for(int j=0;j<n_out;j++) output[j]=node_val[n_in+j];
    }
    // 接続追加の変異(段階4b)。繋がってない2ノードを繋ぐ
    // innov_counter: イノベーション番号を管理する外部カウンタ(参照で受け取る)
    Genome add_connection(int& innov_counter) const {
        Genome c = *this;
        int n_in=(int)sensors.size();
        int n_out=(int)actuators.size();
        // 候補: (from, to) のペアを探す
        // from = 入力 or 中間(出力ノードからは出さない)
        // to   = 中間 or 出力(入力ノードへは入れない)
        std::vector<std::pair<int,int>> candidates;
        for(const auto& a : c.nodes){
            if(a.type==1) continue;  // 出力ノードからは出さない
            for(const auto& b : c.nodes){
                if(b.type==0) continue;  // 入力ノードへは入れない
                if(a.id==b.id) continue; // 自己ループ禁止
                // 既に繋がってるか
                bool exists=false;
                for(const auto& cc : c.conns) if(cc.in_node==a.id && cc.out_node==b.id){ exists=true; break; }
                if(exists) continue;
                // ループにならないか(bからaへ辿れてしまうと逆流=ループ)
                // 簡易チェック: bから出発してaに到達できるなら、a→bはループを作る
                // (今はフィードフォワード維持のため)
                // BFSでbの下流にaがいないか確認
                std::vector<int> stack={b.id}; bool makes_loop=false;
                std::vector<bool> visited(10000,false);
                while(!stack.empty()){
                    int cur=stack.back(); stack.pop_back();
                    if(cur==a.id){ makes_loop=true; break; }
                    if(cur<10000){ if(visited[cur])continue; visited[cur]=true; }
                    for(const auto& cc:c.conns) if(cc.enabled && cc.in_node==cur) stack.push_back(cc.out_node);
                }
                if(makes_loop) continue;
                candidates.push_back({a.id,b.id});
            }
        }
        if(candidates.empty()) return c;  // 追加できる場所がない
        // ランダムに1つ選んで接続を追加
        auto pick = candidates[(int)(dist01(rng)*candidates.size())%candidates.size()];
        float s=1.f/std::sqrt((float)n_in);
        c.conns.push_back({pick.first, pick.second, gauss(rng)*s, true, innov_counter++});
        return c;
    }
    // ノード追加の変異(段階4c)。既存の接続の途中に中間ノードを挿入
    Genome add_node(int& innov_counter, int& next_node_id) const {
        Genome c = *this;
        // 有効な接続を集める
        std::vector<int> enabled_idx;
        for(int i=0;i<(int)c.conns.size();i++) if(c.conns[i].enabled) enabled_idx.push_back(i);
        if(enabled_idx.empty()) return c;
        // ランダムに1本選ぶ
        int ci = enabled_idx[(int)(dist01(rng)*enabled_idx.size())%enabled_idx.size()];
        ConnGene& old = c.conns[ci];
        int a = old.in_node, b = old.out_node;
        float w = old.weight;
        // 元の接続を無効化
        c.conns[ci].enabled = false;
        // 新しい中間ノードを作る
        int newid = next_node_id++;
        c.nodes.push_back({newid, 2});   // type=2(中間)
        // A→N(重み1.0)、N→B(元の重みw)
        c.conns.push_back({a, newid, 1.0f, true, innov_counter++});
        c.conns.push_back({newid, b, w,    true, innov_counter++});
        return c;
    }
    Genome mutate_weights(float rate) const {
        Genome c = *this;
        for(auto& conn : c.conns){
            conn.weight += gauss(rng) * rate;
        }
        return c;
    }
};
// 案B: 入力→出力を全結合、中間ノードなし の初期ゲノムを作る
// sensors と actuators を渡すと、それに応じた脳を作る
Genome make_initial_genome(const std::vector<Sensor>& sensors, const std::vector<Actuator>& actuators){
    Genome g;
    g.sensors = sensors;
    g.actuators = actuators;
    int n_in = (int)sensors.size();
    int n_out = (int)actuators.size();
    // ノード: 入力ノード(id 0..n_in-1) + 出力ノード(id n_in..n_in+n_out-1)
    for(int i=0;i<n_in;i++)  g.nodes.push_back({i, 0});          // 入力
    for(int j=0;j<n_out;j++) g.nodes.push_back({n_in+j, 1});     // 出力
    // 接続: 全入力 → 全出力(全結合)
    int innov=0;
    float s=1.f/std::sqrt((float)n_in);
    for(int i=0;i<n_in;i++)for(int j=0;j<n_out;j++){
        g.conns.push_back({i, n_in+j, gauss(rng)*s, true, innov++});
    }
    g.finalize();
    return g;
}
// ゲノムの総合変異(重み + たまに構造)
Genome mutate_genome(const Genome& g){
    Genome c = g.mutate_weights(MUT_RATE);
    if(dist01(rng) < 0.05f) c = c.add_connection(g_innov_counter);
    if(dist01(rng) < 0.03f) c = c.add_node(g_innov_counter, g_next_node_id);
    c.finalize();
    return c;
}
// 交叉: 母を土台に、父から遺伝子(接続)を混ぜる。整合性は保証しない(流産は別途検査)
Genome crossover(const Genome& mother, const Genome& father){
    Genome child = mother;   // 母のノード構造を土台に
    for(auto& c : child.conns){
        // 母の各接続に、父の同じinnovがあれば、50%で父の重みを採用
        for(const auto& fc : father.conns){
            if(fc.innov == c.innov){
                if(dist01(rng) < 0.5f) c.weight = fc.weight;
                break;
            }
        }
    }
    // 父にしかない接続を、50%で取り込む(← ここで不整合が起きうる = 発生異常の種)
    for(const auto& fc : father.conns){
        bool in_mother=false;
        for(const auto& c : mother.conns) if(c.innov==fc.innov){ in_mother=true; break; }
        if(!in_mother && dist01(rng)<0.5f){
            child.conns.push_back(fc);   // 父の接続を追加(母にそのノードが無ければ不整合→流産)
        }
    }
    child.finalize();
    return child;
}
// ゲノムが整合しているか(全接続の参照ノードが存在するか)
bool is_valid_genome(const Genome& g){
    // ノードidの集合
    int max_id=0; for(const auto& n:g.nodes) if(n.id>max_id) max_id=n.id;
    std::vector<bool> exists(max_id+2, false);
    for(const auto& n:g.nodes) exists[n.id]=true;
    for(const auto& c:g.conns){
        if(!c.enabled) continue;
        if(c.in_node>max_id || c.out_node>max_id) return false;   // 範囲外
        if(!exists[c.in_node] || !exists[c.out_node]) return false; // 参照先が無い=発生異常
    }
    return true;
}
// 2つのゲノムの構造的距離(NEATの互換性距離)
float genome_distance(const Genome& a, const Genome& b){
    const float C1=1.0f;   // 構造の違い(片方にしかない接続)の重み
    const float C3=0.4f;   // 重みの差の重み
    int disjoint=0;              // 片方にしかない接続の数
    float weight_diff=0.f;       // 共通接続の重み差の合計
    int matching=0;              // 共通接続の数
    // aの各接続について
    for(const auto& ca : a.conns){
        bool found=false;
        for(const auto& cb : b.conns){
            if(ca.innov==cb.innov){
                weight_diff += std::fabs(ca.weight - cb.weight);
                matching++;
                found=true;
                break;
            }
        }
        if(!found) disjoint++;   // aにしかない
    }
    // bにしかない接続を数える
    for(const auto& cb : b.conns){
        bool found=false;
        for(const auto& ca : a.conns){ if(ca.innov==cb.innov){ found=true; break; } }
        if(!found) disjoint++;   // bにしかない
    }
    // 正規化: 大きいゲノムのサイズで割る(NEATの流儀)
    int N = std::max((int)a.conns.size(), (int)b.conns.size());
    if(N<1) N=1;
    float avg_wdiff = (matching>0) ? weight_diff/matching : 0.f;
    return C1 * disjoint / N + C3 * avg_wdiff;
}

struct Genes { float size, vision, eat_gain, speed; };
Genes random_genes(){return {BASE_SIZE*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION)),PREY_VISION*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION)),EAT_GAIN*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION)),SPEED*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION))};}
float mut1(float v,float w){float m=v*(1.f+frand(-w,w));return m<0.01f?0.01f:m;}
Genes mutate_genes(const Genes& g,float stress){float w=MUT_RATE_GENE*(1.f+stress*MUT_STRESS_FACTOR);return {mut1(g.size,w),mut1(g.vision,w),mut1(g.eat_gain,w),mut1(g.speed,w)};}

struct Agent {
    float x,y,vx,vy,energy,max_energy; bool alive,is_prey; int age,repro_counter;
    float fear,affinity,stress; Genome genome; Genes genes;
    Agent(){}
    int id;
    int repro_count, food_count;
    int last_eat_frame;
    int mate_ready_counter;
    int cluster;
    bool saw_pred, saw_ally;
};
struct Plant { float x,y; int repro_counter; };
struct Meat { float x,y,energy; };

struct Field {
    int cols,rows; float diff,decay; std::vector<float> grid;
    Field(int c,int r,float df,float dc):cols(c),rows(r),diff(df),decay(dc){grid.assign(c*r,0.f);}
    int idx(int cx,int cy)const{return cy*cols+cx;}
    void emit(int cp,float x,float y,float amt){int cx=(int)x/cp,cy=(int)y/cp;cx=(cx%cols+cols)%cols;cy=(cy%rows+rows)%rows;grid[idx(cx,cy)]+=amt;if(grid[idx(cx,cy)]<0)grid[idx(cx,cy)]=0;}
    float get(int cp,float x,float y)const{int cx=(int)x/cp,cy=(int)y/cp;cx=(cx%cols+cols)%cols;cy=(cy%rows+rows)%rows;return grid[idx(cx,cy)];}
    void update(){
        std::vector<float> ng(cols*rows,0.f); float keep=1.f-diff,share=diff/4.f;
        for(int cy=0;cy<rows;cy++)for(int cx=0;cx<cols;cx++){float v=grid[idx(cx,cy)];ng[idx(cx,cy)]+=v*keep;
            int l=(cx-1+cols)%cols,r=(cx+1)%cols,u=(cy-1+rows)%rows,d=(cy+1)%rows;
            ng[idx(l,cy)]+=v*share;ng[idx(r,cy)]+=v*share;ng[idx(cx,u)]+=v*share;ng[idx(cx,d)]+=v*share;}
        for(auto&v:ng)v*=decay; grid=std::move(ng);
    }
};

int cell_index(int cx,int cy){return cy*GRID_COLS+cx;}
float torus_delta(float a,float b,float size){float d=b-a;if(d>size/2)d-=size;else if(d<-size/2)d+=size;return d;}

struct PreyRec { int fit; Genome genome; Genes genes; };
std::vector<PreyRec> prey_hall;
std::vector<std::pair<int,Genome>> pred_hall;
int generation=1;

// preyのセンサー37個を作る(植物12 + 敵12 + 同種12 + energy1)
std::vector<Sensor> make_prey_sensors(){
    std::vector<Sensor> s;
    for(int i=0;i<12;i++) s.push_back({0, TWO_PI*i/12, PREY_VISION});  // 植物12方向
    for(int i=0;i<12;i++) s.push_back({1, TWO_PI*i/12, PREY_VISION});  // 敵12方向
    for(int i=0;i<12;i++) s.push_back({2, TWO_PI*i/12, PREY_VISION});  // 同種12方向
    s.push_back({4, 0.f, 0.f});  // energy(内部センサー)
    return s;
}
// predatorのセンサー37個(prey12 + 肉12 + 同種12 + energy1)
std::vector<Sensor> make_pred_sensors(){
    std::vector<Sensor> s;
    for(int i=0;i<12;i++) s.push_back({1, PREDATOR_FOV*i/11, PREDATOR_VISION});  // prey12(FOV内)
    for(int i=0;i<12;i++) s.push_back({3, PREDATOR_FOV*i/11, PREDATOR_VISION});  // 肉12
    for(int i=0;i<12;i++) s.push_back({2, PREDATOR_FOV*i/11, PREDATOR_VISION});  // 同種12
    s.push_back({4, 0.f, 0.f});  // energy
    return s;
}
std::vector<Actuator> make_actuators(){
    return {{0},{1}};  // vx, vy
}

Agent make_prey(const Genes& g){
    Agent a;a.x=frand(0,WIDTH);a.y=frand(0,HEIGHT);a.vx=frand(-SPEED,SPEED);a.vy=frand(-SPEED,SPEED);
    a.energy=INIT_ENERGY;a.alive=true;a.is_prey=true;a.age=0;a.repro_count=0; a.food_count=0;a.repro_counter=0;a.fear=0;a.affinity=0;a.stress=0;a.id=next_id++;a.last_eat_frame=-1000;a.mate_ready_counter=0;
    a.genes=g;                                          // ← これを追加
    a.max_energy=g.size*MAXENERGY_COEF;                 // ← これも追加(max_energyの設定)
    a.genome=make_initial_genome(make_prey_sensors(), make_actuators());
    return a;
}
Agent make_predator(){
    Agent a;a.x=frand(0,WIDTH);a.y=frand(0,HEIGHT);a.vx=frand(-SPEED,SPEED);a.vy=frand(-SPEED,SPEED);
    a.energy=INIT_ENERGY;a.alive=true;a.is_prey=false;a.age=0;a.repro_count=0; a.food_count=0;a.repro_counter=0;a.fear=0;a.affinity=0;a.stress=0;a.id=next_id++;a.last_eat_frame=-1000;a.mate_ready_counter=0;
    a.genes={BASE_SIZE,PREDATOR_VISION,EAT_GAIN,SPEED};  // ← predatorのgenesを追加
    a.max_energy=MAX_ENERGY;                             // ← これも追加
    a.genome=make_initial_genome(make_pred_sensors(), make_actuators());
    return a;
}
void trim_prey(){std::sort(prey_hall.begin(),prey_hall.end(),[](auto&a,auto&b){return a.fit>b.fit;});if((int)prey_hall.size()>HALL_SIZE)prey_hall.resize(HALL_SIZE);}
void trim_pred(){std::sort(pred_hall.begin(),pred_hall.end(),[](auto&a,auto&b){return a.first>b.first;});if((int)pred_hall.size()>HALL_SIZE)pred_hall.resize(HALL_SIZE);}

int calc_fit(const Agent& a){
    return (int)(a.repro_count*100 + a.food_count*5 + a.age*0.1f);
}

struct Lineage {
    int id; Genome rep; int birth_tick; int death_tick; int parent; bool alive; int pop; int max_pop;
};
std::vector<Lineage> lineages;
int next_lineage_id=0;

// ===== セーブ/ロード(遺伝的セーブ: 脳と遺伝子を残す) =====
void write_genome(std::ofstream& f, const Genome& g){
    f << "G " << g.sensors.size() << " " << g.actuators.size() << " " << g.nodes.size() << " " << g.conns.size() << "\n";
    for(const auto& s:g.sensors) f << s.target << " " << s.angle << " " << s.range << "\n";
    for(const auto& a:g.actuators) f << a.type << "\n";
    for(const auto& n:g.nodes) f << n.id << " " << n.type << "\n";
    for(const auto& c:g.conns) f << c.in_node << " " << c.out_node << " " << c.weight << " " << (c.enabled?1:0) << " " << c.innov << "\n";
}
Genome read_genome(std::ifstream& f){
    Genome g; std::string tag; int ns=0,na=0,nn=0,nc=0;
    f >> tag >> ns >> na >> nn >> nc;
    if(!f || tag!="G" || ns<0||ns>10000 || na<0||na>1000 || nn<0||nn>200000 || nc<0||nc>2000000)
        throw std::runtime_error("bad genome record");
    g.sensors.resize(ns); for(auto& s:g.sensors) f >> s.target >> s.angle >> s.range;
    g.actuators.resize(na); for(auto& a:g.actuators) f >> a.type;
    g.nodes.resize(nn); for(auto& n:g.nodes) f >> n.id >> n.type;
    g.conns.resize(nc); for(auto& c:g.conns){ int en; f >> c.in_node >> c.out_node >> c.weight >> en >> c.innov; c.enabled=(en!=0); }
    g.finalize();
    return g;
}
void save_state(const std::string& fname, std::vector<Agent>& preys, std::vector<Agent>& predators, int frame){
    std::ofstream f(fname); if(!f){ printf("save failed: %s\n",fname.c_str()); return; }
    f.precision(9);
    f << "ECOSAVE 3\n";
    f << generation << " " << g_innov_counter << " " << g_next_node_id << " " << frame << "\n";
    int npy=0; for(auto&p:preys) if(p.alive)npy++;
    f << "PREYS " << npy << "\n";
    for(auto&p:preys) if(p.alive){ f << p.genes.size << " " << p.genes.vision << " " << p.genes.eat_gain << " " << p.genes.speed << "\n"; write_genome(f,p.genome); }
    int npd=0; for(auto&pd:predators) if(pd.alive)npd++;
    f << "PREDATORS " << npd << "\n";
    for(auto&pd:predators) if(pd.alive){ f << pd.genes.size << " " << pd.genes.vision << " " << pd.genes.eat_gain << " " << pd.genes.speed << "\n"; write_genome(f,pd.genome); }
    f << "PREYHALL " << prey_hall.size() << "\n";
    for(auto&r:prey_hall){ f << r.fit << " " << r.genes.size << " " << r.genes.vision << " " << r.genes.eat_gain << " " << r.genes.speed << "\n"; write_genome(f,r.genome); }
    f << "LINEAGES " << lineages.size() << " " << next_lineage_id << "\n";
    for(auto& L:lineages){
        f << L.id << " " << L.birth_tick << " " << L.death_tick << " " << L.parent << " " << (L.alive?1:0) << " " << L.pop << " " << L.max_pop << "\n";
        write_genome(f, L.rep);
    }
    f << "PREDHALL " << pred_hall.size() << "\n";
    for(auto&r:pred_hall){ f << r.first << "\n"; write_genome(f,r.second); }
    printf("saved: %s (prey %d, pred %d)\n",fname.c_str(),npy,npd);
}
bool load_state(const std::string& fname, std::vector<Agent>& preys, std::vector<Agent>& predators, int& frame){
    std::ifstream f(fname); if(!f){ printf("load failed (no file): %s\n",fname.c_str()); return false; }
    std::string tag; int ver; f >> tag >> ver;
    if(tag!="ECOSAVE"){ printf("load failed (bad format)\n"); return false; }
    try {
    if(ver>3){ printf("load failed (newer version %d)\n",ver); return false; }
    f >> generation >> g_innov_counter >> g_next_node_id;
    if(ver>=2) f >> frame; else frame=0;
    preys.clear(); predators.clear(); prey_hall.clear(); pred_hall.clear();
    std::string sec; int cnt;
    f >> sec >> cnt;
    for(int i=0;i<cnt;i++){ Genes ge; f >> ge.size >> ge.vision >> ge.eat_gain >> ge.speed; Genome g=read_genome(f); Agent a=make_prey(ge); a.genome=g; preys.push_back(std::move(a)); }
    f >> sec >> cnt;
    for(int i=0;i<cnt;i++){ Genes ge; f >> ge.size >> ge.vision >> ge.eat_gain >> ge.speed; Genome g=read_genome(f); Agent a=make_predator(); a.genes=ge; a.genome=g; predators.push_back(std::move(a)); }
    f >> sec >> cnt;
    for(int i=0;i<cnt;i++){ PreyRec r; f >> r.fit >> r.genes.size >> r.genes.vision >> r.genes.eat_gain >> r.genes.speed; r.genome=read_genome(f); prey_hall.push_back(std::move(r)); }
    lineages.clear(); next_lineage_id=0;
    if(ver>=3){
        f >> sec >> cnt >> next_lineage_id;
        for(int i=0;i<cnt;i++){
            Lineage L; int en;
            f >> L.id >> L.birth_tick >> L.death_tick >> L.parent >> en >> L.pop >> L.max_pop;
            L.alive=(en!=0); L.rep=read_genome(f);
            lineages.push_back(L);
        }
    }
    f >> sec >> cnt;
    for(int i=0;i<cnt;i++){ int fit; f >> fit; Genome g=read_genome(f); pred_hall.push_back({fit,g}); }
    printf("loaded: %s (prey %d, pred %d)\n",fname.c_str(),(int)preys.size(),(int)predators.size());
    return true;
    } catch(const std::exception& e){
        printf("load failed (corrupt/incompatible): %s\n", e.what());
        preys.clear(); predators.clear(); lineages.clear();
        return false;
    }
}

sf::Color cluster_color(int i){
    static const sf::Color pal[16]={
        {0,120,255},{60,220,60},{200,80,255},{0,220,220},
        {120,255,160},{255,120,200},{100,140,255},{160,160,160},
        {80,200,120},{140,100,220},{60,180,255},{200,120,255},
        {120,220,180},{255,150,220},{100,255,200},{180,180,255}
    };
    return pal[((i%16)+16)%16];
}

// 生きているpreyを種に分け、血統(履歴)を追跡する。観測のみ、行動には影響しない。
void update_lineages(std::vector<Agent>& preys, float threshold, int tick){
    // 1. 今回のクラスタを作る(その場)
    std::vector<Genome> reps;
    for(auto& p:preys){
        if(!p.alive) continue;
        int best=-1; float bestd=1e9f;
        for(int i=0;i<(int)reps.size();i++){
            float d=genome_distance(p.genome, reps[i]);
            if(d<bestd){ bestd=d; best=i; }
        }
        if(best<0 || bestd>=threshold){ best=(int)reps.size(); reps.push_back(p.genome); }
        p.cluster=best;   // 一時的な今回インデックス
    }
    // 2. 今回のクラスタを、既存の血統に照合
    int oldN=(int)lineages.size();
    std::vector<int> cluster_lineage(reps.size(), -1);
    std::vector<bool> claimed(oldN, false);
    for(int c=0;c<(int)reps.size();c++){
        int best=-1; float bestd=1e9f;
        for(int L=0;L<oldN;L++){
            if(!lineages[L].alive) continue;
            float d=genome_distance(reps[c], lineages[L].rep);
            if(d<bestd){ bestd=d; best=L; }
        }
        if(best>=0 && bestd<threshold && !claimed[best]){
            cluster_lineage[c]=best; claimed[best]=true; lineages[best].rep=reps[c];  // 継続
        } else {
            Lineage nl; nl.id=next_lineage_id++; nl.rep=reps[c];
            nl.birth_tick=tick; nl.death_tick=-1;
            nl.parent=(best>=0)?lineages[best].id:-1; nl.alive=true; nl.pop=0; nl.max_pop=0;  // 誕生(分岐)
            cluster_lineage[c]=(int)lineages.size(); lineages.push_back(nl);
        }
    }
    // 3. 今回どのクラスタにも一致しなかった血統 → 絶滅
    for(int L=0;L<oldN;L++){
        if(lineages[L].alive && !claimed[L]){ lineages[L].alive=false; lineages[L].death_tick=tick; }
    }
    // 4. preyに血統IDを割り当て + 個体数を数える
    for(auto& lin:lineages) lin.pop=0;
    for(auto& p:preys){
        if(!p.alive) continue;
        int Lidx=cluster_lineage[p.cluster];
        p.cluster=lineages[Lidx].id;      // ← p.cluster に"消えない血統ID"が入る
        lineages[Lidx].pop++;
    }
    for(auto& lin:lineages) if(lin.pop>lin.max_pop) lin.max_pop=lin.pop;
}

int main(){
// === 段階5a-1(改)テスト: 交叉 + 流産検査 ===
    // === OpenMP 動作確認 ===
    #ifdef _OPENMP
        printf("OpenMP is ENABLED. max threads = %d\n", omp_get_max_threads());
    #else
        printf("OpenMP is NOT enabled\n");
    #endif
    
    sf::VideoMode desktop = sf::VideoMode::getDesktopMode();
    SCREEN_W = (float)desktop.size.x;
    SCREEN_H = (float)desktop.size.y;
    zoom = SCREEN_H/HEIGHT;   // 新しい画面高さに合わせて初期ズームを再計算
    sf::RenderWindow window(desktop,"Ecosystem C++",sf::State::Fullscreen);
    sf::Font font;
    bool font_ok = font.openFromFile("/System/Library/Fonts/Helvetica.ttc");
    if(!font_ok) font_ok = font.openFromFile("/System/Library/Fonts/Supplemental/Arial.ttf");
    if(!font_ok) font_ok = font.openFromFile("C:\\Windows\\Fonts\\arial.ttf");
    if(!font_ok) font_ok = font.openFromFile("C:\\Windows\\Fonts\\segoeui.ttf");
    window.setFramerateLimit(60);

    std::vector<Agent> preys, predators;
    for(int i=0;i<NUM_PREY;i++)preys.push_back(make_prey(random_genes()));
    for(int i=0;i<NUM_PREDATOR;i++)predators.push_back(make_predator());
    std::vector<Plant> plants;
    for(int i=0;i<NUM_PLANT;i++)plants.push_back({frand(0,WIDTH),frand(0,HEIGHT),0});
    std::vector<Meat> meats;
    Field nutrient(NUT_COLS,NUT_ROWS,NUT_DIFFUSION,NUT_DECAY);
    Field ph_fear(PH_COLS,PH_ROWS,PH_DIFFUSION,PH_DECAY);
    Field ph_aff(PH_COLS,PH_ROWS,PH_DIFFUSION,PH_DECAY);

    std::vector<std::vector<int>> plant_grid(GRID_COLS*GRID_ROWS),prey_grid(GRID_COLS*GRID_ROWS),pred_grid(GRID_COLS*GRID_ROWS),meat_grid(GRID_COLS*GRID_ROWS);
    int pred_cr=(int)std::ceil(PREDATOR_VISION/CELL_SIZE);

    sf::CircleShape pc(3.f), plc(2.f), rc(3.f), mc(2.f);
    plc.setFillColor(sf::Color(0,200,0));
    sf::RectangleShape fieldRect;
    bool show_nutrient=false, show_phero=false, dragging=false;
    float drag_sx=0,drag_sy=0,cam_sx=0,cam_sy=0;
    int frame=0;
    int autosave_timer=0;
    int autosave_slot=0;
    const int AUTOSAVE_INTERVAL=7200;   // 2分ごと(60fps×120秒)
    int selected_id=-1;
    float species_threshold=0.35f;
    bool show_tree=false;
    int ui_mode=0;                    // 0=通常, 1=セーブ名入力, 2=ロード選択
    std::string input_text="";
    bool ui_skip_char=false;
    std::vector<std::string> save_files;

    float measured_fps=60.f;   // 実測FPS(初期値60)
    auto last_time=std::chrono::high_resolution_clock::now();

    int peak_prey=0, peak_pred=0, peak_plant=0, peak_meat=0;

    long prey_starve=0, prey_killed=0, pred_starve=0;
    
    std::vector<int> hist_prey, hist_pred, hist_plant, hist_meat;
    
    auto build=[&](std::vector<std::vector<int>>& g, auto& items){
        for(auto&c:g)c.clear();
        for(int i=0;i<(int)items.size();i++){
            int cx=(int)items[i].x/CELL_SIZE,cy=(int)items[i].y/CELL_SIZE;
            if(cx<0)cx=0;if(cx>=GRID_COLS)cx=GRID_COLS-1;if(cy<0)cy=0;if(cy>=GRID_ROWS)cy=GRID_ROWS-1;
            g[cell_index(cx,cy)].push_back(i);
        }
    };
    bool paused=false;
    bool is_fullscreen=true;   // 起動時は全画面

    char prof_str[256]="";

    while(window.isOpen()){
        while(const std::optional event=window.pollEvent()){
            if(event->is<sf::Event::Closed>())window.close();
            // 文字入力(セーブ名/ロード名)
            if(const auto* te=event->getIf<sf::Event::TextEntered>()){
                if(ui_mode!=0){
                    char32_t u=te->unicode;
                    if(ui_skip_char){ ui_skip_char=false; }                          // モード開始のk/lを捨てる
                    else if(u==8){ if(!input_text.empty()) input_text.pop_back(); }
                    else if(u>=32 && u<127){ if(input_text.size()<40) input_text+=(char)u; }
                }
            }
            if(const auto* k=event->getIf<sf::Event::KeyPressed>()){
                if(ui_mode!=0){
                    // 入力モード中: Enterで確定, Escでキャンセル
                    if(k->code==sf::Keyboard::Key::Enter){
                        if(!input_text.empty()){
                            if(ui_mode==1) save_state(input_text+".txt", preys, predators, frame);
                            else if(ui_mode==2) load_state(input_text+".txt", preys, predators, frame);
                        }
                        ui_mode=0; input_text="";
                    }
                    if(k->code==sf::Keyboard::Key::Escape){ ui_mode=0; input_text=""; }
                } else {
                    if(k->code==sf::Keyboard::Key::G)show_nutrient=!show_nutrient;
                    if(k->code==sf::Keyboard::Key::F)show_phero=!show_phero;
                    if(k->code==sf::Keyboard::Key::T)show_tree=!show_tree;
                    if(k->code==sf::Keyboard::Key::Space)paused=!paused;
                    if(k->code==sf::Keyboard::Key::Escape){
                    is_fullscreen=!is_fullscreen;
                    if(is_fullscreen){
                        window.create(desktop,"Ecosystem C++",sf::State::Fullscreen);
                        SCREEN_W=(float)desktop.size.x; SCREEN_H=(float)desktop.size.y;
                    }else{
                        window.create(sf::VideoMode({1600,900}),"Ecosystem C++",sf::Style::Titlebar|sf::Style::Close);
                        SCREEN_W=1600.f; SCREEN_H=900.f;
                    }
                    window.setFramerateLimit(60);
                    }
                    if(k->code==sf::Keyboard::Key::K){ ui_mode=1; input_text=""; ui_skip_char=true; }
                    if(k->code==sf::Keyboard::Key::L){
                        ui_mode=2; input_text=""; ui_skip_char=true;
                        save_files.clear();
                        try{
                            for(auto& e: std::filesystem::directory_iterator(".")){
                                if(e.is_regular_file() && e.path().extension().string()==".txt")
                                    save_files.push_back(e.path().stem().string());
                            }
                            std::sort(save_files.begin(),save_files.end());
                        }catch(...){}
                    }
                }
            }
            if(const auto* mb=event->getIf<sf::Event::MouseButtonPressed>()){
                if(mb->button==sf::Mouse::Button::Left){
                    float wx=s2w_x(mb->position.x), wy=s2w_y(mb->position.y);
                    selected_id=-1; float best=1e18f;
                    for(auto& p:preys){float dx=wx-p.x,dy=wy-p.y;float d2=dx*dx+dy*dy;if(d2<best){best=d2;selected_id=p.id;}}
                    for(auto& pd:predators){float dx=wx-pd.x,dy=wy-pd.y;float d2=dx*dx+dy*dy;if(d2<best){best=d2;selected_id=pd.id;}}
                    float lim=50.f/zoom; if(best>lim*lim)selected_id=-1;
                }
            }
            if(const auto* mb=event->getIf<sf::Event::MouseButtonReleased>()){
                if(mb->button==sf::Mouse::Button::Left)dragging=false;
            }
            if(const auto* mw=event->getIf<sf::Event::MouseWheelScrolled>()){
                sf::Vector2i mp=sf::Mouse::getPosition(window);
                float mx=(float)mp.x,my=(float)mp.y;
                float wbx=s2w_x(mx),wby=s2w_y(my);
                if(mw->delta>0.f)zoom*=1.25f; else if(mw->delta<0.f)zoom/=1.25f;
                float zmin=SCREEN_H/HEIGHT*0.5f; if(zoom<zmin)zoom=zmin; if(zoom>5.f)zoom=5.f;
                cam_x+=wbx-s2w_x(mx); cam_y+=wby-s2w_y(my);
            }
        }
        if(dragging){
            sf::Vector2i m=sf::Mouse::getPosition(window);
            cam_x=cam_sx-(m.x-drag_sx)/zoom; cam_y=cam_sy-(m.y-drag_sy)/zoom;
        }

        // 実測FPS(なめらかに追従)
        auto cur_time=std::chrono::high_resolution_clock::now();
        float dt=std::chrono::duration<float>(cur_time-last_time).count();
        last_time=cur_time;
        if(dt>0.f) measured_fps=measured_fps*0.95f + (1.f/dt)*0.05f;  // 平滑化
        
        if(!paused){
    
        auto _t6=now();
        build(plant_grid,plants); build(prey_grid,preys); build(pred_grid,predators);
        t_grid += ms(_t6,now());

        std::vector<Plant> new_plants;
        for(auto& pl:plants){
            if(++pl.repro_counter>=PLANT_REPRO_TIME){pl.repro_counter=0;
                if((int)(plants.size()+new_plants.size())<PLANT_MAX)
                    for(int c=0;c<PLANT_CHILDREN;c++){
                        float nx=pl.x+frand(-PLANT_SPREAD,PLANT_SPREAD),ny=pl.y+frand(-PLANT_SPREAD,PLANT_SPREAD);
                        if(nx<0)nx+=WIDTH;if(nx>=WIDTH)nx-=WIDTH;if(ny<0)ny+=HEIGHT;if(ny>=HEIGHT)ny-=HEIGHT;
                        float bonus=std::min(1.f,nutrient.get(NUT_CELL,nx,ny)/NUT_FULL)*(1.f-PLANT_BASE_SURVIVAL);
                        if(dist01(rng)<PLANT_BASE_SURVIVAL+bonus){new_plants.push_back({nx,ny,0});nutrient.emit(NUT_CELL,nx,ny,-NUT_CONSUME);}
                    }
            }
        }
        for(auto&np:new_plants)plants.push_back(np);

        // 植物が少なすぎたら、ランダムな場所に補充(絶滅防止・滑らかに底支え)
        if((int)plants.size() < PLANT_MIN){
            int need = PLANT_MIN - (int)plants.size();
            for(int i=0;i<need;i++) plants.push_back({frand(0,WIDTH),frand(0,HEIGHT),0});
        }

        auto _t0=now();
        #pragma omp parallel for schedule(dynamic, 64)
        for(int idx=0; idx<(int)preys.size(); idx++){
            Agent& p = preys[idx];
            if(!p.alive)continue;
            float vp[NUM_RAYS]={0},vpr[NUM_RAYS]={0},vself[NUM_RAYS]={0};   // vself追加
            float facing=std::atan2(p.vy,p.vx); float pv=p.genes.vision;
            int cr=(int)std::ceil(pv/CELL_SIZE); int scx=(int)p.x/CELL_SIZE,scy=(int)p.y/CELL_SIZE;
            for(int dx=-cr;dx<=cr;dx++)for(int dy=-cr;dy<=cr;dy++){
                int cx=scx+dx,cy=scy+dy; if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                for(int pi:plant_grid[cell_index(cx,cy)]){
                    float ddx=torus_delta(p.x,plants[pi].x,WIDTH),ddy=torus_delta(p.y,plants[pi].y,HEIGHT);
                    float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                    float d=std::sqrt(d2);float rel=std::fmod(std::atan2(ddy,ddx)-facing+TWO_PI*2,TWO_PI);
                    int ri=(int)(rel/TWO_PI*NUM_RAYS+0.5f)%NUM_RAYS;float val=1.f-d/pv;if(val>vp[ri])vp[ri]=val;
                }
                for(int pi:pred_grid[cell_index(cx,cy)]){
                    if(!predators[pi].alive)continue;
                    float ddx=torus_delta(p.x,predators[pi].x,WIDTH),ddy=torus_delta(p.y,predators[pi].y,HEIGHT);
                    float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                    float d=std::sqrt(d2);float rel=std::fmod(std::atan2(ddy,ddx)-facing+TWO_PI*2,TWO_PI);
                    int ri=(int)(rel/TWO_PI*NUM_RAYS+0.5f)%NUM_RAYS;float val=1.f-d/pv;if(val>vpr[ri])vpr[ri]=val;
                }
                for(int pi:prey_grid[cell_index(cx,cy)]){
                    if(!preys[pi].alive || preys[pi].id==p.id)continue;   // 自分は除く
                    float ddx=torus_delta(p.x,preys[pi].x,WIDTH),ddy=torus_delta(p.y,preys[pi].y,HEIGHT);
                    float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                    float d=std::sqrt(d2);float rel=std::fmod(std::atan2(ddy,ddx)-facing+TWO_PI*2,TWO_PI);
                    int ri=(int)(rel/TWO_PI*NUM_RAYS+0.5f)%NUM_RAYS;float val=1.f-d/pv;if(val>vself[ri])vself[ri]=val;
                }
            }
            // Genome用のセンサー値(37個)を作る
            std::vector<float> sv(37, 0.f);
            for(int i=0;i<12;i++){ sv[i]=vp[i]; sv[12+i]=vpr[i]; sv[24+i]=vself[i]; }   // 同種も入れる
            sv[36]=p.energy/p.max_energy;
            // 感情用の記録
            bool sp=false, sa=false;
            for(int i=0;i<12;i++){ if(vpr[i]>0.01f)sp=true; if(vself[i]>0.01f)sa=true; }
            p.saw_pred=sp; p.saw_ally=sa;
            std::vector<float> out;
            p.genome.forward(sv, out);
            p.vx=out[0]*p.genes.speed; p.vy=out[1]*p.genes.speed;
        }
        t_prey_nn += ms(_t0,now());
        auto _t1=now();
        #pragma omp parallel for schedule(dynamic, 32)
        for(int idx=0; idx<(int)predators.size(); idx++){
            Agent& pd = predators[idx];
            if(!pd.alive)continue;
            float vpr[NUM_RAYS]={0};float facing=std::atan2(pd.vy,pd.vx);
            int scx=(int)pd.x/CELL_SIZE,scy=(int)pd.y/CELL_SIZE;
            for(int dx=-pred_cr;dx<=pred_cr;dx++)for(int dy=-pred_cr;dy<=pred_cr;dy++){
                int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                for(int pi:prey_grid[cell_index(cx,cy)]){
                    if(!preys[pi].alive)continue;
                    float ddx=torus_delta(pd.x,preys[pi].x,WIDTH),ddy=torus_delta(pd.y,preys[pi].y,HEIGHT);
                    float d2=ddx*ddx+ddy*ddy;if(d2>PREDATOR_VISION*PREDATOR_VISION||d2<1)continue;
                    float rel=std::fmod(std::atan2(ddy,ddx)-facing+PI*3,TWO_PI)-PI;
                    if(rel<-PREDATOR_FOV/2||rel>PREDATOR_FOV/2)continue;
                    float d=std::sqrt(d2);int ri=(int)((rel+PREDATOR_FOV/2)/PREDATOR_FOV*(NUM_RAYS-1)+0.5f);
                    if(ri<0)ri=0;if(ri>=NUM_RAYS)ri=NUM_RAYS-1;float val=1.f-d/PREDATOR_VISION;if(val>vpr[ri])vpr[ri]=val;
                }
            }
            // Genome用センサー値(37個: prey12 + 肉12 + 同種12 + energy)
            std::vector<float> sv(37, 0.f);
            for(int i=0;i<12;i++){ sv[i]=vpr[i]; /* sv[12+i]=肉(後で), sv[24+i]=同種(後で) */ }
            sv[36]=pd.energy/pd.max_energy;
            std::vector<float> out;
            pd.genome.forward(sv, out);
            pd.vx=out[0]*SPEED; pd.vy=out[1]*SPEED;
        }
        t_pred_nn += ms(_t1,now());

        auto _t2=now();
        for(auto& p:preys){
            if(!p.alive)continue;
            p.x+=p.vx;p.y+=p.vy;
            if(p.x<0)p.x+=WIDTH;if(p.x>=WIDTH)p.x-=WIDTH;if(p.y<0)p.y+=HEIGHT;if(p.y>=HEIGHT)p.y-=HEIGHT;
            p.age++;
            float sp=std::sqrt(p.vx*p.vx+p.vy*p.vy);
            p.energy-=p.genes.size*METABOLISM_COEF+0.5f*p.genes.size*sp*sp*MOVE_COEF+p.genome.cached_active_conns*BRAIN_COST_COEF;
            if(p.energy<=0){p.alive=false;prey_starve++;meats.push_back({p.x,p.y,p.genes.size*MEAT_SIZE_COEF+std::max(0.f,p.energy)});continue;}
            bool pred_sight=p.saw_pred, ally_sight=p.saw_ally;   // NNで計算済みを使う
            if(pred_sight)p.fear+=FEAR_RISE;else p.fear-=FEAR_DECAY;p.fear=clamp01(p.fear);
            if(ally_sight)p.affinity+=AFFINITY_RISE;else p.affinity-=AFFINITY_DECAY;p.affinity=clamp01(p.affinity);
            ph_fear.emit(PH_CELL,p.x,p.y,p.fear*0.1f); ph_aff.emit(PH_CELL,p.x,p.y,p.affinity*0.1f);
            float red=ph_fear.get(PH_CELL,p.x,p.y),blue=ph_aff.get(PH_CELL,p.x,p.y);
            bool R=red>=PHERO_THRESHOLD,B=blue>=PHERO_THRESHOLD;
            if(R&&B)p.stress+=STRESS_MIX;else if(R)p.stress+=STRESS_RISE;else if(B)p.stress-=STRESS_FALL;else p.stress-=STRESS_DECAY;
            if(pred_sight)p.stress+=STRESS_PREDATOR_SIGHT; p.stress=clamp01(p.stress);
            // 交叉条件: ストレス低 かつ 友情満タン が続いたら、カウンタを増やす
            if(p.stress < 0.05f && p.affinity > 0.95f) p.mate_ready_counter++;
            else p.mate_ready_counter = 0;
        }
        for(auto& pd:predators){
            if(!pd.alive)continue;
            pd.x+=pd.vx;pd.y+=pd.vy;
            if(pd.x<0)pd.x+=WIDTH;if(pd.x>=WIDTH)pd.x-=WIDTH;if(pd.y<0)pd.y+=HEIGHT;if(pd.y>=HEIGHT)pd.y-=HEIGHT;
            pd.age++;pd.energy-=METABOLISM;
            if(pd.energy<=0){pd.alive=false;pred_starve++;meats.push_back({pd.x,pd.y,pd.genes.size*MEAT_SIZE_COEF+std::max(0.f,pd.energy)});}
        }
        t_move += ms(_t2,now());

        auto _t3=now();
        for(auto& p:preys){
            if(!p.alive)continue;
            int cxmin=(int)(p.x-EAT_DISTANCE)/CELL_SIZE, cxmax=(int)(p.x+EAT_DISTANCE)/CELL_SIZE;
            int cymin=(int)(p.y-EAT_DISTANCE)/CELL_SIZE, cymax=(int)(p.y+EAT_DISTANCE)/CELL_SIZE;
            for(int cx=cxmin;cx<=cxmax;cx++)for(int cy=cymin;cy<=cymax;cy++){
                if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                for(int pi:plant_grid[cell_index(cx,cy)]){
                    if(plants[pi].x<0)continue;
                    if(frame - p.last_eat_frame < 18)continue;   // 0.3秒(18フレーム)クールタイム
                    float ddx=torus_delta(p.x,plants[pi].x,WIDTH),ddy=torus_delta(p.y,plants[pi].y,HEIGHT);
                    if(ddx*ddx+ddy*ddy<EAT_DISTANCE*EAT_DISTANCE){p.energy+=p.genes.eat_gain; if(p.energy>p.max_energy)p.energy=p.max_energy; plants[pi].x=-1; p.food_count++; p.last_eat_frame=frame;}
                }
            }
        }
        plants.erase(std::remove_if(plants.begin(),plants.end(),[](const Plant&p){return p.x<0;}),plants.end());

        for(auto& pd:predators){
            if(!pd.alive)continue;
            int scx=(int)pd.x/CELL_SIZE,scy=(int)pd.y/CELL_SIZE;
            for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++){
                int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                for(int pi:prey_grid[cell_index(cx,cy)]){
                    if(!preys[pi].alive)continue;
                    float ddx=torus_delta(pd.x,preys[pi].x,WIDTH),ddy=torus_delta(pd.y,preys[pi].y,HEIGHT);
                    if(ddx*ddx+ddy*ddy<EAT_DISTANCE*EAT_DISTANCE){
                        preys[pi].alive=false; prey_killed++; pd.food_count++;prey_hall.push_back({calc_fit(preys[pi]),preys[pi].genome,preys[pi].genes});
                        meats.push_back({preys[pi].x,preys[pi].y,preys[pi].genes.size*MEAT_SIZE_COEF+std::max(0.f,preys[pi].energy)});
                    }
                }
            }
        }
        build(meat_grid,meats);
        for(auto& pd:predators){
            if(!pd.alive)continue;
            int scx=(int)pd.x/CELL_SIZE,scy=(int)pd.y/CELL_SIZE;
            for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++){
                int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                for(int mi:meat_grid[cell_index(cx,cy)]){
                    if(meats[mi].energy<=0)continue;
                    float ddx=torus_delta(pd.x,meats[mi].x,WIDTH),ddy=torus_delta(pd.y,meats[mi].y,HEIGHT);
                    if(ddx*ddx+ddy*ddy<EAT_DISTANCE*EAT_DISTANCE){float bite=std::min(MEAT_EAT_RATE,meats[mi].energy);meats[mi].energy-=bite;pd.energy+=bite;if(pd.energy>pd.max_energy)pd.energy=pd.max_energy;}
                }
            }
        }
        for(auto& m:meats){float rot=std::min(MEAT_ROT_RATE,m.energy);m.energy-=rot;nutrient.emit(NUT_CELL,m.x,m.y,rot);}
        meats.erase(std::remove_if(meats.begin(),meats.end(),[](const Meat&m){return m.energy<=0;}),meats.end());
        t_eat += ms(_t3,now());
        auto _t4=now();
        if(frame%3==0){ nutrient.update();ph_fear.update();ph_aff.update(); }
        t_field += ms(_t4,now());
        if(frame%180==0) update_lineages(preys, species_threshold, frame);

        std::vector<Agent> babies;
        for(auto& p:preys){
            if(!p.alive)continue;
            if(p.energy>p.max_energy*(5.f/6.f)){if(++p.repro_counter>=PREY_REPRO_TIME){p.repro_counter=0;
                float success=1.f-p.stress*STRESS_REPRO_PENALTY;
                if(dist01(rng)<success){float cost=p.genes.size*REPRO_COST_COEF;p.energy-=cost; p.repro_count++;
                    Genome child_genome; bool miscarried=false;
                    // 交叉可能状態(4秒継続)なら、すぐ隣の個体と交叉を試みる
                    const Agent* mate=nullptr;
                    if(p.mate_ready_counter>=240){
                        int scx=(int)p.x/CELL_SIZE,scy=(int)p.y/CELL_SIZE;
                        for(int dx=-1;dx<=1&&!mate;dx++)for(int dy=-1;dy<=1&&!mate;dy++){
                            int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                            for(int pi:prey_grid[cell_index(cx,cy)]){
                                if(!preys[pi].alive||preys[pi].id==p.id)continue;
                                float ddx=torus_delta(p.x,preys[pi].x,WIDTH),ddy=torus_delta(p.y,preys[pi].y,HEIGHT);
                                if(ddx*ddx+ddy*ddy < (EAT_DISTANCE*2)*(EAT_DISTANCE*2)){ mate=&preys[pi]; break; }
                            }
                        }
                    }
                    if(mate){
                        child_genome = crossover(p.genome, mate->genome);
                        child_genome = mutate_genome(child_genome);
                        if(!is_valid_genome(child_genome)) miscarried=true;
                    } else {
                        child_genome = mutate_genome(p.genome);  // 分裂
                    }
                    if(miscarried){
                        // 流産 → 肉ドロップ、子は生まれない
                        meats.push_back({p.x,p.y,p.genes.size*MEAT_SIZE_COEF*0.3f});
                    } else {
                        Agent c=make_prey(mutate_genes(p.genes,p.stress));
                        c.genome=child_genome;
                        c.cluster=p.cluster;   // 母の血統を継ぐ(次のクラスタリングまで正しい色)
                        c.x=p.x+frand(-10,10);c.y=p.y+frand(-10,10);c.energy=cost*0.5f;babies.push_back(std::move(c));
                    }
                }
            }}else p.repro_counter=0;
        }
        for(auto&b:babies)preys.push_back(std::move(b));
        std::vector<Agent> pbabies;
        for(auto& pd:predators){
            if(!pd.alive)continue;
            if(pd.energy>pd.max_energy*(5.f/6.f)){if(++pd.repro_counter>=PREDATOR_REPRO_TIME){pd.repro_counter=0;pd.energy*=0.5f;pd.repro_count++;
                Agent c=make_predator();c.genome=mutate_genome(pd.genome);c.x=pd.x+frand(-10,10);c.y=pd.y+frand(-10,10);c.energy=pd.energy;pbabies.push_back(std::move(c));
            }}else pd.repro_counter=0;
        }
        for(auto&b:pbabies)predators.push_back(std::move(b));

        {
            std::vector<Agent> al;
            for(auto&p:preys){if(p.alive)al.push_back(std::move(p));else prey_hall.push_back({calc_fit(p),p.genome,p.genes});}
            preys=std::move(al);
            std::vector<Agent> ap;
            for(auto&pd:predators){if(pd.alive)ap.push_back(std::move(pd));else pred_hall.push_back({calc_fit(pd),pd.genome});}
            predators=std::move(ap);
        }
        if((int)prey_hall.size()>HALL_SIZE*10)trim_prey();
        if((int)pred_hall.size()>HALL_SIZE*10)trim_pred();

        // 捕食者だけ絶滅 → 殿堂から本番世界に再導入(世界はリセットしない=連続性を保つ)
        if(predators.empty() && !preys.empty() && !pred_hall.empty()){
            for(int i=0;i<NUM_PREDATOR;i++){
                int ix=(int)(dist01(rng)*pred_hall.size()); if(ix>=(int)pred_hall.size())ix=pred_hall.size()-1;
                Agent na=make_predator(); na.genome=mutate_genome(pred_hall[ix].second);
                predators.push_back(std::move(na));
            }
            printf("Predators re-introduced from hall (frame %d)\n", frame);
        }

        if(preys.empty()){
            for(auto&p:preys)prey_hall.push_back({calc_fit(p),p.genome,p.genes});
            for(auto&pd:predators)pred_hall.push_back({calc_fit(pd),pd.genome});
            trim_prey();trim_pred();
            int bf=prey_hall.empty()?0:prey_hall[0].fit,bpf=pred_hall.empty()?0:pred_hall[0].first;
            printf("Gen %d: survived %ds, prey fit %d, pred fit %d\n",generation,frame/60,bf,bpf);
            generation++;preys.clear();predators.clear();meats.clear();
            nutrient.grid.assign(NUT_COLS*NUT_ROWS,0.f);ph_fear.grid.assign(PH_COLS*PH_ROWS,0.f);ph_aff.grid.assign(PH_COLS*PH_ROWS,0.f);
            for(int i=0;i<NUM_PREY;i++){if(!prey_hall.empty()){int ix=(int)(dist01(rng)*prey_hall.size());if(ix>=(int)prey_hall.size())ix=prey_hall.size()-1;Agent na=make_prey(prey_hall[ix].genes);na.genome=mutate_genome(prey_hall[ix].genome);preys.push_back(std::move(na));}else preys.push_back(make_prey(random_genes()));}
            for(int i=0;i<NUM_PREDATOR;i++){if(!pred_hall.empty()){int ix=(int)(dist01(rng)*pred_hall.size());if(ix>=(int)pred_hall.size())ix=pred_hall.size()-1;Agent na=make_predator();na.genome=mutate_genome(pred_hall[ix].second);predators.push_back(std::move(na));}else predators.push_back(make_predator());}
            plants.clear();for(int i=0;i<NUM_PLANT;i++)plants.push_back({frand(0,WIDTH),frand(0,HEIGHT),0});
            frame=0;
        }
    }
        // WASDでカメラ移動
        float pan=15.f/zoom;
        if(ui_mode==0){
            if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::W))cam_y-=pan;
            if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S))cam_y+=pan;
            if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A))cam_x-=pan;
            if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D))cam_x+=pan;
        }
        if(cam_x<0)cam_x=0; if(cam_x>WIDTH)cam_x=WIDTH;
        if(cam_y<0)cam_y=0; if(cam_y>HEIGHT)cam_y=HEIGHT;

        // 描画(カメラ変換 + カリング)
        auto _t5=now();
        window.clear(sf::Color::Black);
        // 世界の範囲を暗いグレーで描く(境界表示)
        {
            sf::RectangleShape worldRect;
            worldRect.setSize({WIDTH*zoom, HEIGHT*zoom});
            worldRect.setPosition({w2s_x(0.f), w2s_y(0.f)});
            worldRect.setFillColor(sf::Color(30,30,30));      // 内側=暗いグレー
            worldRect.setOutlineColor(sf::Color(80,80,80));   // 枠線=明るいグレー
            worldRect.setOutlineThickness(1.f);
            window.draw(worldRect);
        }
        
        if(show_nutrient){
            float sz=NUT_CELL*zoom;
            for(int cy=0;cy<NUT_ROWS;cy++)for(int cx=0;cx<NUT_COLS;cx++){
                float v=nutrient.grid[cy*NUT_COLS+cx]; if(v<=0.1f)continue;
                float wx=cx*NUT_CELL,wy=cy*NUT_CELL; if(!on_screen(wx,wy,sz))continue;
                int g=(int)std::min(255.f,v*80.f);
                int a=(int)std::min(200.f,v*60.f); 
                fieldRect.setSize({sz,sz});fieldRect.setPosition({w2s_x(wx),w2s_y(wy)});fieldRect.setFillColor(sf::Color(0,g,0,a));window.draw(fieldRect);
            }
        }
        if(show_phero){
            float sz=PH_CELL*zoom;
            for(int cy=0;cy<PH_ROWS;cy++)for(int cx=0;cx<PH_COLS;cx++){
                float fv=ph_fear.grid[cy*PH_COLS+cx],av=ph_aff.grid[cy*PH_COLS+cx]; if(fv<=0.01f&&av<=0.01f)continue;
                float mx=std::max(fv,av); if(mx<=0.02f)continue;
                float wx=cx*PH_CELL,wy=cy*PH_CELL; if(!on_screen(wx,wy,sz))continue;
                int r=(int)std::min(255.f,fv*500.f),b=(int)std::min(255.f,av*500.f);
                int a=(int)std::min(200.f,mx*800.f);
                fieldRect.setSize({sz,sz});fieldRect.setPosition({w2s_x(wx),w2s_y(wy)});fieldRect.setFillColor(sf::Color(r,0,b,a));window.draw(fieldRect);
            }
        }
        {
            sf::VertexArray pts(sf::PrimitiveType::Triangles);
            float ps=std::max(1.5f,2.f*zoom);
            for(const auto& pl:plants){
                if(!on_screen(pl.x,pl.y))continue;
                float x=w2s_x(pl.x),y=w2s_y(pl.y);
                sf::Color col(0,200,0);
                // 小さい四角(三角形2つ)
                sf::Vertex a,b,c,d;
                a.position={x-ps,y-ps};a.color=col; b.position={x+ps,y-ps};b.color=col;
                c.position={x+ps,y+ps};c.color=col; d.position={x-ps,y+ps};d.color=col;
                pts.append(a);pts.append(b);pts.append(c);
                pts.append(a);pts.append(c);pts.append(d);
            }
            window.draw(pts);
        }
        mc.setFillColor(sf::Color(230,200,40));
        {
            sf::VertexArray mts(sf::PrimitiveType::Triangles);
            for(const auto& m:meats){
                if(!on_screen(m.x,m.y))continue;
                float r=std::max(1.f,m.energy*MEAT_RADIUS_SCALE*zoom);
                float x=w2s_x(m.x),y=w2s_y(m.y);
                sf::Color col(230,200,40);
                sf::Vertex a,b,c,d;
                a.position={x-r,y-r};a.color=col; b.position={x+r,y-r};b.color=col;
                c.position={x+r,y+r};c.color=col; d.position={x-r,y+r};d.color=col;
                mts.append(a);mts.append(b);mts.append(c);
                mts.append(a);mts.append(c);mts.append(d);
            }
            window.draw(mts);
        }
        // prey を VertexArray でバッチ描画(四角)
        {
            sf::VertexArray pts(sf::PrimitiveType::Triangles);
            for(const auto& p:preys){
                if(!on_screen(p.x,p.y))continue;
                float r=std::max(1.f,p.genes.size*zoom);
                float x=w2s_x(p.x),y=w2s_y(p.y);
                sf::Color col=cluster_color(p.cluster);
                sf::Vertex a,b,c,d;
                a.position={x-r,y-r};a.color=col; b.position={x+r,y-r};b.color=col;
                c.position={x+r,y+r};c.color=col; d.position={x-r,y+r};d.color=col;
                pts.append(a);pts.append(b);pts.append(c);
                pts.append(a);pts.append(c);pts.append(d);
            }
            window.draw(pts);
        }
        // predator を VertexArray でバッチ描画(四角)
        {
            sf::VertexArray pts(sf::PrimitiveType::Triangles);
            for(const auto& pd:predators){
                if(!on_screen(pd.x,pd.y))continue;
                float r=std::max(1.f,3.f*zoom);
                float x=w2s_x(pd.x),y=w2s_y(pd.y);
                sf::Color col(255,50,50);
                sf::Vertex a,b,c,d;
                a.position={x-r,y-r};a.color=col; b.position={x+r,y-r};b.color=col;
                c.position={x+r,y+r};c.color=col; d.position={x-r,y+r};d.color=col;
                pts.append(a);pts.append(b);pts.append(c);
                pts.append(a);pts.append(c);pts.append(d);
            }
            window.draw(pts);
        }
        // ===== NN可視化(選択個体) =====
        if(selected_id>=0){
            const Agent* sel=nullptr;
            for(const auto& p:preys)if(p.id==selected_id){sel=&p;break;}
            if(!sel)for(const auto& pd:predators)if(pd.id==selected_id){sel=&pd;break;}
            if(sel){
                // --- 選択個体の視覚を再計算 ---
                int NIN=sel->is_prey?PREY_IN:PRED_IN;
                std::vector<float> in(NIN,0.f);
                float facing=std::atan2(sel->vy,sel->vx);
                if(sel->is_prey){
                    float vp[NUM_RAYS]={0},vpr[NUM_RAYS]={0}; float pv=sel->genes.vision;
                    int cr=(int)std::ceil(pv/CELL_SIZE); int scx=(int)sel->x/CELL_SIZE,scy=(int)sel->y/CELL_SIZE;
                    for(int dx=-cr;dx<=cr;dx++)for(int dy=-cr;dy<=cr;dy++){
                        int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                        for(int pi:plant_grid[cell_index(cx,cy)]){
                            float ddx=torus_delta(sel->x,plants[pi].x,WIDTH),ddy=torus_delta(sel->y,plants[pi].y,HEIGHT);
                            float d2=ddx*ddx+ddy*ddy;if(d2>pv*pv||d2<1)continue;float d=std::sqrt(d2);
                            float rel=std::fmod(std::atan2(ddy,ddx)-facing+TWO_PI*2,TWO_PI);int ri=(int)(rel/TWO_PI*NUM_RAYS+0.5f)%NUM_RAYS;
                            float val=1.f-d/pv;if(val>vp[ri])vp[ri]=val;
                        }
                        for(int pi:pred_grid[cell_index(cx,cy)]){
                            if(!predators[pi].alive)continue;
                            float ddx=torus_delta(sel->x,predators[pi].x,WIDTH),ddy=torus_delta(sel->y,predators[pi].y,HEIGHT);
                            float d2=ddx*ddx+ddy*ddy;if(d2>pv*pv||d2<1)continue;float d=std::sqrt(d2);
                            float rel=std::fmod(std::atan2(ddy,ddx)-facing+TWO_PI*2,TWO_PI);int ri=(int)(rel/TWO_PI*NUM_RAYS+0.5f)%NUM_RAYS;
                            float val=1.f-d/pv;if(val>vpr[ri])vpr[ri]=val;
                        }
                    }
                    for(int i=0;i<NUM_RAYS;i++){in[i]=vp[i];in[NUM_RAYS+i]=vpr[i];} in[24]=sel->energy/sel->max_energy;
                } else {
                    float vpr[NUM_RAYS]={0};
                    int scx=(int)sel->x/CELL_SIZE,scy=(int)sel->y/CELL_SIZE;
                    for(int dx=-pred_cr;dx<=pred_cr;dx++)for(int dy=-pred_cr;dy<=pred_cr;dy++){
                        int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                        for(int pi:prey_grid[cell_index(cx,cy)]){
                            if(!preys[pi].alive)continue;
                            float ddx=torus_delta(sel->x,preys[pi].x,WIDTH),ddy=torus_delta(sel->y,preys[pi].y,HEIGHT);
                            float d2=ddx*ddx+ddy*ddy;if(d2>PREDATOR_VISION*PREDATOR_VISION||d2<1)continue;
                            float rel=std::fmod(std::atan2(ddy,ddx)-facing+PI*3,TWO_PI)-PI;
                            if(rel<-PREDATOR_FOV/2||rel>PREDATOR_FOV/2)continue;float d=std::sqrt(d2);
                            int ri=(int)((rel+PREDATOR_FOV/2)/PREDATOR_FOV*(NUM_RAYS-1)+0.5f);if(ri<0)ri=0;if(ri>=NUM_RAYS)ri=NUM_RAYS-1;
                            float val=1.f-d/PREDATOR_VISION;if(val>vpr[ri])vpr[ri]=val;
                        }
                    }
                    for(int i=0;i<NUM_RAYS;i++)in[i]=vpr[i]; in[12]=sel->energy/sel->max_energy;
                }
                // ===== 選択個体の視線(レイ)=====
                {
                    float pv = sel->is_prey ? sel->genes.vision : PREDATOR_VISION;
                    float sx=w2s_x(sel->x), sy=w2s_y(sel->y);
                    float rFood[NUM_RAYS]={0}, rEnemy[NUM_RAYS]={0}, rAlly[NUM_RAYS]={0}, rMeat[NUM_RAYS]={0};
                    int scx=(int)sel->x/CELL_SIZE, scy=(int)sel->y/CELL_SIZE;
                    int cr=(int)std::ceil(pv/CELL_SIZE);
                    auto rayIndex=[&](float ddx,float ddy)->int{
                        if(sel->is_prey){
                            float rel=std::fmod(std::atan2(ddy,ddx)-facing+TWO_PI*2,TWO_PI);
                            return (int)(rel/TWO_PI*NUM_RAYS+0.5f)%NUM_RAYS;
                        }else{
                            float rel=std::fmod(std::atan2(ddy,ddx)-facing+PI*3,TWO_PI)-PI;
                            if(rel<-PREDATOR_FOV/2||rel>PREDATOR_FOV/2)return -1;
                            int ri=(int)((rel+PREDATOR_FOV/2)/PREDATOR_FOV*(NUM_RAYS-1)+0.5f);
                            if(ri<0)ri=0; if(ri>=NUM_RAYS)ri=NUM_RAYS-1; return ri;
                        }
                    };
                    for(int dx=-cr;dx<=cr;dx++)for(int dy=-cr;dy<=cr;dy++){
                        int cx=scx+dx,cy=scy+dy; if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                        int ci=cell_index(cx,cy);
                        if(sel->is_prey){
                            for(int pi:plant_grid[ci]){
                                float ddx=torus_delta(sel->x,plants[pi].x,WIDTH),ddy=torus_delta(sel->y,plants[pi].y,HEIGHT);
                                float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                                int ri=rayIndex(ddx,ddy); if(ri<0)continue; float val=1.f-std::sqrt(d2)/pv; if(val>rFood[ri])rFood[ri]=val;
                            }
                            for(int pi:pred_grid[ci]){ if(!predators[pi].alive)continue;
                                float ddx=torus_delta(sel->x,predators[pi].x,WIDTH),ddy=torus_delta(sel->y,predators[pi].y,HEIGHT);
                                float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                                int ri=rayIndex(ddx,ddy); if(ri<0)continue; float val=1.f-std::sqrt(d2)/pv; if(val>rEnemy[ri])rEnemy[ri]=val;
                            }
                            for(int pi:prey_grid[ci]){ if(!preys[pi].alive||preys[pi].id==sel->id)continue;
                                float ddx=torus_delta(sel->x,preys[pi].x,WIDTH),ddy=torus_delta(sel->y,preys[pi].y,HEIGHT);
                                float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                                int ri=rayIndex(ddx,ddy); if(ri<0)continue; float val=1.f-std::sqrt(d2)/pv; if(val>rAlly[ri])rAlly[ri]=val;
                            }
                        } else {
                            for(int pi:prey_grid[ci]){ if(!preys[pi].alive)continue;
                                float ddx=torus_delta(sel->x,preys[pi].x,WIDTH),ddy=torus_delta(sel->y,preys[pi].y,HEIGHT);
                                float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                                int ri=rayIndex(ddx,ddy); if(ri<0)continue; float val=1.f-std::sqrt(d2)/pv; if(val>rFood[ri])rFood[ri]=val;
                            }
                            for(int pi:meat_grid[ci]){ if(meats[pi].energy<=0)continue;
                                float ddx=torus_delta(sel->x,meats[pi].x,WIDTH),ddy=torus_delta(sel->y,meats[pi].y,HEIGHT);
                                float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                                int ri=rayIndex(ddx,ddy); if(ri<0)continue; float val=1.f-std::sqrt(d2)/pv; if(val>rMeat[ri])rMeat[ri]=val;
                            }
                            for(int pi:pred_grid[ci]){ if(!predators[pi].alive||predators[pi].id==sel->id)continue;
                                float ddx=torus_delta(sel->x,predators[pi].x,WIDTH),ddy=torus_delta(sel->y,predators[pi].y,HEIGHT);
                                float d2=ddx*ddx+ddy*ddy; if(d2>pv*pv||d2<1)continue;
                                int ri=rayIndex(ddx,ddy); if(ri<0)continue; float val=1.f-std::sqrt(d2)/pv; if(val>rAlly[ri])rAlly[ri]=val;
                            }
                        }
                    }
                    for(int i=0;i<NUM_RAYS;i++){
                        float ang;
                        if(sel->is_prey) ang = facing + TWO_PI*((float)i/NUM_RAYS);
                        else ang = facing - PREDATOR_FOV/2 + PREDATOR_FOV*((float)i/(NUM_RAYS-1));
                        float mx=0.01f; sf::Color rc2(60,60,60,120);
                        if(rFood[i]>mx){ mx=rFood[i]; rc2=sf::Color(60,255,60,200); }
                        if(rEnemy[i]>mx){ mx=rEnemy[i]; rc2=sf::Color(255,60,60,200); }
                        if(rMeat[i]>mx){ mx=rMeat[i]; rc2=sf::Color(230,200,40,200); }
                        if(rAlly[i]>mx){ mx=rAlly[i]; rc2=sf::Color(60,160,255,200); }
                        float ex=sx+std::cos(ang)*pv*zoom, ey=sy+std::sin(ang)*pv*zoom;
                        sf::VertexArray line(sf::PrimitiveType::Lines,2);
                        line[0].position={sx,sy}; line[0].color=rc2;
                        line[1].position={ex,ey}; line[1].color=rc2;
                        window.draw(line);
                    }
                }

                // ===== ステータスパネル(NNパネルの下)=====
                if(font_ok){
                    float PX=20.f, PY=20.f, PW=SCREEN_W/4.f*0.8f, PH=SCREEN_H/2.f*0.8f;
                    float SX=PX, SY=PY+PH+10.f, SW=PW;
                    // バー描画ヘルパー
                    auto drawBar=[&](float bx,float by,float bw,float val,sf::Color col,const std::string& label){
                        sf::RectangleShape bg({bw,14.f}); bg.setPosition({bx,by}); bg.setFillColor(sf::Color(50,50,50)); window.draw(bg);
                        sf::RectangleShape fg({bw*clamp01(val),14.f}); fg.setPosition({bx,by}); fg.setFillColor(col); window.draw(fg);
                        sf::Text t(font,label,12); t.setPosition({bx+bw+6.f,by}); t.setFillColor(sf::Color::White); window.draw(t);
                    };
                    float bx=SX+60.f, bw=SW-120.f, by=SY, lh=20.f;
                    // ラベル(左)
                    auto lab=[&](float y,const std::string& s){ sf::Text t(font,s,12); t.setPosition({SX,y}); t.setFillColor(sf::Color::White); window.draw(t); };
                    char bb[64];
                    lab(by,"Energy"); snprintf(bb,64,"%.1f / %.1f",sel->energy,sel->max_energy); drawBar(bx,by,bw,sel->energy/sel->max_energy,sf::Color(0,200,80),bb); by+=lh;
                    lab(by,"Stress"); snprintf(bb,64,"%.2f / 1.00",sel->stress); drawBar(bx,by,bw,sel->stress,sf::Color(220,60,60),bb); by+=lh;
                    lab(by,"Fear");   snprintf(bb,64,"%.2f / 1.00",sel->fear); drawBar(bx,by,bw,sel->fear,sf::Color(255,120,0),bb); by+=lh;
                    lab(by,"Affin");  snprintf(bb,64,"%.2f / 1.00",sel->affinity); drawBar(bx,by,bw,sel->affinity,sf::Color(0,120,255),bb); by+=lh;
                    // 数値情報
                    char buf[256];
                    snprintf(buf,sizeof(buf),"Age: %.1fs", sel->age/60.f);
                    { sf::Text t(font,buf,13); t.setPosition({SX,by}); t.setFillColor(sf::Color::White); window.draw(t); by+=lh; }
                    if(sel->is_prey){
                        sf::RectangleShape sq({12.f,12.f}); sq.setPosition({SX,by+1.f}); sq.setFillColor(cluster_color(sel->cluster)); window.draw(sq);
                        char sb2[64]; snprintf(sb2,64,"    Species #%d", sel->cluster);
                        sf::Text t(font,sb2,13); t.setPosition({SX,by}); t.setFillColor(sf::Color::White); window.draw(t); by+=lh;
                    }
                    if(sel->is_prey){
                        snprintf(buf,sizeof(buf),"Size x%.2f  Vision x%.2f", sel->genes.size/BASE_SIZE, sel->genes.vision/PREY_VISION);
                        { sf::Text t(font,buf,13); t.setPosition({SX,by}); t.setFillColor(sf::Color::White); window.draw(t); by+=lh; }
                        snprintf(buf,sizeof(buf),"Eat x%.2f   Speed x%.2f", sel->genes.eat_gain/EAT_GAIN, sel->genes.speed/SPEED);
                        { sf::Text t(font,buf,13); t.setPosition({SX,by}); t.setFillColor(sf::Color::White); window.draw(t); by+=lh; }
                        snprintf(buf,sizeof(buf),"Repro: %d   Food: %d", sel->repro_count, sel->food_count);
                    } else {
                        snprintf(buf,sizeof(buf),"Repro: %d   Kill: %d", sel->repro_count, sel->food_count);
                    }
                    { sf::Text t(font,buf,13); t.setPosition({SX,by}); t.setFillColor(sf::Color::White); window.draw(t); by+=lh; }
                }
                // 選択個体の黄色い輪(シミュレーション画面上)
                {
                    float rr=(sel->is_prey?sel->genes.size:3.f)*zoom;
                    float sx=w2s_x(sel->x),sy=w2s_y(sel->y);
                    sf::CircleShape ring(rr+4.f);
                    ring.setFillColor(sf::Color::Transparent);
                    ring.setOutlineColor(sf::Color(255,255,0));
                    ring.setOutlineThickness(2.f);
                    ring.setPosition({sx-(rr+4.f),sy-(rr+4.f)});
                    window.draw(ring);
                }
                // ===== NEAT可視化: Genomeの脳を描く =====
                {
                    const Genome& gen = sel->genome;
                    // forward して各ノードの活性値を得る
                    std::vector<float> sv(37, 0.f);
                    for(int i=0;i<(int)in.size() && i<37;i++) sv[i]=in[i];
                    // 各ノードの値を計算(forwardと同じロジックで node_val を得る)
                    int max_id=0; for(const auto& n:gen.nodes) if(n.id>max_id) max_id=n.id;
                    std::vector<float> nodeval(max_id+1, 0.f);
                    int n_in=(int)gen.sensors.size();
                    for(int i=0;i<n_in;i++) nodeval[i]=sv[i];
                    // 簡易順伝播(計算順序: forwardと同じトポロジカル順)
                    {
                        std::vector<int> deg(max_id+1,0);
                        for(const auto& c:gen.conns) if(c.enabled) deg[c.out_node]++;
                        std::vector<int> q; for(int i=0;i<n_in;i++) q.push_back(i);
                        size_t qi=0; std::vector<int> order;
                        for(int i=0;i<n_in;i++) order.push_back(i);
                        while(qi<q.size()){
                            int nid=q[qi++];
                            for(const auto& c:gen.conns){ if(!c.enabled)continue; if(c.in_node==nid){ deg[c.out_node]--; if(deg[c.out_node]==0){q.push_back(c.out_node);order.push_back(c.out_node);} } }
                        }
                        for(int nid:order){ if(nid<n_in)continue; float s=0.f; for(const auto& c:gen.conns){if(!c.enabled)continue; if(c.out_node==nid)s+=nodeval[c.in_node]*c.weight;} nodeval[nid]=std::tanh(s); }
                    }

                    // レイアウト
                    float PX=20.f, PY=20.f, PW=SCREEN_W/4.f*0.8f, PH=SCREEN_H/2.f*0.8f;
                    sf::RectangleShape panel({PW,PH}); panel.setPosition({PX,PY});
                    panel.setFillColor(sf::Color(0,0,0,180)); panel.setOutlineColor(sf::Color(80,80,80)); panel.setOutlineThickness(1.f);
                    window.draw(panel);

                    int n_out=(int)gen.actuators.size();
                    // 各ノードの画面座標を決める(入力=左列, 出力=右列, 中間=中央にばらまく)
                    std::vector<sf::Vector2f> pos(max_id+1, {0,0});
                    // 入力ノード(id 0..n_in-1): 左列に縦に並べる
                    for(int i=0;i<n_in;i++){
                        float t=(n_in==1)?0.5f:(float)i/(n_in-1);
                        pos[i]={PX+PW*0.10f, PY+20.f+(PH-40.f)*t};
                    }
                    // 出力ノード(id n_in..n_in+n_out-1): 右列
                    for(int j=0;j<n_out;j++){
                        float center=0.5f, spread=0.15f;
                        float t=(n_out==1)?center : center + spread*(2.f*j/(n_out-1)-1.f);
                        pos[n_in+j]={PX+PW*0.90f, PY+20.f+(PH-40.f)*t};
                    }
                    // 中間ノード(type==2): 中央にばらまく(idベースで散らす)
                    int hidden_count=0; for(const auto& n:gen.nodes) if(n.type==2) hidden_count++;
                    int hi=0;
                    for(const auto& n:gen.nodes){
                        if(n.type==2){
                            float t=(hidden_count==1)?0.5f:(float)hi/(hidden_count-1);
                            float hx=PX+PW*0.35f + (PW*0.30f)*((hi%2==0)?0.f:1.f);  // 2列に散らす
                            pos[n.id]={hx, PY+30.f+(PH-60.f)*t};
                            hi++;
                        }
                    }

                    // 接続(線)を描く
                    sf::VertexArray va(sf::PrimitiveType::Lines);
                    for(const auto& c:gen.conns){
                        if(!c.enabled)continue;
                        float aw=std::fabs(c.weight); if(aw<0.15f)continue;
                        int alpha=(int)std::min(150.f,aw*120.f);
                        sf::Color col = c.weight>0 ? sf::Color(80,160,255,alpha) : sf::Color(255,80,80,alpha);
                        sf::Vertex v1,v2;
                        v1.position=pos[c.in_node]; v1.color=col;
                        v2.position=pos[c.out_node]; v2.color=col;
                        va.append(v1); va.append(v2);
                    }
                    window.draw(va);

                    // ノードを描く(活性値で色分け)
                    for(const auto& n:gen.nodes){
                        float a=nodeval[n.id];
                        int inten=(int)std::min(255.f,std::fabs(a)*255.f);
                        sf::Color col = a>=0 ? sf::Color(40,80+inten*0.7f,255) : sf::Color(255,60,60);
                        float rad = (n.type==2)?4.f:5.f;   // 中間は少し小さく
                        sf::CircleShape node(rad); node.setFillColor(col);
                        node.setPosition({pos[n.id].x-rad, pos[n.id].y-rad});
                        window.draw(node);
                    }

                    // ノード数・接続数を表示
                    if(font_ok){
                        char nb[64]; int active_conns=0; for(const auto&c:gen.conns)if(c.enabled)active_conns++;
                        snprintf(nb,64,"Nodes: %d  Conns: %d  Hidden: %d",(int)gen.nodes.size(),active_conns,hidden_count);
                        sf::Text t(font,nb,12); t.setPosition({PX+6.f,PY+PH-18.f}); t.setFillColor(sf::Color(200,200,200));
                        window.draw(t);
                    }
                }
            }
        }
        // ===== 個体数グラフ(右上に4つ)=====
        if(font_ok){
            float GW=SCREEN_W/4.f*0.8f;
            float GH=200.f;
            float GX=SCREEN_W-GW-20.f;
            float GY0=20.f;
            float PADT=26.f, PADB=6.f, PADL=6.f, PADR=44.f;  // 内側の余白(上=ラベル,右=目盛り数字)
            struct GraphInfo { const std::vector<int>* data; const char* label; sf::Color col; int peak; };
            GraphInfo graphs[4]={
                {&hist_prey,"Prey",sf::Color(0,120,255),peak_prey},
                {&hist_pred,"Predator",sf::Color(255,60,60),peak_pred},
                {&hist_plant,"Plant",sf::Color(0,200,0),peak_plant},
                {&hist_meat,"Meat",sf::Color(230,200,40),peak_meat},
            };
            auto niceMax=[](int mx)->int{
                if(mx<=0)return 1;
                float e=std::floor(std::log10((float)mx));
                float base=std::pow(10.f,e);
                float f=mx/base;
                float nf; if(f<=1)nf=1; else if(f<=2)nf=2; else if(f<=5)nf=5; else nf=10;
                return (int)(nf*base);
            };
            for(int g=0;g<4;g++){
                float GY=GY0+g*(GH+30.f);
                const std::vector<int>& d=*graphs[g].data;
                sf::Color col=graphs[g].col;
                // 背景 + 色付き枠
                sf::RectangleShape bg({GW,GH}); bg.setPosition({GX,GY});
                bg.setFillColor(sf::Color(0,0,0,180));
                bg.setOutlineColor(col); bg.setOutlineThickness(2.f);   // 枠を各色に
                window.draw(bg);
                // プロット領域(枠の内側)
                float plotX=GX+PADL, plotY=GY+PADT;
                float plotW=GW-PADL-PADR, plotH=GH-PADT-PADB;
                int rawmx=1; for(int v:d)if(v>rawmx)rawmx=v;
                int mx=niceMax(rawmx);
                // --- 縦軸目盛り ---
                int ticks[4]={0, mx/3, mx*2/3, mx};
                for(int ti=0;ti<4;ti++){
                    float ratio=ticks[ti]/(float)mx;
                    float ly=plotY+plotH-plotH*ratio;
                    sf::VertexArray gl(sf::PrimitiveType::Lines,2);
                    gl[0].position={plotX,ly}; gl[0].color=sf::Color(60,60,60);
                    gl[1].position={plotX+plotW,ly}; gl[1].color=sf::Color(60,60,60);
                    window.draw(gl);
                    char nb[32]; snprintf(nb,32,"%d",ticks[ti]);
                    sf::Text nt(font,nb,10); nt.setFillColor(sf::Color(150,150,150));
                    nt.setPosition({plotX+plotW+3.f, ly-6.f});
                    window.draw(nt);
                }
                // --- 塗りつぶし(上端=濃い、0=薄い。高さ基準)---
                if(d.size()>=2){
                    sf::VertexArray fill(sf::PrimitiveType::Triangles);
                    float baseY=plotY+plotH;   // 0の位置(下端)
                    float topY=plotY;          // mx(一番上のキリ数字)の位置
                    // 高さ→α を返す(topで濃い、baseで薄い)
                    auto alphaAt=[&](float y)->int{
                        float t=(baseY-y)/(baseY-topY);  // 0(下)〜1(上)
                        if(t<0)t=0; if(t>1)t=1;
                        return (int)(15 + t*(150-15));   // 15〜150
                    };
                    for(int i=0;i+1<(int)d.size();i++){
                        float px0=plotX+plotW*i/(float)(HIST_MAX-1);
                        float px1=plotX+plotW*(i+1)/(float)(HIST_MAX-1);
                        float py0=plotY+plotH-plotH*d[i]/(float)mx;
                        float py1=plotY+plotH-plotH*d[i+1]/(float)mx;
                        sf::Color tl_c=col; tl_c.a=alphaAt(py0);
                        sf::Color tr_c=col; tr_c.a=alphaAt(py1);
                        sf::Color b_c=col;  b_c.a=alphaAt(baseY);   // 下端=一番薄い
                        sf::Vertex tl,tr,br,bl;
                        tl.position={px0,py0}; tl.color=tl_c;
                        tr.position={px1,py1}; tr.color=tr_c;
                        br.position={px1,baseY}; br.color=b_c;
                        bl.position={px0,baseY}; bl.color=b_c;
                        fill.append(tl); fill.append(tr); fill.append(br);
                        fill.append(tl); fill.append(br); fill.append(bl);
                    }
                    window.draw(fill);
                }
                // --- 折れ線(プロット領域内)---
                if(d.size()>=2){
                    sf::VertexArray line(sf::PrimitiveType::LineStrip);
                    for(int i=0;i<(int)d.size();i++){
                        float px=plotX+plotW*i/(float)(HIST_MAX-1);
                        float py=plotY+plotH-plotH*d[i]/(float)mx;
                        sf::Vertex v; v.position={px,py}; v.color=col; line.append(v);
                    }
                    window.draw(line);
                }
                // --- ラベル + 現在値 + max ---
                char buf[80]; int cur=d.empty()?0:d.back();
                snprintf(buf,80,"%s: %d  (max %d)",graphs[g].label,cur,graphs[g].peak);
                sf::Text t(font,buf,13); t.setPosition({GX+6.f,GY+5.f}); t.setFillColor(col);
                window.draw(t);
                // --- 横軸(時間、4分割)---
                float total_sec=HIST_MAX/measured_fps;
                for(int xi=0;xi<4;xi++){
                    float ratio=xi/3.f;
                    float px=plotX+plotW*ratio;
                    float sec_ago=total_sec*(1.f-ratio);
                    char tb[32]; snprintf(tb,32,"%.0fs",sec_ago);
                    sf::Text tt(font,tb,10); tt.setFillColor(sf::Color(150,150,150));
                    tt.setPosition({px-((xi==3)?16.f:0.f),GY+GH+3.f});
                    window.draw(tt);
                }
            }
        }
        if(font_ok){
            float ms_per_frame=1000.f/measured_fps;
            float elapsed_sec=frame/measured_fps;
            // 上段:世代・tick・時間・ms
            char sb[160];
            snprintf(sb,160,"Gen: %d    Tick: %d    Time: %.1fs    %.1f ms/frame",generation,frame,elapsed_sec,ms_per_frame);
            sf::Text st(font,sb,20); st.setFillColor(sf::Color::White);
            sf::FloatRect b1=st.getLocalBounds();
            float tw=b1.size.x, th=b1.size.y;
            // 下段:最大fit(小さく)
            int best_prey_fit=0; for(auto&r:prey_hall)if(r.fit>best_prey_fit)best_prey_fit=r.fit;
            int best_pred_fit=0; for(auto&p:pred_hall)if(p.first>best_pred_fit)best_pred_fit=p.first;

            char fb[128];
            snprintf(fb,128,"Best Fit  -  Prey: %d   Predator: %d",best_prey_fit,best_pred_fit);
            sf::Text ft(font,fb,14); ft.setFillColor(sf::Color(180,180,180));
            sf::FloatRect b2=ft.getLocalBounds();
            float fw=b2.size.x;
            // 枠(2段を囲む)
            float cx=SCREEN_W/2.f;
            float boxW=std::max(tw,fw)+24.f;
            float boxH=th+b2.size.y+30.f;
            float boxY=SCREEN_H-boxH-14.f;
            sf::RectangleShape box({boxW,boxH});
            box.setPosition({cx-boxW/2.f,boxY});
            box.setFillColor(sf::Color(0,0,0,180)); box.setOutlineColor(sf::Color(150,150,150)); box.setOutlineThickness(2.f);
            window.draw(box);
            // 上段テキスト
            st.setPosition({cx-tw/2.f, boxY+8.f-b1.position.y});
            window.draw(st);
            // 下段テキスト
            ft.setPosition({cx-fw/2.f, boxY+8.f+th+6.f-b2.position.y});
            window.draw(ft);
        }
        
        // 脳の複雑さ(平均) prey vs predator
        if(font_ok){
            double pH=0,pC=0,dH=0,dC=0; int npy=0,npd=0;
            int pHmax=0,pCmax=0,dHmax=0,dCmax=0;
            for(const auto& p:preys){ if(!p.alive)continue;
                int h=(int)p.genome.nodes.size()-p.genome.cached_n_in-(int)p.genome.actuators.size();
                int c=p.genome.cached_active_conns;
                pH+=h; pC+=c; if(h>pHmax)pHmax=h; if(c>pCmax)pCmax=c; npy++; }
            for(const auto& pd:predators){ if(!pd.alive)continue;
                int h=(int)pd.genome.nodes.size()-pd.genome.cached_n_in-(int)pd.genome.actuators.size();
                int c=pd.genome.cached_active_conns;
                dH+=h; dC+=c; if(h>dHmax)dHmax=h; if(c>dCmax)dCmax=c; npd++; }
            char bs[240];
            snprintf(bs,240,"Brain | Prey: hidden %.1f (max %d), conns %.0f (max %d)   Predator: hidden %.1f (max %d), conns %.0f (max %d)",
                     npy?pH/npy:0.0, pHmax, npy?pC/npy:0.0, pCmax, npd?dH/npd:0.0, dHmax, npd?dC/npd:0.0, dCmax);
            sf::Text bt(font,bs,14); bt.setFillColor(sf::Color(180,220,255));
            bt.setPosition({12.f, SCREEN_H-46.f});
            window.draw(bt);
        }

        if(font_ok){
            long ptot=prey_starve+prey_killed;
            char ds[220];
            snprintf(ds,220,"Prey deaths: starved %.0f%% / eaten %.0f%% (total %ld)   Pred starved: %ld",
                     ptot?100.0*prey_starve/ptot:0.0, ptot?100.0*prey_killed/ptot:0.0, ptot, pred_starve);
            sf::Text dt(font,ds,14); dt.setFillColor(sf::Color(255,210,150));
            dt.setPosition({12.f, SCREEN_H-66.f});
            window.draw(dt);
        }
        
        if(font_ok){
            int alive_lin=0; for(auto&l:lineages) if(l.alive)alive_lin++;
            char cs[160];
            snprintf(cs,160,"Species alive: %d   lineages ever: %d   threshold %.2f", alive_lin, (int)lineages.size(), species_threshold);
            sf::Text ct(font,cs,14); ct.setFillColor(sf::Color(200,255,200));
            ct.setPosition({12.f, SCREEN_H-86.f});
            window.draw(ct);
        }

        if(font_ok && prof_str[0]!='\0'){
            sf::Text pt(font,prof_str,14); pt.setFillColor(sf::Color::White);
            pt.setPosition({12.f, SCREEN_H-26.f});
            window.draw(pt);
        }
        if(font_ok && ui_mode!=0){
            int nf=(int)save_files.size();
            float bw=700.f;
            float bh=(ui_mode==1)?70.f : std::min(SCREEN_H*0.9f, 58.f+std::max(1,nf)*18.f);
            float bx=SCREEN_W/2.f-bw/2.f, by=SCREEN_H/2.f-bh/2.f;
            sf::RectangleShape box({bw,bh}); box.setPosition({bx,by});
            box.setFillColor(sf::Color(0,0,0,225)); box.setOutlineColor(sf::Color(200,200,200)); box.setOutlineThickness(2.f);
            window.draw(box);
            char pr[160];
            snprintf(pr,160,"%s  %s_", ui_mode==1?"Save as:":"Load (type name, Enter):", input_text.c_str());
            sf::Text t(font,pr,18); t.setFillColor(sf::Color::White); t.setPosition({bx+14.f,by+12.f});
            window.draw(t);
            if(ui_mode==2){
                float ly=by+44.f;
                int maxshow=(int)((bh-50.f)/18.f);
                for(int i=0;i<nf && i<maxshow;i++){
                    sf::Text ft(font,save_files[i],14); ft.setFillColor(sf::Color(180,220,255));
                    ft.setPosition({bx+18.f,ly}); window.draw(ft); ly+=18.f;
                }
                if(nf>maxshow){
                    char more[48]; snprintf(more,48,"...and %d more",nf-maxshow);
                    sf::Text mt(font,more,12); mt.setFillColor(sf::Color(150,150,150));
                    mt.setPosition({bx+18.f,ly}); window.draw(mt);
                }
            }
        }
        // ===== 系統樹(Tキーで表示。上=過去, 下=現在。画面いっぱい)=====
        if(show_tree && font_ok){
            int N=(int)lineages.size();
            float panelX=SCREEN_W*0.02f, panelY=SCREEN_H*0.03f;
            float panelW=SCREEN_W*0.96f, panelH=SCREEN_H*0.94f;
            sf::RectangleShape bg({panelW,panelH}); bg.setPosition({panelX,panelY});
            bg.setFillColor(sf::Color(0,0,0,235)); bg.setOutlineColor(sf::Color(150,150,150)); bg.setOutlineThickness(2.f);
            window.draw(bg);
            float inX=panelX+30.f, inY=panelY+55.f, inW=panelW-60.f, inH=panelH-95.f;
            std::vector<float> sxpos(N,0.f), sy0(N,0.f), sy1(N,0.f);
            if(N>0){
                std::vector<float> lx(N,-1.f); int leafc=0;
                std::function<float(int)> assign=[&](int i)->float{
                    std::vector<int> ch; for(int j=0;j<N;j++) if(lineages[j].parent==lineages[i].id) ch.push_back(j);
                    if(ch.empty()){ lx[i]=(float)(leafc++); return lx[i]; }
                    float s=0; for(int c:ch) s+=assign(c); lx[i]=s/ch.size(); return lx[i];
                };
                for(int i=0;i<N;i++) if(lineages[i].parent<0) assign(i);
                float xspan=(leafc>1)?(float)(leafc-1):1.f;
                float T1=(float)std::max(1,frame);
                auto sx=[&](float slot){ return inX + (leafc>1 ? slot/xspan : 0.5f)*inW; };
                auto ty=[&](int tick){ return inY + ((float)tick/T1)*inH; };
                sf::VertexArray va(sf::PrimitiveType::Lines);
                for(int i=0;i<N;i++){
                    float x=sx(lx[i]);
                    int bt=lineages[i].birth_tick;
                    int dt=lineages[i].alive?frame:lineages[i].death_tick;
                    sxpos[i]=x; sy0[i]=ty(bt); sy1[i]=ty(dt);
                    sf::Color col=cluster_color(lineages[i].id);
                    sf::Vertex a,b; a.position={x,ty(bt)}; a.color=col; b.position={x,ty(dt)}; b.color=col;
                    va.append(a); va.append(b);
                    if(lineages[i].parent>=0){
                        int pj=-1; for(int j=0;j<N;j++) if(lineages[j].id==lineages[i].parent){pj=j;break;}
                        if(pj>=0){
                            float px=sx(lx[pj]); sf::Color pc=col; pc.a=150;
                            sf::Vertex c,d; c.position={px,ty(bt)}; c.color=pc; d.position={x,ty(bt)}; d.color=pc;
                            va.append(c); va.append(d);
                        }
                    }
                }
                window.draw(va);
                sf::VertexArray xm(sf::PrimitiveType::Lines);
                for(int i=0;i<N;i++){
                    if(lineages[i].alive) continue;
                    float x=sxpos[i], y=sy1[i], r=5.f; sf::Color col=cluster_color(lineages[i].id);
                    sf::Vertex a,b,c,d;
                    a.position={x-r,y-r}; a.color=col; b.position={x+r,y+r}; b.color=col;
                    c.position={x-r,y+r}; c.color=col; d.position={x+r,y-r}; d.color=col;
                    xm.append(a);xm.append(b);xm.append(c);xm.append(d);
                }
                window.draw(xm);
            }
            int alive_lin=0; for(auto&l:lineages) if(l.alive)alive_lin++;
            char tt[128]; snprintf(tt,128,"Phylogeny   alive %d / ever %d   (T close)", alive_lin, N);
            sf::Text ttl(font,tt,18); ttl.setFillColor(sf::Color::White); ttl.setPosition({panelX+16.f,panelY+14.f}); window.draw(ttl);
            sf::Text pastT(font,"past",13); pastT.setFillColor(sf::Color(150,150,150)); pastT.setPosition({panelX+6.f,inY-4.f}); window.draw(pastT);
            sf::Text nowT(font,"now",13); nowT.setFillColor(sf::Color(150,150,150)); nowT.setPosition({panelX+6.f,inY+inH-4.f}); window.draw(nowT);
            sf::Vector2i mp=sf::Mouse::getPosition(window);
            int hit=-1; float bestdx=8.f;
            for(int i=0;i<N;i++){
                if(mp.y>=sy0[i]-3 && mp.y<=sy1[i]+3){
                    float dx=std::fabs((float)mp.x - sxpos[i]);
                    if(dx<bestdx){ bestdx=dx; hit=i; }
                }
            }
            if(hit>=0){
                char hb[128];
                snprintf(hb,128,"species #%d   now %d   max %d%s", lineages[hit].id, lineages[hit].pop, lineages[hit].max_pop, lineages[hit].alive?"":"  (extinct)");
                sf::Text ht(font,hb,14); ht.setFillColor(sf::Color::White);
                sf::FloatRect hbnd=ht.getLocalBounds();
                float boxW=hbnd.size.x+14.f, boxH=24.f;
                float bx=(float)mp.x+12.f, by=(float)mp.y+12.f;
                if(bx+boxW > SCREEN_W-10.f) bx=(float)mp.x-12.f-boxW;   // 右端 → 左に出す
                if(by+boxH > SCREEN_H-10.f) by=(float)mp.y-12.f-boxH;   // 下端 → 上に出す
                if(bx<10.f) bx=10.f;
                sf::RectangleShape hbg({boxW,boxH}); hbg.setPosition({bx-4.f,by-2.f});
                hbg.setFillColor(sf::Color(20,20,20,235)); hbg.setOutlineColor(cluster_color(lineages[hit].id)); hbg.setOutlineThickness(1.5f);
                window.draw(hbg);
                ht.setPosition({bx+2.f,by}); window.draw(ht);
            }
        }
        window.display();
        t_draw += ms(_t5,now());
        if(!paused){
            // 個体数履歴を記録
            hist_prey.push_back((int)preys.size());
            hist_pred.push_back((int)predators.size());
            hist_plant.push_back((int)plants.size());
            hist_meat.push_back((int)meats.size());
            if((int)hist_prey.size()>HIST_MAX){
                hist_prey.erase(hist_prey.begin());
                hist_pred.erase(hist_pred.begin());
                hist_plant.erase(hist_plant.begin());
                hist_meat.erase(hist_meat.begin());
            }
            if((int)preys.size()>peak_prey)peak_prey=(int)preys.size();
            if((int)predators.size()>peak_pred)peak_pred=(int)predators.size();
            if((int)plants.size()>peak_plant)peak_plant=(int)plants.size();
            if((int)meats.size()>peak_meat)peak_meat=(int)meats.size();
            frame++;
            if(++autosave_timer >= AUTOSAVE_INTERVAL){
                autosave_timer=0;
                char afn[32]; snprintf(afn,32,"autosave%d.txt",autosave_slot);
                save_state(afn, preys, predators, frame);
                autosave_slot=(autosave_slot+1)%3;   // autosave0→1→2→0 と回す
            }
        }
        prof_frames++;
        if(prof_frames>=60){
            snprintf(prof_str,256,"preyNN=%.2f predNN=%.2f move=%.2f eat=%.2f field=%.2f draw=%.2f",
                     t_prey_nn/prof_frames,t_pred_nn/prof_frames,t_move/prof_frames,t_eat/prof_frames,t_field/prof_frames,t_draw/prof_frames);
            t_prey_nn=t_pred_nn=t_move=t_eat=t_field=t_draw=0; prof_frames=0;
        }
    }
    return 0;
}