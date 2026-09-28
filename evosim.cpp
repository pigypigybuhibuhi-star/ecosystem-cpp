#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <cstdio>
#include <random>
#include <algorithm>
#include <fstream>
#include <string>
#include <unordered_map>

const int GRID_W = 320;
const int GRID_H = 180;
const int SCALE  = 4;
const float PREDATE_EFF   = 0.20f;   // 一噛みで奪う割合
const float PREDATE_LOSS  = 0.30f;   // 捕食時のエネルギー損失(食物連鎖の非効率)
const float INIT_TROPHIC  = 0.01f;   // 最初は植物寄り
int next_cell_id = 0;

const float INIT_ADHESION  = 0.0f;   // 最初は接着しない
const float ADHESION_COST  = 0.003f; // 接着分子の産生コスト(毎フレーム)
const float ADHESION_MIN   = 0.15f;  // これ以下だと結合できない
const float SIG_TOLERANCE  = 0.10f;  // 署名の一致許容範囲
const float BOND_DIST      = 1.5f;   // 目標距離(グリッド単位)
const float BOND_BREAK     = 3.0f;   // これ以上離れたら切れる
const float BOND_STIFF     = 0.3f;   // 距離制約の強さ
const float SHARE_RATE = 0.15f;   // 結合相手とのエネルギー移動率

std::vector<float> solar(GRID_W * GRID_H, 0.f);

std::vector<float> energy(GRID_W * GRID_H, 0.f);   // 場に溜まっているエネルギー
const float INFLOW  = 0.02f;   // 注入係数
const float DISSIPATE = 0.01f; // 散逸率(毎フレーム逃げる割合)

inline int idx(int x, int y){ return y * GRID_W + x; }

const int MAX_BONDS = 6;

struct Cell {
    float x, y;
    float energy;
    bool alive;
    float absorb;
    float trophic;
    float dir;
    float speed;
    float turn;
    int   id;                  // 消えないID(配列が詰められても変わらない)
    float adhesion;            // 遺伝子: 接着分子の産生量
    float signature;           // 遺伝子: 接着署名(一致した相手とだけ結合)
    int   bonds[MAX_BONDS];    // 結合相手のID(-1 = 空)
    int   n_bonds;
};

std::vector<Cell> cells;
std::vector<std::vector<int>> cell_grid(GRID_W * GRID_H);
struct Corpse {
    float x, y;
    float energy;
};
std::vector<Corpse> corpses;
std::vector<std::vector<int>> corpse_grid(GRID_W * GRID_H);

// 履歴(グラフ用)
const int HIST_MAX = 6000;
const int SAMPLE_INTERVAL = 300;   // 5秒ごとに1点(60fps換算)
std::vector<float> h_trophic, h_absorb, h_speed;
std::vector<int>   h_pop, h_corpse;

const float SCAVENGE_EFF = 0.60f;   // 死骸から奪う割合(生体より高い=安く食える)
const float SCAVENGE_LOSS = 0.10f;  // 腐肉食の損失(捕食より小さい)
const float ROT_RATE = 0.001f;       // 腐敗して場に戻る割合
const float INIT_ABSORB      = 0.5f;
const float METABOLISM       = 0.01f;
const float INIT_ENERGY      = 1.0f;
const float DIVIDE_THRESHOLD = 4.0f;   // これを超えたら分裂
const float MUT_STRENGTH     = 0.02f;  // 変異の大きさ
const float MOVE_COST   = 0.004f;  // 速度あたりの移動コスト
const float TURN_COST   = 0.001f;  // 旋回コスト
const float MAX_SPEED   = 0.5f;    // 最大速度(グリッド/フレーム)
const float INIT_SPEED  = 0.0f;    // 最初は動かない

std::mt19937 rng(12345);
std::uniform_real_distribution<float> dist01(0.f, 1.f);

inline float frand(float a, float b){ return a + (b - a) * dist01(rng); }

inline float mutate(float v){
    float nv = v + frand(-MUT_STRENGTH, MUT_STRENGTH);
    if(nv < 0.01f) nv = 0.01f;
    if(nv > 1.00f) nv = 1.00f;
    return nv;
}

void init_solar(){
    for(int y = 0; y < GRID_H; y++){
        float lat = std::fabs((float)y - (GRID_H - 1) / 2.f) / ((GRID_H - 1) / 2.f);
        float v = std::cos(lat * 3.14159265f / 2.f);
        for(int x = 0; x < GRID_W; x++) solar[idx(x,y)] = v;
    }
}

sf::Color value_to_color(float v){
    if(v < 0.f) v = 0.f;
    if(v > 1.f) v = 1.f;
    if(v < 0.5f){
        float t = v / 0.5f;
        return sf::Color((std::uint8_t)(10 + 20*t), (std::uint8_t)(20 + 120*t), (std::uint8_t)(60 + 40*t));
    } else {
        float t = (v - 0.5f) / 0.5f;
        return sf::Color((std::uint8_t)(30 + 225*t), (std::uint8_t)(140 + 100*t), (std::uint8_t)(100 - 60*t));
    }
}

void update_field(){
    for(int i = 0; i < GRID_W * GRID_H; i++){
        energy[i] += solar[i] * INFLOW;
        energy[i] -= energy[i] * DISSIPATE;
        if(energy[i] < 0.f) energy[i] = 0.f;
    }
}


void update_cells(){
    // 空間グリッドを作り直す
    for(auto& v : cell_grid) v.clear();
    for(int i = 0; i < (int)cells.size(); i++){
        if(!cells[i].alive) continue;
        int gx = (int)cells[i].x, gy = (int)cells[i].y;
        if(gx < 0 || gx >= GRID_W || gy < 0 || gy >= GRID_H) continue;
        cell_grid[idx(gx,gy)].push_back(i);
    }
    for(auto& v : corpse_grid) v.clear();
    for(int i = 0; i < (int)corpses.size(); i++){
        int gx = (int)corpses[i].x, gy = (int)corpses[i].y;
        if(gx < 0 || gx >= GRID_W || gy < 0 || gy >= GRID_H) continue;
        corpse_grid[idx(gx,gy)].push_back(i);
    }

    // ID → 配列位置 の対応表(毎フレーム作り直す)
    static std::unordered_map<int,int> id2idx;
    id2idx.clear();
    for(int i = 0; i < (int)cells.size(); i++)
        if(cells[i].alive) id2idx[cells[i].id] = i;

    // --- 既存の結合: 距離制約と、切れる判定 ---
    for(int i = 0; i < (int)cells.size(); i++){
        Cell& a = cells[i];
        if(!a.alive) continue;
        int w = 0;
        for(int b = 0; b < a.n_bonds; b++){
            auto it = id2idx.find(a.bonds[b]);
            if(it == id2idx.end()) continue;          // 相手が死んだ → 結合が消える
            Cell& o = cells[it->second];
            float dx = o.x - a.x, dy = o.y - a.y;
            if(dx >  GRID_W/2.f) dx -= GRID_W;         // 東西ループ
            if(dx < -GRID_W/2.f) dx += GRID_W;
            float d = std::sqrt(dx*dx + dy*dy);
            if(d > BOND_BREAK) continue;               // 離れすぎ → 切れる
            if(d > 0.001f){
                float corr = (d - BOND_DIST) * BOND_STIFF * 0.5f;
                a.x += dx/d * corr;  a.y += dy/d * corr;
                o.x -= dx/d * corr;  o.y -= dy/d * corr;
            }
            // 資源共有: 豊かな方から、乏しい方へ流れる
            float flow = (a.energy - o.energy) * SHARE_RATE * 0.5f;
            a.energy -= flow;
            o.energy += flow;
            a.bonds[w++] = a.bonds[b];                 // 生きてる結合だけ残す
        }
        a.n_bonds = w;
    }

    // --- 新しい結合: 同じマスにいて、署名が合えばくっつく ---
    for(int i = 0; i < (int)cells.size(); i++){
        Cell& a = cells[i];
        if(!a.alive || a.adhesion < ADHESION_MIN || a.n_bonds >= MAX_BONDS) continue;
        int gx = (int)a.x, gy = (int)a.y;
        if(gx < 0 || gx >= GRID_W || gy < 0 || gy >= GRID_H) continue;
        for(int j : cell_grid[idx(gx,gy)]){
            if(j == i) continue;
            Cell& o = cells[j];
            if(!o.alive || o.adhesion < ADHESION_MIN || o.n_bonds >= MAX_BONDS) continue;
            if(std::fabs(a.signature - o.signature) > SIG_TOLERANCE) continue;
            bool already = false;
            for(int b = 0; b < a.n_bonds; b++) if(a.bonds[b] == o.id) already = true;
            if(already) continue;
            a.bonds[a.n_bonds++] = o.id;
            o.bonds[o.n_bonds++] = a.id;
            if(a.n_bonds >= MAX_BONDS) break;
        }
    }

    std::vector<Cell> babies;

    for(int i = 0; i < (int)cells.size(); i++){
        Cell& c = cells[i];
        if(!c.alive) continue;
        int gx = (int)c.x, gy = (int)c.y;
        if(gx < 0 || gx >= GRID_W || gy < 0 || gy >= GRID_H) continue;

        // 場から吸う
        float take = energy[idx(gx,gy)] * c.absorb * (1.f - c.trophic * c.trophic);
        energy[idx(gx,gy)] -= take;
        c.energy += take;

        if(c.trophic > 0.02f){
            // 死骸を漁る(抵抗しないので安く食える = 谷の橋)
            const std::vector<int>& dead = corpse_grid[idx(gx,gy)];
            if(!dead.empty()){
                int pick = dead[(int)(frand(0.f, (float)dead.size() - 0.001f))];
                float bite = corpses[pick].energy * SCAVENGE_EFF * c.trophic * (1.f - c.adhesion);
                corpses[pick].energy -= bite;
                c.energy += bite * (1.f - SCAVENGE_LOSS);
            }
            // 生きた相手を襲う
            const std::vector<int>& here = cell_grid[idx(gx,gy)];
            if(here.size() > 1){
                int pick = here[(int)(frand(0.f, (float)here.size() - 0.001f))];
                // 結合相手は、接触面が結合に占有されているので包み込めない
                bool bound_to = false;
                for(int b = 0; b < c.n_bonds; b++)
                    if(pick >= 0 && pick < (int)cells.size() && c.bonds[b] == cells[pick].id) bound_to = true;
                if(pick != i && !bound_to && cells[pick].alive){
                    // 接着分子が表面を占有するぶん、捕食に使える面が減る
                    float bite = cells[pick].energy * PREDATE_EFF * c.trophic * (1.f - c.adhesion);
                    cells[pick].energy -= bite;
                    c.energy += bite * (1.f - PREDATE_LOSS);
                    if(cells[pick].energy <= 0.f){
                        cells[pick].alive = false;
                        corpses.push_back({cells[pick].x, cells[pick].y, cells[pick].energy + 0.5f});
                    }
                }
            }
        }
        // 移動(コストを払って動く)
        if(c.speed > 0.001f){
            c.dir += frand(-c.turn, c.turn);
            float sp = c.speed * MAX_SPEED;
            c.x += std::cos(c.dir) * sp;
            c.y += std::sin(c.dir) * sp;
            if(c.x < 0.f)      c.x += GRID_W;
            if(c.x >= GRID_W)  c.x -= GRID_W;
            if(c.y < 0.f)      { c.y = 0.f; c.dir = -c.dir; }
            if(c.y >= GRID_H)  { c.y = GRID_H - 1.f; c.dir = -c.dir; }
            c.energy -= sp * MOVE_COST + c.turn * TURN_COST;
        }
        c.energy -= METABOLISM + c.adhesion * ADHESION_COST;
        if(c.energy <= 0.f){
            c.alive = false;
            corpses.push_back({c.x, c.y, 0.3f});   // 餓死体は痩せている
            continue;
        }

        if(c.energy >= DIVIDE_THRESHOLD){
            c.energy *= 0.5f;
            Cell child = c;
            child.x = c.x + frand(-2.f, 2.f);
            child.y = c.y + frand(-2.f, 2.f);
            if(child.x < 0.f)      child.x += GRID_W;
            if(child.x >= GRID_W)  child.x -= GRID_W;
            if(child.y < 0.f)      child.y = 0.f;
            if(child.y >= GRID_H)  child.y = GRID_H - 1.f;
            child.absorb  = mutate(c.absorb);
            child.trophic = mutate(c.trophic);
            child.speed = mutate(c.speed);
            child.turn  = mutate(c.turn);
            child.adhesion  = mutate(c.adhesion);
            child.signature = mutate(c.signature);
            child.id        = next_cell_id++;
            child.n_bonds   = 0;
            for(int b = 0; b < MAX_BONDS; b++) child.bonds[b] = -1;
            child.dir   = frand(0.f, 6.2832f);
            babies.push_back(child);
        }
    }

    for(auto& b : babies) cells.push_back(b);

    cells.erase(std::remove_if(cells.begin(), cells.end(),
        [](const Cell& c){ return !c.alive; }), cells.end());

    // 死骸が腐って、場のエネルギーに戻る(物質循環)
    for(auto& cp : corpses){
        float rot = cp.energy * ROT_RATE;
        cp.energy -= rot;
        int gx = (int)cp.x, gy = (int)cp.y;
        if(gx >= 0 && gx < GRID_W && gy >= 0 && gy < GRID_H)
            energy[idx(gx,gy)] += rot;
    }
    corpses.erase(std::remove_if(corpses.begin(), corpses.end(),
        [](const Corpse& c){ return c.energy < 0.01f; }), corpses.end());
}

void save_state(const std::string& fname, int frame){
    std::ofstream f(fname);
    if(!f){ printf("save failed: %s\n", fname.c_str()); return; }
    f.precision(9);
    f << "EVOSAVE 1\n";
    f << frame << "\n";

    f << "FIELD " << (GRID_W*GRID_H) << "\n";
    for(int i = 0; i < GRID_W*GRID_H; i++) f << energy[i] << "\n";

    int n = 0; for(const auto& c : cells) if(c.alive) n++;
    f << "CELLS " << n << "\n";
    for(const auto& c : cells){
        if(!c.alive) continue;
        f << c.x << " " << c.y << " " << c.energy << " "
          << c.absorb << " " << c.trophic << " "
          << c.dir << " " << c.speed << " " << c.turn << "\n";
    }

    f << "CORPSES " << corpses.size() << "\n";
    for(const auto& cp : corpses)
        f << cp.x << " " << cp.y << " " << cp.energy << "\n";

    printf("saved: %s (cells %d, corpses %d)\n", fname.c_str(), n, (int)corpses.size());
}

bool load_state(const std::string& fname, int& frame){
    std::ifstream f(fname);
    if(!f){ printf("load failed (no file): %s\n", fname.c_str()); return false; }
    std::string tag; int ver;
    f >> tag >> ver;
    if(tag != "EVOSAVE" || ver > 1){ printf("load failed (bad format)\n"); return false; }
    try {
        f >> frame;

        std::string sec; int cnt;
        f >> sec >> cnt;
        if(cnt != GRID_W*GRID_H) throw std::runtime_error("grid size mismatch");
        for(int i = 0; i < cnt; i++) f >> energy[i];

        f >> sec >> cnt;
        if(cnt < 0 || cnt > 5000000) throw std::runtime_error("bad cell count");
        cells.clear();
        for(int i = 0; i < cnt; i++){
            Cell c;
            f >> c.x >> c.y >> c.energy >> c.absorb >> c.trophic
              >> c.dir >> c.speed >> c.turn;
            c.alive = true;
            cells.push_back(c);
                           
        }

        f >> sec >> cnt;
        if(cnt < 0 || cnt > 5000000) throw std::runtime_error("bad corpse count");
        corpses.clear();
        for(int i = 0; i < cnt; i++){
            Corpse cp;
            f >> cp.x >> cp.y >> cp.energy;
            corpses.push_back(cp);
        }

        printf("loaded: %s (cells %d, corpses %d)\n",
               fname.c_str(), (int)cells.size(), (int)corpses.size());
        return true;
    } catch(const std::exception& e){
        printf("load failed (corrupt): %s\n", e.what());
        return false;
    }
}

int main(){
    init_solar();
    int frame = 0;
    int autosave_slot = 0;
    {
        Cell c;
        c.x = GRID_W/2.f; c.y = GRID_H/2.f;
        c.energy = INIT_ENERGY; c.alive = true;
        c.absorb = INIT_ABSORB; c.trophic = INIT_TROPHIC;
        c.dir = 0.f; c.speed = INIT_SPEED; c.turn = 0.1f;
        c.id = next_cell_id++;
        c.adhesion = INIT_ADHESION; c.signature = 0.5f;
        c.n_bonds = 0;
        for(int b = 0; b < MAX_BONDS; b++) c.bonds[b] = -1;
        cells.push_back(c);
    }

    sf::RenderWindow window(
        sf::VideoMode({(unsigned)(GRID_W*SCALE), (unsigned)(GRID_H*SCALE)}),
        "Evo Ecosystem - Phase 2");
    window.setFramerateLimit(60);

    sf::Font font;
    bool font_ok = font.openFromFile("/System/Library/Fonts/Helvetica.ttc");
    if(!font_ok) font_ok = font.openFromFile("/System/Library/Fonts/Supplemental/Arial.ttf");
    if(!font_ok) font_ok = font.openFromFile("C:\\Windows\\Fonts\\arial.ttf");

    sf::Image image({(unsigned)GRID_W, (unsigned)GRID_H}, sf::Color::Black);
    sf::Texture texture;
    if(!texture.loadFromImage(image)){ printf("texture load failed\n"); return 1; }
    sf::Sprite sprite(texture);
    sprite.setScale({(float)SCALE, (float)SCALE});



    while(window.isOpen()){
        while(const std::optional event = window.pollEvent()){
            if(event->is<sf::Event::Closed>()) window.close();
            if(const auto* k = event->getIf<sf::Event::KeyPressed>()){
                if(k->code == sf::Keyboard::Key::Escape) window.close();
                if(k->code == sf::Keyboard::Key::K) save_state("save.txt", frame);
                if(k->code == sf::Keyboard::Key::L) load_state("save.txt", frame);
            }
        }

        update_field();
        update_cells();            // ← ここ
        
        frame++;
        if(frame % SAMPLE_INTERVAL == 0){
            double sa=0.0, st=0.0, ss=0.0;
            for(const auto& c : cells){ sa+=c.absorb; st+=c.trophic; ss+=c.speed; }
            int n = (int)cells.size();
            h_absorb.push_back(n ? (float)(sa/n) : 0.f);
            h_trophic.push_back(n ? (float)(st/n) : 0.f);
            h_speed.push_back(n ? (float)(ss/n) : 0.f);
            h_pop.push_back(n);
            h_corpse.push_back((int)corpses.size());
            if((int)h_pop.size() > HIST_MAX){
                h_absorb.erase(h_absorb.begin());
                h_trophic.erase(h_trophic.begin());
                h_speed.erase(h_speed.begin());
                h_pop.erase(h_pop.begin());
                h_corpse.erase(h_corpse.begin());
            }
        }
        if(frame % 60 == 0){
            double sa=0.0, st=0.0, ss=0.0, sad=0.0;
            int nmove=0, nb=0, nbonded=0;
            for(const auto& c : cells){
                sa += c.absorb; st += c.trophic; ss += c.speed; sad += c.adhesion;
                if(c.speed > 0.2f) nmove++;
                nb += c.n_bonds;
                if(c.n_bonds > 0) nbonded++;
            }
            int n = (int)cells.size();
            printf("t:%d pop:%d absorb:%.3f trophic:%.3f speed:%.3f adhesion:%.3f bonds:%d bonded:%d corpse:%d\n",
                   frame, n,
                   n ? sa/n : 0.0, n ? st/n : 0.0, n ? ss/n : 0.0, n ? sad/n : 0.0,
                   nb/2, nbonded, (int)corpses.size());
        }
        if(frame % 7200 == 0){
            char fn[32]; snprintf(fn, 32, "autosave%d.txt", autosave_slot);
            save_state(fn, frame);
            autosave_slot = (autosave_slot + 1) % 3;
        }


        for(int y = 0; y < GRID_H; y++)
            for(int x = 0; x < GRID_W; x++)
                image.setPixel({(unsigned)x, (unsigned)y}, value_to_color(energy[idx(x,y)]));

        texture.update(image);



        

        window.clear();
        window.draw(sprite);
        // 死骸をバッチ描画(暗い赤褐色)
        static sf::VertexArray cps(sf::PrimitiveType::Triangles);
        cps.clear();
        {
            const float r = SCALE * 0.8f;
            sf::Color col(120, 40, 30);
            for(const auto& cp : corpses){
                float px = cp.x * SCALE, py = cp.y * SCALE;
                sf::Vertex v0,v1,v2,v3;
                v0.position={px-r,py-r}; v0.color=col;
                v1.position={px+r,py-r}; v1.color=col;
                v2.position={px+r,py+r}; v2.color=col;
                v3.position={px-r,py+r}; v3.color=col;
                cps.append(v0); cps.append(v1); cps.append(v2);
                cps.append(v0); cps.append(v2); cps.append(v3);
            }
        }
        window.draw(cps);

        // 細胞をバッチ描画(緑=独立栄養, 赤=従属栄養)
        static sf::VertexArray pts(sf::PrimitiveType::Triangles);
        pts.clear();
        {
            const float r = SCALE * 1.2f;
            for(const auto& c : cells){
                if(!c.alive) continue;
                sf::Color col(
                    (std::uint8_t)(60 + 195 * c.trophic),
                    (std::uint8_t)(220 - 160 * c.trophic),
                    (std::uint8_t)(80));
                float px = c.x * SCALE, py = c.y * SCALE;
                sf::Vertex v0,v1,v2,v3;
                v0.position={px-r,py-r}; v0.color=col;
                v1.position={px+r,py-r}; v1.color=col;
                v2.position={px+r,py+r}; v2.color=col;
                v3.position={px-r,py+r}; v3.color=col;
                pts.append(v0); pts.append(v1); pts.append(v2);
                pts.append(v0); pts.append(v2); pts.append(v3);
            }
        }
        window.draw(pts);
        {
            static std::unordered_map<int,int> vmap;
            vmap.clear();
            for(int i = 0; i < (int)cells.size(); i++) vmap[cells[i].id] = i;
            sf::VertexArray bl(sf::PrimitiveType::Lines);
            for(const auto& a : cells){
                for(int b = 0; b < a.n_bonds; b++){
                    auto it = vmap.find(a.bonds[b]);
                    if(it == vmap.end()) continue;
                    const Cell& o = cells[it->second];
                    if(std::fabs(o.x - a.x) > GRID_W/2.f) continue;
                    sf::Vertex v0, v1;
                    v0.position = {a.x*SCALE, a.y*SCALE}; v0.color = sf::Color(255,255,255,150);
                    v1.position = {o.x*SCALE, o.y*SCALE}; v1.color = sf::Color(255,255,255,150);
                    bl.append(v0); bl.append(v1);
                }
            }
            window.draw(bl);
        }
        // ===== 履歴グラフ =====
        if(h_pop.size() > 1){
            float GW = GRID_W * SCALE * 0.30f;
            float GH = GRID_H * SCALE * 0.28f;
            float GX = GRID_W * SCALE - GW - 12.f;
            float GY = 12.f;

            sf::RectangleShape bg({GW, GH});
            bg.setPosition({GX, GY});
            bg.setFillColor(sf::Color(0, 0, 0, 190));
            bg.setOutlineColor(sf::Color(120, 120, 120));
            bg.setOutlineThickness(1.f);
            window.draw(bg);

            int N = (int)h_pop.size();
            float pad = 6.f;
            float px0 = GX + pad, py0 = GY + pad;
            float pw = GW - pad*2, ph = GH - pad*2;

            // 0..1 の形質を3本
            auto drawTrait = [&](const std::vector<float>& v, sf::Color col){
                sf::VertexArray ln(sf::PrimitiveType::LineStrip);
                for(int i = 0; i < N; i++){
                    sf::Vertex vt;
                    vt.position = { px0 + pw * i / (float)(N-1),
                                    py0 + ph * (1.f - std::min(1.f, v[i])) };
                    vt.color = col;
                    ln.append(vt);
                }
                window.draw(ln);
            };
            drawTrait(h_trophic, sf::Color(255, 90, 70));    // 赤: 従属栄養度
            drawTrait(h_absorb,  sf::Color(90, 230, 110));   // 緑: 吸収効率
            drawTrait(h_speed,   sf::Color(110, 170, 255));  // 青: 速度

            // 個体数(ピークで正規化)
            int peak = 1;
            for(int p : h_pop) if(p > peak) peak = p;
            sf::VertexArray lp(sf::PrimitiveType::LineStrip);
            for(int i = 0; i < N; i++){
                sf::Vertex vt;
                vt.position = { px0 + pw * i / (float)(N-1),
                                py0 + ph * (1.f - h_pop[i] / (float)peak) };
                vt.color = sf::Color(200, 200, 200, 140);
                lp.append(vt);
            }
            window.draw(lp);

            if(font_ok){
                char buf[160];
                snprintf(buf, 160, "trophic %.3f  absorb %.3f  speed %.3f",
                         h_trophic.back(), h_absorb.back(), h_speed.back());
                sf::Text t1(font, buf, 11);
                t1.setFillColor(sf::Color(220,220,220));
                t1.setPosition({GX, GY + GH + 3.f});
                window.draw(t1);

                snprintf(buf, 160, "pop %d (peak %d)   corpse %d   t=%d",
                         h_pop.back(), peak, h_corpse.back(), frame);
                sf::Text t2(font, buf, 11);
                t2.setFillColor(sf::Color(180,180,180));
                t2.setPosition({GX, GY + GH + 17.f});
                window.draw(t2);
            }
        }
        window.display();
    }
    return 0;
}
