#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <algorithm>
#include <functional>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <ctime>
#include <string>
#include <stdexcept>
#include <unordered_map>
#include <queue>
#ifdef __APPLE__
  #include <IOKit/pwr_mgt/IOPMLib.h>
#elif defined(_WIN32)
  #include <windows.h>
#endif
#include <thread>
#include <mutex>
#include <condition_variable>
#include <functional>
#include "src/params.h"
#include "src/state.h"
#include "src/lineage.h"
#include "src/pool.h"
#include "src/terrain.h"
#include "src/water.h"
#include "src/air.h"
#include "src/pathogen_odor.h"
#include "src/slide.h"
#include "src/quake.h"
#include "src/microbe.h"
#include "src/weather_fire.h"
#include "src/plant.h"
#include "src/brain.h"
#include "src/field.h"
#include "src/animal.h"
#include "src/save.h"
#include "src/session.h"
#include "src/draw.h"
#include "src/assay.h"
#include "src/census.h"
int main(int argc,char** argv){
    // 起動引数: --seed N  --stop-at N  --record FILE  --replay FILE --shuffle-seed N  --headless  --quake T M
    std::string rec_path, replay_path; unsigned shuffle_seed=1;
    std::string assay_path, assay_save, assay_out; int assay_max=1000;
    for(int i=1;i<argc;i++){
        std::string a=argv[i];
        if(a=="--seed"    && i+1<argc){ unsigned sd=(unsigned)std::stoul(argv[++i]); rng.seed(sd); wx_rng.seed(sd^0x2545F491u); lt_rng.seed(sd^0x9E3779B9u); qk_rng.seed(sd^0x85EBCA6Bu); }
        if(a=="--stop-at" && i+1<argc) g_stop_at=std::stol(argv[++i]);
        if(a=="--record"  && i+1<argc) rec_path=argv[++i];
        if(a=="--replay"  && i+1<argc) replay_path=argv[++i];
        if(a=="--shuffle-seed" && i+1<argc) shuffle_seed=(unsigned)std::stoul(argv[++i]);
        if(a=="--headless") g_headless=true;
        if(a=="--dump-brains" && i+1<argc) g_dump_every=std::stoi(argv[++i]);
        if(a=="--dump-max"    && i+1<argc) g_dump_max=std::max(1,std::stoi(argv[++i]));
        if(a=="--assay"       && i+1<argc) assay_path=argv[++i];
        if(a=="--assay-save"  && i+1<argc) assay_save=argv[++i];
        if(a=="--assay-max"   && i+1<argc) assay_max=std::max(1,std::stoi(argv[++i]));
        if(a=="--assay-out"   && i+1<argc) assay_out=argv[++i];
        if(a=="--track"       && i+1<argc) g_track_n=std::max(0,std::stoi(argv[++i]));
        if(a=="--track-from"  && i+1<argc) g_track_from=std::stoi(argv[++i]);
        if(a=="--trigger-log" && i+1<argc) g_trig_w=std::max(0,std::stoi(argv[++i]));
        if(a=="--trigger-max" && i+1<argc) g_trig_max=std::max(0,std::stoi(argv[++i]));
        if(a=="--positive-control" && i+1<argc) g_pc_n=std::max(0,std::stoi(argv[++i]));
        if(a=="--no-recurrent"     && i+1<argc) g_lesion_n=std::max(0,std::stoi(argv[++i]));
        if(a=="--assay-series" && i+1<argc) g_series_path=argv[++i];
        if(a=="--census"      && i+1<argc) g_census_every=std::max(0,std::stoi(argv[++i]));
        if(a=="--quake"       && i+2<argc){ g_quake_force_t=std::stol(argv[++i]); g_quake_force_m=std::stof(argv[++i]); }
    }
    field_open();
    census_open();
    // 試験器は世界を作らずに動き、終わったら抜ける
    if(!assay_path.empty()||!assay_save.empty()) return run_assay(assay_path,assay_save,assay_max,assay_out);
    if(!rec_path.empty()){
        g_rec=fopen(rec_path.c_str(),"w");
        if(!g_rec){ fprintf(stderr,"record: cannot open %s\n",rec_path.c_str()); return 1; }
        fprintf(g_rec,"# evosim environment events v1: tick kind x y strength\n");
        fprintf(g_rec,"# fire strength = fire intensity at ignition; slide strength = slope; storm strength = peak rain rate mm/h (x y may lie upwind outside the plot); lightning strength = 1 per cloud-to-ground flash; quake strength = magnitude (x y = epicenter on the fault outside the plot)\n");
    }
    if(!replay_path.empty()){
        if(!ev_load_replay(replay_path.c_str(),shuffle_seed)) return 1;
        g_replay=true;
        printf("replay: %s shuffle-seed %u (fire/slide/storm/quake physics off)\n",replay_path.c_str(),shuffle_seed);
    }
    prevent_sleep_begin();
    for(int k=0;k<N_ODOR;k++) odor[k].assign(GRID_W*GRID_H,0.f);
    odor_buf.assign(GRID_W*GRID_H,0.f);
    for(int k=0;k<N_IMMUNE;k++){
        pstrain[k].assign(GRID_W*GRID_H,0.5f);
        astrain[k].assign(GRID_W*GRID_H,0.5f);
    }
    init_terrain();
    compute_currents();      // 吹送流を求める
    compute_rivers();
    {   // 地形の要約。土地の形が報告書で確かめられるように
        long land=0, steep=0, riv=0; double ssum=0, smax=0, hmax=0, dmax=0;
        for(int i=0;i<GRID_W*GRID_H;i++){
            float zm=(elev[i]-SEA_LEVEL)*ELEV_M;
            if(sea[i]){ if(-zm>dmax) dmax=-zm; continue; }
            land++; ssum+=slope[i]; if(slope[i]>smax) smax=slope[i];
            if(slope[i]>=SLIDE_SLOPE) steep++;
            if(river[i]) riv++;
            if(zm>hmax) hmax=zm;
        }
        printf("terrain: lat %.1fN | land %.1f%% | top %.0f m  sea depth %.0f m | slope mean %.3f max %.3f | >=30deg %.2f%% | river %.2f%%\n",
               LATITUDE_DEG,100.0*land/(GRID_W*GRID_H),hmax,dmax,land?ssum/land:0.0,smax,
               land?100.0*steep/land:0.0,land?100.0*riv/land:0.0);
    }
    update_climate(0);
    {   // 起動時の土壌水分は条件。反復で近づけず、釣り合いの値(緩和の目標)を最初から置く
        // 降水は陸からの再循環で土壌水分に依る。仮の値で降水を決め、その降水の目標を置き、
        // 置いた土壌水分で降水を決め直す。残る食い違いは再循環の分だけで、ログに出す
        for(int i=0;i<GRID_W*GRID_H;i++) soil_water[i]=sea[i]?FIELD_CAP:SOIL_W_INIT;
        compute_precip();
        for(int i=0;i<GRID_W*GRID_H;i++) if(!sea[i]) soil_water[i]=water_target(i,1.f);  // t=0 は雨季と乾季の境(wet=1)
        compute_precip();
        double sw=0, gap=0, gmax=0; long land=0;
        for(int i=0;i<GRID_W*GRID_H;i++){
            if(sea[i]) continue;
            double g=std::fabs(water_target(i,1.f)-soil_water[i]);
            sw+=soil_water[i]; gap+=g; if(g>gmax) gmax=g; land++;
        }
        printf("water init: land soil mean %.4f | |target-soil| mean %.5f max %.5f\n",
               land?sw/land:0.0,land?gap/land:0.0,gmax);
    }
    classify_climate();
    // 区画は惑星の一部で、まわりの陸と海から種が入ってくるのが始まりの条件(D-2026-10-24-1 の作者の答え (b))。
    // 最初の PLANT_INIT_LAND 株は陸に、陸に合う遺伝子で置き、残りを浅い海に置く。
    // 陸の株は、種子の発芽と同じく土の水が足りるマスにだけ置く
    int placed=0,guard=0;
    while(placed<PLANT_INIT_N && guard++<400000){
        Plant p;
        p.x=frand(0.f,(float)GRID_W); p.y=frand(0.f,(float)GRID_H);
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        bool on_land=placed<PLANT_INIT_LAND;
        int ci0=idx(gx,gy);
        if(on_land){
            if(sea[ci0]) continue;
            if(soil_water[ci0]<0.25f) continue;
        } else {
            if(!sea[ci0]) continue;
            if(depthn(ci0)>0.35f) continue;         // 海の株は浅い海に
        }
        sea_plant_genes(p);
        // 陸の株は、まわりの陸の温帯の草本の性質で入る: 水に頼らず(aquatic 0)、根を張り(浮く強さ 0.3)、
        // 生育の下限の気温を下回ると休眠に入る
        if(on_land){ p.aquatic=0.0f; p.salt_tol=0.0f; p.buoyancy=0.3f; p.dorm_temp=PLANT_LAND_DORM; }
        // 条件(D-2026-10-24-1 の作者の答え (c)): 最初の海の植物の半分は海底に根を張る(浮く強さ 0.3。水面の漂流を受けない)。
        // 浮く強さは遺伝子なので、どちらの生き方が残るかは進化に任せる。乱数は引かない
        else if(placed%2==1) p.buoyancy=0.3f;
        plants.push_back(p); placed++;
    }
    // 分解者を撒く。死骸だけを食い、炭素と窒素を還す第三の界
    for(int i=0;i<MICROBE_INIT;i++){
        Microbe m;
        // 陸に撒く。海底は届かない
        int tries=0;
        do { m.x=frand(0.f,(float)GRID_W); m.y=frand(0.f,(float)GRID_H); }
        while(sea[idx((int)m.x,(int)m.y)] && ++tries<50);
        m.energy=0.4f; m.alive=true; m.age=0;
        m.rate=0.3f; m.cold_tol=0.2f;
        m.id=next_microbe_id++;
        microbes.push_back(m);
    }
    sf::VideoMode desktop=sf::VideoMode::getDesktopMode();
    SCREEN_W=(float)desktop.size.x; SCREEN_H=(float)desktop.size.y;
    sf::RenderWindow window;   // --headless では開かない
    if(!g_headless){
        window.create(desktop,"Phase 3",sf::State::Fullscreen);
        window.setFramerateLimit(60);
    }
    bool fullscreen=true;

    sf::Font font;
    // 日本語の字形を持つフォントを順に試す
    bool font_ok=false;
    {
        const char* cands[]={
            "/System/Library/Fonts/ヒラギノ角ゴシック W3.ttc",
            "/System/Library/Fonts/Hiragino Sans GB.ttc",
            "/System/Library/Fonts/Supplemental/Arial Unicode.ttf",
            "/Library/Fonts/Arial Unicode.ttf",
            "C:/Windows/Fonts/meiryo.ttc",
            "C:/Windows/Fonts/msgothic.ttc",
            "C:/Windows/Fonts/YuGothM.ttc",
            "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
            "/System/Library/Fonts/Helvetica.ttc"
        };
        for(const char* p:cands){
            if(font.openFromFile(p)){ font_ok=true; printf("font: %s\n",p); break; }
        }
        if(!font_ok) printf("font: none found\n");
    }

    sf::Image image({(unsigned)GRID_W,(unsigned)GRID_H},sf::Color::Black);
    sf::Texture texture;
    if(!g_headless && !texture.loadFromImage(image)){ printf("texture failed\n"); return 1; }
    sf::Sprite sprite(texture);

    float zoom=std::min(SCREEN_W/(float)GRID_W,SCREEN_H/(float)GRID_H);
    float cam_x=GRID_W/2.f, cam_y=GRID_H/2.f;
    float cam_vx=0.f, cam_vy=0.f;   // カメラの速度。慣性で滑らかに動く
    bool show_help=false;           // 左上の「?」で操作一覧を開く
    // 東西も南北と同じく区画の縁で終わる
    auto w2sx=[&](float wx){ return (wx-cam_x)*zoom+SCREEN_W/2.f; };
    auto w2sy=[&](float wy){ return (wy-cam_y)*zoom+SCREEN_H/2.f; };
    auto s2wx=[&](float sx){
        return (sx-SCREEN_W/2.f)/zoom+cam_x;
    };
    auto s2wy=[&](float sy){ return (sy-SCREEN_H/2.f)/zoom+cam_y; };
    auto clamp_cam=[&](){
        float hw=SCREEN_W/(2.f*zoom);            // 東西も壁
        if(hw>=GRID_W*0.5f) cam_x=GRID_W*0.5f;
        else { if(cam_x<hw) cam_x=hw; if(cam_x>GRID_W-hw) cam_x=GRID_W-hw; }
        float hh=SCREEN_H/(2.f*zoom);            // 南北は壁
        if(hh>=GRID_H*0.5f) cam_y=GRID_H*0.5f;
        else { if(cam_y<hh) cam_y=hh; if(cam_y>GRID_H-hh) cam_y=GRID_H-hh; }
    };
    auto zoom_at=[&](float f,float mxs,float mys){
        float wx=(mxs-SCREEN_W/2.f)/zoom+cam_x;   // 折り返さない生の座標
        float wy=(mys-SCREEN_H/2.f)/zoom+cam_y;
        zoom*=f;
        float zmin=std::min(SCREEN_W/(float)GRID_W,SCREEN_H/(float)GRID_H)*0.5f;
        if(zoom<zmin) zoom=zmin;
        if(zoom>80.f) zoom=80.f;
        cam_x = wx - (mxs-SCREEN_W/2.f)/zoom;
        cam_y = wy - (mys-SCREEN_H/2.f)/zoom;
        clamp_cam();
    };
    int   nrep  = 0;      // 横方向に何枚並べるか
    float wrapW = 0.f;    // 世界一周ぶんの画面幅

    int frame=0, autosave_slot=0, view_mode=0;
    bool released=false, paused=false, show_graph=true;
    int sel_animal=-1, sel_plant=-1, sel_lineage=-1;
    int ui_mode=0; std::string input_text; bool ui_skip_char=false;
    std::vector<std::string> save_files;
    std::vector<std::time_t> save_times;
    bool sort_by_time=false;
    bool show_tree=false, color_species=false;
    float tree_zoom=1.f, tree_ox=0.f, tree_oy=0.f;
    float tree_vx=0.f, tree_vy=0.f;   // 系統樹の移動も慣性で滑らかにする
    std::string tree_search;
    int filter_mode=0;
    // --- 下部バーと拡大グラフ ---
    struct Hit { float x,y,w,h; int id; };
    std::vector<Hit> ui_hits;              // クリックできる領域
    int   hover_save=-1;                   // セーブ一覧でカーソルが乗っている行
    int   del_confirm=-1;                  // 削除の確認待ちの行
    int   graph_big=-1;                    // 拡大中のグラフ(-1で無し)
    int   graph_scale=0;                   // 0:1千 1:1万 2:10万 3:全期間
    int   graph_off=0;                     // 時間軸の移動量
    float mouse_x=0.f, mouse_y=0.f;
    // WASD は毎フレームの polling だと、1フレームが長くなったとき取りこぼす。
    // 押下と離上をイベントで受けて、状態として保持する
    bool key_w=false,key_a=false,key_s=false,key_d=false;
    std::unordered_map<std::string,std::string> save_info;  // ホバー用の概要

    // 系統樹の当たり判定用(描画時に更新)
    struct TreeHit { int li; float x,y0,y1; };
    std::vector<TreeHit> tree_hits;
    struct FilterHit { int f; float x,y,w,h; };
    std::vector<FilterHit> filter_hits;

    bool running=true;
    while(running && (g_headless || window.isOpen())){
        // 1フレームの実時間を測る。t>=50000 の平均 fps を夜間の比較に使う
        { static auto last=std::chrono::steady_clock::now();
          auto now=std::chrono::steady_clock::now();
          double dt=std::chrono::duration<double>(now-last).count();
          last=now;
          if(frame>=50000){ g_fps_time+=dt; g_fps_n++; } }
        while(const std::optional event = g_headless ? std::optional<sf::Event>{} : window.pollEvent()){
            if(event->is<sf::Event::Closed>()) window.close();

            if(const auto* te=event->getIf<sf::Event::TextEntered>()){
                char32_t u=te->unicode;
                if(ui_mode!=0){
                    if(ui_skip_char){ ui_skip_char=false; }
                    else if(u==8){ if(!input_text.empty()) input_text.pop_back(); }
                    else if(u>=32&&u<127){ if(input_text.size()<40) input_text+=(char)u; }
                } else if(show_tree){
                    if(u==8){ if(!tree_search.empty()) tree_search.pop_back(); }
                    else if(u>='0'&&u<='9'){ if(tree_search.size()<6) tree_search+=(char)u; }
                }
            }

            if(const auto* kr=event->getIf<sf::Event::KeyReleased>()){
                if(kr->code==sf::Keyboard::Key::W) key_w=false;
                if(kr->code==sf::Keyboard::Key::A) key_a=false;
                if(kr->code==sf::Keyboard::Key::S) key_s=false;
                if(kr->code==sf::Keyboard::Key::D) key_d=false;
            }
            if(const auto* k=event->getIf<sf::Event::KeyPressed>()){
                if(k->code==sf::Keyboard::Key::W) key_w=true;
                if(k->code==sf::Keyboard::Key::A) key_a=true;
                if(k->code==sf::Keyboard::Key::S) key_s=true;
                if(k->code==sf::Keyboard::Key::D) key_d=true;
                if(ui_mode!=0){
                    if(k->code==sf::Keyboard::Key::Enter){
                        if(!input_text.empty()){
                            if(ui_mode==1) save_state(input_text+".txt",frame);
                            else {
                                std::string nm=input_text;
                                size_t p=nm.find("   (ver");
                                if(p!=std::string::npos) nm=nm.substr(0,p);
                                load_state(nm+".txt",frame);
                            }
                        }
                        ui_mode=0; input_text.clear();
                    }
                    if(k->code==sf::Keyboard::Key::Escape){
                        ui_mode=0; input_text.clear(); del_confirm=-1; }
                    if(ui_mode==2 && k->code==sf::Keyboard::Key::Tab){
                        sort_by_time=!sort_by_time;
                        std::vector<std::pair<std::string,std::time_t>> fs;
                        for(size_t i=0;i<save_files.size();i++)
                            fs.push_back({save_files[i],save_times[i]});
                        if(sort_by_time)
                            std::sort(fs.begin(),fs.end(),[](auto&a,auto&b){
                                return a.second>b.second; });
                        else
                            std::sort(fs.begin(),fs.end(),[](auto&a,auto&b){
                                return a.first<b.first; });
                        save_files.clear(); save_times.clear();
                        for(auto& f:fs){ save_files.push_back(f.first);
                                         save_times.push_back(f.second); }
                    }
                } else if(show_tree){
                    if(k->code==sf::Keyboard::Key::T){ show_tree=false; tree_search.clear(); }
                    if(k->code==sf::Keyboard::Key::Z) tree_zoom=std::max(0.6f,tree_zoom/1.2f);
                    if(k->code==sf::Keyboard::Key::X){
                        tree_zoom=std::max(1.f,tree_zoom/1.2f);
                        if(tree_zoom<=1.01f){ tree_ox=0.f; tree_oy=0.f; }
                    }
                    if(k->code==sf::Keyboard::Key::F){
                        for(int f=0;f<FILTER_N;f++) filter_on[f]=false;
                        filter_on[0]=true;
                    }
                    if(k->code==sf::Keyboard::Key::Escape){ show_tree=false; tree_search.clear(); }
                } else {
                    if(k->code==sf::Keyboard::Key::Escape){
                        fullscreen=!fullscreen;
                        if(fullscreen){ window.create(desktop,"Phase 3",sf::State::Fullscreen);
                                        SCREEN_W=(float)desktop.size.x; SCREEN_H=(float)desktop.size.y; }
                        else { window.create(sf::VideoMode({1600,900}),"Phase 3");
                               SCREEN_W=1600.f; SCREEN_H=900.f; }
                        window.setFramerateLimit(60);
                    }
                    if(k->code==sf::Keyboard::Key::Space) paused=!paused;
                    if(k->code==sf::Keyboard::Key::G) show_graph=!show_graph;
                    if(k->code==sf::Keyboard::Key::C) color_species=!color_species;
                    if(k->code==sf::Keyboard::Key::R) release_animals(40,frame);
                    if(k->code==sf::Keyboard::Key::P){
                        // 海に光合成プランクトンを放つ。二遺伝子の谷を渡せないので、
                        // 対岸から始めて定着するかを見る
                        int n=0;
                        for(int t=0;t<3000&&n<600;t++){
                            int gx=(int)frand(0.f,(float)GRID_W);
                            int gy=(int)frand(0.f,(float)GRID_H);
                            int ci=idx(gx,gy);
                            if(!sea[ci]||ice[ci]) continue;
                            if(temper[ci]<4.f) continue;
                            Plant p{};
                            p.x=gx+frand(0.f,1.f); p.y=gy+frand(0.f,1.f);
                            p.energy=PLANT_INIT_E; p.height=0.05f; p.eff_height=0.05f;
                            p.alive=true; p.age=0;
                            p.absorb=0.95f; p.max_height=0.10f; p.shade_tol=0.05f;
                            p.tough=0.f; p.clonal=0.f; p.disperse=0.5f;
                            p.cold_tol=0.10f; p.aquatic=1.0f; p.salt_tol=1.0f; p.drought_tol=0.f;
                            p.dorm_temp=-20.f; p.dormant=false;
                            p.store_cap=0.f; p.water_store=0.f;
                            p.fire_tol=0.f; p.climb=0.f; p.sexual=0.1f;
                            p.brood=0.5f; p.n_fix=0.3f; p.n_store=1.0f;
                            p.lifespan=0.05f; p.parasite=0.f; p.carnivory=0.f;
                            p.fire_seed=0.f; p.repro_alloc=0.6f; p.seed_bank=0.f;
                            p.defense=0.f; p.seed_wait=0; p.buoyancy=1.0f;
                            p.mut_rate=0.3f;        p.mother=-1; p.father=-1; p.gen=0;
                            for(int k2=0;k2<N_IMMUNE;k2++) p.im[k2]=frand(0.f,1.f);
                            p.id=next_plant_id++; p.genet=p.id; p.lin=-1;
                            plants.push_back(p); n++;
                            g_bud.c_immig+=(p.energy+p.height*GROW_COST)*CO2_PER_ENERGY;   // 外から入れた浮遊植物の炭素と窒素
                            g_bud.n_immig+=p.n_store+p.height*N_PER_GROWTH;
                        }
                        printf("=== released %d plankton at t=%d ===\n",n,frame);
                    }
                    if(k->code==sf::Keyboard::Key::N){ show_night=!show_night; g_last_view=-1; }
                    if(k->code==sf::Keyboard::Key::T){ show_tree=true; sel_lineage=-1; }
                    if(k->code==sf::Keyboard::Key::Z) zoom_at(0.80f,SCREEN_W/2.f,SCREEN_H/2.f);
                    if(k->code==sf::Keyboard::Key::X) zoom_at(1.25f,SCREEN_W/2.f,SCREEN_H/2.f);
                    if(k->code==sf::Keyboard::Key::K){ ui_mode=1; input_text.clear(); ui_skip_char=true; }
                    if(k->code==sf::Keyboard::Key::L){
                        ui_mode=2; input_text.clear(); ui_skip_char=true;
                        save_files.clear(); save_times.clear();
                        try{
                            std::vector<std::pair<std::string,std::time_t>> fs;
                            for(auto& e:std::filesystem::directory_iterator("."))
                                if(e.is_regular_file()&&e.path().extension()==".txt"){
                                    // 版が古いファイルは読めない。一覧で分かるようにする
                                    std::ifstream tf(e.path());
                                    std::string tg; int vv=0;
                                    if(tf) tf >> tg >> vv;
                                    std::string nm=e.path().stem().string();
                                    if(tg!="P3SAVE"||vv<10) nm+="   (ver "+std::to_string(vv)+" - too old)";
                                    fs.push_back({nm,file_time_to_t(e.last_write_time())});
                                }
                            if(sort_by_time)
                                std::sort(fs.begin(),fs.end(),[](auto&a,auto&b){
                                    return a.second>b.second; });   // 新しい順
                            else
                                std::sort(fs.begin(),fs.end(),[](auto&a,auto&b){
                                    return a.first<b.first; });
                            for(auto& f:fs){ save_files.push_back(f.first);
                                             save_times.push_back(f.second); }
                        }catch(...){}
                    }
                }
            }
            if(const auto* mm=event->getIf<sf::Event::MouseMoved>()){
                mouse_x=(float)mm->position.x; mouse_y=(float)mm->position.y;
            }
            if(const auto* mw=event->getIf<sf::Event::MouseWheelScrolled>()){
                if(graph_big>=0){
                    float wx=(float)mw->position.x, wy=(float)mw->position.y;
                    bool inside=false;
                    for(const auto& h:ui_hits)
                        if(h.id==220&&wx>=h.x&&wx<=h.x+h.w&&wy>=h.y&&wy<=h.y+h.h)
                            inside=true;
                    if(inside){
                        graph_off+=(int)(mw->delta*12.f);
                        continue;
                    }
                }
                if(show_tree){
                    if(mw->delta>0.f) tree_zoom=std::min(12.f,tree_zoom*1.2f);
                    else if(mw->delta<0.f){
                        tree_zoom=std::max(1.f,tree_zoom/1.2f);
                        if(tree_zoom<=1.01f){ tree_ox=0.f; tree_oy=0.f; }
                    }
                } else {
                    sf::Vector2i mp=sf::Mouse::getPosition(window);
                    if(mw->delta>0.f)      zoom_at(1.25f,(float)mp.x,(float)mp.y);
                    else if(mw->delta<0.f) zoom_at(0.80f,(float)mp.x,(float)mp.y);
                }
            }

            if(const auto* mb=event->getIf<sf::Event::MouseButtonPressed>()){
                // 下部バーのボタン
                {
                    float bx=(float)mb->position.x, by=(float)mb->position.y;
                    bool consumed=false;
                    for(const auto& h:ui_hits){
                        if(bx<h.x||bx>h.x+h.w||by<h.y||by>h.y+h.h) continue;
                        if(h.id==100){ view_mode=(view_mode+10)%11; g_last_view=-1; }
                        if(h.id==101){ view_mode=(view_mode+1)%11;  g_last_view=-1; }
                        if(h.id==102){ show_tree=true; }
                        if(h.id==103){ show_tree=false; tree_search.clear(); }
                        if(h.id==104){ show_help=!show_help; }
                        if(h.id>=500&&h.id<600){
                            int row=h.id-500;
                            if(row<(int)save_files.size()){
                                input_text=save_files[row];
                                del_confirm=-1;
                            }
                        }
                        if(h.id>=300&&h.id<400){
                            int row=h.id-300;
                            if(del_confirm==row){
                                // 二度目のクリック。実行する
                                std::string p=save_files[row]+".txt";
                                if(std::remove(p.c_str())==0)
                                    printf("deleted: %s\n",p.c_str());
                                else
                                    printf("delete failed: %s\n",p.c_str());
                                save_info.erase(save_files[row]);
                                save_files.erase(save_files.begin()+row);
                                if(row<(int)save_times.size())
                                    save_times.erase(save_times.begin()+row);
                                del_confirm=-1;
                            } else del_confirm=row;   // 一度目。確認へ
                        }
                        if(h.id>=200&&h.id<210){ graph_scale=h.id-200; graph_off=0; }
                        if(h.id==210){ graph_big=-1; }
                        if(h.id>=400&&h.id<410){ graph_big=h.id-400; graph_off=0; }
                        consumed=true; break;
                    }
                    if(consumed) continue;
                    if(ui_mode==2) del_confirm=-1;   // 別の場所を押したら取り消す
                }
                if(mb->button==sf::Mouse::Button::Left && ui_mode==0){
                    if(show_tree){
                        float bx=(float)mb->position.x, by=(float)mb->position.y;
                        // まずフィルタのチェックボックスを見る
                        bool hit_filter=false;
                        for(const auto& f:filter_hits){
                            if(bx>=f.x&&bx<=f.x+f.w&&by>=f.y&&by<=f.y+f.h){
                                filter_on[f.f]=!filter_on[f.f];
                                hit_filter=true; break;
                            }
                        }
                        if(hit_filter) continue;
                        int best=-1; float bd=14.f;
                        for(const auto& h:tree_hits){
                            if(by<h.y0-6.f||by>h.y1+6.f) continue;
                            float d=std::fabs(bx-h.x);
                            if(d<bd){ bd=d; best=h.li; }
                        }
                        sel_lineage=best;
                    } else {
                        float wx=s2wx((float)mb->position.x), wy=s2wy((float)mb->position.y);
                        sel_animal=-1; sel_plant=-1;
                        // 動物と植物を同じ条件で比べ、近い方を選ぶ。
                        // 動物を先に確定させると、密度の関係で植物に届かない
                        float ba=9e9f, bp=9e9f; int ia=-1, ip=-1;
                        for(const auto& a:animals){
                            if(!a.alive) continue;
                            float dx=a.x-wx,dy=a.y-wy,d=dx*dx+dy*dy;
                            if(d<ba){ ba=d; ia=a.id; }
                        }
                        for(const auto& p:plants){
                            if(!p.alive||p.seed_wait>0) continue;   // 種子は見えない
                            float dx=p.x-wx,dy=p.y-wy,d=dx*dx+dy*dy;
                            if(d<bp){ bp=d; ip=p.id; }
                        }
                        // 拾える範囲は拡大率で決める。引いているときほど広く
                        float pick=18.f/std::max(1.f,zoom); pick*=pick;
                        if(pick<1.5f) pick=1.5f;
                        if(ba<=bp && ba<pick){
                            sel_animal=ia;
                            if(life_target!=ia){ life_target=ia; life_log.clear(); }
                        } else if(bp<pick) sel_plant=ip;
                    }
                }
            }
        }

        if(ui_mode==0 && show_tree){
            float tp=14.f;
            float ttx=0.f, tty=0.f;
            if(key_w) tty+=tp;
            if(key_s) tty-=tp;
            if(key_a) ttx+=tp;
            if(key_d) ttx-=tp;
            tree_vx += (ttx-tree_vx)*0.55f;
            tree_vy += (tty-tree_vy)*0.55f;
            if(ttx==0.f) tree_vx*=0.65f;
            if(tty==0.f) tree_vy*=0.65f;
            tree_ox += tree_vx;
            tree_oy += tree_vy;
            float lx=SCREEN_W*0.5f*(tree_zoom-1.f)+40.f;
            float ly=SCREEN_H*0.5f*(tree_zoom-1.f)+40.f;
            if(tree_ox> lx) tree_ox= lx;  if(tree_ox<-lx) tree_ox=-lx;
            if(tree_oy> ly) tree_oy= ly;  if(tree_oy<-ly) tree_oy=-ly;
        } else if(ui_mode==0){
            // 速度として持つ。キーを離しても滑らかに減速する
            float pan=12.f/zoom*6.f;
            const float CAM_ACCEL=0.55f, CAM_FRICTION=0.65f;
            if(ui_mode==0){
                float tx=0.f, ty=0.f;
                if(key_w) ty-=pan;
                if(key_s) ty+=pan;
                if(key_a) tx-=pan;
                if(key_d) tx+=pan;
                cam_vx += (tx-cam_vx)*CAM_ACCEL;
                cam_vy += (ty-cam_vy)*CAM_ACCEL;
                if(tx==0.f) cam_vx*=CAM_FRICTION;   // 離したら摩擦で止まる
                if(ty==0.f) cam_vy*=CAM_FRICTION;
                if(std::fabs(cam_vx)<0.01f) cam_vx=0.f;
                if(std::fabs(cam_vy)<0.01f) cam_vy=0.f;
                cam_x += cam_vx;
                cam_y += cam_vy;
            }
            clamp_cam();
        }

        if(!paused){
            g_frame=frame;
            prof_begin(); update_climate(frame); update_storms();    prof_end(0);
            prof_begin();
            if(frame%PRECIP_INTERVAL==0){ compute_precip(); classify_climate(); }
                                                                      prof_end(1);
            prof_begin(); update_water(frame); update_groundwater(frame); prof_end(2);
            prof_begin(); update_animals();                           prof_end(3);
            prof_begin();
            if(frame%PATH_INTERVAL==0) update_pathogens();
            update_plants();
            prof_end(4);
            { // 微生物は個体数が桁違いに多い。独立して測る
              auto ta=std::chrono::steady_clock::now();
              update_microbes();
              auto tb=std::chrono::steady_clock::now();
              static double am=0,af=0,ao=0,as=0; static int cm=0;
              am+=std::chrono::duration<double,std::milli>(tb-ta).count();
              auto tc=std::chrono::steady_clock::now();
              update_fire(frame);
              auto td=std::chrono::steady_clock::now();
              af+=std::chrono::duration<double,std::milli>(td-tc).count();
              auto te=std::chrono::steady_clock::now();
              if(frame%ODOR_INTERVAL==0) update_odor();
              auto tf2=std::chrono::steady_clock::now();
              ao+=std::chrono::duration<double,std::milli>(tf2-te).count();
              auto tg=std::chrono::steady_clock::now();
              update_quakes();
              update_slides();
              if(g_replay) ev_replay_advance();
              if(g_rec && frame%600==0) fflush(g_rec);
              auto th2=std::chrono::steady_clock::now();
              as+=std::chrono::duration<double,std::milli>(th2-tg).count();
              if(++cm>=60){
                  printf("  [prof] 微生物 %.1f  火 %.1f  匂い %.1f  崖崩れ %.1f ms (個体 %zu)\n",
                         am/cm,af/cm,ao/cm,as/cm,microbes.size());
                  am=af=ao=as=0; cm=0;
              }
            }
            prof_begin();
            if(frame%GAS_INTERVAL==0) update_gases();   // 拡散は4tickに一度で足りる
            prof_end(5);
            frame++;
            if(g_dump_every>0 && frame%g_dump_every==0) dump_brains(frame);
            if(g_stop_at>0 && frame>=g_stop_at){
                printf("=== stop at t=%d  mean fps (t>=50000) %.2f over %ld frames ===\n",
                       frame, g_fps_time>0? g_fps_n/g_fps_time : 0.0, g_fps_n);
                if(g_replay) printf("replay: done fire %ld slide %ld storm %ld lightning %ld quake %ld of %zu\n",
                                    g_replay_done[EV_FIRE],g_replay_done[EV_SLIDE],g_replay_done[EV_STORM],g_replay_done[EV_LIGHTNING],g_replay_done[EV_QUAKE],g_replay_ev.size());
                if(g_rec){ fclose(g_rec); g_rec=nullptr; }
                field_close();
                census_close();
                fflush(stdout);
                running=false;
                window.close();
            }

            if((int)plants.size() >peak_plant)  peak_plant =(int)plants.size();
            if((int)animals.size()>peak_animal) peak_animal=(int)animals.size();
            if((int)corpses.size()>peak_corpse) peak_corpse=(int)corpses.size();

            // 植物の数だけを待つと、その間に抑制なしで増え続ける。
            // 時間でも切って、暴走する前に消費者を入れる
            if(!released && ((int)plants.size()>ANIMAL_RELEASE || frame>=RELEASE_TICK)){
                int put=0,g2=0, id0=next_animal_id;
                while(put<ANIMAL_INIT_N && g2++<200000){
                    Animal a;
                    a.x=frand(0.f,(float)GRID_W); a.y=frand(0.f,(float)GRID_H);
                    int gx=(int)a.x, gy=(int)a.y;
                    if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
                    if(!sea[idx(gx,gy)]) continue;
        if(depthn(idx(gx,gy))>0.12f) continue;
        if(shelf[idx(gx,gy)]<0.5f) continue;   // 岸のそばに置く
                    a.dir=frand(0.f,6.2832f); a.energy=3.0f; a.alive=true;
                    a.age=0; a.eaten=0;
                    a.speed=0.5f; a.turn=0.2f; a.reach=0.7f;
                    a.size=0.5f; a.diet=0.0f; a.jaw=1.0f; a.cold_tol=0.3f;
                    a.body=0.10f; a.aquatic=1.0f; a.salt_tol=1.0f; a.resp=0.3f; a.detritus=0.0f;
                    a.sexual=0.2f; a.arboreal=0.0f; a.perch=0.0f; a.aerobic=1.0f;
                    a.mut_rate=0.3f; a.bite_size=0.5f; a.digest=0.5f;                    a.mother=-1; a.father=-1; a.gen=0;
                    a.lesion=false;
                    for(int k=0;k<N_REC;k++) a.rec[k]=0.f;
                    for(int k=0;k<N_IMMUNE;k++) a.im[k]=frand(0.f,1.f);
                    a.brood=0.1f; a.lifespan=0.2f;
                    a.toxin=0.f; a.tox_resist=0.f; a.burrow=0.f;
                    a.nocturnal=0.f; a.store_fat=0.2f; a.repro_alloc=0.5f;
                    a.fat=0.f; a.depth=0.f;
                    a.id=next_animal_id++; a.lin=-1;
                    for(int k=0;k<N_W;k++) a.w[k]=0.f;
                    // 動く神経網を種として置く。
                    // 全てゼロだと、器官が使われる前に維持費で淘汰される
                    {
                        // 第1層の3つだけを働かせて始める。深さと幅は進化が決める
                        for(int k=0;k<N_NODE;k++) a.ngain[k]=0.05f;
                        a.ngain[0]=1.0f; a.ngain[1]=1.0f; a.ngain[2]=1.0f;
                        // 節0: 何か見えたか  節1: 常にオン  節2: 空腹か
                        for(int k=0;k<N_SENSOR;k++) a.w[0*N_IN+k]=1.4f;
                        a.w[1*N_IN+(N_SENSOR+4)]=1.2f;
                        a.w[2*N_IN+(N_SENSOR+0)]=-1.5f;
                        a.w[2*N_IN+(N_SENSOR+4)]=0.8f;
                        int o=W_L0+W_LL;
                        a.w[o+0*MAX_NODE+0]=-1.2f;   // 見えたら旋回を減らす
                        a.w[o+1*MAX_NODE+0]=+1.8f;   // 見えたら前へ出る
                        a.w[o+0*MAX_NODE+1]=+0.6f;   // 何もなければ旋回する
                        a.w[o+1*MAX_NODE+1]=-0.3f;
                        a.w[o+1*MAX_NODE+2]=+0.7f;   // 空腹なら動く
                        // 採食と攻撃は、最初は常に行う。
                        // そこから「いつ食うか」を進化が削っていく
                        int od=W_L0+W_LL+W_OUT;
                        a.w[od+2*N_IN+(N_SENSOR+4)]=+1.5f;   // eat ← バイアス
                        a.w[od+3*N_IN+(N_SENSOR+4)]=+0.8f;   // attack ← バイアス
                        a.w[od+2*N_IN+(N_SENSOR+0)]=-1.0f;   // 満腹なら食わない
                    }
                    // 何に向いた目を持つかは、はじめから散らしておく。
                    // 壁から始めると変異では谷を渡れない
                    a.sen[0]={0.0f,6.0f,0.02f};                  // 前方に植物を見る目
                    a.sen[1]={frand(-0.9f,0.9f),5.0f,frand(0.f,1.f)};
                    for(int k=2;k<N_SENSOR;k++)
                        a.sen[k]={frand(-3.1416f,3.1416f),
                                  (dist01(rng)<0.55f)?frand(2.f,9.f):0.f,
                                  frand(0.f,1.f)};
                    animals.push_back(a); put++;
                    g_bud.c_immig+=(a.energy+a.fat+a.body*CORPSE_RETURN)*CO2_PER_ENERGY;   // 外から放した動物が持ち込む炭素
                }
                released=true;
                printf("=== animals released at t=%d ===\n",frame);
                // 対照の動物(測定の道具 4)。放した順に先頭から、乱数を使わずに選ぶ。
                // 陽性対照は試験器と同じ手で組んだ神経網と感覚器に置き換え、陰性対照は再帰の入力を切る
                if(g_pc_n>0||g_lesion_n>0){
                    Animal pc=assay_positive_control();
                    int npc=0, nls=0;
                    for(auto& a:animals){
                        if(a.id<id0) continue;
                        if(npc<g_pc_n){
                            for(int k=0;k<N_W;k++) a.w[k]=pc.w[k];
                            for(int k=0;k<N_NODE;k++) a.ngain[k]=pc.ngain[k];
                            for(int k=0;k<N_SENSOR;k++) a.sen[k]=pc.sen[k];
                            npc++;
                        } else if(nls<g_lesion_n){ a.lesion=true; nls++; }
                    }
                    printf("controls: positive %d  no-recurrent %d\n",npc,nls);
                }
            }

            prof_begin();
            if(frame%LINEAGE_INTERVAL==0){
                update_lineages(0,frame);
                update_lineages(1,frame);
                if(!animals.empty()) remember_animals();
            }
            prof_end(6);

            if(frame%SAMPLE_INTERVAL==0){
                double a1=0,a2=0,a3=0,a4=0;
                for(const auto& p:plants){ a1+=p.max_height;a2+=p.shade_tol;a3+=p.tough;a4+=p.cold_tol; }
                int np=(int)plants.size();
                h_maxH .push_back(np?(float)(a1/np):0.f);
                h_shade.push_back(np?(float)(a2/np):0.f);
                h_tough.push_back(np?(float)(a3/np):0.f);
                h_coldP.push_back(np?(float)(a4/np):0.f);
                h_plant.push_back(np);
                double b1=0,b2=0,b3=0,b4=0;
                for(const auto& a:animals){ b1+=a.reach;b2+=a.diet;b3+=a.jaw;b4+=a.cold_tol; }
                int na=(int)animals.size();
                h_reach.push_back(na?(float)(b1/na):0.f);
                h_diet .push_back(na?(float)(b2/na):0.f);
                h_jaw  .push_back(na?(float)(b3/na):0.f);
                h_coldA.push_back(na?(float)(b4/na):0.f);
                h_animal.push_back(na);
                h_corpse.push_back((int)corpses.size());
                double c1=0,c2=0,c3=0; int land=0;
                for(int i=0;i<GRID_W*GRID_H;i++){
                    c1+=co2[i]; c2+=o2[i];
                    if(!sea[i]){ c3+=temper[i]; land++; }
                }
                h_co2.push_back((float)(c1/(GRID_W*GRID_H)));
            s_plant.push((float)plants.size());
            s_animal.push((float)animals.size());
            s_co2.push(h_co2.empty()?0.f:h_co2.back());
            s_o2 .push(h_o2.empty()?0.f:h_o2.back());
            if(!h_maxH.empty())  s_maxh.push(h_maxH.back());
            if(!h_reach.empty()) s_reach.push(h_reach.back());
            if(!h_temp.empty())  s_temp.push(h_temp.back());
            if(!h_soil.empty())  s_soil.push(h_soil.back());
            if(!h_nut.empty())   s_nut.push(h_nut.back());
            if(!h_shade.empty()) s_shade.push(h_shade.back());
            if(!h_diet.empty())  s_diet.push(h_diet.back());
            s_corpse.push((float)corpses.size());
            s_micro.push((float)microbes.size());
            if(!h_tough.empty()) s_tough.push(h_tough.back());
            if(!h_jaw.empty())   s_jaw.push(h_jaw.back());
            { double sd=0; int nd=0;
              for(const auto& p:plants){ if(!p.alive||p.seed_wait>0) continue;
                  sd+=p.defense; nd++; }
              s_defense.push(nd?(float)(sd/nd):0.f); }
            if(!h_coldP.empty()) s_coldP.push(h_coldP.back());
            if(!h_coldA.empty()) s_coldA.push(h_coldA.back());
            if(!h_coldP.empty()) s_coldP.push(h_coldP.back());
            if(!h_coldA.empty()) s_coldA.push(h_coldA.back());
            { double sw=0,sn=0,sp=0,so=0,ss=0; int nl=0,ni=0,nse=0;
              for(int i=0;i<GRID_W*GRID_H;i++){
                  if(sea[i]){ ss+=sst_anom[i]; nse++; if(ice[i]) ni++; }
                  else { sw+=soil_water[i]; sn+=soil_n[i]; nl++; }
                  sp+=pload[i]; so+=odor[0][i]+odor[1][i]+odor[2][i];
              }
              h_soil.push_back(nl?(float)(sw/nl):0.f);
              h_nut .push_back(nl?(float)(sn/nl):0.f);
              h_ice .push_back((float)ni/(GRID_W*GRID_H));
              h_path.push_back((float)(sp/(GRID_W*GRID_H)));
              h_odor.push_back((float)(so/(GRID_W*GRID_H)));
              h_sst .push_back(nse?(float)(ss/nse):0.f);
              auto cap=[](std::vector<float>& v){
                  while(v.size()>h_co2.size()) v.erase(v.begin()); };
              cap(h_soil);cap(h_nut);cap(h_ice);cap(h_path);cap(h_odor);cap(h_sst);
            }
                h_o2 .push_back((float)(c2/(GRID_W*GRID_H)));
                h_temp.push_back(land?(float)(c3/land):0.f);
                if((int)h_plant.size()>HIST_MAX){
                    h_maxH.erase(h_maxH.begin());   h_shade.erase(h_shade.begin());
                    h_tough.erase(h_tough.begin()); h_coldP.erase(h_coldP.begin());
                    h_plant.erase(h_plant.begin()); h_reach.erase(h_reach.begin());
                    h_diet.erase(h_diet.begin());   h_jaw.erase(h_jaw.begin());
                    h_coldA.erase(h_coldA.begin()); h_animal.erase(h_animal.begin());
                    h_corpse.erase(h_corpse.begin());
                    h_co2.erase(h_co2.begin());     h_o2.erase(h_o2.begin());
                    h_temp.erase(h_temp.begin());
                }
            }

            if(frame%AUTOSAVE_EVERY==0){
                char fn[32]; snprintf(fn,32,"autosave%d.txt",autosave_slot);
                save_state(fn,frame); autosave_slot=(autosave_slot+1)%3;
            }

            if(g_census_f && frame%g_census_every==0) census_write(frame);

            if(frame%PRINT_EVERY==0){
                int n=(int)plants.size();
                double sa=0,sh=0,st=0,stg=0,scl=0,sdp=0,scd=0;
                for(const auto& p:plants){ sa+=p.absorb;sh+=p.max_height;st+=p.shade_tol;
                    stg+=p.tough;scl+=p.clonal;sdp+=p.disperse;scd+=p.cold_tol; }
                scan_forms();
                double saq=0,sdr=0,sdt=0,sst=0; int ndorm=0;
                for(const auto& p:plants){ saq+=p.aquatic; sdr+=p.drought_tol;
                    sdt+=p.dorm_temp; sst+=p.store_cap; if(p.dormant) ndorm++; }
                double sw=0; int land=0; int kc[9]={0,0,0,0,0,0,0,0,0};
                for(int i=0;i<GRID_W*GRID_H;i++) if(!sea[i]){
                    sw+=soil_water[i]; land++; kc[climate[i]]++; }
                int n_seed=0; for(const auto& p:plants) if(p.seed_wait>0) n_seed++;
                // 合計に種子が含まれているので、割る数も全数のままにする
                // (待機中の数は waiting として別に出している)
                printf("t:%d P:%d absorb:%.3f maxH:%.3f shade:%.3f tough:%.3f clonal:%.3f disp:%.3f cold:%.3f aq:%.3f dry:%.3f | soil:%.3f\n",
                       frame,n,n?sa/n:0,n?sh/n:0,n?st/n:0,n?stg/n:0,n?scl/n:0,n?sdp/n:0,
                       n?scd/n:0,n?saq/n:0,n?sdr/n:0,land?sw/land:0);
                if(frame%BUDGET_EVERY==0){
                    // 在庫。生体の炭素はエネルギー×CO2_PER_ENERGY(光合成の換算と同じ)で数える。
                    // 植物の体(高さ)は建設に GROW_COST を払ったぶんの炭素を持つとみなす
                    double ep=0,ea=0,ec=0,em=0,np=0,nc=0,ns=0,nl=0,wl=0,wp=0,wg=0,sl=0,sg=0;
                    // 植物の窒素 = 貯え + 体(高さ)を作るのに使った N_PER_GROWTH ぶん
                    for(const auto& p:plants){ ep+=p.energy+p.height*GROW_COST; np+=p.n_store+p.height*N_PER_GROWTH; wp+=p.water_store; }
                    // 動物の体の組織は CORPSE_RETURN ぶんの炭素を持つ
                    for(const auto& a:animals) ea+=a.energy+a.fat+a.body*CORPSE_RETURN;
                    for(const auto& c:corpses){ ec+=c.energy; nc+=c.n; }
                    for(const auto& m:microbes) if(m.alive) em+=m.energy;
                    for(int i=0;i<GRID_W*GRID_H;i++){
                        if(sea[i]) ns+=soil_n[i]; else { nl+=soil_n[i]; wl+=soil_water[i]; wg+=gw_h[i]*GW_W; sl+=salt_soil[i]; sg+=salt_gw[i]; } }
                    double catm=grid_sum(co2)+grid_sum(co2_res), oatm=grid_sum(o2)+grid_sum(o2_res);   // 持ち越しの端数も大気の在庫
                    double cbio=(ep+ea+ec+em)*CO2_PER_ENERGY;
                    double ctot=catm+cbio, ntot=np+nc+ns+nl;
                    for(int i=0;i<GRID_W*GRID_H;i++){ g_bud.w_transp+=w_transp_cell[i]; w_transp_cell[i]=0.0; }
                    double wtot=wl+wp+wg;
                    double stot=sl+sg;
                    static bool have=false; static double pc=0,po=0,pn=0,pw=0,pcb=0,pst=0;
                    if(have){
                        double cext=g_bud.c_relax+g_bud.c_edge+g_bud.c_trans+g_bud.c_immig;
                        double next=g_bud.n_sink+g_bud.n_runoff+g_bud.n_upwell+g_bud.n_clamp+g_bud.n_fix+g_bud.n_fire+g_bud.n_carn+g_bud.n_immig;
                        printf("      budget C: atm(dev) %.1f bio %.1f (P %.1f A %.1f dead %.1f mic %.1f) | dC %+.2f = ext %+.2f (relax %+.2f edge %+.2f trans %+.3f immig %+.2f) + internal %+.2f\n",
                               catm,cbio,ep*CO2_PER_ENERGY,ea*CO2_PER_ENERGY,ec*CO2_PER_ENERGY,em*CO2_PER_ENERGY,
                               ctot-pc,cext,g_bud.c_relax,g_bud.c_edge,g_bud.c_trans,g_bud.c_immig,ctot-pc-cext);
                        // O2 は単独では保存しない。光合成で有機物の炭素 1 につき O2 を 1 出し、呼吸・分解・燃焼で 1:1 で使う。
                        // 内部の O2 の増減は生体の炭素の増減と等しいはずで、その差(redox)が破れ。
                        // 外から入った生き物の炭素は、この土地の光合成で作ったものではないので除く
                        double oint=oatm-po-g_bud.o_relax-g_bud.o_edge-g_bud.o_trans;
                        printf("      budget O2: atm(dev) %.1f | dO2 %+.2f = relax %+.2f edge %+.2f trans %+.3f + internal %+.2f (dCbio %+.2f, redox %+.2f)\n",
                               oatm,oatm-po,g_bud.o_relax,g_bud.o_edge,g_bud.o_trans,oint,cbio-pcb-g_bud.c_immig,oint-(cbio-pcb-g_bud.c_immig));
                        printf("      budget N: soil land %.1f sea %.1f plant %.1f dead %.1f | dN %+.2f = sink %+.2f runoff %+.2f upwell %+.2f clamp %+.3f fix %+.2f fire %+.3f carn %+.3f immig %+.3f + internal %+.2f\n",
                               nl,ns,np,nc,ntot-pn,g_bud.n_sink,g_bud.n_runoff,g_bud.n_upwell,g_bud.n_clamp,g_bud.n_fix,g_bud.n_fire,g_bud.n_carn,g_bud.n_immig,ntot-pn-next);
                        double wext=g_bud.w_relax+g_bud.w_transp+g_bud.w_sea+g_bud.w_storm
                                   +g_bud.w_gw_sea+g_bud.w_gw_edge+g_bud.w_gw_out+g_bud.w_gw_clamp;
                        printf("      budget W: soil land %.1f plant %.1f ground %.1f | dW %+.2f = relax %+.2f transp %+.2f sea %+.3f storm %+.2f gwSea %+.2f gwEdge %+.2f gwOut %+.2f gwClamp %+.3f + internal %+.2f\n",
                               wl,wp,wg,wtot-pw,g_bud.w_relax,g_bud.w_transp,g_bud.w_sea,g_bud.w_storm,
                               g_bud.w_gw_sea,g_bud.w_gw_edge,g_bud.w_gw_out,g_bud.w_gw_clamp,wtot-pw-wext);
                        // 塩の単位は g/kg × 土の水の単位。土の水 1.0 は 100 mm なので、1 単位は 1 m² あたり 100 g(マス 1 つで 10 kg)
                        double sext=g_bud.s_sea+g_bud.s_edge+g_bud.s_out+g_bud.s_clamp;
                        printf("      budget S: soil %.1f ground %.1f | dS %+.2f = sea %+.2f edge %+.2f out %+.2f clamp %+.3f + internal %+.2f\n",
                               sl,sg,stot-pst,g_bud.s_sea,g_bud.s_edge,g_bud.s_out,g_bud.s_clamp,stot-pst-sext);
                    }
                    have=true; pc=ctot; po=oatm; pn=ntot; pw=wtot; pcb=cbio; pst=stot;
                    g_bud=Budget();
                }
                { double sfi=0,scl2=0,sx=0;
                  for(const auto& p:plants){ sfi+=p.fire_tol; scl2+=p.climb; sx+=p.sexual; }
                  double asx=0; for(const auto& a:animals) asx+=a.sexual;
                  // 燃料の湿りの陸の平均と、燃えうる(消炎含水率 FIRE_MX を下回り、燃料が下限を超える)陸の割合
                  double fm=0; long nl=0, nburnable=0;
                  for(int i=0;i<GRID_W*GRID_H;i++){ if(sea[i]||fuel_moist[i]<0.f) continue;
                      fm+=fuel_moist[i]; nl++; if(fuel_moist[i]<FIRE_MX && fuel[i]>=FIRE_FUEL_MIN) nburnable++; }
                  printf("      fire:%.0f cells (%ld ignitions) | fuelM %.3f burnable %.1f%% | fireTol:%.3f climb:%.3f | sexual P:%.3f A:%.3f\n",
                         fire_area,fire_events,nl?fm/nl:0.0,nl?100.0*nburnable/nl:0.0,n?sfi/n:0.0,n?scl2/n:0.0,
                         n?sx/n:0.0,(int)animals.size()?asx/animals.size():0.0); }
                { double sab=0,spe=0; int nup=0;
                  for(const auto& a:animals){ sab+=a.arboreal; spe+=a.perch;
                                              if(a.perch>0.3f) nup++; }
                  int na2=(int)animals.size();
                  printf("      arboreal:%.3f  in canopy:%.1f%%  mean perch:%.2f\n",
                         na2?sab/na2:0.0, na2?100.0*nup/na2:0.0, na2?spe/na2:0.0); }
                { double srs=0,sae=0;
                  for(const auto& a:animals){ srs+=a.resp; sae+=a.aerobic; }
                  int na3=(int)animals.size();
                  double o2m2=0,tot2=0;
                  for(int i=0;i<GRID_W*GRID_H;i++){ o2m2+=o2_at(i); tot2+=o2_at(i)+co2_at(i)+N2_COL; }
                  printf("      O2 partial:%.1f%% (Earth 21%%) | resp:%.3f  aerobic:%.0f%%\n",
                         tot2>0?100.0*o2m2/tot2:0.0, na3?srs/na3:0.0,
                         na3?100.0*sae/na3:0.0); }
                { double sbr=0,snf=0,sls=0,sag=0;
                  for(const auto& p:plants){ sbr+=p.brood; snf+=p.n_fix;
                                             sls+=p.lifespan; sag+=p.age; }
                  double abr=0,als=0,aag=0;
                  for(const auto& a:animals){ abr+=a.brood; als+=a.lifespan; aag+=a.age; }
                  int na4=(int)animals.size();
                  double sn=0; int ld2=0;
                  for(int i=0;i<GRID_W*GRID_H;i++) if(!sea[i]){ sn+=soil_n[i]; ld2++; }
                  double spa=0,sca=0,sfs=0,sra=0;
                  for(const auto& p:plants){ spa+=p.parasite; sca+=p.carnivory;
                                             sfs+=p.fire_seed; sra+=p.repro_alloc; }
                  double atx=0,atr=0,abu=0,anc=0,afa=0,ara=0;
                  for(const auto& a:animals){ atx+=a.toxin; atr+=a.tox_resist;
                      abu+=a.burrow; anc+=a.nocturnal; afa+=a.fat; ara+=a.repro_alloc; }
                  int na5=(int)animals.size();
                  if(dbg_n>0){
                    // tick当たりに揃える。収入は日照時のみ、支出はPLANT_STRIDE毎
                    double g_t=dbg_gain/(double)PRINT_EVERY/std::max(1,n);
                    double c_t=dbg_cost/(double)PRINT_EVERY*PLANT_STRIDE/std::max(1,n);
                    printf("      BALANCE/tick  gain %.4f  cost %.4f  net %+.4f | light %.3f tf %.2f wf %.2f co2 %.2f\n",
                           g_t, c_t, g_t-c_t,
                           dbg_light/dbg_n, dbg_tf/dbg_n, dbg_wf/dbg_n, dbg_co2/dbg_n);
                  }
                  dbg_gain=dbg_cost=0; dbg_n=0;
                  dbg_light=dbg_tf=dbg_wf=dbg_co2=0;
                  { double sdf=0,sbk=0; int nseed=0, nice=0;
                    for(const auto& p:plants){ sdf+=p.defense; sbk+=p.seed_bank;
                                               if(p.seed_wait>0) nseed++; }
                    for(int i=0;i<GRID_W*GRID_H;i++) if(ice[i]) nice++;
                    printf("      defense %.3f (tough %.3f) | seedBank %.3f waiting %d | ice %d cells | slides %ld | quakes %ld\n",
                           n?sdf/n:0.0, n?stg/n:0.0, n?sbk/n:0.0, nseed, nice,
                           slide_events, quake_events); }
                  { int nl=0, nsat=0, ngs=0, nss=0; float cold=0.f; double sh=0, smax=0;
                    for(int i=0;i<GRID_W*GRID_H;i++){ if(sea[i]) continue; nl++; cold=std::min(cold,storm_dT[i]);
                        sh+=gw_h[i]; if(gw_h[i]>=GW_DEPTH_M*0.999f) nsat++;
                        if(gw_h[i]>0.f && salt_gw[i]/(gw_h[i]*GW_W)>1.f) ngs++;   // 1 g/kg を超えると淡水ではない(汽水)
                        float ss=soil_salinity(i); if(ss>1.f) nss++; smax=std::max(smax,(double)ss); }
                    printf("      groundwater: saturated thickness land mean %.3f m | water table at surface %d cells\n",
                           nl?sh/nl:0.0, nsat);
                    { double stp=0, sta=0;
                      for(const auto& p:plants) stp+=p.salt_tol;
                      for(const auto& a:animals) sta+=a.salt_tol;
                      printf("      salt: brackish groundwater %d cells | saline root zone %d cells (max %.1f g/kg) | salt_tol P:%.3f A:%.3f\n",
                           ngs, nss, smax, plants.empty()?0.0:stp/plants.size(), animals.empty()?0.0:sta/animals.size()); }
                    printf("      weather: storms %ld (active %zu) | lightning %ld | rain land mean %.2f mm since last | cold pool now %.1f C | gust max %.1f m/s since last\n",
                           storm_events, storms.size(), lightning_events, nl?storm_rain_mm/nl:0.0, cold, gust_max);
                    storm_rain_mm=0; gust_max=0.f; }
                         printf("      plant: parasite %.3f carnivory %.3f fireSeed %.3f reproAlloc %.3f\n",
                                 n?spa/n:0.0,
                                 n?sca/n:0.0,
                                 n?sfs/n:0.0,
                                 n?sra/n:0.0);
                  { int sm[N_ODOR]={0,0,0,0}, sv=0, nkin=0;
                    for(const auto& a:animals)
                      for(int k=0;k<N_SENSOR;k++){
                        if(a.sen[k].range<SENSOR_MIN) continue;
                        int bd=(int)(a.sen[k].tune*N_TUNE);
                        if(bd>N_TUNE-1) bd=N_TUNE-1;
                        if(bd>=7&&bd<=10) sm[bd-7]++;
                        else if(bd==11) nkin++;
                        else sv++;
                      }
                    double bw=0, bd_=0; int na6=(int)animals.size();
                    int lw[MAX_LAYER]={0};
                    for(const auto& a:animals){
                        int depth=0;
                        for(int L=0;L<MAX_LAYER;L++){
                            int act=0;
                            for(int h=0;h<MAX_NODE;h++)
                                if(a.ngain[L*MAX_NODE+h]>NODE_ON) act++;
                            lw[L]+=act;
                            if(act) depth=L+1;
                            bw+=act;
                        }
                        bd_+=depth;
                    }
                    printf("      odor field: bio %.3f rot %.3f wound %.3f smoke %.3f | sensors sight:%d smell:%d/%d/%d/%d kin:%d | brain %.1f/%d\n",
                           odor_total[0],odor_total[1],odor_total[2],odor_total[3],
                           sv,sm[0],sm[1],sm[2],sm[3],nkin,
                           na6?bw/na6:0.0,N_NODE);
                    printf("      brain: depth %.2f  width %.1f/%.1f/%.1f (max %d layers x %d)\n",
                           na6?bd_/na6:0.0,
                           na6?(double)lw[0]/na6:0.0, na6?(double)lw[1]/na6:0.0,
                           na6?(double)lw[2]/na6:0.0, MAX_LAYER, MAX_NODE); }
                { int pd=0,pn=0;
                  for(int i=0;i<GRID_W*GRID_H;i++){
                      if(sea[i]) continue;
                      if(daylit[i]>0.5f) pd++; else pn++; }
                  // 極夜・白夜の判定。太陽が一日中沈まない/昇らない緯度
                  float circ=90.f-AXIAL_TILT;
                  float dd=sun_dec*180.f/3.14159265f;
                  printf("      sun dec %+.1f deg | land in day %d / night %d | polar day: %s\n",
                         dd, pd, pn,
                         dd>0.5f?"north":(dd<-0.5f?"south":"none")); }
                  { double sbs=0,sdg=0; int na7=(int)animals.size();
                    for(const auto& a:animals){ sbs+=a.bite_size; sdg+=a.digest; }
                    printf("      feeding: bite %.3f  digest %.3f\n",
                           na7?sbs/na7:0.0, na7?sdg/na7:0.0); }
                  printf("      animal: toxin %.3f resist %.3f burrow %.3f nocturn %.3f fat %.2f repro %.2f | day %.2f\n",
                         na5?atx/na5:0.0,na5?atr/na5:0.0,na5?abu/na5:0.0,
                         na5?anc/na5:0.0,na5?afa/na5:0.0,na5?ara/na5:0.0,day_light);
                  printf("      brood P:%.2f A:%.2f | lifespan P:%.2f A:%.2f | age P:%.0f A:%.0f | nFix:%.3f soilN:%.2f\n",
                         n?sbr/n:0.0, na4?abr/na4:0.0,
                         n?sls/n:0.0, na4?als/na4:0.0,
                         n?sag/n:0.0, na4?aag/na4:0.0,
                         n?snf/n:0.0, ld2?sn/ld2:0.0);
                  { double ns=0,nl=0; int cs=0,cl=0;
                    for(int i=0;i<GRID_W*GRID_H;i++){
                        if(sea[i]){ ns+=soil_n[i]; cs++; } else { nl+=soil_n[i]; cl++; } }
                { double amin=1e9,amax=-1e9,aabs=0; int nc=0;
                  float wmax=0.f;
                  for(int i=0;i<GRID_W*GRID_H;i++){
                      if(!sea[i]) continue;
                      double v=sst_anom[i];
                      if(v<amin)amin=v; if(v>amax)amax=v;
                      aabs+=std::fabs(v); nc++;
                      float sp=std::sqrt(cur_u[i]*cur_u[i]+cur_v[i]*cur_v[i]);
                      if(sp>wmax) wmax=sp;
                  }
                  printf("      current: sst %.2f .. %+.2f (mean|%.2f|) maxflow %.2f\n",
                         nc?amin:0.0, nc?amax:0.0, nc?aabs/nc:0.0, wmax); }
                    printf("      soilN  sea %.3f  land %.3f\n",
                           cs?ns/cs:0.0, cl?nl/cl:0.0); }
                }
                { double lp=0,la=0; int cp2=0,ca2=0;
                  for(int i=0;i<GRID_W*GRID_H;i++){
                      if(pload[i]>0.02f){ lp+=pload[i]; cp2++; }
                      if(aload[i]>0.02f){ la+=aload[i]; ca2++; } }
                  // 免疫型の多様さ: 集団の分散を見る
                  double vim=0;
                  if(n>1){
                      for(int k=0;k<N_IMMUNE;k++){
                          double m=0; for(const auto& p:plants) m+=p.im[k]; m/=n;
                          double v=0; for(const auto& p:plants){ double d=p.im[k]-m; v+=d*d; }
                          vim+=v/n;
                      }
                      vim=std::sqrt(vim/N_IMMUNE);
                  }
                  printf("      pathogen P:%.2f(%d cells) A:%.2f(%d cells) | immune sd:%.3f\n",
                         cp2?lp/cp2:0.0,cp2, ca2?la/ca2:0.0,ca2, vim); }
                printf("      dormT:%.1f dormant:%.0f%% store:%.2f\n",
                       n?sdt/n:0.0, n?100.0*ndorm/n:0.0, n?sst/n:0.0);
                if(land) printf("      climate: EF%.0f%% ET%.0f%% BW%.0f%% BS%.0f%% D%.0f%% C%.0f%% Aw%.0f%% Af%.0f%%\n",
                       100.0*kc[1]/land,100.0*kc[2]/land,100.0*kc[3]/land,100.0*kc[4]/land,
                       100.0*kc[5]/land,100.0*kc[6]/land,100.0*kc[7]/land,100.0*kc[8]/land);
                int an=(int)animals.size();
                double p1=0,p2=0,p3=0,p4=0,p5=0,p6=0;
                for(const auto& a:animals){ p1+=a.speed;p2+=a.reach;p3+=a.body;
                    p4+=a.diet;p5+=a.jaw;p6+=a.cold_tol; }
                double aaq=0,adt=0,asr=0,asn=0; int cp0=0, sp_t[3]={0,0,0};
                for(const auto& a:animals){
                    aaq+=a.aquatic; adt+=a.detritus;
                    for(int k=0;k<N_SENSOR;k++) if(a.sen[k].range>=SENSOR_MIN){
                        asn+=1; asr+=a.sen[k].range;
                        sp_t[a.sen[k].tune<0.34f?0:(a.sen[k].tune<0.67f?1:2)]++;
                    }
                }
                for(const auto& c:corpses) if(c.origin==0) cp0++;
                double c1=0,c2=0;
                for(int i=0;i<GRID_W*GRID_H;i++){ c1+=co2[i]; c2+=o2[i]; }
                int alp=0,ala=0;
                for(const auto& L:lineages) if(L.alive){ if(L.kind==0) alp++; else ala++; }
                printf("      A:%d speed:%.3f reach:%.3f body:%.3f diet:%.3f detr:%.3f jaw:%.3f cold:%.3f aq:%.3f | CO2dev:%.4f O2dev:%.4f corpse:%d(plant %d) | sp P%d A%d\n",
                       an,an?p1/an:0,an?p2/an:0,an?p3/an:0,an?p4/an:0,an?adt/an:0,
                       an?p5/an:0,an?p6/an:0,an?aaq/an:0,
                       c1/(GRID_W*GRID_H),c2/(GRID_W*GRID_H),
                       (int)corpses.size(),cp0,alp,ala);
                printf("      sensors/indiv:%.2f avg range:%.1f | tuned to plant:%d animal:%d corpse:%d\n",
                       an?asn/an:0.0, asn?asr/asn:0.0, sp_t[0],sp_t[1],sp_t[2]);
                // 繁殖した親と集団全体を比べる。正なら押し上げ、負なら削られている
                printf("      [debug] repro events since last: A=%ld P=%ld\n",
                       sg_a_nr,sg_p_nr);
                if(sg_a_nr>20){
                    for(auto& a:animals){ if(!a.alive) continue;
                        double v[N_SG_A]; sg_vec_a(a,v);
                        for(int k=0;k<N_SG_A;k++) sg_a_pop[k]+=v[k];
                        sg_a_np++; }
                    printf("      selection A:");
                    for(int k=0;k<N_SG_A;k++){
                        double d=sg_a_rep[k]/sg_a_nr - sg_a_pop[k]/std::max(1L,sg_a_np);
                        printf(" %s%+.4f",SG_A_NAME[k],d);
                    }
                    printf("\n");
                    for(int k=0;k<N_SG_A;k++){ sg_a_rep[k]=0; sg_a_pop[k]=0; }
                    sg_a_nr=0; sg_a_np=0;
                }
                if(sg_p_nr>20){
                    for(auto& p:plants){ if(!p.alive||p.seed_wait>0) continue;
                        double v[N_SG_P]; sg_vec_p(p,v);
                        for(int k=0;k<N_SG_P;k++) sg_p_pop[k]+=v[k];
                        sg_p_np++; }
                    printf("      selection P:");
                    for(int k=0;k<N_SG_P;k++){
                        double d=sg_p_rep[k]/sg_p_nr - sg_p_pop[k]/std::max(1L,sg_p_np);
                        printf(" %s%+.4f",SG_P_NAME[k],d);
                    }
                    printf("\n");
                    for(int k=0;k<N_SG_P;k++){ sg_p_rep[k]=0; sg_p_pop[k]=0; }
                    sg_p_nr=0; sg_p_np=0;
                }
                { double gp=0,ga=0; long sex_p=0,sex_a=0;
                  for(const auto& p:plants){ if(!p.alive) continue;
                      gp+=p.gen; if(p.father>=0) sex_p++; }
                  for(const auto& a:animals){ if(!a.alive) continue;
                      ga+=a.gen; if(a.father>=0) sex_a++; }
                  int np2=n, na2=(int)animals.size();
                  printf("      genealogy: gen P%.0f A%.0f | 有性由来 P%.0f%% A%.0f%%\n",
                         np2?gp/np2:0.0, na2?ga/na2:0.0,
                         np2?100.0*sex_p/np2:0.0, na2?100.0*sex_a/na2:0.0); }
                { double mr2=0,mc2=0; int nm=(int)microbes.size();
                  for(const auto& m:microbes){ mr2+=m.rate; mc2+=m.cold_tol; }
                  printf("      microbe: %d  rate %.3f  cold %.3f  deaths %ld  capped %.1f%%  colonized %ld\n",
                         nm, nm?mr2/nm:0.0, nm?mc2/nm:0.0, d_micro, g_mic_upd?100.0*g_mic_cap/g_mic_upd:0.0, g_mic_col);
                  g_mic_upd=g_mic_cap=g_mic_col=0; }
                compute_averages();   // レーダー図に重ねる集団平均
                printf("      deaths: starve %ld  cold %ld  preyed %ld  slide %ld\n",
                       d_starve,d_cold,d_preyed,d_slide);
            }
        }

        // ===== 背景 =====
        // ここから window.display() までは描画だけ。--headless では丸ごと飛ばす
        if(!g_headless){
        ui_hits.clear();   // 当たり判定はフレームの先頭で作り直す
        prof_begin();
        if(frame%2==0 || view_mode!=g_last_view || paused){
        g_last_view=view_mode;
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y); sf::Color c;
            if(view_mode==1){
                float t=clamp01((temper[i]+30.f)/70.f);
                c=sf::Color((std::uint8_t)(255*t),
                            (std::uint8_t)(60+80*(1-std::fabs(t-0.5f)*2)),
                            (std::uint8_t)(255*(1-t)));
                if(sea[i]) c=sf::Color((std::uint8_t)(c.r*0.55f),(std::uint8_t)(c.g*0.55f),
                                       (std::uint8_t)(c.b*0.55f));
            } else if(view_mode==2){
                float ca=clamp01(co2_at(i)/(CO2_COL*1.5f));
                float cb=clamp01(o2_at(i)/(O2_COL*1.5f));
                c=sf::Color((std::uint8_t)(40+215*ca),(std::uint8_t)(40+215*cb),
                            (std::uint8_t)(50+50*std::min(ca,cb)));
            } else if(view_mode==3){
                float d=clamp01(corpse_density[i]/6.f);
                if(sea[i]) c=sf::Color((std::uint8_t)(12+90*d),(std::uint8_t)(26+60*d),
                                       (std::uint8_t)(50+20*d));
                else       c=sf::Color((std::uint8_t)(22+170*d),(std::uint8_t)(20+105*d),
                                       (std::uint8_t)(18+50*d));
            } else if(view_mode==4){
                if(sea[i]) c=sf::Color(16,34,66);
                else switch(climate[i]){
                    case 1: c=sf::Color(206,224,244); break;  // EF 氷雪
                    case 2: c=sf::Color(158,140,172); break;  // ET ツンドラ
                    case 3: c=sf::Color(232,128,56);  break;  // BW 砂漠
                    case 4: c=sf::Color(236,202,110); break;  // BS ステップ
                    case 5: c=sf::Color(118,150,202); break;  // D 亜寒帯
                    case 6: c=sf::Color(108,200,108); break;  // C 温帯
                    case 7: c=sf::Color(202,220,92);  break;  // Aw サバナ
                    default:c=sf::Color(28,130,90);   break;  // Af 熱帯雨林
                }
            } else if(view_mode==10){
                if(!sea[i]){
                    float h=clamp01((elev[i]-SEA_LEVEL)/(1.f-SEA_LEVEL));
                    c=sf::Color((std::uint8_t)(40+50*h),(std::uint8_t)(40+50*h),
                                (std::uint8_t)(38+46*h));
                } else {
                    float a=clamp01((sst_anom[i]+CURRENT_HEAT)/(2.f*CURRENT_HEAT));
                    c=sf::Color((std::uint8_t)(20+230*a),(std::uint8_t)(40+80*(1-std::fabs(a-0.5f)*2)),
                                (std::uint8_t)(20+230*(1-a)));
                }
            } else if(view_mode==9){
                // 傾斜(m/m)。崖崩れの起きうる 30° 以上を赤で示す
                if(sea[i]) c=sf::Color(16,34,66);
                else {
                    float t=clamp01(slope[i]/REPOSE_SLOPE);
                    if(slope[i]>=SLIDE_SLOPE) c=sf::Color(230,70,50);
                    else c=sf::Color((std::uint8_t)(40+160*t),(std::uint8_t)(60+150*t),
                                     (std::uint8_t)(40+60*t));
                }
            } else if(view_mode==8){
                float nv=clamp01(soil_n[i]/1.2f);
                if(sea[i]){
                    if(upwell[i]) c=sf::Color(40,(std::uint8_t)(90+150*nv),160);
                    else c=sf::Color(12,(std::uint8_t)(24+120*nv),(std::uint8_t)(50+40*nv));
                } else c=sf::Color((std::uint8_t)(40+60*nv),
                                   (std::uint8_t)(50+180*nv),
                                   (std::uint8_t)(30+50*nv));
            } else if(view_mode==7){
                float d=daylit[i];
                float t=clamp01(solar[i]/SOLAR_GAIN);
                if(sea[i])
                    c=sf::Color((std::uint8_t)(10+40*d),(std::uint8_t)(20+90*d),
                                (std::uint8_t)(40+110*d));
                else
                    c=sf::Color((std::uint8_t)(18+237*t),(std::uint8_t)(16+200*t),
                                (std::uint8_t)(30+90*d));
            } else if(view_mode==6){
                float o0=clamp01(odor[0][i]/ODOR_SCALE);
                float o1=clamp01(odor[1][i]/ODOR_SCALE);
                float o2=clamp01(odor[2][i]/ODOR_SCALE);
                float o3=clamp01(odor[3][i]/ODOR_SCALE);
                c=sf::Color((std::uint8_t)clampf(20+220*(o1+o3*0.7f),0.f,255.f),
                            (std::uint8_t)clampf(18+200*(o0+o2*0.8f),0.f,255.f),
                            (std::uint8_t)clampf(26+180*(o0*0.5f+o3),0.f,255.f));
            } else if(view_mode==5){
                float lp=clamp01(pload[i]/PATH_LOAD_MAX);
                float la=clamp01(aload[i]/PATH_LOAD_MAX);
                if(sea[i]) c=sf::Color(14,26,52);
                else c=sf::Color((std::uint8_t)(30+200*lp),(std::uint8_t)(24+60*la),
                                 (std::uint8_t)(40+160*la));
            } else {
                if(sea[i]&&ice[i]){
                    float sl=show_night?(0.35f+0.65f*daylit[i]):0.94f;
                    c=sf::Color((std::uint8_t)(214*sl),(std::uint8_t)(228*sl),
                                (std::uint8_t)(242*sl));
                } else if(sea[i]){
                    // 深さを非線形に強調する
                    float d=std::pow(clamp01((SEA_LEVEL-elev[i])/SEA_LEVEL),0.55f);
                    float sl=show_night?(0.32f+0.68f*daylit[i]):0.92f;
                    c=sf::Color((std::uint8_t)((8+10*(1-d))*sl),
                                (std::uint8_t)((25+35*(1-d))*sl),
                                (std::uint8_t)((60+55*(1-d))*sl));
                } else if(fire[i]>0.05f){
                    float f=std::min(1.f,fire[i]);
                    c=sf::Color((std::uint8_t)(200+55*f),(std::uint8_t)(90+90*f),
                                (std::uint8_t)(30*f));
                } else if(river[i]){
                    c=sf::Color(55,120,180);
                } else {
                    // ケッペンの気候区分ごとの植生の色
                    float r,g,b;
                    switch(climate[i]){
                        case 1: r=236;g=241;b=246; break;  // 氷雪
                        case 2: r=146;g=143;b=124; break;  // ツンドラ
                        case 3: r=204;g=178;b=122; break;  // 砂漠
                        case 4: r=172;g=160;b=98;  break;  // ステップ
                        case 5: r=58; g=96; b=78;  break;  // 亜寒帯
                        case 6: r=64; g=122;b=56;  break;  // 温帯
                        case 7: r=118;g=134;b=52;  break;  // サバナ
                        default:r=24; g=92; b=38;  break;  // 熱帯雨林
                    }
                    // 高い所ほど明るく塗り、起伏を見えるようにする
                    {
                        float m=0.25f*clamp01((elev[i]-SEA_LEVEL)/(1.f-SEA_LEVEL));
                        r=r*(1.f-m)+180*m; g=g*(1.f-m)+170*m; b=b*(1.f-m)+150*m;
                    }
                    float s=show_night?(0.30f+0.70f*daylit[i]):0.92f;
                    c=sf::Color((std::uint8_t)clampf(r*s,0.f,255.f),
                                (std::uint8_t)clampf(g*s,0.f,255.f),
                                (std::uint8_t)clampf(b*s,0.f,255.f));
                }
            }
            // 河川はどの地図でも見えるようにする
            if(!sea[i] && river[i] && view_mode!=9){
                float w=clamp01(flow_acc[i]/(RIVER_THRESHOLD*20.f));
                c=sf::Color((std::uint8_t)(c.r*0.25f+40*(1.f-w)),
                            (std::uint8_t)(c.g*0.25f+110+60*w),
                            (std::uint8_t)(c.b*0.25f+150+80*w));
            }
            image.setPixel({(unsigned)x,(unsigned)y},c);
        }
        texture.update(image);
        }
        prof_end(7);

        prof_begin();

        wrapW = GRID_W*zoom;
        nrep  = 0;                               // 区画は一つだけ描く(東西は開いた境界で、複製は無い)
        sprite.setScale({zoom,zoom});
        window.clear();
        {
            float bx0=w2sx(0.f), by=w2sy(0.f);
            for(int k=-nrep;k<=nrep;k++){
                float bx=bx0+k*wrapW;
                if(bx>SCREEN_W+2.f||bx+wrapW<-2.f) continue;
                sprite.setPosition({bx,by});
                window.draw(sprite);
            }
        }

        float ps=std::max(1.f,zoom*0.5f);
        static std::vector<sf::Vertex> vbuf;
        auto quad=[&](float sx,float sy,float r,sf::Color col){
            sf::Vertex v0,v1,v2,v3;
            v0.position={sx-r,sy-r}; v0.color=col;
            v1.position={sx+r,sy-r}; v1.color=col;
            v2.position={sx+r,sy+r}; v2.color=col;
            v3.position={sx-r,sy+r}; v3.color=col;
            vbuf.push_back(v0);vbuf.push_back(v1);vbuf.push_back(v2);
            vbuf.push_back(v0);vbuf.push_back(v2);vbuf.push_back(v3);
        };
        // 画面に入る複製の範囲を直接求める。空振りするループを回さない
        auto draw_wrapped=[&](float sx0,float sy,float r,sf::Color col){
            if(sy<-20.f||sy>SCREEN_H+20.f) return;
            if(sx0<-20.f||sx0>SCREEN_W+20.f) return;
            quad(sx0,sy,r,col);
        };
        // 引きの画面では点で描く。頂点数が6分の1になる
        auto point_wrapped=[&](float sx0,float sy,sf::Color col){
            if(sy<-20.f||sy>SCREEN_H+20.f) return;
            if(sx0<-20.f||sx0>SCREEN_W+20.f) return;
            sf::Vertex v; v.position={sx0,sy}; v.color=col;
            vbuf.push_back(v);
        };
        auto flush_points=[&](){
            if(!vbuf.empty())
                window.draw(vbuf.data(),vbuf.size(),sf::PrimitiveType::Points);
            vbuf.clear();
        };
        auto flush=[&](){
            if(!vbuf.empty())
                window.draw(vbuf.data(),vbuf.size(),sf::PrimitiveType::Triangles);
            vbuf.clear();
        };
        static std::vector<int> dorder;
        dorder.clear();
        dorder.reserve(plants.size());
        // 画面の外にある株は、はじめから並べない
        float vy0=cam_y-SCREEN_H/(2.f*zoom)-2.f;
        float vy1=cam_y+SCREEN_H/(2.f*zoom)+2.f;
        // 引いた画面では1マスが数ピクセルしかない。
        // 同じマスの株は重なって見えないので、代表だけを描く
        if(zoom<3.f){
            static std::vector<int> seen;
            seen.assign(GRID_W*GRID_H,-1);
            for(int i=0;i<(int)plants.size();i++){
                const Plant& p=plants[i];
                if(p.seed_wait>0) continue;
                if(p.y<vy0||p.y>vy1) continue;
                int gx=(int)p.x, gy=(int)p.y;
                if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
                int ci=idx(gx,gy);
                // 同じマスでは、最も背の高い株を代表にする
                if(seen[ci]<0||p.height>plants[seen[ci]].height) seen[ci]=i;
            }
            for(int ci=0;ci<GRID_W*GRID_H;ci++)
                if(seen[ci]>=0) dorder.push_back(seen[ci]);
        } else {
            for(int i=0;i<(int)plants.size();i++){
                if(plants[i].seed_wait>0) continue;
                if(plants[i].y<vy0||plants[i].y>vy1) continue;
                dorder.push_back(i);
            }
        }
        // 1セルが数ピクセルしかない引きの画面では、重なりの順序は見えない
        if(zoom>=4.f)
            std::sort(dorder.begin(),dorder.end(),[](int a,int b){
                return plants[a].height<plants[b].height; });

        bool hide_life = (view_mode==4);   // 気候を見るときは生物で覆わない
        vbuf.clear();
        if(!hide_life){
        vbuf.reserve(dorder.size()*6+64);
        for(int i:dorder){
            const Plant& p=plants[i];
            float sx=w2sx(p.x), sy=w2sy(p.y);
            if(sy<-20||sy>SCREEN_H+20) continue;
            float hn=std::min(1.f,p.height/3.f);
            sf::Color col((std::uint8_t)(30+80*hn),(std::uint8_t)(90+165*hn),
                          (std::uint8_t)(40+60*p.shade_tol));
            if(p.cold_tol>0.6f){ col.r=(std::uint8_t)std::min(255.f,(float)col.r+70.f);
                                 col.b=(std::uint8_t)std::min(255.f,(float)col.b+70.f); }
            if(p.aquatic>0.5f){ col.g=(std::uint8_t)(col.g*0.6f);
                                col.b=(std::uint8_t)std::min(255.f,(float)col.b+110.f); }
            if(color_species && p.lin>=0) col=lineage_color(p.lin);
            if(zoom<3.f) point_wrapped(sx,sy,col);
            else         draw_wrapped(sx,sy,ps*(0.8f+1.2f*hn),col);
        }
        if(zoom<3.f) flush_points(); else flush();
        }

        if(!hide_life){
        vbuf.reserve(animals.size()*6+64);
        for(const auto& a:animals){
            float sx=w2sx(a.x), sy=w2sy(a.y);
            if(sy<-20||sy>SCREEN_H+20) continue;
            sf::Color col((std::uint8_t)(120+135*a.diet),(std::uint8_t)(180-120*a.diet),
                          (std::uint8_t)(70+120*a.cold_tol));
            if(a.aquatic>0.5f){ col.r=(std::uint8_t)(col.r*0.5f);
                                col.b=(std::uint8_t)std::min(255.f,(float)col.b+130.f); }
            if(a.perch>0.3f){                   // 樹上の個体は明るく浮かせて描く
                col.r=(std::uint8_t)std::min(255.f,(float)col.r+55.f);
                col.g=(std::uint8_t)std::min(255.f,(float)col.g+55.f);
                sy-=std::min(8.f,a.perch*zoom*0.35f);
            }
            if(color_species && a.lin>=0) col=lineage_color(a.lin);
            if(zoom<3.f) point_wrapped(sx,sy,col);
            else         draw_wrapped(sx,sy,ps*(0.6f+0.8f*std::min(1.f,a.body/2.f)),col);
        }
        if(zoom<3.f) flush_points(); else flush();
        }

        // ===== 個体パネル =====
        if(font_ok && !show_tree && (sel_animal>=0||sel_plant>=0)){
            const Animal* sa=nullptr; const Plant* sp2=nullptr;
            if(sel_animal>=0) for(const auto& a:animals) if(a.id==sel_animal){ sa=&a; break; }
            if(sel_plant >=0) for(const auto& p:plants ) if(p.id==sel_plant ){ sp2=&p; break; }
            if(sa||sp2){
                float mx=sa?w2sx(sa->x):w2sx(sp2->x);
                float my=sa?w2sy(sa->y):w2sy(sp2->y);
                sf::CircleShape ring(15.f); ring.setOrigin({15.f,15.f});
                ring.setPosition({mx,my}); ring.setFillColor(sf::Color::Transparent);
                ring.setOutlineColor(sf::Color::Yellow); ring.setOutlineThickness(2.f);
                window.draw(ring);
                float PX=18.f,PY=18.f,PW=400.f,PH=sa?620.f:660.f;
                sf::RectangleShape bg({PW,PH}); bg.setPosition({PX,PY});
                bg.setFillColor(sf::Color(0,0,0,212));
                bg.setOutlineColor(sf::Color(150,150,150)); bg.setOutlineThickness(2.f);
                window.draw(bg);
                float bx=PX+80.f,bw=PW-104.f,by=PY+12.f,lh=20.f;
                auto lab=[&](float y,const std::string& s){
                    sf::Text t(font,jp(s.c_str()),12); t.setPosition({PX+10.f,y});
                    t.setFillColor(sf::Color::White); window.draw(t); };
                auto bar=[&](float y,float v,sf::Color col,const std::string& tx){
                    sf::RectangleShape b0({bw,13.f}); b0.setPosition({bx,y});
                    b0.setFillColor(sf::Color(50,50,50)); window.draw(b0);
                    sf::RectangleShape b1({bw*clamp01(v),13.f}); b1.setPosition({bx,y});
                    b1.setFillColor(col); window.draw(b1);
                    sf::Text t(font,tx,11); t.setPosition({bx+4.f,y});
                    t.setFillColor(sf::Color::White); window.draw(t); };
                char b[120];
                if(sa){
                    int ci=idx(std::max(0,std::min(GRID_W-1,(int)sa->x)),
                               std::max(0,std::min(GRID_H-1,(int)sa->y)));
                    int num=-1; for(const auto& L:lineages) if(L.id==sa->lin){ num=L.num; break; }
                    if(num>0) snprintf(b,120,"動物 #%d   (種 #A%d)",sa->id,num);
                    else      snprintf(b,120,"動物 #%d",sa->id);
                    sf::Text ttl(font,jp(b),14); ttl.setPosition({PX+10.f,by});
                    ttl.setFillColor(sf::Color(255,200,120)); window.draw(ttl); by+=lh+2;
                    lab(by,"体力"); snprintf(b,120,"%.2f / %.1f",sa->energy,A_DIVIDE_TH);
                    bar(by,sa->energy/A_DIVIDE_TH,sf::Color(0,200,80),b); by+=lh+6;
                    {
                        static const char* AN[13]={"肉食","体格","届く高さ","速度",
                                                   "顎","耐寒","水生","分解",
                                                   "感覚","樹上","呼吸","産子数","寿命"};
                        float sr=0.f;
                        for(int k=0;k<N_SENSOR;k++)
                            if(sa->sen[k].range>=SENSOR_MIN) sr+=sa->sen[k].range;
                        float av[13]={ sa->diet, sa->body/4.f, sa->reach/5.f, sa->speed,
                                       sa->jaw/2.f, sa->cold_tol, sa->aquatic,
                                       sa->detritus, sr/(SENSOR_MAX*N_SENSOR),
                                       sa->arboreal, sa->resp, sa->brood, sa->lifespan };
                        float ar[13]={ sa->diet, sa->body, sa->reach, sa->speed,
                                       sa->jaw, sa->cold_tol, sa->aquatic,
                                       sa->detritus, sr, sa->arboreal,
                                       sa->resp, sa->brood, sa->lifespan };
                        draw_radar(window,font,PX+PW*0.5f,by+104.f,70.f,
                                   AN,av,ar,13,sf::Color(255,190,110));
                        by+=248.f;   // 13軸はラベルが外へ張り出すので余白を広く取る
                        sf::Text nt(font,jp("神経網"),12);
                        nt.setFillColor(sf::Color(200,210,230));
                        nt.setPosition({PX+10.f,by}); window.draw(nt);
                        by+=16.f;
                        draw_nn(window,font,*sa,PX+6.f,by,PW-12.f,116.f);
                        by+=124.f;
                    }
                    { sf::Text t(font,jp(animal_summary(*sa).c_str()),13);
                      t.setFillColor(sf::Color(255,215,140));
                      t.setPosition({PX+10.f,by}); window.draw(t); by+=lh; }
                    snprintf(b,120,"齢 %d   捕食回数 %d",sa->age,sa->eaten);
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color::White); window.draw(t); by+=lh; }
                    snprintf(b,120,"気温 %.1f度  適温 %.1f  致死 %.1f",
                             temper[ci],opt_temp_a(*sa),kill_temp_a(*sa));
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(180,220,255)); window.draw(t); by+=lh; }
                    { int nv=0; float rs=0.f;
                      for(int k=0;k<N_SENSOR;k++)
                          if(sa->sen[k].range>=SENSOR_MIN){ nv++; rs+=sa->sen[k].range; }
                      snprintf(b,120,"感覚器 %d本  射程 %.1f  樹上 %.2f  酸素 %.0f%%",
                               nv,rs,sa->perch,sa->aerobic*100.f); }
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(200,220,255)); window.draw(t); by+=16.f; }
                    for(int k=0;k<N_SENSOR;k++){
                        const Animal::Sensor& S=sa->sen[k];
                        if(S.range<SENSOR_MIN) continue;
                        static const char* TL[N_TUNE]={
                            "植物","自分より小さい動物","自分より大きい動物","死骸",
                            "炎","林冠の高さ","病原体",
                            "匂い:生体","匂い:腐敗","匂い:食害","匂い:煙",
                            "同種","水際"};
                        int bd2=(int)(S.tune*N_TUNE); if(bd2>N_TUNE-1) bd2=N_TUNE-1;
                        snprintf(b,120,"  [%d] %+.0f度  射程%.1f  %s",
                                 k,S.angle*57.3f,S.range,TL[bd2]);
                        sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                        t.setFillColor(sf::Color(170,190,210)); window.draw(t); by+=16.f;
                    }
                    snprintf(b,120," ");
                    { sf::Text t(font,b,12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(170,170,170)); window.draw(t); }
                } else {
                    int ci=idx(std::max(0,std::min(GRID_W-1,(int)sp2->x)),
                               std::max(0,std::min(GRID_H-1,(int)sp2->y)));
                    int num=-1; for(const auto& L:lineages) if(L.id==sp2->lin){ num=L.num; break; }
                    if(num>0) snprintf(b,120,"植物 #%d   (種 #P%d)",sp2->id,num);
                    else      snprintf(b,120,"植物 #%d",sp2->id);
                    sf::Text ttl(font,jp(b),14); ttl.setPosition({PX+10.f,by});
                    ttl.setFillColor(sf::Color(140,235,140)); window.draw(ttl); by+=lh+2;
                    lab(by,"体力");  snprintf(b,120,"%.2f / %.1f",sp2->energy,PLANT_DIVIDE_TH);
                    bar(by,sp2->energy/PLANT_DIVIDE_TH,sf::Color(0,200,80),b); by+=lh;
                    lab(by,"樹高");  snprintf(b,120,"%.2f / %.2f",sp2->height,sp2->max_height);
                    bar(by,sp2->height/std::max(0.01f,sp2->max_height),
                        sf::Color(90,230,110),b); by+=lh+6;
                    {
                        static const char* PN[13]={
                            "最大樹高","吸収","耐陰","防御","地下茎","散布","耐寒",
                            "水生","耐乾","つる","耐火","貯水","寿命"};
                        float pv[13]={ sp2->max_height/5.f, sp2->absorb, sp2->shade_tol,
                                       sp2->tough, sp2->clonal, sp2->disperse,
                                       sp2->cold_tol, sp2->aquatic, sp2->drought_tol,
                                       sp2->climb, sp2->fire_tol,
                                       sp2->store_cap/STORE_MAX, sp2->lifespan };
                        float pr[13]={ sp2->max_height, sp2->absorb, sp2->shade_tol,
                                       sp2->tough, sp2->clonal, sp2->disperse,
                                       sp2->cold_tol, sp2->aquatic, sp2->drought_tol,
                                       sp2->climb, sp2->fire_tol,
                                       sp2->store_cap, sp2->lifespan };
                        draw_radar(window,font,PX+PW*0.5f,by+104.f,70.f,
                                   PN,pv,pr,13,sf::Color(140,235,140));
                        by+=248.f;   // 13軸はラベルが外へ張り出すので余白を広く取る
                        // 同じ型の、世界で最小・最大の個体と並べて姿を見せる
                        int fm=plant_form(*sp2);
                        float lo=form_min[fm], hi=form_max[fm];
                        char fb[80];
                        snprintf(fb,80,"%s   %d individuals   brood %.2f  nFix %.2f",
                                 FORM_NAME[fm],form_cnt[fm],sp2->brood,sp2->n_fix);
                        sf::Text ft(font,jp(fb),12);
                        ft.setFillColor(sf::Color(180,235,180));
                        ft.setPosition({PX+10.f,by}); window.draw(ft);
                        by+=18.f;
                        float baseY=by+96.f, HMAX=88.f;
                        sf::RectangleShape gl({PW-24.f,1.f});
                        gl.setPosition({PX+12.f,baseY});
                        gl.setFillColor(sf::Color(90,90,96)); window.draw(gl);
                        float span=std::max(0.01f,hi-lo);
                        draw_plant_form(window,fm,PX+PW*0.22f,baseY,
                                        HMAX*std::max(0.12f,lo/std::max(0.01f,hi)),
                                        sf::Color(90,130,90));
                        draw_plant_form(window,fm,PX+PW*0.50f,baseY,
                                        HMAX*std::max(0.12f,sp2->max_height/std::max(0.01f,hi)),
                                        sf::Color(150,245,150));
                        draw_plant_form(window,fm,PX+PW*0.78f,baseY,HMAX,
                                        sf::Color(90,130,90));
                        by=baseY+8.f;
                        snprintf(fb,80,"min %.2f",lo);
                        { sf::Text t(font,fb,10); t.setFillColor(sf::Color(140,160,140));
                          t.setPosition({PX+PW*0.22f-22.f,by}); window.draw(t); }
                        snprintf(fb,80,"%.2f",sp2->max_height);
                        { sf::Text t(font,fb,11); t.setFillColor(sf::Color(180,250,180));
                          t.setPosition({PX+PW*0.50f-14.f,by}); window.draw(t); }
                        snprintf(fb,80,"max %.2f",hi);
                        { sf::Text t(font,fb,10); t.setFillColor(sf::Color(140,160,140));
                          t.setPosition({PX+PW*0.78f-22.f,by}); window.draw(t); }
                        by+=18.f;
                        float rel=(sp2->max_height-lo)/span;
                        sf::RectangleShape rb0({PW-40.f,8.f});
                        rb0.setPosition({PX+20.f,by}); rb0.setFillColor(sf::Color(52,52,56));
                        window.draw(rb0);
                        sf::RectangleShape rb1({(PW-40.f)*clamp01(rel),8.f});
                        rb1.setPosition({PX+20.f,by}); rb1.setFillColor(sf::Color(150,245,150));
                        window.draw(rb1);
                        by+=18.f;
                    }
                    snprintf(b,120,"齢 %d   光の飽和点 %.2f",sp2->age,saturation(*sp2));
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color::White); window.draw(t); by+=lh; }
                    snprintf(b,120,"気温 %.1f度  適温 %.1f  致死 %.1f",
                             temper[ci],opt_temp_p(*sp2),kill_temp_p(*sp2));
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(180,220,255)); window.draw(t); by+=lh; }
                    snprintf(b,120,"標高 %.2f  CO2 %.2f",elev[ci],co2[ci]);
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(170,170,170)); window.draw(t); }
                }
            } else { sel_animal=-1; sel_plant=-1; }
        }

        // ===== グラフ =====
        if(show_graph && !show_tree && h_plant.size()>1){
            int N=(int)h_plant.size();
            float GW=SCREEN_W*0.26f, GH=SCREEN_H*0.17f;
            float GX=SCREEN_W-GW-16.f, gy=16.f;
            const float GAP=44.f;
            float ph=GH-34;
            char b[240];
            auto panel=[&](const char* title,int gid)->float{
                sf::RectangleShape bg({GW,GH}); bg.setPosition({GX,gy});
                bg.setFillColor(sf::Color(0,0,0,195));
                bg.setOutlineColor(sf::Color(120,120,120)); bg.setOutlineThickness(1.f);
                window.draw(bg);
                if(font_ok){ sf::Text t(font,jp(title),12); t.setPosition({GX+6.f,gy+4.f});
                             t.setFillColor(sf::Color(200,200,200)); window.draw(t); }
                bool hov=(mouse_x>=GX&&mouse_x<=GX+GW&&mouse_y>=gy&&mouse_y<=gy+GH);
                if(hov){
                    sf::RectangleShape hl({GW,GH}); hl.setPosition({GX,gy});
                    hl.setFillColor(sf::Color(0,0,0,0));
                    hl.setOutlineColor(sf::Color(180,220,240)); hl.setOutlineThickness(2.f);
                    window.draw(hl);
                }
                ui_hits.push_back({GX,gy,GW,GH,400+gid});
                return gy+22.f; };
            auto lineF=[&](float py0,const std::vector<float>& v,float lo,float hi,sf::Color col){
                sf::VertexArray ln(sf::PrimitiveType::LineStrip);
                float px0=GX+6, pw=GW-12;
                // 指定の上限を超えたら、実際の最大に合わせて伸ばす。振り切れない
                for(int i=0;i<N&&i<(int)v.size();i++) if(v[i]>hi) hi=v[i]*1.08f;
                for(int i=0;i<N&&i<(int)v.size();i++){
                    float t=(v[i]-lo)/std::max(1e-6f,hi-lo);
                    sf::Vertex vt;
                    vt.position={px0+pw*i/(float)(N-1),py0+ph*(1.f-clamp01(t))};
                    vt.color=col; ln.append(vt);
                }
                window.draw(ln); };
            auto lineI=[&](float py0,const std::vector<int>& v,int pk,sf::Color col){
                sf::VertexArray ln(sf::PrimitiveType::LineStrip);
                float px0=GX+6, pw=GW-12;
                for(int i=0;i<N&&i<(int)v.size();i++) if(v[i]>pk) pk=v[i];
                for(int i=0;i<N&&i<(int)v.size();i++){
                    sf::Vertex vt;
                    vt.position={px0+pw*i/(float)(N-1),
                                 py0+ph*(1.f-clamp01((float)v[i]/std::max(1,pk)))};
                    vt.color=col; ln.append(vt);
                }
                window.draw(ln); };
            auto legend=[&](const char* s,sf::Color col){
                if(!font_ok) return;
                sf::Text t(font,jp(s),11); t.setFillColor(col);
                sf::FloatRect b2=t.getLocalBounds();
                sf::RectangleShape bp({b2.size.x+10.f,16.f});
                bp.setPosition({GX-3.f,gy+GH+1.f});
                bp.setFillColor(sf::Color(0,0,0,170));
                window.draw(bp);
                t.setPosition({GX,gy+GH+2.f}); window.draw(t);
            };

            if(view_mode==0){
                float py=panel("個体数",0);
                lineI(py,h_plant ,peak_plant ,sf::Color(90,230,110));
                lineI(py,h_animal,peak_animal,sf::Color(255,150,70));
                lineI(py,h_corpse,peak_corpse,sf::Color(180,90,80));
                snprintf(b,240,"植物 %d   動物 %d   死骸 %d",
                         h_plant.back(),h_animal.back(),h_corpse.back());
                legend(b,sf::Color(200,220,200)); gy+=GH+GAP;
                py=panel("植物の形質",1);
                lineF(py,h_maxH ,0.f,3.f,sf::Color(90,230,110));
                lineF(py,h_shade,0.f,1.f,sf::Color(70,200,230));
                lineF(py,h_tough,0.f,1.f,sf::Color(200,200,90));
                snprintf(b,240,"最大樹高 %.2f   耐陰性 %.2f   防御 %.2f",
                         h_maxH.back(),h_shade.back(),h_tough.back());
                legend(b,sf::Color(200,220,200)); gy+=GH+GAP;
                py=panel("動物の形質",2);
                lineF(py,h_reach,0.f,3.f,sf::Color(255,150,70));
                lineF(py,h_diet ,0.f,1.f,sf::Color(255,80,80));
                lineF(py,h_jaw  ,0.f,2.f,sf::Color(200,120,255));
                snprintf(b,240,"届く高さ %.2f   肉食性 %.2f   顎 %.2f",
                         h_reach.back(),h_diet.back(),h_jaw.back());
                legend(b,sf::Color(200,220,200));
            } else if(view_mode==1){
                float py=panel("陸の平均気温",3);
                lineF(py,h_temp,-30.f,40.f,sf::Color(255,140,90));
                snprintf(b,240,"気温 %.1f度   (目盛 -30〜40)",h_temp.back());
                legend(b,sf::Color(255,180,140)); gy+=GH+GAP;
                py=panel("耐寒性",4);
                lineF(py,h_coldP,0.f,1.f,sf::Color(150,220,255));
                lineF(py,h_coldA,0.f,1.f,sf::Color(255,255,255));
                snprintf(b,240,"植物 %.2f   動物 %.2f",h_coldP.back(),h_coldA.back());
                legend(b,sf::Color(200,220,255)); gy+=GH+GAP;
                py=panel("個体数",0);
                lineI(py,h_plant ,peak_plant ,sf::Color(90,230,110));
                lineI(py,h_animal,peak_animal,sf::Color(255,150,70));
                snprintf(b,240,"植物 %d   動物 %d",h_plant.back(),h_animal.back());
                legend(b,sf::Color(200,220,200));
            } else {
                float py=panel("大気",5);
                lineF(py,h_co2,-1.f,1.f,sf::Color(230,90,200));   // 柱の量からのずれのマスの平均(気体の単位)
                lineF(py,h_o2 ,-1.f,1.f,sf::Color(70,180,240));
                // 今見ている地図に対応するグラフを重ねる
                switch(view_mode){
                case 4: lineF(py,h_soil,0.f,1.f,sf::Color(90,190,255));  break;
                case 5: lineF(py,h_path,0.f,2.f,sf::Color(230,90,120));  break;
                case 6: lineF(py,h_odor,0.f,4.f,sf::Color(200,150,255)); break;
                case 7: lineF(py,h_ice ,0.f,0.5f,sf::Color(220,235,250));break;
                case 8: lineF(py,h_nut ,0.f,1.0f,sf::Color(150,230,120));break;
                case 10:lineF(py,h_sst ,-9.f,9.f,sf::Color(255,140,90)); break;
                }
                snprintf(b,240,"CO2 %.3f   O2 %.3f",h_co2.back(),h_o2.back());
                legend(b,sf::Color(220,180,220)); gy+=GH+GAP;
                py=panel("生産者と消費者",6);
                lineI(py,h_plant ,peak_plant ,sf::Color(90,230,110));
                lineI(py,h_animal,peak_animal,sf::Color(255,150,70));
                snprintf(b,240,"植物 %d   動物 %d   死骸 %d",
                         h_plant.back(),h_animal.back(),h_corpse.back());
                legend(b,sf::Color(200,220,200));
            }
            if(font_ok){
                const char* vm[11]={"地形","気温","大気","腐植",
                                    "気候区分","病原体","匂い","昼夜","栄養",
                                    "傾斜","海流"};
                int span=N*SAMPLE_INTERVAL;
                snprintf(b,240,"span -%.1fk .. now   t %d   zoom %.1fx   [V] %s%s",
                         span/1000.f,frame,zoom,vm[view_mode],paused?"   [PAUSED]":"");
                sf::Text t(font,b,11); t.setFillColor(sf::Color(190,190,190));
                sf::FloatRect b3=t.getLocalBounds();
                sf::RectangleShape bp2({b3.size.x+10.f,16.f});
                bp2.setPosition({GX-3.f,gy+GH+15.f});
                bp2.setFillColor(sf::Color(0,0,0,170));
                window.draw(bp2);
                t.setPosition({GX,gy+GH+16.f}); window.draw(t);
            }

            // プロファイラ: 所要時間の長い順に並べる
            if(font_ok){
                int ord[PROF_N]; for(int k=0;k<PROF_N;k++) ord[k]=k;
                std::sort(ord,ord+PROF_N,[](int a,int b){ return prof[a]>prof[b]; });
                double sum=0; for(int k=0;k<PROF_N;k++) sum+=prof[k];
                float py=gy+GH+36.f;
                sf::RectangleShape bp({GW,18.f*(PROF_N+1)+8.f});
                bp.setPosition({GX-3.f,py-3.f});
                bp.setFillColor(sf::Color(0,0,0,185));
                bp.setOutlineColor(sf::Color(110,110,110)); bp.setOutlineThickness(1.f);
                window.draw(bp);
                char pb[96];
                snprintf(pb,96,"profile   total %.1f ms  (%.0f fps)",sum,sum>0?1000.0/sum:0.0);
                sf::Text ph(font,pb,12); ph.setFillColor(sf::Color(220,220,140));
                ph.setPosition({GX+2.f,py}); window.draw(ph);
                py+=18.f;
                for(int k=0;k<PROF_N;k++){
                    int j=ord[k];
                    float frac=(sum>0)?(float)(prof[j]/sum):0.f;
                    sf::RectangleShape bar({(GW-110.f)*frac,12.f});
                    bar.setPosition({GX+96.f,py+1.f});
                    bar.setFillColor(frac>0.35f?sf::Color(230,90,70)
                                    :frac>0.15f?sf::Color(230,180,70)
                                               :sf::Color(90,180,120));
                    window.draw(bar);
                    snprintf(pb,96,"%-8s %5.2f",PROF_NAME[j],prof[j]);
                    sf::Text pt(font,pb,12);
                    pt.setFillColor(sf::Color(200,200,200));
                    pt.setPosition({GX+2.f,py}); window.draw(pt);
                    py+=18.f;
                }
            }
        }

        // ===== 系統樹 =====
        if(show_tree && font_ok){
            sf::RectangleShape bg({SCREEN_W,SCREEN_H});
            bg.setFillColor(sf::Color(0,0,0,240)); window.draw(bg);
            tree_hits.clear();

            int minb=frame;
            for(const auto& L:lineages) if(L.birth_tick<minb) minb=L.birth_tick;
            float span=(float)std::max(1,frame-minb);
            float topY=84.f, botY=SCREEN_H-80.f;
            auto ty=[&](int t){
                float base=topY+(botY-topY)*((float)(t-minb)/span);
                return (base-SCREEN_H*0.5f)*tree_zoom+SCREEN_H*0.5f+tree_oy; };

            int snum=-1;
            if(!tree_search.empty()){ int s=0; for(char ch:tree_search) s=s*10+(ch-'0'); snum=s; }
            sf::Vector2i mp=sf::Mouse::getPosition(window);
            int hover=-1; float hbest=12.f;

            static std::unordered_map<int,float> cache_lx[2];
            static std::unordered_map<int,int>   cache_id2i[2];
            static std::vector<int>              cache_ix[2];
            static int cache_leaf[2]={1,1};
            static int cache_n[2]={-1,-1};

            for(int kind=0;kind<2;kind++){
                float L0=(kind==0)?SCREEN_W*0.04f:SCREEN_W*0.53f;
                float LW=SCREEN_W*0.43f;

                int total_k=0;
                for(const auto& L:lineages) if(L.kind==kind) total_k++;
                if(total_k==0) continue;

                // 系統が増えた時だけ、配置を組み直す
                if(cache_n[kind]!=total_k){
                    cache_n[kind]=total_k;
                    auto& cix=cache_ix[kind]; cix.clear();
                    for(int i=0;i<(int)lineages.size();i++)
                        if(lineages[i].kind==kind) cix.push_back(i);
                    auto& cid=cache_id2i[kind]; cid.clear();
                    for(int i:cix) cid[lineages[i].id]=i;
                    std::unordered_map<int,std::vector<int>> ch;   // 子リストを一度で作る
                    for(int j:cix) ch[lineages[j].parent].push_back(j);
                    auto& clx=cache_lx[kind]; clx.clear();
                    int leaf=0;
                    std::function<float(int,int)> assign=[&](int i,int depth)->float{
                        if(depth>1500){ clx[i]=(float)(leaf++); return clx[i]; }
                        auto it=ch.find(lineages[i].id);
                        if(it==ch.end()||it->second.empty()){ clx[i]=(float)(leaf++); return clx[i]; }
                        float s=0; for(int c:it->second) s+=assign(c,depth+1);
                        clx[i]=s/it->second.size(); return clx[i]; };
                    for(int i:cix)
                        if(lineages[i].parent<0||cid.find(lineages[i].parent)==cid.end())
                            assign(i,0);
                    cache_leaf[kind]=leaf;
                }

                const auto& ix=cache_ix[kind];
                const auto& id2i=cache_id2i[kind];
                const auto& lx=cache_lx[kind];
                int leaf=cache_leaf[kind];
                float spanx=(leaf>1)?(float)(leaf-1):1.f;
                auto tx=[&](float slot){
                    float base=L0+(leaf>1?slot/spanx:0.5f)*LW;
                    return (base-SCREEN_W*0.5f)*tree_zoom+SCREEN_W*0.5f+tree_ox; };

                sf::VertexArray va(sf::PrimitiveType::Lines);
                // 画面に入る縦の範囲を先に求め、外れる系統は map を引く前に捨てる
                float vt0=-40.f, vt1=SCREEN_H+40.f;
                for(int i:ix){
                    const Lineage& L=lineages[i];
                    int bt=L.birth_tick, dt=L.alive?frame:L.death_tick;
                    float yb=ty(bt), yd=ty(dt);
                    if(yd<vt0||yb>vt1) continue;
                    auto lit=lx.find(i);
                    if(lit==lx.end()) continue;
                    float x=tx(lit->second);
                    if(x<-40.f||x>SCREEN_W+40.f) continue;
                    bool match=lineage_match_multi(L);
                    sf::Color col=lineage_color(L.id);
                    if(!match) col=sf::Color(70,70,70,120);
                    sf::Vertex a,b2;
                    a.position={x,ty(bt)}; a.color=col;
                    b2.position={x,ty(dt)}; b2.color=col;
                    va.append(a); va.append(b2);
                    auto pit=id2i.find(L.parent);
                    auto plit=(pit!=id2i.end())?lx.find(pit->second):lx.end();
                    if(plit!=lx.end()){
                        float px=tx(plit->second);
                        sf::Color pc=col; pc.a=match?150:70;
                        sf::Vertex c,d;
                        c.position={px,ty(bt)}; c.color=pc;
                        d.position={x, ty(bt)}; d.color=pc;
                        va.append(c); va.append(d);
                    }
                    tree_hits.push_back({i,x,ty(bt),ty(dt)});
                    float mdx=std::fabs((float)mp.x-x);
                    if(mp.y>=ty(bt)-5&&mp.y<=ty(dt)+5&&mdx<hbest){ hbest=mdx; hover=i; }
                    if(snum>=0&&L.num==snum){
                        static int anim=0; anim++;
                        float pulse=0.55f+0.45f*std::sin(anim*0.25f);
                        float rr=7.f+5.f*pulse;
                        for(float yy:{ty(bt),ty(dt)}){
                            sf::CircleShape ring(rr); ring.setOrigin({rr,rr});
                            ring.setPosition({x,yy});
                            ring.setFillColor(sf::Color::Transparent);
                            ring.setOutlineColor(sf::Color(255,255,90,(int)(pulse*255)));
                            ring.setOutlineThickness(2.5f); window.draw(ring);
                        }
                    }
                    if(sel_lineage==i){
                        sf::RectangleShape hl({7.f,std::max(3.f,ty(dt)-ty(bt))});
                        hl.setPosition({x-3.5f,ty(bt)});
                        hl.setFillColor(sf::Color(255,255,255,70)); window.draw(hl);
                    }
                }
                window.draw(va);

                sf::VertexArray xm(sf::PrimitiveType::Lines);
                for(int i:ix){
                    auto lit2=lx.find(i);
                    if(lineages[i].alive||lit2==lx.end()) continue;
                    if(!lineage_match_multi(lineages[i])) continue;
                    float x=tx(lit2->second), y=ty(lineages[i].death_tick), r=5.f;
                    if(x<-40.f||x>SCREEN_W+40.f||y<-40.f||y>SCREEN_H+40.f) continue;
                    sf::Color col=lineage_color(lineages[i].id);
                    sf::Vertex a,b2,c,d;
                    a.position={x-r,y-r}; a.color=col; b2.position={x+r,y+r}; b2.color=col;
                    c.position={x-r,y+r}; c.color=col; d.position={x+r,y-r}; d.color=col;
                    xm.append(a);xm.append(b2);xm.append(c);xm.append(d);
                }
                window.draw(xm);

                int al=0,ev=0,mt=0;
                for(int i:ix){ ev++; if(lineages[i].alive) al++;
                               if(lineage_match_multi(lineages[i])) mt++; }
                char hb[160];
                snprintf(hb,160,"%s   現存 %d / 累計 %d   該当 %d",
                         kind==0?"植物  (#P)":"動物  (#A)",al,ev,mt);
                sf::Text ht(font,jp(hb),16);
                ht.setFillColor(kind==0?sf::Color(140,235,140):sf::Color(255,190,120));
                ht.setPosition({L0,50.f}); window.draw(ht);
            }

            char tb[260];
            // 閉じるボタン。右上の角
            {
                float cw=34.f, cxp=SCREEN_W-20.f-cw, cyp=16.f;
                bool hov=(mouse_x>=cxp&&mouse_x<=cxp+cw&&mouse_y>=cyp&&mouse_y<=cyp+cw);
                sf::RectangleShape cb3({cw,cw});
                cb3.setPosition({cxp,cyp});
                cb3.setFillColor(hov?sf::Color(110,40,40,235):sf::Color(40,22,22,220));
                cb3.setOutlineColor(hov?sf::Color(255,160,160):sf::Color(170,110,110));
                cb3.setOutlineThickness(1.f);
                window.draw(cb3);
                sf::Text ct3(font,jp("×"),20);
                ct3.setFillColor(hov?sf::Color::White:sf::Color(230,180,180));
                sf::FloatRect cr3=ct3.getLocalBounds();
                ct3.setPosition({cxp+cw*0.5f-cr3.size.x*0.5f,cyp+3.f});
                window.draw(ct3);
                ui_hits.push_back({cxp,cyp,cw,cw,103});
            }
            snprintf(tb,260,"系統樹   t %d   倍率 %.1f   [T 閉じる] [WASD 移動] [Z/X 拡大縮小] [F フィルタ解除] [数字 検索] [クリック 種の詳細]",
                     frame,tree_zoom);
            sf::Text tt(font,jp(tb),15); tt.setFillColor(sf::Color::White);
            tt.setPosition({SCREEN_W*0.04f,18.f}); window.draw(tt);

            {
                bool found=false;
                if(snum>=0) for(const auto& L:lineages) if(L.num==snum){ found=true; break; }
                char sb2[140];
                snprintf(sb2,140,"検索 #: %s_%s",
                         tree_search.c_str(),
                         tree_search.empty()?"  (数字を入力)":(found?"  (見つかった)":"  (該当なし)"));

                // フィルタ: クリックで切り替えるチェックボックス
                filter_hits.clear();
                {
                    float fx=SCREEN_W*0.04f, fy=SCREEN_H-28.f;
                    for(int f=0;f<FILTER_N;f++){
                        sf::Text lt(font,jp(FILTER_NAME[f]),12);
                        float tw=lt.getLocalBounds().size.x;
                        float bw=tw+26.f;
                        if(fx+bw>SCREEN_W-20.f){ fx=SCREEN_W*0.04f; fy-=20.f; }
                        sf::RectangleShape bg2({bw,17.f});
                        bg2.setPosition({fx,fy});
                        bg2.setFillColor(filter_on[f]?sf::Color(40,70,50,230)
                                                     :sf::Color(20,20,20,200));
                        bg2.setOutlineColor(filter_on[f]?sf::Color(120,230,140)
                                                        :sf::Color(90,90,90));
                        bg2.setOutlineThickness(1.f);
                        window.draw(bg2);
                        sf::Text ck(font,filter_on[f]?"x":" ",12);
                        ck.setFillColor(sf::Color(150,255,170));
                        ck.setPosition({fx+5.f,fy+1.f});
                        window.draw(ck);
                        lt.setFillColor(filter_on[f]?sf::Color(230,255,235)
                                                    :sf::Color(150,150,150));
                        lt.setPosition({fx+18.f,fy+1.f});
                        window.draw(lt);
                        filter_hits.push_back({f,fx,fy,bw,17.f});
                        fx += bw+5.f;
                    }
                }
                sf::Text st(font,jp(sb2),14);
                st.setFillColor(found?sf::Color(255,255,90):sf::Color(180,180,180));
                st.setPosition({SCREEN_W*0.04f,SCREEN_H-70.f}); window.draw(st);
            }

            if(hover>=0){
                const Lineage& L=lineages[hover];
                char hb[200];
                snprintf(hb,200,"#%s%d   現在 %d   最大 %d   誕生 t%d%s",
                         L.kind==0?"P":"A",L.num,L.pop,L.max_pop,L.birth_tick,
                         L.alive?"":"   (絶滅)");
                sf::Text t(font,jp(hb),14); t.setFillColor(sf::Color::White);
                sf::FloatRect bn=t.getLocalBounds();
                float bx=(float)mp.x+14.f, by=(float)mp.y+14.f;
                if(bx+bn.size.x+16.f>SCREEN_W) bx=(float)mp.x-14.f-bn.size.x-16.f;
                if(by+26.f>SCREEN_H) by=(float)mp.y-30.f;
                sf::RectangleShape bgb({bn.size.x+16.f,25.f});
                bgb.setPosition({bx-4.f,by-3.f});
                bgb.setFillColor(sf::Color(20,20,20,240));
                bgb.setOutlineColor(lineage_color(L.id)); bgb.setOutlineThickness(1.5f);
                window.draw(bgb);
                t.setPosition({bx+2.f,by}); window.draw(t);
            }

            // 選択した種の平均ステータス
            if(sel_lineage>=0 && sel_lineage<(int)lineages.size()){
                const Lineage& L=lineages[sel_lineage];
                float PX=18.f,PY=110.f,PW=400.f,PH=430.f;
                sf::RectangleShape bg2({PW,PH}); bg2.setPosition({PX,PY});
                bg2.setFillColor(sf::Color(0,0,0,225));
                bg2.setOutlineColor(lineage_color(L.id)); bg2.setOutlineThickness(2.f);
                window.draw(bg2);
                float bx=PX+84.f,bw=PW-108.f,by=PY+12.f,lh=20.f;
                auto lab=[&](float y,const std::string& s){
                    sf::Text t(font,jp(s.c_str()),12); t.setPosition({PX+10.f,y});
                    t.setFillColor(sf::Color::White); window.draw(t); };
                auto bar=[&](float y,float v,sf::Color col,const std::string& tx){
                    sf::RectangleShape b0({bw,13.f}); b0.setPosition({bx,y});
                    b0.setFillColor(sf::Color(50,50,50)); window.draw(b0);
                    sf::RectangleShape b1({bw*clamp01(v),13.f}); b1.setPosition({bx,y});
                    b1.setFillColor(col); window.draw(b1);
                    sf::Text t(font,tx,11); t.setPosition({bx+4.f,y});
                    t.setFillColor(sf::Color::White); window.draw(t); };
                char b[130];
                // 生きていれば現在の平均、絶滅していれば最後に見えていた姿
                const float* vm_ = view_mean(L);
                snprintf(b,130,"種 #%s%d  (平均)",L.kind==0?"P":"A",L.num);
                sf::Text ttl(font,jp(b),15); ttl.setPosition({PX+10.f,by});
                ttl.setFillColor(lineage_color(L.id)); window.draw(ttl); by+=lh+4;
                snprintf(b,130,"個体数 %d   最大 %d   誕生 t%d%s",
                         L.alive?L.pop:L.final_pop,L.max_pop,L.birth_tick,
                         L.alive?"":"  絶滅");
                { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                  t.setFillColor(sf::Color(200,200,200)); window.draw(t); by+=lh; }
                if(L.kind==0){
                    static const char* PN[9]={"最大樹高","吸収","耐陰","防御",
                                              "地下茎","散布","耐寒","水生","耐乾"};
                    float pv[9]={ vm_[1]/5.f, vm_[0], vm_[2], vm_[3],
                                  vm_[4], vm_[5], vm_[6], vm_[8], vm_[9] };
                    float pr[9]={ vm_[1], vm_[0], vm_[2], vm_[3],
                                  vm_[4], vm_[5], vm_[6], vm_[8], vm_[9] };
                    draw_radar(window,font,PX+PW*0.5f,by+98.f,68.f,
                               PN,pv,pr,9,lineage_color(L.id));
                    by+=228.f;
                    snprintf(b,130,"樹高 %.2f   水生 %.2f",vm_[7],vm_[8]);
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(190,200,190)); window.draw(t); }
                } else {
                    static const char* AN[9]={"肉食","体格","届く高さ","速度",
                                              "顎","耐寒","水生","感覚","樹上"};
                    float av[9]={ vm_[4], vm_[3]/4.f, vm_[2]/5.f, vm_[0],
                                  vm_[5]/2.f, vm_[6], vm_[10],
                                  vm_[11]/(SENSOR_MAX*N_SENSOR), vm_[12] };
                    float ar[9]={ vm_[4], vm_[3], vm_[2], vm_[0],
                                  vm_[5], vm_[6], vm_[10], vm_[11], vm_[12] };
                    draw_radar(window,font,PX+PW*0.5f,by+98.f,68.f,
                               AN,av,ar,9,lineage_color(L.id));
                    by+=228.f;
                    snprintf(b,130,"体格 %.2f   感覚器の射程 %.1f",vm_[7],vm_[11]);
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(200,200,190)); window.draw(t); by+=18.f; }
                    snprintf(b,130,"神経網の規模 %+.2f  層の数 %.1f",vm_[8],vm_[9]);
                    { sf::Text t(font,jp(b),12); t.setPosition({PX+10.f,by});
                      t.setFillColor(sf::Color(170,170,170)); window.draw(t); }
                }
            }
        }

        // ===== セーブ/ロードUI =====
        if(font_ok && ui_mode!=0){
            int nf=(int)save_files.size();
            float bw=700.f;
            float bh=(ui_mode==1)?76.f:std::min(SCREEN_H*0.85f,60.f+std::max(1,nf)*18.f);
            float bx=SCREEN_W/2.f-bw/2.f, by=SCREEN_H/2.f-bh/2.f;
            sf::RectangleShape box({bw,bh}); box.setPosition({bx,by});
            box.setFillColor(sf::Color(0,0,0,228));
            box.setOutlineColor(sf::Color(200,200,200)); box.setOutlineThickness(2.f);
            window.draw(box);
            char pr[140];
            if(ui_mode==1)
                snprintf(pr,140,"保存名:  %s_",input_text.c_str());
            else
                snprintf(pr,140,"読み込み (名前を入力してEnter)   [Tab] 並び: %s   %s_",
                         sort_by_time?"新しい順":"名前順", input_text.c_str());
            sf::Text t(font,jp(pr),18); t.setFillColor(sf::Color::White);
            t.setPosition({bx+16.f,by+14.f}); window.draw(t);
            if(ui_mode==2){
                float ly=by+46.f; int maxshow=(int)((bh-52.f)/18.f);
                hover_save=-1;
                for(int i=0;i<nf&&i<maxshow;i++){
                    float rowy=ly-2.f, rowh=18.f;
                    bool hov=(mouse_x>=bx+12.f&&mouse_x<=bx+bw-12.f&&
                              mouse_y>=rowy&&mouse_y<=rowy+rowh);
                    if(hov){
                        hover_save=i;
                        sf::RectangleShape hb({bw-24.f,rowh});
                        hb.setPosition({bx+12.f,rowy});
                        hb.setFillColor(sf::Color(40,70,100,150));
                        window.draw(hb);
                    }
                    // 削除ボタン
                    float dxp=bx+bw-46.f, dw=20.f;
                    bool dhov=(mouse_x>=dxp&&mouse_x<=dxp+dw&&
                               mouse_y>=rowy&&mouse_y<=rowy+rowh);
                    bool pend=(del_confirm==i);
                    sf::RectangleShape db({dw,rowh});
                    db.setPosition({dxp,rowy});
                    db.setFillColor(pend?sf::Color(180,40,40,230)
                                   :(dhov?sf::Color(120,50,50,200):sf::Color(0,0,0,0)));
                    window.draw(db);
                    sf::Text dt(font,jp("×"),13);
                    dt.setFillColor(pend?sf::Color::White
                                   :(dhov?sf::Color(255,180,180):sf::Color(140,110,110)));
                    sf::FloatRect dr=dt.getLocalBounds();
                    dt.setPosition({dxp+dw*0.5f-dr.size.x*0.5f,ly});
                    window.draw(dt);
                    ui_hits.push_back({dxp,rowy,dw,rowh,300+i});
                    // 行をクリックしたら、その名前を入力欄に書き込む
                    ui_hits.push_back({bx+12.f,rowy,bw-70.f,rowh,500+i});
                    sf::Text ft(font,jp(save_files[i].c_str()),13);
                    ft.setFillColor(pend?sf::Color(255,190,190):sf::Color(180,220,255));
                    ft.setPosition({bx+20.f,ly}); window.draw(ft);
                    if(i<(int)save_times.size()){
                        std::time_t tt=save_times[i];
                        char ts[32]; std::tm lt{};
                        #ifdef _WIN32
                          localtime_s(&lt,&tt);
                        #else
                          localtime_r(&tt,&lt);
                        #endif
                        std::strftime(ts,32,"%m/%d %H:%M",&lt);
                        sf::Text tt2(font,ts,12);
                        tt2.setFillColor(sf::Color(140,150,165));
                        tt2.setPosition({bx+bw-130.f,ly}); window.draw(tt2);
                    }
                    ly+=18.f;
                }
                // 削除の確認
                if(del_confirm>=0&&del_confirm<nf){
                    char cb[120];
                    snprintf(cb,120,"「%s」を削除します。もう一度 × を押してください",
                             save_files[del_confirm].c_str());
                    sf::Text ct(font,jp(cb),13);
                    ct.setFillColor(sf::Color(255,170,170));
                    ct.setPosition({bx+20.f,by+bh-24.f});
                    window.draw(ct);
                }
                // ホバーした一件の概要
                if(hover_save>=0&&hover_save<nf){
                    const std::string& nm=save_files[hover_save];
                    auto it=save_info.find(nm);
                    if(it==save_info.end()){
                        save_info[nm]=save_summary(nm+".txt");
                        it=save_info.find(nm);
                    }
                    std::string txt=it->second;
                    if(hover_save<(int)save_times.size()){
                        std::tm lt{}; std::time_t tt=save_times[hover_save];
                        #ifdef _WIN32
                          localtime_s(&lt,&tt);
                        #else
                          localtime_r(&tt,&lt);
                        #endif
                        char ts[40]; std::strftime(ts,40,"   %m/%d %H:%M",&lt);
                        txt+=ts;
                    }
                    sf::Text tt2(font,jp(txt.c_str()),13);
                    sf::FloatRect tr=tt2.getLocalBounds();
                    float tipx=std::min(mouse_x+14.f,SCREEN_W-tr.size.x-26.f);
                    float tipy=mouse_y+18.f;
                    sf::RectangleShape tbg({tr.size.x+20.f,26.f});
                    tbg.setPosition({tipx-10.f,tipy-4.f});
                    tbg.setFillColor(sf::Color(10,16,22,240));
                    tbg.setOutlineColor(sf::Color(130,150,170));
                    tbg.setOutlineThickness(1.f);
                    window.draw(tbg);
                    tt2.setFillColor(sf::Color(225,235,245));
                    tt2.setPosition({tipx,tipy}); window.draw(tt2);
                }
                if(nf>maxshow){
                    char m[48]; snprintf(m,48,"...他 %d 件",nf-maxshow);
                    sf::Text mt(font,jp(m),12); mt.setFillColor(sf::Color(150,150,150));
                    mt.setPosition({bx+20.f,ly}); window.draw(mt);
                }
            }
        }

        prof_end(8);
        // 気候マップの凡例
        if(view_mode==4 && font_ok){
            struct KL { const char* n; sf::Color c; };
            const KL kl[8]={
                {"Af 熱帯雨林",sf::Color(28,130,90)},  {"Aw サバナ",sf::Color(202,220,92)},
                {"BW 砂漠",   sf::Color(232,128,56)},  {"BS ステップ",sf::Color(236,202,110)},
                {"C  温帯",   sf::Color(108,200,108)}, {"D  亜寒帯", sf::Color(118,150,202)},
                {"ET ツンドラ",sf::Color(158,140,172)},{"EF 氷雪",   sf::Color(206,224,244)}
            };
            // 左下はヘルプのボタンが使う。その上へ逃がす
            float lx=18.f, ly=SCREEN_H-16.f-34.f-12.f-(8*20.f+6.f);
            sf::RectangleShape bg3({186.f,8*20.f+6.f});
            bg3.setPosition({lx-6.f,ly-4.f});
            bg3.setFillColor(sf::Color(0,0,0,195));
            bg3.setOutlineColor(sf::Color(120,120,120)); bg3.setOutlineThickness(1.f);
            window.draw(bg3);
            for(int k=0;k<8;k++){
                sf::RectangleShape sw({14.f,14.f});
                sw.setPosition({lx,ly+k*20.f}); sw.setFillColor(kl[k].c);
                window.draw(sw);
                sf::Text t(font,jp(kl[k].n),12);
                t.setFillColor(sf::Color(225,225,225));
                t.setPosition({lx+20.f,ly+k*20.f}); window.draw(t);
            }
        }
        // ===== 拡大グラフ =====
        if(graph_big>=0 && font_ok && !show_tree){
            // 背景を薄く覆う
            sf::RectangleShape vlq({(float)SCREEN_W,(float)SCREEN_H});
            vlq.setFillColor(sf::Color(0,0,0,170)); window.draw(vlq);
            // 横に広すぎると読みにくい。縦に三段積むので、幅は抑えて高さを取る
            float PW2=std::min(SCREEN_W*0.52f,760.f), PH2=SCREEN_H*0.80f;
            float PX2=(SCREEN_W-PW2)*0.5f, PY2=(SCREEN_H-PH2)*0.5f;
            sf::RectangleShape bg2({PW2,PH2});
            bg2.setPosition({PX2,PY2});
            bg2.setFillColor(sf::Color(8,14,12,245));
            bg2.setOutlineColor(sf::Color(130,150,160)); bg2.setOutlineThickness(1.f);
            window.draw(bg2);

            struct GDef { const char* name; Series* a; Series* b; Series* c;
                          const char* la; const char* lb; const char* lc;
                          float lo,hi; };
            GDef gd[7]={
                {"個体数",     &s_plant,&s_animal,&s_corpse,
                               "植物","動物","死がい",            0.f,-1.f},
                {"植物の形質", &s_maxh, &s_shade, &s_tough,
                               "最大樹高","耐陰性","防御",         0.f,-1.f},
                {"動物の形質", &s_reach,&s_diet,  &s_jaw,
                               "届く高さ","肉食性","顎",           0.f,-1.f},
                {"陸の平均気温",&s_temp, nullptr,  nullptr,
                               "気温",nullptr,nullptr,           -30.f,40.f},
                {"耐寒性",     &s_coldP,&s_coldA, nullptr,
                               "植物","動物",nullptr,              0.f,1.f},
                {"大気",       &s_co2,  &s_o2,    nullptr,
                               "CO2","O2",nullptr,                0.f,-1.f},
                {"生産者と消費者",&s_plant,&s_animal,&s_micro,
                               "植物","動物","分解者",             0.f,-1.f},
            };
            int gi=std::max(0,std::min(graph_big,6));
            GDef& g=gd[gi];

            sf::Text ttl(font,jp(g.name),20);
            ttl.setFillColor(sf::Color::White);
            ttl.setPosition({PX2+20.f,PY2+14.f}); window.draw(ttl);

            // 尺の切り替え
            float sx2=PX2+PW2-20.f;
            for(int s=N_SCALE-1;s>=0;s--){
                sf::Text st2(font,jp(SCALE_NAME[s]),14);
                sf::FloatRect sr=st2.getLocalBounds();
                float bw3=sr.size.x+22.f, bh3=26.f;
                sx2-=bw3+6.f;
                bool on=(graph_scale==s);
                sf::RectangleShape sb2({bw3,bh3});
                sb2.setPosition({sx2,PY2+14.f});
                sb2.setFillColor(on?sf::Color(40,90,60):sf::Color(28,34,38));
                sb2.setOutlineColor(on?sf::Color(140,230,160):sf::Color(100,110,120));
                sb2.setOutlineThickness(1.f);
                window.draw(sb2);
                st2.setFillColor(on?sf::Color(200,250,210):sf::Color(170,180,190));
                st2.setPosition({sx2+11.f,PY2+17.f}); window.draw(st2);
                ui_hits.push_back({sx2,PY2+14.f,bw3,bh3,200+s});
            }
            // 閉じる
            {
                float cw=28.f, cxp=PX2+PW2-20.f-cw, cyp=PY2+PH2-20.f-cw;
                sf::RectangleShape cb2({cw,cw});
                cb2.setPosition({cxp,cyp});
                cb2.setFillColor(sf::Color(60,30,30,220));
                cb2.setOutlineColor(sf::Color(190,120,120)); cb2.setOutlineThickness(1.f);
                window.draw(cb2);
                sf::Text ct2(font,jp("×"),18); ct2.setFillColor(sf::Color(255,200,200));
                sf::FloatRect cr=ct2.getLocalBounds();
                ct2.setPosition({cxp+cw*0.5f-cr.size.x*0.5f,cyp+2.f});
                window.draw(ct2);
                ui_hits.push_back({cxp,cyp,cw,cw,210});
            }

            // 折れ線
            // 二本を重ねると尺度が違って読めない。上下に分けて、各々に目盛を付ける
            const std::vector<float>& va=g.a->v[graph_scale];
            const std::vector<float>* vb=g.b? &g.b->v[graph_scale] : nullptr;
            const std::vector<float>* vc=g.c? &g.c->v[graph_scale] : nullptr;
            int npane=1+(vb?1:0)+(vc?1:0);
            float ax=PX2+62.f, aw2=PW2-82.f;
            float atop=PY2+60.f, atot=PH2-140.f;
            float agap=18.f;
            float ah=(atot-agap*(npane-1))/npane;
            int N2=(int)va.size();
            int win=std::min(N2,300);
            int maxoff=std::max(0,N2-win);
            if(graph_off>maxoff) graph_off=maxoff;
            if(graph_off<0) graph_off=0;
            int st3=N2-win-graph_off; if(st3<0) st3=0;
            char lb2[64];
            // 一段ぶんを描く。系列ごとに自分の目盛を持つので、尺度が違っても読める
            auto pane=[&](const std::vector<float>& v,sf::Color col,
                          const char* name,float py,float ph2){
                sf::RectangleShape fr3({aw2,ph2});
                fr3.setPosition({ax,py});
                fr3.setFillColor(sf::Color(0,0,0,0));
                fr3.setOutlineColor(sf::Color(70,80,90)); fr3.setOutlineThickness(1.f);
                window.draw(fr3);
                float lo=g.lo, hi=g.hi;
                if(hi<0.f){                        // 自動目盛
                    hi=1.f;
                    for(int i=st3;i<st3+win&&i<(int)v.size();i++)
                        if(v[i]>hi) hi=v[i];
                    hi*=1.1f;
                }
                for(int k=0;k<=4;k++){
                    float t=k/4.f, yy=py+ph2*(1.f-t);
                    sf::RectangleShape gl2({aw2,1.f});
                    gl2.setPosition({ax,yy}); gl2.setFillColor(sf::Color(50,58,66));
                    window.draw(gl2);
                    snprintf(lb2,64,"%.1f",lo+(hi-lo)*t);
                    sf::Text t2(font,jp(lb2),11);
                    t2.setFillColor(sf::Color(150,160,170));
                    t2.setPosition({PX2+14.f,yy-8.f}); window.draw(t2);
                }
                if((int)v.size()>=2){
                    sf::VertexArray ln(sf::PrimitiveType::LineStrip);
                    for(int i=st3;i<st3+win&&i<(int)v.size();i++){
                        float t=(v[i]-lo)/std::max(1e-6f,hi-lo);
                        sf::Vertex vt;
                        vt.position={ax+aw2*(i-st3)/(float)std::max(1,win-1),
                                     py+ph2*(1.f-clamp01(t))};
                        vt.color=col; ln.append(vt);
                    }
                    window.draw(ln);
                }
                sf::Text nt2(font,jp(name),14);
                nt2.setFillColor(col);
                nt2.setPosition({ax+8.f,py+4.f}); window.draw(nt2);
            };
            pane(va,sf::Color(120,235,140),g.la,atop,ah);
            if(vb) pane(*vb,sf::Color(255,170,90),g.lb,atop+ah+agap,ah);
            if(vc) pane(*vc,sf::Color(200,140,255),g.lc,atop+(ah+agap)*2,ah);
            float ay=atop;   // 以降の案内表示の基準
            snprintf(lb2,64,"スクロールで時間移動   %d / %d",graph_off,maxoff);
            { sf::Text t2(font,jp(lb2),12); t2.setFillColor(sf::Color(140,150,160));
              t2.setPosition({ax,ay+ah+32.f}); window.draw(t2); }
            ui_hits.push_back({ax,ay,aw2,ah,220});   // スクロール領域
        }
        // ===== 操作の案内 =====
        if(font_ok && !show_tree){
            // 個体パネルが左上に出るので、案内は左下に置く
            const float HB=34.f, HM=16.f;
            float hx=HM, hy=SCREEN_H-HM-HB;
            bool hov=(mouse_x>=hx&&mouse_x<=hx+HB&&mouse_y>=hy&&mouse_y<=hy+HB);
            sf::RectangleShape hb({HB,HB});
            hb.setPosition({hx,hy});
            hb.setFillColor(show_help?sf::Color(30,64,44,235)
                           :(hov?sf::Color(40,50,58,225):sf::Color(0,0,0,195)));
            hb.setOutlineColor(show_help||hov?sf::Color(140,225,160)
                                             :sf::Color(120,135,150));
            hb.setOutlineThickness(1.f);
            window.draw(hb);
            sf::Text ht(font,jp("?"),20);
            ht.setFillColor(show_help||hov?sf::Color(200,250,210):sf::Color(190,205,215));
            sf::FloatRect hr=ht.getLocalBounds();
            ht.setPosition({hx+HB*0.5f-hr.size.x*0.5f,hy+3.f});
            window.draw(ht);
            ui_hits.push_back({hx,hy,HB,HB,104});

            if(show_help){
                struct HL { const char* k; const char* d; };
                static const HL hl[]={
                    {"緑の点","植物"},
                    {"橙の点","動物(赤いほど肉食)"},
                    {"紫の点","分解者"},
                    {"クリック","生き物を選ぶ"},
                    {"W A S D","地図を動かす"},
                    {"Z / X","縮小・拡大"},
                    {"‹ ›","地図の切りかえ"},
                    {"N","昼夜の陰影を消す"},
                    {"C","種ごとに色分け"},
                    {"G","グラフの表示"},
                    {"Space","一時停止"},
                    {"K / L","保存・読み込み"},
                    {"R","動物を再投入"},
                    {"P","プランクトンを放つ"},
                };
                int nh=(int)(sizeof(hl)/sizeof(hl[0]));
                float pw2=252.f, rowh=22.f;
                float ph2=nh*rowh+20.f;
                // 気候の凡例と重ならないよう、ボタンの右へ開く
                float px2=hx+HB+8.f, py2=SCREEN_H-16.f-ph2;
                sf::RectangleShape pb({pw2,ph2});
                pb.setPosition({px2,py2});
                pb.setFillColor(sf::Color(0,0,0,225));
                pb.setOutlineColor(sf::Color(120,135,150));
                pb.setOutlineThickness(1.f);
                window.draw(pb);
                float ly2=py2+10.f;
                for(int i=0;i<nh;i++){
                    sf::Text kt(font,jp(hl[i].k),13);
                    kt.setFillColor(sf::Color(140,225,160));
                    kt.setPosition({px2+12.f,ly2}); window.draw(kt);
                    sf::Text dt2(font,jp(hl[i].d),13);
                    dt2.setFillColor(sf::Color(205,215,225));
                    dt2.setPosition({px2+96.f,ly2}); window.draw(dt2);
                    ly2+=rowh;
                }
            }
        }

        // ===== 下部バー =====
        if(font_ok && !show_tree){
            const float BH=72.f, BM=16.f, TW=96.f;
            float bw2=std::min(SCREEN_W*0.46f,620.f);   // 中身に合わせた幅
            float by0=SCREEN_H-BM-BH;
            float bx0=(SCREEN_W-bw2-TW-12.f)*0.5f;      // バーとボタンをまとめて中央へ
            // 本体
            sf::RectangleShape bar({bw2,BH});
            bar.setPosition({bx0,by0});
            bar.setFillColor(sf::Color(0,0,0,205));
            bar.setOutlineColor(sf::Color(110,120,130)); bar.setOutlineThickness(1.f);
            window.draw(bar);
            // 上段と下段を分ける線
            sf::RectangleShape sep({bw2-24.f,1.f});
            sep.setPosition({bx0+12.f,by0+BH*0.5f});
            sep.setFillColor(sf::Color(90,100,110)); window.draw(sep);

            char sb[96];
            // 左: 季節。両半球を併記する
            float dd=sun_dec*180.f/3.14159265f;
            const char* nseas = dd> 8.f?"夏":(dd<-8.f?"冬":(sun_dec>=0.f?"春":"秋"));
            const char* sseas = dd> 8.f?"冬":(dd<-8.f?"夏":(sun_dec>=0.f?"秋":"春"));
            snprintf(sb,96,"北 %s / 南 %s",nseas,sseas);
            { sf::Text t(font,jp(sb),13); t.setFillColor(sf::Color(200,215,225));
              t.setPosition({bx0+14.f,by0+8.f}); window.draw(t); }
            // 中央: tick
            snprintf(sb,96,"t = %d",frame);
            { sf::Text t(font,jp(sb),17); t.setFillColor(sf::Color::White);
              sf::FloatRect r=t.getLocalBounds();
              t.setPosition({bx0+bw2*0.5f-r.size.x*0.5f,by0+7.f}); window.draw(t); }
            // 右: 描画速度
            { static float ms_avg=16.f;
              double tot=0; for(int k=0;k<PROF_N;k++) tot+=prof[k];
              ms_avg=ms_avg*0.9f+(float)tot*0.1f;
              snprintf(sb,96,"%.1f ms/frame",ms_avg);
              sf::Text t(font,jp(sb),13); t.setFillColor(sf::Color(200,215,225));
              sf::FloatRect r=t.getLocalBounds();
              t.setPosition({bx0+bw2-14.f-r.size.x,by0+8.f}); window.draw(t); }

            // 下段: マップの切り替え
            {
                const char* vmn[11]={"地形","気温","大気","腐植","気候区分",
                                     "病原体","匂い","昼夜","栄養","傾斜","海流"};
                snprintf(sb,96,"%d / 11   %s",view_mode+1,vmn[view_mode]);
                sf::Text t(font,jp(sb),16); t.setFillColor(sf::Color(230,238,245));
                sf::FloatRect r=t.getLocalBounds();
                float cx=bx0+bw2*0.5f, ty=by0+BH*0.5f+8.f;
                t.setPosition({cx-r.size.x*0.5f,ty}); window.draw(t);
                // 前後の矢印
                float aw=28.f, ay=by0+BH*0.5f+4.f;
                float lxp=cx-r.size.x*0.5f-aw-14.f, rxp=cx+r.size.x*0.5f+14.f;
                for(int d=0;d<2;d++){
                    float axp=d?rxp:lxp;
                    bool hov=(mouse_x>=axp&&mouse_x<=axp+aw&&
                              mouse_y>=ay&&mouse_y<=ay+aw);
                    sf::RectangleShape ab({aw,aw});
                    ab.setPosition({axp,ay});
                    ab.setFillColor(hov?sf::Color(70,90,105):sf::Color(38,46,54));
                    ab.setOutlineColor(sf::Color(120,135,150)); ab.setOutlineThickness(1.f);
                    window.draw(ab);
                    sf::Text at(font,jp(d?"›":"‹"),20);
                    at.setFillColor(sf::Color::White);
                    sf::FloatRect ar=at.getLocalBounds();
                    at.setPosition({axp+aw*0.5f-ar.size.x*0.5f,ay+2.f});
                    window.draw(at);
                    ui_hits.push_back({axp,ay,aw,aw,d?101:100});
                }
            }

            // 系統樹ボタン。バーの外、右下の角に下端を揃えて置く
            {
                float tx=bx0+bw2+12.f, ty=by0, th=BH;
                bool hov=(mouse_x>=tx&&mouse_x<=tx+TW&&mouse_y>=ty&&mouse_y<=ty+th);
                sf::RectangleShape tb({TW,th});
                tb.setPosition({tx,ty});
                tb.setFillColor(hov?sf::Color(30,64,44,235):sf::Color(0,0,0,205));
                tb.setOutlineColor(hov?sf::Color(140,225,160):sf::Color(110,120,130));
                tb.setOutlineThickness(1.f);
                window.draw(tb);
                // 枝の記号
                float cx=tx+TW*0.5f, cy=ty+22.f;
                sf::Color lc=hov?sf::Color(160,240,175):sf::Color(150,175,190);
                auto ln=[&](float x1,float y1,float x2,float y2){
                    sf::VertexArray v(sf::PrimitiveType::Lines);
                    sf::Vertex a,b; a.position={x1,y1}; b.position={x2,y2};
                    a.color=lc; b.color=lc; v.append(a); v.append(b); window.draw(v);
                };
                ln(cx,cy-10.f,cx,cy+12.f);
                ln(cx-14.f,cy+2.f,cx+14.f,cy+2.f);
                ln(cx-14.f,cy+2.f,cx-14.f,cy+12.f);
                ln(cx+14.f,cy+2.f,cx+14.f,cy+12.f);
                sf::Text t(font,jp("系統樹"),13);
                t.setFillColor(hov?sf::Color(200,250,210):sf::Color(190,205,215));
                sf::FloatRect r=t.getLocalBounds();
                t.setPosition({cx-r.size.x*0.5f,ty+th-26.f});
                window.draw(t);
                ui_hits.push_back({tx,ty,TW,th,102});
            }
        }
        window.display();
        }   // !g_headless
    }
    prevent_sleep_end();
    return 0;
}

//起動コマンドは以下の通り
//clang++ -std=c++17 phase3.cpp -o phase3 $(pkg-config --cflags --libs sfml-graphics) -framework IOKit -framework CoreFoundation && ./phase3