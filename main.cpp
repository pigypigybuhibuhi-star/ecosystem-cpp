#include <SFML/Graphics.hpp>
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <chrono>
// 計測用
double t_prey_nn=0,t_pred_nn=0,t_move=0,t_eat=0,t_field=0,t_draw=0,t_grid=0;
int prof_frames=0;
auto now=[](){ return std::chrono::high_resolution_clock::now(); };
auto ms=[](auto a,auto b){ return std::chrono::duration<double,std::milli>(b-a).count(); };
// 画面と世界を分離
const float SCREEN_W=2240.f, SCREEN_H=1200.f;
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

struct Brain {
    int n_in; std::vector<float> W1,b1,W2,b2,W3,b3;
    Brain(int in=PREY_IN):n_in(in){
        W1.resize(n_in*N_H1);b1.resize(N_H1,0.f);W2.resize(N_H1*N_H2);b2.resize(N_H2,0.f);W3.resize(N_H2*N_OUT);b3.resize(N_OUT,0.f);
        float s1=1.f/std::sqrt((float)n_in),s2=1.f/std::sqrt((float)N_H1),s3=1.f/std::sqrt((float)N_H2);
        for(auto&w:W1)w=gauss(rng)*s1;for(auto&w:W2)w=gauss(rng)*s2;for(auto&w:W3)w=gauss(rng)*s3;
    }
    void forward(const float* in,float* h1,float* h2,float* out)const{
        for(int j=0;j<N_H1;j++){float s=b1[j];for(int i=0;i<n_in;i++)s+=in[i]*W1[i*N_H1+j];h1[j]=std::tanh(s);}
        for(int j=0;j<N_H2;j++){float s=b2[j];for(int i=0;i<N_H1;i++)s+=h1[i]*W2[i*N_H2+j];h2[j]=std::tanh(s);}
        for(int j=0;j<N_OUT;j++){float s=b3[j];for(int i=0;i<N_H2;i++)s+=h2[i]*W3[i*N_OUT+j];out[j]=std::tanh(s);}
    }
    Brain mutated(float rate)const{
        Brain c=*this;
        for(auto&w:c.W1)w+=gauss(rng)*rate;for(auto&w:c.b1)w+=gauss(rng)*rate;
        for(auto&w:c.W2)w+=gauss(rng)*rate;for(auto&w:c.b2)w+=gauss(rng)*rate;
        for(auto&w:c.W3)w+=gauss(rng)*rate;for(auto&w:c.b3)w+=gauss(rng)*rate;
        return c;
    }
};
struct Genes { float size, vision, eat_gain, speed; };
Genes random_genes(){return {BASE_SIZE*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION)),PREY_VISION*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION)),EAT_GAIN*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION)),SPEED*(1.f+frand(-SIZE_VARIATION,SIZE_VARIATION))};}
float mut1(float v,float w){float m=v*(1.f+frand(-w,w));return m<0.01f?0.01f:m;}
Genes mutate_genes(const Genes& g,float stress){float w=MUT_RATE_GENE*(1.f+stress*MUT_STRESS_FACTOR);return {mut1(g.size,w),mut1(g.vision,w),mut1(g.eat_gain,w),mut1(g.speed,w)};}

struct Agent {
    float x,y,vx,vy,energy,max_energy; bool alive,is_prey; int age,repro_counter;
    float fear,affinity,stress; Brain brain; Genes genes;
    Agent(int in):brain(in){}
    int id;
    int repro_count, food_count;
    int last_eat_frame;
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

struct PreyRec { int fit; Brain brain; Genes genes; };
std::vector<PreyRec> prey_hall;
std::vector<std::pair<int,Brain>> pred_hall;
int generation=1;

Agent make_prey(const Brain& b,const Genes& g){
    Agent a(PREY_IN);a.x=frand(0,WIDTH);a.y=frand(0,HEIGHT);a.vx=frand(-SPEED,SPEED);a.vy=frand(-SPEED,SPEED);
    a.energy=INIT_ENERGY;a.alive=true;a.is_prey=true;a.age=0;a.repro_count=0; a.food_count=0;a.repro_counter=0;a.fear=0;a.affinity=0;a.stress=0;a.id=next_id++;a.last_eat_frame=-1000;
    a.brain=b;a.genes=g;a.max_energy=g.size*MAXENERGY_COEF;return a;
}
Agent make_predator(const Brain& b){
    Agent a(PRED_IN);a.x=frand(0,WIDTH);a.y=frand(0,HEIGHT);a.vx=frand(-SPEED,SPEED);a.vy=frand(-SPEED,SPEED);
    a.energy=INIT_ENERGY;a.alive=true;a.is_prey=false;a.age=0;a.repro_count=0; a.food_count=0;a.repro_counter=0;a.fear=0;a.affinity=0;a.stress=0;a.id=next_id++;a.last_eat_frame=-1000;
    a.brain=b;a.genes={BASE_SIZE,PREDATOR_VISION,EAT_GAIN,SPEED};a.max_energy=MAX_ENERGY;return a;
}
void trim_prey(){std::sort(prey_hall.begin(),prey_hall.end(),[](auto&a,auto&b){return a.fit>b.fit;});if((int)prey_hall.size()>HALL_SIZE)prey_hall.resize(HALL_SIZE);}
void trim_pred(){std::sort(pred_hall.begin(),pred_hall.end(),[](auto&a,auto&b){return a.first>b.first;});if((int)pred_hall.size()>HALL_SIZE)pred_hall.resize(HALL_SIZE);}

int calc_fit(const Agent& a){
    return (int)(a.repro_count*100 + a.food_count*5 + a.age*0.1f);
}

int main(){
    sf::RenderWindow window(sf::VideoMode({2240,1200}),"Ecosystem C++",sf::Style::Titlebar|sf::Style::Close);
    sf::Font font;
    bool font_ok = font.openFromFile("/System/Library/Fonts/Helvetica.ttc");
    if(!font_ok) font_ok = font.openFromFile("/System/Library/Fonts/Supplemental/Arial.ttf");
    window.setFramerateLimit(60);

    std::vector<Agent> preys, predators;
    for(int i=0;i<NUM_PREY;i++)preys.push_back(make_prey(Brain(PREY_IN),random_genes()));
    for(int i=0;i<NUM_PREDATOR;i++)predators.push_back(make_predator(Brain(PRED_IN)));
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
    int selected_id=-1;

    float measured_fps=60.f;   // 実測FPS(初期値60)
    auto last_time=std::chrono::high_resolution_clock::now();

    int peak_prey=0, peak_pred=0, peak_plant=0, peak_meat=0;
    
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

    char prof_str[256]="";

    while(window.isOpen()){
        while(const std::optional event=window.pollEvent()){
            if(event->is<sf::Event::Closed>())window.close();
            if(const auto* k=event->getIf<sf::Event::KeyPressed>()){
                if(k->code==sf::Keyboard::Key::G)show_nutrient=!show_nutrient;
                if(k->code==sf::Keyboard::Key::F)show_phero=!show_phero;
                if(k->code==sf::Keyboard::Key::Space)paused=!paused;
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

        auto _t0=now();
        for(auto& p:preys){
            if(!p.alive)continue;
            float vp[NUM_RAYS]={0},vpr[NUM_RAYS]={0};
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
            }
            float in[PREY_IN]; for(int i=0;i<NUM_RAYS;i++){in[i]=vp[i];in[NUM_RAYS+i]=vpr[i];} in[24]=p.energy/p.max_energy;
            float h1[N_H1],h2[N_H2],out[N_OUT]; p.brain.forward(in,h1,h2,out);
            p.vx=out[0]*p.genes.speed;p.vy=out[1]*p.genes.speed;
        }
        t_prey_nn += ms(_t0,now());
        auto _t1=now();
        for(auto& pd:predators){
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
            float in[PRED_IN];for(int i=0;i<NUM_RAYS;i++)in[i]=vpr[i];in[12]=pd.energy/pd.max_energy;
            float h1[N_H1],h2[N_H2],out[N_OUT];pd.brain.forward(in,h1,h2,out);
            pd.vx=out[0]*SPEED;pd.vy=out[1]*SPEED;
        }
        t_pred_nn += ms(_t1,now());

        auto _t2=now();
        for(auto& p:preys){
            if(!p.alive)continue;
            p.x+=p.vx;p.y+=p.vy;
            if(p.x<0)p.x+=WIDTH;if(p.x>=WIDTH)p.x-=WIDTH;if(p.y<0)p.y+=HEIGHT;if(p.y>=HEIGHT)p.y-=HEIGHT;
            p.age++;
            float sp=std::sqrt(p.vx*p.vx+p.vy*p.vy);
            p.energy-=p.genes.size*METABOLISM_COEF+0.5f*p.genes.size*sp*sp*MOVE_COEF;
            if(p.energy<=0){p.alive=false;meats.push_back({p.x,p.y,p.genes.size*MEAT_SIZE_COEF+std::max(0.f,p.energy)});continue;}
            bool pred_sight=false,ally_sight=false; float pv=p.genes.vision;
            int scx=(int)p.x/CELL_SIZE,scy=(int)p.y/CELL_SIZE;int cr=(int)std::ceil(pv/CELL_SIZE);
            for(int dx=-cr;dx<=cr&&!(pred_sight&&ally_sight);dx++)for(int dy=-cr;dy<=cr;dy++){
                int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
                for(int pi:pred_grid[cell_index(cx,cy)]){if(!predators[pi].alive)continue;
                    float ddx=torus_delta(p.x,predators[pi].x,WIDTH),ddy=torus_delta(p.y,predators[pi].y,HEIGHT);
                    if(ddx*ddx+ddy*ddy<pv*pv){pred_sight=true;break;}}
                for(int pi:prey_grid[cell_index(cx,cy)]){
                    float ddx=torus_delta(p.x,preys[pi].x,WIDTH),ddy=torus_delta(p.y,preys[pi].y,HEIGHT);
                    float d2=ddx*ddx+ddy*ddy;if(d2>1&&d2<pv*pv){ally_sight=true;break;}}
            }
            if(pred_sight)p.fear+=FEAR_RISE;else p.fear-=FEAR_DECAY;p.fear=clamp01(p.fear);
            if(ally_sight)p.affinity+=AFFINITY_RISE;else p.affinity-=AFFINITY_DECAY;p.affinity=clamp01(p.affinity);
            ph_fear.emit(PH_CELL,p.x,p.y,p.fear*0.1f); ph_aff.emit(PH_CELL,p.x,p.y,p.affinity*0.1f);
            float red=ph_fear.get(PH_CELL,p.x,p.y),blue=ph_aff.get(PH_CELL,p.x,p.y);
            bool R=red>=PHERO_THRESHOLD,B=blue>=PHERO_THRESHOLD;
            if(R&&B)p.stress+=STRESS_MIX;else if(R)p.stress+=STRESS_RISE;else if(B)p.stress-=STRESS_FALL;else p.stress-=STRESS_DECAY;
            if(pred_sight)p.stress+=STRESS_PREDATOR_SIGHT; p.stress=clamp01(p.stress);
        }
        for(auto& pd:predators){
            if(!pd.alive)continue;
            pd.x+=pd.vx;pd.y+=pd.vy;
            if(pd.x<0)pd.x+=WIDTH;if(pd.x>=WIDTH)pd.x-=WIDTH;if(pd.y<0)pd.y+=HEIGHT;if(pd.y>=HEIGHT)pd.y-=HEIGHT;
            pd.age++;pd.energy-=METABOLISM;
            if(pd.energy<=0){pd.alive=false;meats.push_back({pd.x,pd.y,pd.genes.size*MEAT_SIZE_COEF+std::max(0.f,pd.energy)});}
        }
        t_move += ms(_t2,now());

        auto _t3=now();
        for(auto& p:preys){
            if(!p.alive)continue;
            int scx=(int)p.x/CELL_SIZE,scy=(int)p.y/CELL_SIZE;
            for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++){
                int cx=scx+dx,cy=scy+dy;if(cx<0||cx>=GRID_COLS||cy<0||cy>=GRID_ROWS)continue;
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
                        preys[pi].alive=false; pd.food_count++;prey_hall.push_back({calc_fit(preys[pi]),preys[pi].brain,preys[pi].genes});
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

        std::vector<Agent> babies;
        for(auto& p:preys){
            if(!p.alive)continue;
            if(p.energy>p.max_energy*(5.f/6.f)){if(++p.repro_counter>=PREY_REPRO_TIME){p.repro_counter=0;
                float success=1.f-p.stress*STRESS_REPRO_PENALTY;
                if(dist01(rng)<success){float cost=p.genes.size*REPRO_COST_COEF;p.energy-=cost; p.repro_count++;
                    Agent c=make_prey(p.brain.mutated(MUT_RATE),mutate_genes(p.genes,p.stress));
                    c.x=p.x+frand(-10,10);c.y=p.y+frand(-10,10);c.energy=cost*0.5f;babies.push_back(std::move(c));}
            }}else p.repro_counter=0;
        }
        for(auto&b:babies)preys.push_back(std::move(b));
        std::vector<Agent> pbabies;
        for(auto& pd:predators){
            if(!pd.alive)continue;
            if(pd.energy>pd.max_energy*(5.f/6.f)){if(++pd.repro_counter>=PREDATOR_REPRO_TIME){pd.repro_counter=0;pd.energy*=0.5f;pd.repro_count++;
                Agent c=make_predator(pd.brain.mutated(MUT_RATE));c.x=pd.x+frand(-10,10);c.y=pd.y+frand(-10,10);c.energy=pd.energy;pbabies.push_back(std::move(c));
            }}else pd.repro_counter=0;
        }
        for(auto&b:pbabies)predators.push_back(std::move(b));

        {
            std::vector<Agent> al;
            for(auto&p:preys){if(p.alive)al.push_back(std::move(p));else prey_hall.push_back({calc_fit(p),p.brain,p.genes});}
            preys=std::move(al);
            std::vector<Agent> ap;
            for(auto&pd:predators){if(pd.alive)ap.push_back(std::move(pd));else pred_hall.push_back({calc_fit(pd),pd.brain});}
            predators=std::move(ap);
        }
        if((int)prey_hall.size()>HALL_SIZE*10)trim_prey();
        if((int)pred_hall.size()>HALL_SIZE*10)trim_pred();

        if(preys.empty()||predators.empty()){
            for(auto&p:preys)prey_hall.push_back({calc_fit(p),p.brain,p.genes});
            for(auto&pd:predators)pred_hall.push_back({calc_fit(pd),pd.brain});
            trim_prey();trim_pred();
            int bf=prey_hall.empty()?0:prey_hall[0].fit,bpf=pred_hall.empty()?0:pred_hall[0].first;
            printf("Gen %d: survived %ds, prey fit %d, pred fit %d\n",generation,frame/60,bf,bpf);
            generation++;preys.clear();predators.clear();meats.clear();
            nutrient.grid.assign(NUT_COLS*NUT_ROWS,0.f);ph_fear.grid.assign(PH_COLS*PH_ROWS,0.f);ph_aff.grid.assign(PH_COLS*PH_ROWS,0.f);
            for(int i=0;i<NUM_PREY;i++){if(!prey_hall.empty()){int ix=(int)(dist01(rng)*prey_hall.size());if(ix>=(int)prey_hall.size())ix=prey_hall.size()-1;preys.push_back(make_prey(prey_hall[ix].brain.mutated(MUT_RATE),prey_hall[ix].genes));}else preys.push_back(make_prey(Brain(PREY_IN),random_genes()));}
            for(int i=0;i<NUM_PREDATOR;i++){if(!pred_hall.empty()){int ix=(int)(dist01(rng)*pred_hall.size());if(ix>=(int)pred_hall.size())ix=pred_hall.size()-1;predators.push_back(make_predator(pred_hall[ix].second.mutated(MUT_RATE)));}else predators.push_back(make_predator(Brain(PRED_IN)));}
            plants.clear();for(int i=0;i<NUM_PLANT;i++)plants.push_back({frand(0,WIDTH),frand(0,HEIGHT),0});
            frame=0;
        }
    }
        // WASDでカメラ移動
        float pan=15.f/zoom;   // ズームに応じて移動量調整
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::W))cam_y-=pan;
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::S))cam_y+=pan;
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::A))cam_x-=pan;
        if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::D))cam_x+=pan;
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
        pc.setFillColor(sf::Color(0,100,255));
        for(const auto& p:preys){if(!on_screen(p.x,p.y))continue;float r=std::max(1.f,p.genes.size*zoom);pc.setRadius(r);pc.setPosition({w2s_x(p.x)-r,w2s_y(p.y)-r});window.draw(pc);}
        rc.setFillColor(sf::Color(255,50,50));
        for(const auto& pd:predators){if(!on_screen(pd.x,pd.y))continue;float r=std::max(1.f,3.f*zoom);rc.setRadius(r);rc.setPosition({w2s_x(pd.x)-r,w2s_y(pd.y)-r});window.draw(rc);}
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
                    for(int i=0;i<NUM_RAYS;i++){
                        float ang;
                        if(sel->is_prey) ang = facing + TWO_PI*((float)i/NUM_RAYS);
                        else ang = facing - PREDATOR_FOV/2 + PREDATOR_FOV*((float)i/(NUM_RAYS-1));
                        // このレイに映ってるもの(prey視覚:植物in[i], 捕食者in[12+i])
                        sf::Color rc2(60,60,60,120);
                        if(sel->is_prey){
                            float plantv=in[i], predv=in[NUM_RAYS+i];
                            if(predv>plantv && predv>0.01f) rc2=sf::Color(255,60,60,200);
                            else if(plantv>0.01f) rc2=sf::Color(60,255,60,200);
                        } else {
                            if(in[i]>0.01f) rc2=sf::Color(255,120,60,200);
                        }
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
                // --- forward ---
                float h1[N_H1],h2[N_H2],out[N_OUT];
                sel->brain.forward(in.data(),h1,h2,out);

                // --- レイアウト ---
                float PX=20.f, PY=20.f, PW=SCREEN_W/4.f*0.8f, PH=SCREEN_H/2.f*0.8f;
                // 背景パネル
                sf::RectangleShape panel({PW,PH}); panel.setPosition({PX,PY});
                panel.setFillColor(sf::Color(0,0,0,180)); panel.setOutlineColor(sf::Color(80,80,80)); panel.setOutlineThickness(1.f);
                window.draw(panel);
                // 各層のノードのx座標
                float colx[4]={PX+PW*0.10f, PX+PW*0.37f, PX+PW*0.63f, PX+PW*0.90f};
                int counts[4]={NIN,N_H1,N_H2,N_OUT};
                // 各層ノードのy座標を計算する関数
                auto nodeY=[&](int layer,int i)->float{
                    int n=counts[layer]; float top=PY+20.f, bot=PY+PH-20.f;
                    if(n==1)return (top+bot)/2;
                    return top+(bot-top)*i/(n-1);
                };
                // 活性値の配列ポインタ
                const float* acts[4]={in.data(),h1,h2,out};

                // --- 線(重み)を先に描く ---
                auto drawLines=[&](int L,const std::vector<float>& W,int nin,int nout){
                    // W: nin*nout, W[i*nout+j]
                    sf::VertexArray va(sf::PrimitiveType::Lines);
                    for(int i=0;i<nin;i++)for(int j=0;j<nout;j++){
                        float w=W[i*nout+j]; float aw=std::fabs(w);
                        if(aw<0.15f)continue;   // 弱い線は省く(間引き)
                        int alpha=(int)std::min(120.f,aw*100.f);
                        sf::Color c = w>0 ? sf::Color(80,160,255,alpha) : sf::Color(255,80,80,alpha);
                        sf::Vertex v1,v2;
                        v1.position={colx[L],nodeY(L,i)}; v1.color=c;
                        v2.position={colx[L+1],nodeY(L+1,j)}; v2.color=c;
                        va.append(v1); va.append(v2);
                    }
                    window.draw(va);
                };
                drawLines(0,sel->brain.W1,NIN,N_H1);
                drawLines(1,sel->brain.W2,N_H1,N_H2);
                drawLines(2,sel->brain.W3,N_H2,N_OUT);

                // --- ノードを描く ---
                for(int L=0;L<4;L++){
                    for(int i=0;i<counts[L];i++){
                        float a=acts[L][i]; // -1..1
                        int inten=(int)std::min(255.f,std::fabs(a)*255.f);
                        sf::Color c = a>=0 ? sf::Color(60,120+inten/2,255,255) : sf::Color(255,80,80,255);
                        if(a>=0) c=sf::Color(40,80+inten*0.7f,255);
                        else c=sf::Color(255,60,60);
                        float rad=std::max(2.f,PW*0.012f);
                        sf::CircleShape node(rad); node.setFillColor(c);
                        node.setPosition({colx[L]-rad,nodeY(L,i)-rad});
                        window.draw(node);
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
       
        if(font_ok && prof_str[0]!='\0'){
            sf::Text pt(font,prof_str,14); pt.setFillColor(sf::Color::White);
            pt.setPosition({12.f, SCREEN_H-26.f});
            window.draw(pt);
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