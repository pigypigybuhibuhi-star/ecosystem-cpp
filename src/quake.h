// 地震の揺れ(phase3.cpp が順に #include する)
// 区画の外の断層で起きる地震。震源の位置と断層の性質は条件、起き方は G-R 則と改良大森則、
// 揺れは距離減衰式、斜面の崩れはニューマーク法で決める。崩れそのものは slide_at に渡す
struct Mainshock { long tick; float mag, x; };   // 余震を生む本震(起きた tick、規模、震央の x(マス))
std::vector<Mainshock> mainshocks;
std::mt19937 qk_rng(12344u^0x85EBCA6Bu);         // 地震だけの乱数。地震が起きても生き物の乱数列は動かない
std::uniform_real_distribution<double> qk01(0.0,1.0);
long  quake_events=0;
long  g_quake_force_t=-1; float g_quake_force_m=0.f;   // --quake T M: tick T に規模 M の本震を区画の中央の北に起こす(条件)
const float QK_FAULT_Y  = -(float)(QK_FAULT_DIST_KM*1e3/CELL_M);              // 断層の y(マス)。北の端の外
const float QK_FAULT_X0 = GRID_W*0.5f-(float)(QK_FAULT_LEN_KM*1e3/CELL_M)*0.5f;
const float QK_FAULT_X1 = GRID_W*0.5f+(float)(QK_FAULT_LEN_KM*1e3/CELL_M)*0.5f;
// 打ち切った G-R 則から規模を引く(逆関数法)
inline double qk_mag(double b){
    double span=1.0-std::pow(10.0,-b*(QK_MMAX-QK_MMIN));
    return QK_MMIN-std::log10(1.0-qk01(qk_rng)*span)/b;
}
// 破壊の範囲(断層に沿った x、マス)
inline void qk_rupture(float mag,float ex,double& xa,double& xb){
    double L=std::pow(10.0,-2.44+0.59*mag)*1e3/CELL_M;
    xa=std::max((double)QK_FAULT_X0,ex-0.5*L); xb=std::min((double)QK_FAULT_X1,ex+0.5*L);
}
// 区画の各マスの揺れを求め、限界加速度を超えた斜面をニューマーク法の確率で崩す。区画の中の最大加速度(g)を返す
// 揺れはその tick の地面の加速度(g)として shake_g に、起きた tick を shake_t に残し、動物の感覚器(揺れの帯)が読む。
// 動物は ANIMAL_STRIDE tick ごとに交互に処理されるので、揺れは起きた後の ANIMAL_STRIDE tick の間だけ読める(shake_felt)。
// こうすると、どの位相の個体も一度ずつ感じ取る
std::vector<float> shake_g(GRID_W*GRID_H,0.f);
std::vector<int>   shake_t(GRID_W*GRID_H,-(1<<30));
inline float shake_felt(int ci){
    int age=g_frame-shake_t[ci];
    return (age>0 && age<=ANIMAL_STRIDE)? shake_g[ci] : 0.f;
}
float quake_shake(float mag,float ex,bool slides=true){
    double xa,xb; qk_rupture(mag,ex,xa,xb);
    float amax=0.f;
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        double x=ci%GRID_W+0.5, y=ci/GRID_W+0.5;
        double dx= x<xa? xa-x : (x>xb? x-xb : 0.0), dy=y-QK_FAULT_Y;
        double d=std::sqrt(dx*dx+dy*dy)*CELL_M*1e-3, r=std::sqrt(d*d+7.3*7.3);
        double A=std::pow(10.0,-1.02+0.249*mag-std::log10(r)-0.00255*r);
        if(A>amax) amax=(float)A;
        if(g_frame-shake_t[ci]>ANIMAL_STRIDE-1) shake_g[ci]=0.f;   // 前の揺れが全ての位相に読まれ終えていれば消す
        if((float)A>shake_g[ci]) shake_g[ci]=(float)A;
        shake_t[ci]=g_frame;
        if(!slides) continue;
        double s=slope[ci];
        if(sea[ci]||s<=0.0) continue;
        double sa=s/std::sqrt(1.0+s*s);
        double fs=slope_fs(ci), ac=(fs-1.0)*sa;
        if(A<=ac) continue;                         // 限界加速度に届かない斜面は動かない(乱数を引かない)
        double P=0.335;
        if(ac>0.0){
            double q=ac/A, dn=std::pow(10.0,0.215)*std::pow(1.0-q,2.341)*std::pow(q,-1.438);
            P=0.335*(1.0-std::exp(-0.048*std::pow(dn,1.565)));
        }
        if(qk01(qk_rng)<P) slide_at(ci);
    }
    return amax;
}
void quake_at(float mag,float ex,bool main){
    quake_events++;
    ev_record(EV_QUAKE,(int)std::floor(ex),(int)std::floor(QK_FAULT_Y),mag);
    if(main) mainshocks.push_back({(long)g_frame,mag,ex});
    long s0=slide_events;
    float a=quake_shake(mag,ex);
    printf("  quake t=%d M%.2f %s x=%.0f | pga max %.3f g | slides %ld\n",
           g_frame,mag,main?"main":"after",ex,a,slide_events-s0);
}
void update_quakes(){
    if(g_replay){                      // 再生中は記録の各件だけ。再生した揺れからは崩れを起こさない(崩れは記録の各件で再生する)
        // 揺れそのもの(感覚の物理)は通常通り起こし、動物が感じ取れるようにする
        ev_replay_due(EV_QUAKE,[](const EnvEvent& e){ quake_events++; quake_shake(e.strength,e.x+0.5f,false); },false);
        return;
    }
    if(g_frame==g_quake_force_t) quake_at(g_quake_force_m,GRID_W*0.5f,true);
    // 本震。断層に沿って一様に起きる
    {
        double e=QK_RATE*DT;
        int n=(int)e; double f=e-n;
        if(f>0.0 && qk01(qk_rng)<f) n++;
        for(int k=0;k<n;k++){
            float m=(float)qk_mag(QK_B);
            float ex=QK_FAULT_X0+(QK_FAULT_X1-QK_FAULT_X0)*(float)qk01(qk_rng);
            quake_at(m,ex,true);
        }
    }
    // 余震。この tick の区間で改良大森則の率を積分した期待値
    if(mainshocks.empty()) return;
    static std::vector<double> w; w.assign(mainshocks.size(),0.0);
    double e=0.0;
    for(size_t i=0;i<mainshocks.size();i++){
        const Mainshock& s=mainshocks[i];
        double t1=(g_frame-s.tick)*DT/DAY_MIN, t0=t1-DT/DAY_MIN;   // tick の差 × DT で分に直してから日に
        if(t0<0.0) continue;           // 本震と同じ tick には余震を置かない
        double I=(std::pow(t0+QK_RJ_C,1.0-QK_RJ_P)-std::pow(t1+QK_RJ_C,1.0-QK_RJ_P))/(QK_RJ_P-1.0);
        w[i]=std::pow(10.0,QK_RJ_A+QK_RJ_B*(s.mag-QK_MMIN))*I;
        e+=w[i];
    }
    int n=(int)e; double f=e-n;
    if(f>0.0 && qk01(qk_rng)<f) n++;
    for(int k=0;k<n;k++){
        double u=qk01(qk_rng)*e, c=0.0; size_t j=0;
        for(;j+1<mainshocks.size();j++){ c+=w[j]; if(u<c) break; }
        double xa,xb; qk_rupture(mainshocks[j].mag,mainshocks[j].x,xa,xb);
        float m=(float)qk_mag(QK_RJ_B);
        float ex=(float)(xa+(xb-xa)*qk01(qk_rng));
        quake_at(m,ex,false);
    }
}
