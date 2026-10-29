// 描画(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 神経網の図。入力がいくつあり、どう繋がっているかを見せる
// SFML は const char* を Latin-1 として読む。UTF-8 から明示的に変換する
inline sf::String jp(const char* s){
    std::string t(s);
    return sf::String::fromUtf8(t.begin(), t.end());
}

void draw_nn(sf::RenderWindow& win,const sf::Font& font,const Animal& a,
             float x0,float y0,float w,float h)
{
    float lx=x0+58.f, mx2=x0+w*0.56f, rx=x0+w-30.f;
    float ystep=h/(float)std::max(1,N_IN);
    // 結線。太さと色で重みを表す
    auto link=[&](float ax,float ay,float bx,float by,float wv){
        if(std::fabs(wv)<0.10f) return;
        int th=1+(int)std::min(3.f,std::fabs(wv));
        sf::Color c = wv>0 ? sf::Color(110,210,255,(int)std::min(200.f,70+55*std::fabs(wv)))
                           : sf::Color(255,120,110,(int)std::min(200.f,70+55*std::fabs(wv)));
        for(int q=0;q<th;q++){
            sf::VertexArray ln(sf::PrimitiveType::Lines);
            sf::Vertex v0,v1;
            v0.position={ax,ay+q*0.9f}; v0.color=c;
            v1.position={bx,by+q*0.9f}; v1.color=c;
            ln.append(v0); ln.append(v1); win.draw(ln);
        }
    };
    // 働いている層だけを描く
    int actL[MAX_LAYER]={0}, ndep=0;
    for(int L=0;L<MAX_LAYER;L++){
        for(int hh=0;hh<MAX_NODE;hh++)
            if(a.ngain[L*MAX_NODE+hh]>NODE_ON) actL[L]++;
        if(actL[L]) ndep++;
    }
    float colw=(rx-lx)/(float)(ndep+1);
    float px[MAX_LAYER]; int col=0;
    for(int L=0;L<MAX_LAYER;L++){
        px[L]=lx+colw*(col+1);
        if(actL[L]) col++;
    }
    float nstep=h/(float)MAX_NODE;
    int prevL=-1;
    for(int L=0;L<MAX_LAYER;L++){
        if(!actL[L]) continue;
        for(int hh=0;hh<MAX_NODE;hh++){
            float g=a.ngain[L*MAX_NODE+hh];
            if(g<=NODE_ON) continue;
            float hy=y0+nstep*(hh+0.5f);
            if(prevL<0){
                for(int k=0;k<N_IN;k++)
                    link(lx+4.f,y0+ystep*(k+0.5f),px[L]-4.f,hy,a.w[hh*N_IN+k]);
            } else {
                int off=W_L0+(L-1)*MAX_NODE*MAX_NODE+hh*MAX_NODE;
                for(int j=0;j<MAX_NODE;j++){
                    if(a.ngain[prevL*MAX_NODE+j]<=NODE_ON) continue;
                    link(px[prevL]+4.f,y0+nstep*(j+0.5f),px[L]-4.f,hy,a.w[off+j]);
                }
            }
            sf::CircleShape d(3.f+2.f*g); d.setOrigin({3.f+2.f*g,3.f+2.f*g});
            d.setPosition({px[L],hy});
            d.setFillColor(sf::Color(170,150,220)); win.draw(d);
        }
        prevL=L;
    }
    // 最終層 → 出力
    if(prevL>=0){
        for(int hh=0;hh<MAX_NODE;hh++){
            if(a.ngain[prevL*MAX_NODE+hh]<=NODE_ON) continue;
            for(int o=0;o<2;o++)
                link(px[prevL]+4.f,y0+nstep*(hh+0.5f),rx-4.f,
                     y0+h*(o==0?0.32f:0.68f),a.w[W_L0+W_LL+o*MAX_NODE+hh]);
        }
    }
    // 入力 → 出力の直結
    for(int k=0;k<N_IN;k++)
        for(int o=0;o<2;o++)
            link(lx+4.f,y0+ystep*(k+0.5f),rx-4.f,y0+h*(o==0?0.32f:0.68f),
                 a.w[W_L0+W_LL+W_OUT+o*N_IN+k]);
    // 入力の節
    for(int k=0;k<N_IN;k++){
        float iy=y0+ystep*(k+0.5f);
        bool active=true; const char* nm=""; sf::Color nc(150,150,158);
        char buf[32];
        if(k<N_SENSOR){
            const Animal::Sensor& S=a.sen[k];
            active=(S.range>=SENSOR_MIN);
            if(active){
                        static const char* TN[N_TUNE]={
                            "植物","小型","大型","死骸","炎","林冠","病原",
                            "匂:生体","匂:腐敗","匂:食害","匂:煙",
                            "同種","水"};
                        int bd=(int)(S.tune*N_TUNE); if(bd>N_TUNE-1) bd=N_TUNE-1;
                        snprintf(buf,32,"%s %+.0f r%.1f",TN[bd],S.angle*57.3f,S.range);
                        nc = (bd>=7) ? sf::Color(220,170,255)
                           : (bd==0) ? sf::Color(120,230,120)
                           : (bd<=2) ? sf::Color(255,160,90)
                                     : sf::Color(190,140,110);
            } else snprintf(buf,32,"(none)");
            nm=buf;
        } else {
            static const char* IN[5]={"体力","気温","水分","体格","常時"};
            nm=IN[k-N_SENSOR]; nc=sf::Color(190,190,200);
        }
        sf::CircleShape d(active?5.f:3.f); d.setOrigin({active?5.f:3.f,active?5.f:3.f});
        d.setPosition({lx,iy});
        d.setFillColor(active?nc:sf::Color(70,70,76)); win.draw(d);
        sf::Text t(font,jp(nm),10);
        t.setFillColor(active?sf::Color(210,210,218):sf::Color(110,110,118));
        sf::FloatRect bb=t.getLocalBounds();
        t.setPosition({lx-10.f-bb.size.x,iy-6.f}); win.draw(t);
    }
    // 出力の節
    const char* ON[2]={"旋回","前進"};
    for(int o=0;o<2;o++){
        float oy=y0+h*(o==0?0.32f:0.68f);
        sf::CircleShape d(6.f); d.setOrigin({6.f,6.f});
        d.setPosition({rx,oy}); d.setFillColor(sf::Color(255,220,140)); win.draw(d);
        sf::Text t(font,jp(ON[o]),10);
        t.setFillColor(sf::Color(220,220,200));
        t.setPosition({rx+10.f,oy-6.f}); win.draw(t);
    }
}

// 型ごとの姿を手続き的に描く。h は画面上の高さ
void draw_plant_form(sf::RenderWindow& win,int form,float x,float baseY,
                     float h,sf::Color col)
{
    if(h<3.f) h=3.f;
    switch(form){
    case 1: {   // 寄生植物: 宿主に吸器を差し込む
        sf::RectangleShape host({std::max(2.f,h*0.09f),h*0.9f});
        host.setPosition({x+h*0.16f,baseY-h*0.9f});
        host.setFillColor(sf::Color(78,104,72)); win.draw(host);
        sf::VertexArray tw(sf::PrimitiveType::LineStrip);
        for(int t=0;t<=16;t++){
            float u=t/16.f;
            sf::Vertex v;
            v.position={x-h*0.10f+std::sin(u*7.f)*h*0.09f, baseY-u*h*0.75f};
            v.color=col; tw.append(v);
        }
        win.draw(tw);
        for(int k=1;k<=3;k++){
            sf::VertexArray hp(sf::PrimitiveType::Lines);
            sf::Vertex v0,v1;
            float yy=baseY-h*0.22f*k;
            v0.position={x-h*0.06f,yy}; v0.color=col;
            v1.position={x+h*0.18f,yy}; v1.color=sf::Color(214,120,120);
            hp.append(v0); hp.append(v1); win.draw(hp);
        }
        break; }
    case 2: {   // 食虫植物: 漏斗状の捕虫器
        sf::RectangleShape st({std::max(1.5f,h*0.045f),h*0.55f});
        st.setPosition({x-h*0.022f,baseY-h*0.55f});
        st.setFillColor(col); win.draw(st);
        sf::ConvexShape pit; pit.setPointCount(4);
        pit.setPoint(0,{x-h*0.20f,baseY-h*0.95f});
        pit.setPoint(1,{x+h*0.20f,baseY-h*0.95f});
        pit.setPoint(2,{x+h*0.09f,baseY-h*0.52f});
        pit.setPoint(3,{x-h*0.09f,baseY-h*0.52f});
        pit.setFillColor(sf::Color(158,74,96)); win.draw(pit);
        sf::RectangleShape lid({h*0.42f,h*0.05f});
        lid.setPosition({x-h*0.21f,baseY-h*1.02f});
        lid.setFillColor(sf::Color(190,96,110)); win.draw(lid);
        break; }
    case 3: {   // 着生植物: 幹を持たず、高い位置に張りつく
        sf::RectangleShape host({std::max(2.f,h*0.10f),h});
        host.setPosition({x-h*0.05f,baseY-h});
        host.setFillColor(sf::Color(88,72,56)); win.draw(host);
        for(int s=-1;s<=1;s++){
            sf::CircleShape lf(h*0.13f); lf.setOrigin({h*0.13f,h*0.13f});
            lf.setScale({1.7f,0.8f});
            lf.setPosition({x+s*h*0.19f,baseY-h*0.74f});
            lf.setFillColor(col); win.draw(lf);
        }
        break; }
    case 0: {   // プランクトン / 海藻: 揺れる葉状体
        for(int s=-1;s<=1;s++){
            sf::VertexArray f(sf::PrimitiveType::LineStrip);
            for(int t=0;t<=10;t++){
                float u=t/10.f;
                sf::Vertex v;
                v.position={x+s*h*0.14f+std::sin(u*4.2f+s)*h*0.13f, baseY-u*h};
                v.color=sf::Color(col.r,col.g,col.b,(std::uint8_t)(255-u*60));
                f.append(v);
            }
            win.draw(f);
        }
        break; }
    case 4: {   // つる: 支えに巻きついて登る
        // 支え(枯れた幹の見立て)
        sf::RectangleShape sup({std::max(1.5f,h*0.05f),h});
        sup.setPosition({x-h*0.025f,baseY-h});
        sup.setFillColor(sf::Color(96,82,64)); win.draw(sup);
        // 巻きつく茎
        sf::VertexArray vine(sf::PrimitiveType::LineStrip);
        for(int t=0;t<=26;t++){
            float u=t/26.f;
            sf::Vertex v;
            v.position={x+std::sin(u*10.5f)*h*0.13f, baseY-u*h};
            v.color=col; vine.append(v);
        }
        win.draw(vine);
        // 葉
        for(int k=1;k<=5;k++){
            float u=k/5.5f;
            float lx2=x+std::sin(u*10.5f)*h*0.13f;
            sf::CircleShape lf(h*0.075f); lf.setOrigin({h*0.075f,h*0.075f});
            lf.setScale({1.5f,0.7f});
            lf.setPosition({lx2+((k%2)?h*0.11f:-h*0.11f),baseY-u*h});
            lf.setFillColor(col); win.draw(lf);
        }
        break; }
    case 5: {   // 竹型: 節を持つ稈が並ぶ
        for(int s=-1;s<=1;s++){
            float sx=x+s*h*0.15f, sh=h*(s==0?1.f:0.78f);
            sf::RectangleShape st({std::max(1.5f,h*0.045f),sh});
            st.setPosition({sx-h*0.022f,baseY-sh});
            st.setFillColor(col); win.draw(st);
            for(int nd=1;nd<5;nd++){
                sf::RectangleShape n({h*0.10f,1.6f});
                n.setPosition({sx-h*0.05f,baseY-sh*nd/5.f});
                n.setFillColor(sf::Color(col.r,col.g,col.b,190)); win.draw(n);
            }
        }
        break; }
    case 6: {   // 樹木: 幹と樹冠
        sf::RectangleShape tr({std::max(2.f,h*0.10f),h*0.62f});
        tr.setPosition({x-h*0.05f,baseY-h*0.62f});
        tr.setFillColor(sf::Color(108,82,54)); win.draw(tr);
        sf::CircleShape cr(h*0.30f); cr.setOrigin({h*0.30f,h*0.30f});
        cr.setPosition({x,baseY-h*0.70f}); cr.setFillColor(col); win.draw(cr);
        sf::CircleShape c2(h*0.20f); c2.setOrigin({h*0.20f,h*0.20f});
        c2.setPosition({x-h*0.20f,baseY-h*0.55f});
        c2.setFillColor(sf::Color(col.r,col.g,col.b,200)); win.draw(c2);
        break; }
    case 9: {   // 多肉(現在は分類に使われていない)
        sf::CircleShape b(h*0.34f); b.setOrigin({h*0.34f,h*0.34f});
        b.setScale({0.85f,1.25f});
        b.setPosition({x,baseY-h*0.42f}); b.setFillColor(col); win.draw(b);
        sf::VertexArray sp(sf::PrimitiveType::Lines);
        for(int k=0;k<7;k++){
            float a=-1.5708f+(k-3)*0.34f;
            sf::Vertex v0,v1;
            v0.position={x+std::cos(a)*h*0.26f,baseY-h*0.42f+std::sin(a)*h*0.40f};
            v1.position={x+std::cos(a)*h*0.40f,baseY-h*0.42f+std::sin(a)*h*0.56f};
            v0.color=v1.color=sf::Color(226,226,206,170);
            sp.append(v0); sp.append(v1);
        }
        win.draw(sp);
        break; }
    default: {  // 草本: 細い葉が広がる
        for(int s=-2;s<=2;s++){
            sf::VertexArray bl(sf::PrimitiveType::LineStrip);
            for(int t=0;t<=6;t++){
                float u=t/6.f;
                sf::Vertex v;
                v.position={x+s*h*0.11f*u*1.5f, baseY-u*h*(1.f-std::abs(s)*0.13f)};
                v.color=col; bl.append(v);
            }
            win.draw(bl);
        }
        break; }
    }
}

// n角形のレーダー図。値はすべて0..1に正規化して渡す
void draw_radar(sf::RenderWindow& win,const sf::Font& font,
                float cx,float cy,float R,
                const char* const* names,const float* vals,
                const float* raw,int n,sf::Color col,
                const float* avg=nullptr)
{
    const float TAU=6.2831853f, TOP=-1.5707963f;
    auto vtx=[&](float a,float r){ return sf::Vector2f{cx+std::cos(a)*r,cy+std::sin(a)*r}; };
    // 世界平均を灰色で重ねる。どこが尖っているかが一目で分かる
    if(avg){
        sf::VertexArray ag(sf::PrimitiveType::LineStrip);
        for(int k=0;k<=n;k++){
            sf::Vertex v;
            v.position=vtx(TOP+TAU*(k%n)/n,R*clamp01(avg[k%n]));
            v.color=sf::Color(140,140,150,170); ag.append(v);
        }
        win.draw(ag);
    }

    // 目盛りの環
    for(int g=1;g<=4;g++){
        float rr=R*g/4.f;
        sf::VertexArray ring(sf::PrimitiveType::LineStrip);
        for(int k=0;k<=n;k++){
            sf::Vertex v; v.position=vtx(TOP+TAU*(k%n)/n,rr);
            v.color=sf::Color(86,86,94,g==4?210:105);
            ring.append(v);
        }
        win.draw(ring);
    }
    // 軸
    sf::VertexArray ax(sf::PrimitiveType::Lines);
    for(int k=0;k<n;k++){
        sf::Vertex v0,v1;
        v0.position={cx,cy};                  v0.color=sf::Color(72,72,80,150);
        v1.position=vtx(TOP+TAU*k/n,R);       v1.color=sf::Color(72,72,80,150);
        ax.append(v0); ax.append(v1);
    }
    win.draw(ax);
    // 中身の塗り
    sf::VertexArray fill(sf::PrimitiveType::Triangles);
    for(int k=0;k<n;k++){
        float a0=TOP+TAU*k/n, a1=TOP+TAU*((k+1)%n)/n;
        float v0=clamp01(vals[k]), v1=clamp01(vals[(k+1)%n]);
        sf::Vertex c,p0,p1;
        c.position={cx,cy};        c.color =sf::Color(col.r,col.g,col.b,55);
        p0.position=vtx(a0,R*v0);  p0.color=sf::Color(col.r,col.g,col.b,70);
        p1.position=vtx(a1,R*v1);  p1.color=sf::Color(col.r,col.g,col.b,70);
        fill.append(c); fill.append(p0); fill.append(p1);
    }
    win.draw(fill);
    // 輪郭
    sf::VertexArray out(sf::PrimitiveType::LineStrip);
    for(int k=0;k<=n;k++){
        sf::Vertex v;
        v.position=vtx(TOP+TAU*(k%n)/n,R*clamp01(vals[k%n]));
        v.color=col; out.append(v);
    }
    win.draw(out);
    // 頂点の印と名前
    for(int k=0;k<n;k++){
        float a=TOP+TAU*k/n;
        sf::CircleShape dot(2.6f); dot.setOrigin({2.6f,2.6f});
        dot.setPosition(vtx(a,R*clamp01(vals[k])));
        dot.setFillColor(col); win.draw(dot);
        char lb[48];
        if(raw) snprintf(lb,48,"%s %.2f",names[k],raw[k]);
        else    snprintf(lb,48,"%s",names[k]);
        // 軸が13本あると中央揃えでは隣と重なる。
        // 左半分は右端を、右半分は左端を基準にして外側へ逃がす
        sf::Text t(font,jp(lb),13);
        t.setFillColor(sf::Color(222,222,230));
        sf::FloatRect bb=t.getLocalBounds();
        sf::Vector2f lp=vtx(a,R+26.f);
        float ox=lp.x-cx;
        float tx = (ox<-R*0.3f) ? lp.x-bb.size.x
                 : (ox> R*0.3f) ? lp.x
                                : lp.x-bb.size.x*0.5f;
        t.setPosition({tx,lp.y-8.f});
        win.draw(t);
    }
}

// セーブの冒頭だけを読んで概要を作る。一覧を開くたびに全件読むのは重いので、
// カーソルが乗った一件だけを読み、結果は記憶する
std::string save_summary(const std::string& fn){
    std::ifstream f(fn);
    if(!f) return "読めない";
    std::string tag, sec;
    f >> tag;
    if(tag!="P3SAVE2") return "古い形式 — 読み込めない";
    int fr=0,a=0,b=0,c=0,d=0,e=0,gw=0,gh=0;
    f >> sec;
    if(sec!="HEADER") return "壊れている";
    f >> fr >> a >> b >> c >> d >> e >> gw >> gh;
    // 各セクションの件数だけを拾う
    long np=-1,na=-1,nl=-1;
    std::string w; long n;
    while(f >> w){
        if(w=="PLANTS"){ f>>n; np=n; }
        else if(w=="ANIMALS"){ f>>n; na=n; }
        else if(w=="LINEAGES"){ f>>n; nl=n; break; }
        else if(w=="SCHEMA"){ std::string k; f>>k>>n; for(long i=0;i<n;i++) f>>k; }
        else if(w=="CELLS"){ f>>n; break; }   // 本体は読まない
        else f.ignore(1<<20,'\n');
    }
    char b2[160];
    if(np<0) snprintf(b2,160,"t = %d   (件数は本体を読まないと不明)",fr);
    else snprintf(b2,160,"t = %d   植物 %ld   動物 %ld",fr,np,na);
    return std::string(b2);
}

