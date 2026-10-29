// 神経網の計算(brain_forward)(phase3.cpp から分けた。phase3.cpp が順に #include する)
// 神経網の一回の計算。世界の中の動物と、試験器(--assay)が同じ計算を使う。
// in を読んで out を書き、働いたニューロンと結線の数を返す。
// 働いているニューロンだけを通す。すべて無効なら入力から出力への直結だけが残り、単層として振る舞う
static inline void brain_forward(const Animal& a,const float* in,float* out,int& n_node,int& n_syn){
    float bus[MAX_NODE]={0}, nxt[MAX_NODE];
    n_node=0; n_syn=0;
    int   last_w=0;                     // 最後に働いた層の幅
    for(int L=0;L<MAX_LAYER;L++){
        int act=0;
        for(int h=0;h<MAX_NODE;h++) if(a.ngain[L*MAX_NODE+h]>NODE_ON) act++;
        if(!act) continue;              // この層は存在しない
        for(int h=0;h<MAX_NODE;h++){
            float g=a.ngain[L*MAX_NODE+h];
            if(g<=NODE_ON){ nxt[h]=0.f; continue; }
            float s=0.f;
            if(last_w==0){              // 第1層。入力から受ける
                const float* wr=&a.w[h*N_IN];
                for(int k=0;k<N_IN;k++) s+=wr[k]*in[k];
                n_syn+=N_IN;
            } else {                    // 前の層から受ける
                int off=W_L0+(L-1)*MAX_NODE*MAX_NODE+h*MAX_NODE;
                if(off+MAX_NODE>W_L0+W_LL) off=W_L0+h*MAX_NODE;
                const float* wr=&a.w[off];
                for(int j=0;j<MAX_NODE;j++) s+=wr[j]*bus[j];
                n_syn+=last_w;
            }
            nxt[h]=std::tanh(s*g);
            n_node++;
        }
        for(int h=0;h<MAX_NODE;h++) bus[h]=nxt[h];
        last_w=act;
    }
    for(int o=0;o<N_OUT;o++) out[o]=0.f;
    if(last_w>0){
        const float* wo=&a.w[W_L0+W_LL];
        for(int o=0;o<N_OUT;o++)
            for(int h=0;h<MAX_NODE;h++) out[o]+=wo[o*MAX_NODE+h]*bus[h];
        n_syn+=N_OUT*last_w;
    }
    {   // 入力から出力への直結
        const float* wd=&a.w[W_L0+W_LL+W_OUT];
        for(int o=0;o<N_OUT;o++)
            for(int k=0;k<N_IN;k++) out[o]+=wd[o*N_IN+k]*in[k];
    }
    for(int o=0;o<N_OUT;o++) out[o]=std::tanh(out[o]);
}

