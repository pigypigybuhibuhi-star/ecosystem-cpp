// スレッドの溜め置き(parallel_for)(phase3.cpp から分けた。phase3.cpp が順に #include する)
// スレッドを毎回作ると、生成だけで数ミリ秒かかる。
// 起動時に作って寝かせておき、仕事があるときだけ起こす
class Pool {
    std::vector<std::thread> th;
    std::mutex mx;
    std::condition_variable cv, done_cv;
    std::function<void(int,int)> job;
    std::vector<std::pair<int,int>> ranges;
    int next_task=0, running=0, epoch=0;
    bool quit=false;
public:
    int n;
    Pool(){
        n=(int)std::thread::hardware_concurrency();
        if(n<2) n=1; if(n>8) n=8;
        for(int k=0;k<n;k++) th.emplace_back([this]{ worker(); });
    }
    ~Pool(){
        { std::lock_guard<std::mutex> lk(mx); quit=true; }
        cv.notify_all();
        for(auto& t:th) t.join();
    }
    void worker(){
        int my_epoch=0;
        for(;;){
            std::unique_lock<std::mutex> lk(mx);
            cv.wait(lk,[&]{ return quit || epoch!=my_epoch; });
            if(quit) return;
            my_epoch=epoch;
            for(;;){
                if(next_task>=(int)ranges.size()) break;
                auto r=ranges[next_task++];
                lk.unlock();
                job(r.first,r.second);
                lk.lock();
            }
            if(--running==0) done_cv.notify_one();
        }
    }
    template<class F> void run(int a,int b,F f){
        if(n<2||b-a<32){ f(a,b); return; }
        {
            std::lock_guard<std::mutex> lk(mx);
            job=f;
            ranges.clear();
            int span=(b-a+n-1)/n;
            for(int k=0;k<n;k++){
                int s=a+k*span, e=std::min(b,s+span);
                if(s<e) ranges.push_back({s,e});
            }
            next_task=0; running=n; epoch++;
        }
        cv.notify_all();
        std::unique_lock<std::mutex> lk(mx);
        done_cv.wait(lk,[&]{ return running==0; });
    }
};
Pool g_pool;
struct PoolReport { PoolReport(){ printf("thread pool: %d workers\n",g_pool.n); } };
PoolReport g_pool_report;
template<class F> void parallel_for(int a,int b,F f){ g_pool.run(a,b,f); }

