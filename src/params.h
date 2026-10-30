// パラメータ(空間・時間・地形・気候・水・大気・火・天気・植物・動物・匂い・記録の定数)(phase3.cpp から分けた。phase3.cpp が順に #include する)
// ==========================================================
//  パラメータ
// ==========================================================
// --- 空間(作者が決定) ---
// 惑星ではなく一区画の土地。1マス = 10 m、900×500 マス = 9 km × 5 km
// 1 tick = 1 分。1 日 = 1440 tick、1 年 = 365 日 = 525600 tick(移行の段階 4)
const float CELL_M          = 10.0f;   // 1マスの一辺(メートル)
const int   GRID_W          = 900;     // 9 km
const int   GRID_H          = 500;     // 5 km

// --- 時間(移行の段階 3: 時間の下地、段階 4: 1 tick = 1 分) ---
// 速さの定数は「1 分あたりの物理量 × DT」の形で持つ。DT は 1 tick の長さ(分)。
// 段階 3 では DT を今までの 1 tick(1 日 = 14 tick、約 103 分)にとり、挙動を変えなかった。段階 4 で DT = 1 分にした。
// 今までの値は「1 回の更新あたり」の量だったので、次の換算で 1 分あたりに直して持つ。
//   PM(v,n)   n tick ごとの 1 回の更新で v だけ増減していた量 → 1 分あたりの量 v/(n·DT0)
//   PMF(f,n)  n tick ごとの 1 回の更新で割合 f だけ減っていた(起きていた) → 1 分あたりの率 −ln(1−f)/(n·DT0)
//   LIN(r,n)  1 分あたり r の量を、n tick ぶん = r·n·DT
//   FRAC(k,n) 1 分あたりの率 k で、n tick の間に減る(起きる)割合 = 1−exp(−k·n·DT)
//   TICKS(m)  m 分を tick に直す
// 減衰と確率は指数の形で換算する(掛け算で DT 倍すると、DT が大きいとき 1 を超えるため)。
// PM / PMF の中の値は、段階 3 より前の値をそのまま残したもの。段階 5 で実際の 1 分あたりの物理量に置き換える。
// 換算は翻訳時に倍精度で行う。DT = DT0 のとき、float の値は今までとビット単位で同じになる
constexpr double c_log1p(double y){            // log(1+y) = 2·atanh(y/(2+y))。|y| < 1
    double z=y/(2.0+y), z2=z*z, t=z, sum=0.0;
    for(int k=0;k<20000;k++){
        double d=t/(2*k+1); sum+=d; t*=z2;
        if((d<0?-d:d) <= 1e-18*(sum<0?-sum:sum)) break;
    }
    return 2.0*sum;
}
constexpr double c_expm1(double x){            // exp(x)−1 のテイラー級数。小さい x でも桁を落とさない
    double t=x, sum=0.0;
    for(int k=2;k<1000;k++){
        sum+=t; t*=x/k;
        if((t<0?-t:t) <= 1e-18*(sum<0?-sum:sum)) break;
    }
    return sum;
}
constexpr double DAY_MIN  = 1440.0;            // 1 日(分)
constexpr double DT0      = DAY_MIN/14.0;      // 段階 3 より前の 1 tick(分)。換算にだけ使う
constexpr double DT       = 1.0;               // 1 tick(分)
constexpr double YEAR_MIN = 365.0*DAY_MIN;     // 1 年(分)= 365 日 = 525600 分
constexpr double PM  (double v,double n=1){ return v/(n*DT0); }
constexpr double PMF (double f,double n=1){ return -c_log1p(-f)/(n*DT0); }
constexpr float  LIN (double r,double n=1){ return (float)(r*n*DT); }
constexpr float  FRAC(double k,double n=1){ return (float)(-c_expm1(-k*n*DT)); }
constexpr long   TICKS(double m){ return (long)(m/DT+0.5); }
// 回数で決まる速さ(一度に齧る回数、死骸を漁る回数、繁殖の判定)。
// 段階 3 より前に 1 回の更新で n0 回だった行為は、1 分あたり n0/(n·DT0) 回の率で起きる。
// 1 回の更新(n tick = n·DT 分)の回数は 率 × n·DT = n0·(DT/DT0)。端数は確率で決める
constexpr double STEP_RATIO = DT/DT0;          // 段階 3 より前の 1 tick に対する、今の 1 tick の長さ
constexpr double c_sqrt(double x){             // ニュートン法。x >= 0
    if(x<=0.0) return 0.0;
    double r=x>1.0?x:1.0;
    for(int k=0;k<200;k++){ double nr=0.5*(r+x/r); if(nr==r) break; r=nr; }
    return r;
}
// 向きのゆらぎはランダムウォークなので、角度の広がりは時間の平方根に比例する
constexpr float  TURN_SCALE = (float)c_sqrt(STEP_RATIO);
// 更新の間隔。計算の都合で何 tick ごとに処理するかで、物理量ではないので tick で持つ。
// 1 回の更新は間隔ぶんの時間を進める(上の n)
const int   PLANT_STRIDE    = 3;       // 植物の体の更新(光合成は毎 tick)
const int   ANIMAL_STRIDE   = 2;       // 動物も2tickに分けて処理する
const int   SEED_STRIDE     = 8;       // 土中の種子はさらに間引く
const int   MICROBE_STRIDE  = 4;       // 個体数が桁違いに多い。4tickに分けて処理する
const int   PATH_INTERVAL   = 3;       // 病原体を何tickごとに更新するか
const int   ODOR_INTERVAL   = 4;       // 匂いを何tickごとに更新するか
const int   GAS_INTERVAL    = 4;       // 大気の移流と拡散を何tickごとに更新するか
const int   PARA_INTERVAL   = 4;       // 寄生の吸い取りを何tickごとに行うか
const int   GENET_INTERVAL  = 10;      // 地下茎の融通を何tickごとに行うか

// --- 地形 ---
// 惑星ではなく一区画の土地。南の端が海に沈み、北へ丘陵が上がる海岸の土地とする。
// 区画の外の地殻変動(隆起の量と向き)は条件として与え、区画の中の形は
// 隆起・川の侵食(ストリームパワー則)・斜面の崩れ(安息角)から出す
const float SEA_LEVEL       = 0.42f;   // 標高の尺度(0〜1)での海面
const float ELEV_M          = 400.f/(1.f-SEA_LEVEL);   // 標高の尺度 1 が何 m か(約 690 m)。表し方の単位で、起伏は決めない
const float SEA_DEPTH_M     = 60.f;    // 条件: 区画の南の端の海の深さ(m)。沈降する帯の初期の面
const float COAST_FRAC      = 0.25f;   // 条件: 南から何割の帯が沈降して海の下にあるか
const float REPOSE_SLOPE    = 0.70f;   // 安息角 35° の勾配(tan 35°)。これより急な斜面は崩れて残らない
// 地形を作るときだけ使う。すべて実際の単位(m・年)で書き、起伏と傾斜は
// 隆起・川の侵食・斜面拡散の釣り合いから出す。地球の傾斜に合わせて写し直さない
const float TERRAIN_NOISE_M = 2.0f;    // 初期の面の起伏(標準偏差、m)。海岸線を曲げる
const float UPLIFT_MM_YR    = 0.3f;    // 条件: 北の端の隆起の速さ(mm/年)。北ほど大きい(傾動する地塊)。
                                       // 中程度に活発な海岸の山地の値(0.1〜1 mm/年)
const float UPLIFT_VARY     = 0.40f;   // 条件: 隆起の場所ごとのばらつき
// ストリームパワー則 dz/dt = U − K A^m S(n=1、m=0.5)。Cordonnier et al. 2016 の陰的解法
const float EROSION_M       = 0.5f;    // 流域面積の指数
const float EROSION_K       = 1.0e-5f; // 侵食係数(/年、A は m²)。m=0.5 の岩盤の川の値は 1e-6〜1e-4
                                       // (Stock & Montgomery 1999、Harel et al. 2016)の中ほど
const float EROSION_DT_YR   = 10000.f; // 時間刻み(年)。川は陰的解法なので大きく取れる
// 斜面拡散(土壌の匍行)の係数。土に覆われた斜面で 1e-3〜1e-2 m²/年。
// オレゴン海岸山脈で 0.0032 m²/年(Roering et al. 1999)
const float HILL_DIFFUSE_M2_YR = 0.0032f;
const int   REPOSE_PASS     = 20;      // 1 回の侵食のあとに崩れを繰り返す最大の回数(崩れが止まれば打ち切る)
const float FILL_EPS_M      = 0.005f;  // 窪地を埋めた面に付ける 1 マスあたりの最小の傾き(m)
// 川の始まり。流域がこれより小さい斜面は川に削られず、拡散だけで形が決まる
// (Montgomery & Dietrich 1988。谷頭の集水面積はおよそ 0.01〜0.1 km²)
const float CHANNEL_AREA    = 100.0f;    // 何マス(1 マス 100 m²)から川が削るか
const int   EROSION_STEPS   = 300;     // 条件: 地塊が隆起を始めてからの年数(300 × 1 万年 = 300 万年)。
                                       // この長さではまだ釣り合いに届かない(起伏は隆起の年数で決まる)

// --- 気候(惑星規模の条件) ---
// 緯度・季節ごとの日射・卓越風・背景の気温は外から与える。区画は 9 km × 5 km で
// 緯度の差は 0.05° に満たないので、緯度と太陽時は区画の中で一様とする
const float LATITUDE_DEG    = 40.0f;   // 条件: 区画の緯度(北緯)
const float LAT_N           = LATITUDE_DEG/90.f;
// 背景の気温は、緯度の関数として与える(条件)。区画の緯度で一度だけ評価する
const float TEMP_EQUATOR    = 30.0f;
const float TEMP_POLE       = -15.0f;
const float LAPSE_RATE      = 6.5e-3f*ELEV_M;  // 気温減率 6.5 °C/km を標高の尺度あたりに直したもの
const int   SEASON_PERIOD   = (int)TICKS(YEAR_MIN);   // 1 年(tick)
const float SEASON_AMP      = 24.0f;   // 年較差。緯度を掛けて使う
const float SEA_SEASON      = 0.30f;   // 海の年較差は陸の3割ほど
// 海面水温は気温と別の曲線を描く。凍結点より下がれないため、極でも−数度で止まる
const float SST_EQUATOR     = 29.0f;
const float SST_POLE        = -3.5f;
const float CURRENT_HEAT    = 9.0f;    // 水温偏差の表示の幅(度)
const float ICE_LAG         = FRAC(PMF(0.055f));  // 氷の成長と融解は季節に遅れる [率 0.00055 /分]
const int   COAST_RANGE     = 7;       // 陸の冷気が海へ及ぶ距離
const float COAST_CHILL     = 0.55f;   // 沿岸の海が陸の気温に引かれる強さ

// --- 水循環 ---
const float EVAP_COEF       = 0.050f;  // 海面温度あたりの蒸発量
const float MOIST_CAP       = 8.00f;   // 大気が運べる水蒸気の上限
const float LAND_RECYCLE    = 0.040f;  // 陸からの蒸発散(水分の再循環)
const float AIR_IN_RH       = 0.75f;   // 風上の端から入る空気の相対湿度(条件。地表の年平均は海 0.8・陸 0.7 ほど、Dai 2006)
const float RAIN_BASE       = 0.024f;  // 対流による降水率
const float OROGRAPHIC      = 10.0f;   // 地形性降雨(斜面を上る空気は雨を落とす)
const int   ORO_SPAN        = 6;       // 斜面を測る幅。広いほど雨陰がなだらかになる
const int   PRECIP_BLUR     = 6;       // 降水場のぼかし半径
const int   PRECIP_INTERVAL = 300;     // 降水場を計算し直す間隔
// ケッペンの乾燥限界(20T+140)にならい、限界は気温に比例させる
const float ARID_COEF       = 0.0010f; // 乾燥限界の係数
const float ARID_SOFT       = 1.20f;   // 大きいほど乾湿の差がなだらかになる
const float FIELD_CAP       = 1.00f;   // 土壌が保持できる上限
const float WATER_RELAX     = FRAC(PMF(0.020f));  // 土壌水分が目標へ近づく速さ [率 0.000196 /分]
const float RIVER_THRESHOLD = 260.f;   // 侵食後の地形に合わせて下げる
const float RIVER_SUPPLY    = 0.60f;   // 河川が供給する水
const float SOIL_W_INIT     = 0.35f;   // 起動時に降水を決めるための仮の陸の土壌水分(条件)
const float DROUGHT_COST    = LIN(PM(0.055f));  // 耐乾装置の維持費 [1分あたり 0.000535]
const float STORE_COST      = LIN(PM(0.030f));  // 貯水組織の維持費 [1分あたり 0.000292]
const float STORE_FILL      = LIN(PM(0.020f,PLANT_STRIDE),PLANT_STRIDE);  // 湿っているとき貯める速さ [1分あたり 6.48e-05]
const float STORE_MAX       = 3.00f;   // 遺伝子が届く貯水量の上限

// --- 大気 ---
const float CO2_PER_ENERGY  = 0.020f;  // エネルギー 1 を酸化して出る CO2 の量(気体の単位)
const float O2_PER_METAB    = 0.020f;
// 大気の柱 [条件]: 各マスの上の空気の全量を、地上の気圧 101325 Pa ÷ 重力 9.81 m/s² ÷ 空気のモル質量 0.028965 kg/mol
// = 3.57e5 mol/m² とし、体積の割合 CO2 420 ppm・O2 20.95%・N2 と Ar 79.01%(地球の今の値。惑星規模の条件)で分ける。
// 気体の単位はエネルギー 1(= 6 mg C、D-2026-10-24-2)× CO2_PER_ENERGY なので、1 単位 = 0.3 g C = 0.025 mol。
// co2 と o2 の場は、この柱の量からのずれだけを持つ(柱は約 6e5 と 3e8 単位で、float では生き物の小さな出入りが消えるため)
const double AIR_MOL_M2     = 101325.0/9.81/0.028965;            // [3.57e5 mol/m²]
const double GAS_UNIT_MOL   = 6e-3/12.011/CO2_PER_ENERGY;        // 気体の 1 単位のモル数 [0.025 mol]
const double AIR_COL_UNITS  = AIR_MOL_M2*CELL_M*CELL_M/GAS_UNIT_MOL;   // 1 マスの柱の空気の全量 [1.43e9 単位]
const float CO2_COL         = (float)(AIR_COL_UNITS*420e-6);    // [6.0e5 単位]
const float O2_COL          = (float)(AIR_COL_UNITS*0.2095);    // [3.0e8 単位]
const float N2_COL          = (float)(AIR_COL_UNITS*0.7901);    // N2 と Ar [1.13e9 単位]。生き物は出し入れしない
const float SOIL_N_BASE     = 0.35f;   // 土壌の窒素
const float N_PER_GROWTH    = 0.055f;  // 体を作るのに要する窒素
const float N_FIX_COST      = LIN(PM(0.006f));  // 純益0.048に対して現実的な負担に [1分あたり 5.83e-05]
const float N_FIX_RATE      = LIN(PM(0.0180f));  // 共生菌は効率がよい [1分あたり 0.000175]
const float N_STORE_MAX     = 1.5f;    // 植物が体内に貯められる窒素の上限
const float N_SINK          = FRAC(PMF(0.022f));  // 沈降。外洋を貧栄養にするが、枯渇させない [率 0.000216 /分]
const float N_RUNOFF        = LIN(PM(0.045f));  // 河川が沿岸へ運ぶ栄養 [1分あたり 0.000437]
const float N_UPWELL        = LIN(PM(0.040f));  // 湧昇。深層水が栄養を持ち上げる [1分あたり 0.000389]
const int   SHELF_RANGE     = 6;

// --- 海氷 ---
// 海氷の形成には、継続した冷却が要る。瞬間の水温より低い閾値で表す
const float FREEZE_TEMP     = -5.0f;
const float THAW_TEMP       = -1.5f;
const float ICE_LIGHT       = 0.08f;   // 氷が通す光

const float PLANT_RESP_CO2  = 1.00f;   // 植物の呼吸(維持コストのうちCO2を出す割合)。維持の仕事は糖の酸化でまかなうので全量
const float RESP_EFFICIENCY = 0.60f;   // 呼吸効率の遺伝子が持つ効き
const float RESP_COST       = LIN(PM(0.032f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 高効率な呼吸器の維持費 [1分あたり 0.000156]
const float DISSOLVED_O2    = 0.45f;   // 鰓は水から効率よく酸素を取り出す
const float ALT_O2_DECAY    = ELEV_M/8400.f;  // 標高による空気の薄まり。大気のスケールハイト 8.4 km
const float GAS_DIFFUSE     = FRAC(PMF(0.55f,GAS_INTERVAL),GAS_INTERVAL);  // 4tickに一度呼ぶぶん強める [率 0.00194 /分]
const float GAS_RELAX       = FRAC(PMF(0.0060f,GAS_INTERVAL),GAS_INTERVAL);  // [率 1.46e-05 /分]
const float WIND_STRENGTH   = LIN(PM(0.35f,GAS_INTERVAL),GAS_INTERVAL);  // 拡散に負けない強さにする [1分あたり 0.000851]
// 条件: 地上 10 m の風の速さ。北半球の中緯度の陸の観測点の平均はおよそ 3 m/s(Vautard et al. 2010)。
// 向きは卓越風の帯(wind[] の符号)に従う。大気の移流と水面の漂流はこの速さから決める
const float WIND_SURF_MS    = 3.0f;
const float GAS_SHIFT       = LIN(WIND_SURF_MS*60.0/CELL_M,GAS_INTERVAL);  // 1 回の更新で空気が進むマス数 [1分あたり 18 マス]
const float O2_PER_CO2      = 1.00f;

// --- 火 ---
// 酸素が17%未満では燃え広がらず、35%を超えると消火が困難になる(火の窓)
// 地球の21%を1.0とした比。17%未満では広がらず、35%超で消火困難
const float FIRE_O2_MIN     = 0.81f;   // 17% / 21%
const float FIRE_SPREAD     = FRAC(PMF(0.20f));  // 隣へ移る強さ [率 0.00217 /分]
const float FIRE_BURN       = FRAC(PMF(0.45f));  // 一tickで燃える燃料の割合 [率 0.00581 /分]
const float FIRE_DECAY      = FRAC(PMF(0.38f));  // 火勢の減衰。燃え尽きたら消える [率 0.00465 /分]
// 燃料の湿り(乾燥質量あたりの水の比)。土の水分ではなく、大気の相対湿度と気温から決まる平衡含水率へ時定数で追従する
const float FIRE_MX         = 0.25f;   // 消炎含水率。Anderson (1982) の 13 の標準燃料モデルの中央値(12〜40%)
const float FUEL_M_LAG      = FRAC(1.0/60.0);  // 1 時間燃料(細い枯れ草・小枝)の時定数 60 分(Fosberg & Deeming 1971)
const float LAND_RH_MEAN    = 0.70f;   // 条件: 陸の地表の日平均の相対湿度(Dai 2006)。水蒸気圧は一日の中で変わらないとし、日変化は気温から出す
const float FIRE_O2_USE     = CO2_PER_ENERGY;  // 燃焼が消費する酸素。CH2O + O2 → CO2 + H2O で CO2 と 1:1
const float FIRE_TOL_COST   = LIN(PM(0.075f));  // 厚い樹皮は高くつく [1分あたり 0.000729]
const float FIRE_DAMAGE     = LIN(PM(9.0f,PLANT_STRIDE),PLANT_STRIDE);  // 炎そのものの破壊力 [1分あたり 0.0292]
const float FIRE_FUEL_MIN   = 1.5f;

// --- 撹乱 ---
const float SLIDE_SLOPE     = 0.58f;   // これ以上の勾配(30°)で崩れうる。表層崩壊の多くは 30° を超える斜面で起きる
// 崩れる厚さ。表層崩壊は土層と基盤の境で滑り、厚さは 0.5〜2 m(小さな崩壊の中央値はおよそ 1 m、Larsen et al. 2010)
const float SLIDE_DEPTH_M   = 1.0f;
// 崩れの起こり方: 無限斜面の安全率 FS = c'/(γ z sinα cosα) + (1 − m γw/γ) tanφ'/tanα が 1 を下回ったマスが崩れる
// (Skempton & DeLory 1957、Montgomery & Dietrich 1994)。m は崩れる厚さ z のうち地下水で飽和した割合で、
// 地下水面が地表から z より浅く上がると間隙水圧が摩擦を減らす。z は SLIDE_DEPTH_M、内部摩擦角は安息角 REPOSE_SLOPE(35°)
// (粘着の無い粒の安息角は内部摩擦角に等しい)。粘着力は土そのもの 2 kPa と根の補強 3 kPa の和 5 kPa
// (崩積土 0〜5 kPa、草地と低木の根 1〜5 kPa。Montgomery & Dietrich 1994、Schmidt et al. 2001)。単位体積重量 17 kN/m³。
// この値で区画の最も急な斜面(勾配の大きさ約 1)でも乾いた状態の FS は 1 を超え、飽和すると勾配約 0.76(37°)より急な斜面が崩れる。
// 以前は勾配に比例する一定の率(5.71e-10 /分、隆起と釣り合う削られる速さから)で、雨や地下水と無関係だった(D-2026-10-26-2 の答え (a))
const float SLOPE_C_KPA     = 5.0f;
const float SLOPE_GAMMA     = 17.0f;   // kN/m³
const float WATER_GAMMA     = 9.81f;   // kN/m³
// 崩土の到達。見通し角 H/L は体積で決まる: log10(H/L) = −0.085 log10(V) + 0.047(Corominas 1996、全種類の当てはめ)。
// 1 マスの崩壊 V = 厚さ × 10 m × 10 m で H/L ≈ 0.75
const float SLIDE_HL        = std::pow(10.f,-0.085f*std::log10(SLIDE_DEPTH_M*CELL_M*CELL_M)+0.047f);
const int   SLIDE_RUN_MAX   = 200;     // 到達を追う歩数の安全の上限(計算のためだけ)

// --- 地震の揺れ ---
// 条件: 区画の外の断層。北の端(y=0)から 10 km 北を東西に走る長さ 60 km の断層で、区画の中央の北にある。
// 北の端を頂にする傾動地塊(UPLIFT_MM_YR)を持ち上げる断層とみなし、すべりの速さは隆起の速さと同じとする。
// 地震発生層の幅 15 km、剛性率 3e10 Pa(上部地殻の標準の値)
const double QK_FAULT_DIST_KM = 10.0;
const double QK_FAULT_LEN_KM  = 60.0;
const double QK_SEIS_W_KM     = 15.0;
const double QK_RIGIDITY      = 3.0e10;
// 規模の分布: グーテンベルク・リヒター則 log10 N(≥M) = a − bM、b = 1(Gutenberg & Richter 1944)。
// 記録する最小の規模は M4(これより小さい揺れは区画の斜面の限界加速度に届かない)。
// 最大の規模は断層の長さから: M = 5.08 + 1.16 log10(地表の断層の長さ km)(Wells & Coppersmith 1994)で M7.14
const double QK_B             = 1.0;
const double QK_MMIN          = 4.0;
const double QK_MMAX          = 5.08+1.16*std::log10(QK_FAULT_LEN_KM);
// a はモーメントの釣り合いから決める: 断層が蓄えるモーメントの速さ μ L W s(Brune 1968)を、
// 打ち切った G-R 則の地震が解放するモーメント b/(1.5−b) × 10^(a + 9.1 + (1.5−b) Mmax)(M0 = 10^(1.5M + 9.1) N·m、
// Hanks & Kanamori 1979。Youngs & Coppersmith 1985)と等しく置く。a ≈ 2.94、M4 以上の本震は年 0.086 回
const double QK_M0_RATE       = QK_RIGIDITY*QK_FAULT_LEN_KM*1e3*QK_SEIS_W_KM*1e3*UPLIFT_MM_YR*1e-3;   // N·m/年
const double QK_A             = std::log10(QK_M0_RATE*(1.5-QK_B)/QK_B)-9.1-(1.5-QK_B)*QK_MMAX;
const double QK_RATE          = std::pow(10.0,QK_A-QK_B*QK_MMIN)/YEAR_MIN;   // M4 以上の本震 [1分あたり 1.6e-07]
// 余震: 本震の規模 Mm の後、時刻 t(日)に規模 M 以上が起きる率 10^(a' + b'(Mm − M)) (t + c)^−p(改良大森則、
// Reasenberg & Jones 1989 の一般的な値 a' = −1.67、b' = 0.91、p = 1.08、c = 0.05 日)。余震はさらに余震を生まない
// (この式は二次の余震を含めた列全体への当てはめなので)
const double QK_RJ_A          = -1.67;
const double QK_RJ_B          = 0.91;
const double QK_RJ_P          = 1.08;
const double QK_RJ_C          = 0.05;
// 断層の破壊の長さ(地下): log10 L(km) = −2.44 + 0.59 M(Wells & Coppersmith 1994、全種類)
// 揺れ: 最大加速度(g)log10 A = −1.02 + 0.249 M − log10 r − 0.00255 r、r = √(d² + 7.3²) km、
// d は破壊の地表への投影までの距離(Joyner & Boore 1981)。ばらつきは入れず中央の値を使う
// 斜面の限界加速度: ac = (FS − 1) g sinα(Newmark 1965)。FS は崩れの節の無限斜面の安全率(地下水の間隙水圧を含む)を使う
// 崩れ: ニューマークの変位 log10 Dn(cm) = 0.215 + log10[(1 − ac/A)^2.341 (ac/A)^−1.438](Jibson 2007)、
// 崩れる確率 P = 0.335 [1 − exp(−0.048 Dn^1.565)](Jibson et al. 2000、1994 年ノースリッジ地震の当てはめ)

// --- 天気: 対流の雨雲 ---
// 雨雲は日射で暖まった地表の空気が持ち上がって生じ、中層の風に流され、一生を終えて消える。
// 区画は雨雲一つほどの広さしかないので、区画の外(風上)で生じて流れ込む雲も同じ規則で起こす
// 条件: 雨雲を流す中層の風(指向流)の速さ。中緯度の 500 hPa の風はおよそ 10 m/s。向きは卓越風の帯に従う
const float STORM_STEER_MS  = 10.0f;
const float STORM_SHIFT     = LIN(STORM_STEER_MS*60.0/CELL_M);  // 1 分に雨雲が進むマス数 [1分あたり 60 マス]
// 単一の対流セルの一生(積雲期・成熟期・消散期)はおよそ 1 時間(Byers & Braham 1949)
const float STORM_LIFE      = 60.0f;   // 分
// 雨の降る範囲の半径。対流セルの径は 2〜10 km
const float STORM_R_M       = 2500.f;
// 冷気外出流(冷気プール)の半径。降水域より広がり、雨の届く前に気温が下がる(突風前線)
const float STORM_COLD_R_M  = 5000.f;
// 冷気プールによる気温の低下。中緯度の対流の通過で中央値 5〜6 °C(Engerer et al. 2008)
const float STORM_COLD_DT   = 6.0f;
// 突風前線: 冷気プールは周りより重い空気の塊で、密度流として外へ広がる。広がる速さは c = k·√(g'·h)、
// g' = g·ΔT/T(冷えた分の浮力の不足)。k = 0.7(Benjamin 1968 の深い流れの値。観測のフルード数 0.7〜1.0)、
// 冷気の厚さ h は条件で 1000 m(観測で 0.5〜2 km。Charba 1974)。冷気 6 °C で約 10 m/s(観測の突風前線は 5〜15 m/s)。
// 風は冷気の気圧の高まり(静水圧)の勾配に沿って外へ吹く。冷気をガウス形としたので、強さは勾配の形
// (r/σ)·exp((1 − r²/σ²)/2)(r = σ で 1)に従う。3σ より外は無視する
const float GUST_FROUDE     = 0.7f;
const float GUST_DEPTH_M    = 1000.f;
// 成熟期の中心の降水強度。対流性の雨は 10〜100 mm/h。雲ごとに指数分布で、平均は 30 mm/h を
// 20 °C との飽和水蒸気圧の比で増減させる(クラウジウス・クラペイロン、Lenderink & van Meijgaard 2008)
const float STORM_I_MEAN    = 30.0f;   // mm/h
const float STORM_I_MAX     = 150.0f;  // 観測される最大級の降水強度(計算の安全のため)
// 土の水 1.0 が何 mm にあたるか。根の届く層の有効水分量(50〜200 mm)の中ほど
const float SOIL_W_MM       = 100.0f;
// 地下水(D-2026-10-26-2 の答え (a))。岩盤の上の風化した土と砂礫の層を、一定の厚さの帯水層として置く(条件)。
// 斜面の風化層の厚さは 1〜5 m(Heimsath et al. 1997)。比産出率(重力で抜ける間隙の割合)は砂質の土で 0.1〜0.3(Johnson 1967)
const float GW_DEPTH_M      = 3.0f;
const float GW_SY           = 0.2f;
const float GW_W            = GW_SY*1000.f/SOIL_W_MM;   // 飽和した厚さ 1 m が土の水の単位でいくつか
// 地表の浸透の上限(飽和透水係数)。ローム 13 mm/h(Rawls et al. 1982)。これを超える強さの雨は地表を流れ去る(ホートン型の地表流)
const float GW_KS_MM        = 13.0f/60.0f;            // 1 分あたり mm
// 帯水層の水平の透水係数。シルト質の砂 1e-5 m/s(Freeze & Cherry 1979、1e-7〜1e-3 m/s)
const float GW_K_M          = 1.0e-5f*60.0f;          // 1 分あたり m(ダルシーの式の K)
// 地下の流れを何 tick ごとに解くか。陽解法の安定の限界は Δx²·Sy/(4·K·h) ≈ 1100 分なので、60 分は十分小さい
const int   GW_INTERVAL     = 60;
// 塩(D-2026-10-28-1 の作者の答え (b))。海の水の塩分は外洋の基準の値 35 g/kg(Millero et al. 2008、条件)。
// 塩は水に溶けて運ばれるだけで、生まれも消えもしない。土の水が蒸発・蒸散で減っても塩は残り、濃くなる。
// 濃さが NaCl の飽和(水 1 kg に 360 g、溶液で 264 g/kg)を超えた分は塩の結晶として土に残る(量は在庫のまま。濃さの上限)
const float SEA_SALT        = 35.0f;    // g/kg
const float HALITE_SAT      = 264.0f;   // g/kg
// 条件: 対流の発生の頻度(1 km² あたり 1 分あたり、日射の正弦 1・飽和水蒸気圧の比 1 のとき)。
// 中緯度の内陸で雷の日は年 30〜50 日(Changnon 1988)。一地点が年におよそ 30 回雨雲の下に入る値にした
const double STORM_INIT     = 1.3e-6;
// 雷: 対地放電は雨雲の雨と同じ対流から生じる。1 回の対地放電あたりの降水量(rain yield)は
// 大陸の対流でおよそ 1e8 kg(Petersen & Rutledge 1998。大陸 1e7〜1e8、海洋 1e9〜1e10)。
// 雲がその tick に降らせた雨の質量 ÷ この値が、その tick の対地放電の回数の期待値になる
const double LIGHT_RAIN_YIELD_KG = 1.0e8;
// 長い連続電流(40 ms 以上)を伴う対地放電の割合。燃料に着火しうるのはこの放電(Fuquay et al. 1972、約 2 割)。
// 当たった燃料が燃え出す確率は、この割合に湿りによる燃焼の弱まり(Rothermel の η_M)と酸素の比を掛ける
const float LIGHT_LCC       = 0.20f;

// --- 変異 ---
const float MUT_STRENGTH    = 0.02f;
const float MUT_RATE_MAX    = 3.0f;    // 変異率の上限(基準の何倍まで)
const float MUT_RATE_COST   = LIN(PM(0.004f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 修復機構を緩めることの代償 [1分あたり 1.94e-05]

// --- 植物 ---
const float GROW_RATE       = LIN(PM(0.02f));  // [1分あたり 0.000194]
const float GROW_COST       = 3.00f;
const float GROW_RESP       = 0.50f;   // 建設の呼吸。体にする炭素とは別に、その半分を糖の酸化で払う(建設に要るブドウ糖は体の約 1.5 倍、Penning de Vries 1974)
// 体(茎と根)の維持の呼吸。高さ 1 の体は建設に GROW_COST の炭素を持つので、体の炭素に比例する。
// 植物の組織の維持呼吸の実測は 20℃ で体の炭素あたり約 0.5〜3% /日(McCree 1974、Ryan 1991、Amthor 2000)。
// 葉の暗呼吸は RESP_COEF の側が持つので、ここは茎と根の側の低め 1% /日 = 6.94e-6 /分をとる。
// 以前は高さ 1 あたり 5.83e-4 /分(体の炭素あたり約 28% /日)で、実測の数十倍だった
const float HEIGHT_COST     = LIN(GROW_COST*0.01/DAY_MIN);  // 1 分あたり、高さあたり [2.08e-05 /分]
const float RESP_COEF       = LIN(PM(0.120f));  // 真の昼夜で取り込みが半減したぶんを補正 [1分あたり 0.00117]
const float TOUGH_COST      = LIN(PM(0.070f));  // 防御を「今」張るための費用 [1分あたり 0.000681]
const float INDUCE_RATE     = FRAC(PMF(0.22f,PLANT_STRIDE),PLANT_STRIDE);  // 食害の匂いに応じて立ち上がる速さ [率 0.000805 /分]
const float INDUCE_DECAY    = FRAC(PMF(0.020f,PLANT_STRIDE),PLANT_STRIDE);  // 攻撃がなければ緩む [率 6.55e-05 /分]
const float CLONAL_COST     = LIN(PM(0.010f));  // 地下茎の維持費 [1分あたり 9.72e-05]
const float CLIMB_COST      = LIN(PM(0.022f));  // つるの維持費(巻きひげと気根) [1分あたり 0.000214]
const float CLIMB_GAIN      = 0.85f;   // 支えの高さをどれだけ利用できるか
const float CLIMB_PENALTY   = 0.45f;
const float COLD_COST_P     = LIN(PM(0.060f));  // [1分あたり 0.000583]
const float DORM_UPKEEP     = 0.12f;   // 休眠中の維持費の割合
const float DORM_HARDEN     = 30.0f;   // 休眠がもたらす耐凍性
const float DORM_HYST       = 3.0f;    // 覚醒に要する余裕(ばたつき防止)
const float GENET_SHARE_K   = (float)PMF(0.20f);  // 地下茎で同じ株の平均へ近づく 1 分あたりの率の基準(clonal を掛ける)[率 0.00217 /分]
const float SEED_FAR        = 12.0f;
const float PLANT_DIVIDE_TH = 8.0f;
const float PLANT_INIT_E    = 6.0f;    // 立ち上がりの緩衝
const int   PLANT_INIT_N    = 400;
const int   PLANT_INIT_LAND = 200;     // 始まりの条件: そのうち陸に置く株(D-2026-10-24-1 の作者の答え (b))
const float PLANT_LAND_DORM = 5.0f;    // 始まりの条件: 陸の株の休眠の気温 [℃](温帯の草本の生育の基準温度。McMaster & Wilhelm 1997)
const int   PLANT_MAX       = 260000;  // これを超えると繁殖を止める。観察のための上限
const int   ANIMAL_MAX      = 60000;
const float SEX_COST        = LIN(PM(0.015f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 生殖器官の維持費 [1分あたり 7.29e-05]
const float SEX_SEARCH      = 3.0f;    // 相手を探す範囲(マス)

// --- 病原体 ---
// 病原体はその土地で最も多い免疫型に適応する。
// 稀な型が有利になり、組み換えで新しい型を作れる有性生殖に価値が生まれる
const int   N_IMMUNE        = 3;       // 免疫の座位
const float PATH_ADAPT      = FRAC(PMF(0.030f,PATH_INTERVAL),PATH_INTERVAL);  // 病原体が宿主の型へ寄る速さ [率 9.87e-05 /分]
const float PATH_GROW       = LIN(PM(0.020f,PATH_INTERVAL),PATH_INTERVAL);  // 宿主密度に応じた増殖 [1分あたり 6.48e-05]
const float PATH_DECAY      = FRAC(PMF(0.015f,PATH_INTERVAL),PATH_INTERVAL);  // 宿主がいなければ衰える [率 4.9e-05 /分]
const float PATH_LOAD_MAX   = 2.50f;
const float PATH_VIRULENCE  = LIN(PM(0.055f));  // 型が一致したときの被害 [1分あたり 0.000535]
const float PATH_WIDTH      = 0.42f;   // 一致とみなす幅(狭いほど型が効く)

// --- 樹上生活 ---
const float ARBOREAL_COST   = LIN(PM(0.028f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 登攀に適した体の維持費 [1分あたり 0.000136]
const float CLIMB_UP_COST   = 0.55f;   // 高さ1を登るのに要するエネルギー
const float PERCH_GAP       = 0.45f;   // 枝から枝へ渡れる高さの差
const float BRANCH_SUPPORT  = 1.20f;   // 枝が支えられる体重の目安
const float FALL_DAMAGE     = 0.9f;    // 落下で失うエネルギー(高さあたり)
const float TEMP_WIDTH_P    = 18.0f;

// --- 動物 ---
const float A_MAX_SPEED     = LIN(PM(0.20f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 0.000972]
// 基礎代謝(安静時の呼吸)。組織の炭素 1 あたり、1 分あたりの量。
// 地面を歩く小さな昆虫の安静時の O2 消費は 20〜25℃ で約 0.3 ml O2 /g(生重)/時(Chown & Nicolson 2004、Addo-Bediako et al. 2002)。
// 呼吸商 1 で 0.3 ml CO2 = 1.34e-5 mol C = 1.6e-4 g C /g/時。生重 1 g の炭素は約 0.15 g(乾燥 30%、その半分が炭素)なので、
// 組織の炭素あたり 1.07e-3 /時 = 1.8e-5 /分(約 2.6% /日)。代謝は体の量に比例して増える
const float A_SMR           = LIN(1.8e-5,ANIMAL_STRIDE);
// 形質の維持費(顎・射程・耐寒・感覚器・脳など)は、維持する組織の代謝の上乗せ分。組織の量に比例する。
// 以前は 1 個体あたりの定数で、旧い基礎代謝 A_METAB(1 個体あたり 1.94e-4 /分)と並べて置かれていた。
// 形質どうしの費用の比はそのまま残し、旧い基礎代謝に対する比を、今の基礎代謝(組織の炭素あたり)への上乗せの比として使う
const float A_TRAIT_REF     = LIN(PM(0.040f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 形質の費用を測る基準(旧い A_METAB) [1分あたり 0.000194]
const float A_MOVE_COST     = 0.060f;
const float A_REACH_COST    = LIN(PM(0.018f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 8.75e-05]
const float A_SIZE_COST     = LIN(PM(0.020f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 9.72e-05]
const float JAW_COST        = LIN(PM(0.012f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 5.83e-05]
const float COLD_COST_A     = LIN(PM(0.050f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 0.000243]
const float A_BITE          = 0.45f;   // 一口で植物から取る割合。時間は齧る回数の側が持つ(段階 3.5)
const float GRAZE_DEFOLIATE = 0.18f;   // 採食が葉を奪う量。局所的な枯渇を生む
const float A_DIVIDE_TH     = 6.0f;
const float A_EAT_LOSS      = 0.50f;   // 同化効率50%。若い葉や種子は消化しやすい
const float PRED_RATIO      = 1.15f;
const float PRED_EFF        = FRAC(PMF(0.85f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [率 0.00922 /分]
const float SCAV_EFF        = 0.85f;   // 一回漁って死骸から取る割合。死骸は抵抗しないので取り分が大きい。時間は漁る回数の側が持つ(段階 3.5)
const float MEAT_LOSS       = 0.30f;
const float DETRITUS_LOSS   = 0.80f;   // 分解物は栄養価が中間
const float DETRITUS_COST   = LIN(PM(0.008f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 分解酵素の維持費 [1分あたり 3.89e-05]
// 分解者。動かず、死骸だけを食い、炭素と窒素を土に返す
// 微生物の維持の呼吸。体の炭素に比例する(維持係数。Pirt 1965)。土壌の微生物の基礎呼吸の実測(qCO2)は
// 体の炭素あたり約 0.6〜2e-3 /時(Anderson & Domsch 1985、1993)。その低い側 1e-3 /時 = 1.7e-5 /分を、酵素の係数 (0.5+rate) が 1 のときの値とする。
// 以前は 1 個体あたりの定数 3.89e-06 /分(体の炭素 1.0 で実測の約 1/4)で、体の量によらなかった
const float MICROBE_COST    = LIN(1.7e-5);  // 1 分あたり、体の炭素あたり [1.7e-05 /分]
// 炭素利用効率(分解した炭素のうち体になる割合。残りは呼吸で CO2)。0.6 は熱力学の上限に近い。
// 土や落葉の微生物の実測は 0.3〜0.55(Manzoni et al. 2012、Sinsabaugh et al. 2013)。その中ほど
const float MICROBE_EFF     = 0.40f;
const float MICROBE_DIV     = 1.40f;   // 分裂の閾値を上げ、個体数を抑える
const float MICROBE_REF     = 0.75f*MICROBE_DIV;   // 分裂の周期の平均の体の炭素。低温の維持費の換算にだけ使う
const float MICROBE_RETURN  = 0.55f;   // 分解して大気と土に返す割合
// 胞子が飛ぶ距離の平均(マス)。地表から放たれた胞子は、落ちる速さで沈みながら風に運ばれる(Gregory 1973)。
// 落葉の上の子実体が胞子を放つ高さ 5 cm、そこでの風 0.8 m/s(対数分布、u* 0.2 m/s・粗度 1 cm)、
// 直径 5〜10 µm の胞子の終端速度 1 mm/s(ストークス)から、平均の距離 = 高さ × 風 / 終端速度 = 40 m。
// 乱流で遠くへ行くものもあるので、平均をこの値とする指数分布にとる。
// 以前は 0〜2 マス(0〜20 m)の一様で、微生物の広がりが分裂ごとに 20 m に縛られていた
const float SPORE_H         = 0.05f;   // 放つ高さ [m]
const float SPORE_WIND      = 0.80f;   // 放つ高さの風 [m/s]
const float SPORE_VSET      = 1.0e-3f; // 終端速度 [m/s]
const float MICROBE_DISP    = SPORE_H*SPORE_WIND/SPORE_VSET/CELL_M;   // [4 マス = 40 m]
const int   MICROBE_INIT    = 3000;
const float BODY_GROW_RATE  = LIN(PM(0.010f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 4.86e-05]
const float BODY_GROW_COST  = 4.00f;
const float CORPSE_RETURN   = 2.50f;   // 体 1 あたりの組織が持つエネルギー(炭素)。死骸に残るのもこの量
// 動物がエネルギーを使えば、その炭素は酸化されて CO2 になる(O2 は 1:1)。
// 体にある分(正のエネルギー)までしか出せない
inline void co2_add(int c,float d);
inline void o2_add (int c,float d);
inline void a_burn(float* co2g,float* o2g,int cell,float e){
    if(e<=0.f) return;
    float c=e*CO2_PER_ENERGY; o2_add(cell,-std::min(c,O2_COL+o2g[cell])); co2_add(cell,c);
}
// 負のエネルギーは、体に無い燃料で払った借り(その炭素はまだ CO2 に出ていない)。
// 後から体に入った燃料(食べ物・脂肪)は、借りの分だけ酸化されて CO2 になる。返す量を返す
inline float a_repay(float energy,float gain){ return std::min(gain,std::max(0.f,-energy)); }
const float MATURE_FRAC     = 0.80f;
const int   MAX_BROOD       = 12;      // 一度に産める上限
const float BROOD_TAX       = 0.22f;   // 多産による取り分の目減り
// 生き物によらない死骸の減り(光分解・風化)。分解の主役は微生物と動物なので、これは日に当たる落葉の光分解の速さにとる。
// 乾燥地で生物の分解を除いた光分解だけの減りは k ≈ 0.1〜0.5 /年(Austin & Vivanco 2006 など)。その中ほどの 0.3 /年
const float ROT_RATE        = FRAC(0.3/YEAR_MIN);  // [率 5.7e-07 /分]
const int   N_SENSOR        = 5;       // 感覚器の枠。使うかどうかは進化が決める
const int   RAY_STEPS       = 3;       // レイの段数。多いほど重い
const float SENSOR_MAX      = 14.0f;   // 射程の上限
const float SENSOR_MIN      = 0.8f;
const int   N_TUNE          = 13;      // 視覚7 + 嗅覚4 + 同種 + 水
const float SENSOR_COST     = LIN(PM(0.0010f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 射程あたりの維持費 [1分あたり 4.86e-06]
// 入力: 感覚器 + エネルギー・気温・水・体格・バイアス
const int   N_REC           = 2;       // 再帰結合。何を覚えるかは進化が決める
const int   N_IN            = N_SENSOR+5+N_REC;
const int   N_OUT           = 4+N_REC; // 旋回・前進・採食・攻撃 + 再帰
const float ACT_THRESH      = 0.15f;   // これを超えた出力だけが実行される
// 神経網の形も進化する。枠を固定し、各ニューロンに有効・無効の遺伝子を持たせる。
// 無効なニューロンは何も通さず、維持費もかからない。
// 深さと幅が、選択によって決まる
const int   MAX_LAYER       = 3;       // 隠れ層の最大の深さ
const int   MAX_NODE        = 8;       // 一層あたりの最大の幅
const int   N_NODE          = MAX_LAYER*MAX_NODE;
const float NODE_ON         = 0.35f;   // これを超えたニューロンが働く
const int   W_L0            = MAX_NODE*N_IN;              // 入力 → 第1層
const int   W_LL            = (MAX_LAYER-1)*MAX_NODE*MAX_NODE; // 層 → 層
const int   W_OUT           = N_OUT*MAX_NODE;
const int   W_DIR           = N_OUT*N_IN;
const int   N_W             = W_L0+W_LL+W_OUT+W_DIR;
const float NODE_COST       = LIN(PM(0.0022f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 軍拡で代謝が膨らむと、脳が真っ先に切られる [1分あたり 1.07e-05]
const float SYN_COST        = LIN(PM(0.00010f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // [1分あたり 4.86e-07]
const int   ANIMAL_RELEASE  = 20000;
const int   ANIMAL_INIT_N   = 500;

// --- 水域と地形 ---
const float WATER_ATTEN     = 3.0f;    // 水深による光の減衰
const float BUOY_COST       = LIN(PM(0.030f));  // 浮くための組織の維持費 [1分あたり 0.000292]
// 水面に浮く物は、地上 10 m の風の約 3% の速さで風下へ流される(吹送流と波のストークス・ドリフトを合わせた表層の漂流。
// Wu 1983、油の漂流の目安 3%)。浮く強さ(buoyancy)は水面から出ている割合として、この速さに掛ける
const double SURF_DRIFT_FRAC = 0.03;
const double SURF_DRIFT_CPM = SURF_DRIFT_FRAC*WIND_SURF_MS*60.0/CELL_M;   // 水面の漂流の速さ [1分あたり 0.54 マス]
const float BUOY_DRIFT      = LIN(SURF_DRIFT_CPM,PLANT_STRIDE);  // 1 回の更新で流されるマス数 [1分あたり 0.54 マス]
const float HABITAT_COST_P  = LIN(PM(0.35f));  // 生息環境の不適合(植物) [1分あたり 0.0034]
const float HABITAT_COST_A  = LIN(PM(0.10f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 生息環境の不適合(動物) [1分あたり 0.000486]
const float SLOPE_COEF      = 0.50f;   // 勾配による減速。45°(勾配 1)で速さが 2/3。小さな昆虫は体重あたりの登坂の仕事が小さい
const float WATER_SPEED     = 1.30f;   // 水生が水中で得る速度
const float SHALLOW_EASE    = 0.30f;   // 浅瀬の緩さ(小さいほど岸が入りやすい)

// --- 匂い ---
// 4本のチャンネル。視覚と違い、時間に残り、遮蔽を回り込む
const int   N_ODOR          = 4;   // 0=生体 1=腐敗 2=食害 3=災害
const char* ODOR_NAME[4]    = {"生体","腐敗","食害","災害"};
const float ODOR_DECAY      = FRAC(PMF(0.055f,ODOR_INTERVAL),ODOR_INTERVAL);  // 災害(煙)の散逸。寿命 約 5 日 [率 0.000137 /分]
// 生体・腐敗・食害の匂いは揮発性の有機物で、大気の中で OH ラジカル(昼)とオゾン(昼夜)に酸化されて消える。
// 1 分あたりの率 = (k_OH·[OH] + k_O3·[O3])·60。速度定数は 298 K の値(Atkinson & Arey 2003)。
// 代表の物質: 生体 = α-ピネン(植物の放出の主なモノテルペン)、腐敗 = 二硫化ジメチル(死骸の匂いの主成分)、
// 食害 = (Z)-3-ヘキセノール(傷ついた葉の匂い)。煙は粒子で、酸化ではなく沈着で消えるので ODOR_DECAY のまま
// 条件: [OH] は太陽高度の正弦に比例し、真上の太陽で 4e6 個/cm³(中緯度で一日の平均が約 1e6。Prinn et al. 2005)。
// [O3] は背景の 30 ppb = 7.4e11 個/cm³(地表の空気 2.46e19 個/cm³)で、昼夜一定とする
const double ODOR_OH_PEAK   = 4.0e6;
const double ODOR_O3        = 30e-9*2.46e19;
const double ODOR_K_OH[3]   = {5.3e-11, 2.3e-10, 1.08e-10};   // cm³/個/s
const double ODOR_K_O3[3]   = {8.7e-17, 0.0,     6.4e-17};    // 二硫化ジメチルとオゾンの反応は無視できるほど遅い
// 匂いを運ぶ風。地表を歩く虫の匂いは、地表近くの空気の層で運ばれる。
// 地表近くの風は対数分布 u(z) = u*/κ·ln(z/z0) に従う(接地境界層の相似則)。u* は地上 10 m の風(WIND_SURF_MS、条件)から決める。
// 粗度 z0 = 0.01 m(草の短い開けた土地。Wieringa 1993)、匂いの層は地表から 0.1 m までとし、その層の平均の速さで運ぶ
const double ODOR_Z0        = 0.01;   // 地表の粗度(m)
const double ODOR_LAYER_M   = 0.10;   // 匂いを運ぶ地表近くの層の厚さ(m)
const double ODOR_USTAR     = 0.40*WIND_SURF_MS/std::log(10.0/ODOR_Z0);   // 摩擦速度(m/s)。カルマン定数 0.40
const double ODOR_WIND_MS   = ODOR_USTAR/0.40*(ODOR_LAYER_M*std::log(ODOR_LAYER_M/ODOR_Z0)-ODOR_LAYER_M+ODOR_Z0)/(ODOR_LAYER_M-ODOR_Z0);  // 層の平均 約 0.68 m/s
const float ODOR_SHIFT      = LIN(ODOR_WIND_MS*60.0/CELL_M,ODOR_INTERVAL);  // 1 回の更新で匂いが進むマス数 [1分あたり 約 4.1 マス]
// 地上 10 m の突風 1 m/s が、1 回の更新で匂いを運ぶマス数。同じ対数分布で匂いの層の平均へ写す
const float GUST_ODOR_CELLS = LIN(ODOR_WIND_MS/WIND_SURF_MS*60.0/CELL_M,ODOR_INTERVAL);
// 横への広がりは層の中の乱流の拡散。渦拡散係数は混合距離の理論で K = κ·u*·z(Prandtl。接地境界層の中立の場合)とし、
// 層の上端 z = ODOR_LAYER_M で評価する(約 0.007 m²/s)。周りの 3×3 の平均へ寄せる差分では、平均 − 中心 ≈ Δx²∇²/3 なので、
// 1 回の更新で寄せる割合 = 3·K·Δt/Δx²(陽解法。1 を超えると不安定になるが、この値は 0.05 程度)
const double ODOR_K_M2S     = 0.40*ODOR_USTAR*ODOR_LAYER_M;
const float ODOR_DIFFUSE    = (float)(3.0*ODOR_K_M2S*60.0*ODOR_INTERVAL*DT/((double)CELL_M*CELL_M));  // [1 回の更新で 約 0.050]
const float ODOR_SCALE      = 2.5f;   // 感覚器が飽和する濃さ
// 放出量
const float ODOR_PLANT      = LIN(PM(0.010f,ODOR_INTERVAL),ODOR_INTERVAL);  // 植物の体臭(高さと毒に比例) [1分あたり 2.43e-05]
const float ODOR_ANIMAL     = LIN(PM(0.035f,ODOR_INTERVAL),ODOR_INTERVAL);  // 動物の体臭(体格に比例) [1分あたり 8.51e-05]
const float ODOR_ESTRUS     = LIN(PM(0.055f,ODOR_INTERVAL),ODOR_INTERVAL);  // 発情 [1分あたり 0.000134]
const float ODOR_SICK       = LIN(PM(0.045f,ODOR_INTERVAL),ODOR_INTERVAL);  // 病んだ個体 [1分あたり 0.000109]
const float ODOR_ROT        = LIN(PM(0.008f,ODOR_INTERVAL),ODOR_INTERVAL);  // 死骸 [1分あたり 1.94e-05]
const float ODOR_WOUND      = 0.90f;  // 食害を受けた植物
const float ODOR_SMOKE      = LIN(PM(0.45f,ODOR_INTERVAL),ODOR_INTERVAL);  // 煙 [1分あたり 0.00109]

// --- 昼夜 ---
// 太陽高度 sin(h) = sinφ·sinδ + cosφ·cosδ·cosH を各マスで解く。
// 昼夜の境界、極夜、白夜がすべてこの一式から出る
const int   DAY_PERIOD      = (int)TICKS(DAY_MIN);    // 自転の周期(tick)。1 日 = DAY_MIN 分
const float AXIAL_TILT      = 23.4f;   // 地軸の傾き(度)。極夜・白夜の範囲を決める
const float SOLAR_GAIN      = LIN(PM(2.30f));  // 正午に飽和で捨てるぶんを補う [1分あたり 0.0224]
const float NIGHT_RESP      = 0.35f;   // 暗所では光化学系が止まり、呼吸が落ちる
const float DIURNAL_LAND    = 9.0f;    // 陸の日較差(度)
const float DIURNAL_SEA     = 1.2f;    // 海は熱容量が大きく、ほとんど振れない

// --- 植物の新しい戦略 ---
const float PARASITE_COST   = LIN(PM(0.030f));  // 吸器の維持費 [1分あたり 0.000292]
const float PARASITE_DRAIN_K = (float)PMF(0.055f);  // 宿主から奪う 1 分あたりの率の基準(貫く強さを掛ける)[率 0.00055 /分]
const float CARNIVORY_COST  = LIN(PM(0.045f));  // 捕虫器の維持費 [1分あたり 0.000437]
const float CARNIVORY_CATCH = FRAC(PMF(0.14f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 捕獲率 [率 0.000733 /分]
const float CARNIVORY_N     = 0.09f;   // 獲物から得る窒素
const float FIRE_SEED_COST  = LIN(PM(0.012f));  // 休眠種子の維持費 [1分あたり 0.000117]
const int   BURN_WINDOW     = (int)TICKS(900*DT0);  // [64.3 日]
const float BANK_COST       = LIN(PM(0.0035f));  // 土中でも呼吸し、食われ、朽ちる [1分あたり 3.4e-05]
const float BANK_PREDATION  = FRAC(PMF(0.0012f),SEED_STRIDE);  // 種子食者に見つかる確率(SEED_STRIDE tick ごとの 1 回の判定あたり)[率 1.17e-05 /分]
const int   BANK_MAX        = (int)TICKS(1400*DT0);  // 待てる上限 [100 日]

// --- 動物の新しい戦略 ---
const float TOXIN_COST      = LIN(PM(0.035f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 毒の合成費 [1分あたり 0.00017]
const float TOXIN_HARM      = 0.55f;   // 毒が捕食者に与える損害
const float TOXRES_COST     = LIN(PM(0.022f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 解毒の維持費 [1分あたり 0.000107]
const float BURROW_COST     = LIN(PM(0.020f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 掘る体の維持費 [1分あたり 9.72e-05]
const float BURROW_BUFFER   = 0.75f;   // 地中が和らげる気温の割合
const float FAT_COST        = LIN(PM(0.012f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 脂肪を担ぐ維持費 [1分あたり 5.83e-05]
// 脂肪の目標は体の維持費の何日分か(D-2026-10-24-3 の作者の答え)。日数は遺伝子 store_fat(0〜1)× FAT_DAYS_MAX。
// 昆虫の絶食に耐える日数は数日(Drosophila、Rion & Kawecki 2007)から数か月(休眠する種、Hahn & Denlinger 2011)まで幅があるので、
// 上限を 60 日に置き、どこを選ぶかは進化に任せる [条件の幅]
constexpr double FAT_DAYS_MAX = 60.0;
// 糖(グリコーゲンとトレハロース)の短い貯え。絶食した昆虫は先に糖を約 1 日で使い切り、そのあと脂肪を使う
// (Arrese & Soulages 2010)。エネルギーがこの日数分を超えた余りを脂肪に回し、下回れば脂肪を取り崩す [法則]
constexpr double GLYCO_DAYS   = 1.0;
// 1 回の更新の維持費を 1 日分に直す倍率(動物は ANIMAL_STRIDE tick ごとに更新する)
const float UPD_PER_DAY     = (float)(DAY_MIN/(ANIMAL_STRIDE*DT));
const float FAT_RATE        = FRAC(PMF(0.12f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 余剰を脂肪に回す速さ [率 0.000621 /分]
const float NOCTURNAL_COST  = LIN(PM(0.014f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 夜目の維持費 [1分あたり 6.81e-05]

// --- 系統 ---
const float SPECIES_THRESHOLD = 0.85f;
const int   REP_N             = 14;   // Lineage の rep / mean の長さ
const int   LINEAGE_INTERVAL  = 300;
const int   LINEAGE_SAMPLE    = 5000;
const int   LINEAGE_GRACE     = 3;      // 連続でこの回数見失って初めて絶滅とする
// --- コードの中に直に書かれていた速さ(段階 3 で名前を付けた。値は同じ) ---
const float MICROBE_DECOMP_K = (float)PMF(0.10f);  // 微生物が死骸を分解する 1 分あたりの率の基準(rate·act を掛ける)[率 1.02e-03 /分]
// 間引いて MICROBE_STRIDE tick に 1 回しか処理しないぶんは、その間に分解される割合 1−exp(−k·STRIDE·DT) で取り戻す。
// 以前はこれに加えて BOOST=2(「間引いたぶん」)を掛けていたが、間引きの換算は STRIDE で済んでいるので二重だった
// 低温に耐える酵素(不凍タンパク質・低温で働く酵素)の維持費。維持の呼吸は維持する細胞の量に比例する(Pirt 1965)。
// 以前は 1 個体あたりの定数 1.46e-05 /分で、体の炭素が小さい細胞ほど割高だった(体の炭素 0.1 で基礎の維持費の約 9 倍)。
// 分裂の周期の平均の体の炭素(分裂直後 0.7 から閾値 1.4 まで、平均 0.75·MICROBE_DIV = 1.05)で以前と同じ値になる比をとる
const float MICROBE_COLD_COST = LIN(PM(0.0015f));  // 体の炭素 MICROBE_REF の細胞の値 [1分あたり 1.46e-05]
// 微生物の入れ替わり。土の細菌や菌は原生生物・線虫に食われ、ウイルスに溶かされて、1 日あたり数 % ずつ死ぬ
// (土の微生物バイオマスの回転時間は数日〜数十日。Rousk & Bååth 2011、Kuzyakov & Blagodatskaya 2015)。その中ほどの 0.05 /日
// 食った側の呼吸と、溶けた細胞を他の微生物が分解するぶんは、低温で死んだ細胞と同じく全量を CO2 として大気へ戻す
// 入れ替わるのは活動している細胞だけ。基質の無い所で休眠している細胞や胞子は、食われにくく溶けにくく、
// 何か月〜何年も生き残る(Blagodatskaya & Kuzyakov 2013、Lennon & Jones 2011)。基質のあるマスにいるときだけ死ぬ
const float MICROBE_TURNOVER = FRAC(0.05/DAY_MIN,MICROBE_STRIDE);  // 1 回の更新(STRIDE tick)で死ぬ確率 1−exp(−r·STRIDE·DT) [率 3.47e-05 /分]
// 1 個体が 1 分に取り込める炭素の上限(体の炭素あたり)。取り込みは細胞の輸送と酵素の量で飽和する(Monod の Vmax)。
// 土の微生物の最大比増殖速度は、基質を加えた呼吸の解析で 0.1〜0.4 /時(Blagodatsky et al. 2000、Panikov 1995)。
// 下端の 0.1 /時を炭素利用効率 MICROBE_EFF で割り、取り込みの上限とする。以前は上限が無く、マスの死骸が多いほど際限なく取れた
const float MICROBE_VMAX = (float)(0.1/60.0/0.40);  // 1 分あたり、体の炭素あたりの取り込みの上限 [4.17e-03 /分]
// 死骸への胞子の定着。空気には菌の胞子が 1 m³ あたり 10^3〜10^4 個ある(Lacey 1996)。終端速度 1 mm/s で落ちると、
// 1 マス(100 m²)に 1 分あたり 6000 個以上が降るので、どの死骸にも胞子はすぐ届く。定着の遅れを決めるのは、
// 届いた胞子が死骸を食って、ここで 1 個体として数える群落(体の炭素 0.4)まで育つ時間のほう。
// 群落 0.4 は約 2.4 mg C(動物の尺度: 体 1 = 生重 100 mg、組織 1 あたり CORPSE_RETURN 2.5)、直径 10 µm の胞子は約 7.5e-8 mg C。
// その比 3e7 の対数 17.3 を、最大比増殖速度の下端 0.1 /時(MICROBE_VMAX と同じ根拠)で割った 173 時間(約 7 日)を、定着の平均の時間とする。
// 以前は、分裂のときにしか広がらず、最初の 3000 個体から離れた死骸(r6 の測定で死骸のあるマスの 99.9%)は分解されなかった。
// 実際の土と落葉は、どこでも数日のうちに微生物に取りつかれる(Voříšková & Baldrian 2013)
const float MICROBE_COLONIZE_LAG = (float)(17.3/0.1*60.0);   // [分] 173 時間
const float MICROBE_COLONIZE = FRAC(1.0/MICROBE_COLONIZE_LAG,MICROBE_STRIDE);   // 1 回の確かめ(STRIDE tick)で定着する確率 [率 9.63e-05 /分]
const float MICROBE_COLONY_E = 0.40f;   // 1 個体として数え始める群落の体の炭素(撒く初めの個体と同じ)
const float STORE_DRAW      = FRAC(PMF(0.02f));   // 土の水が足りないとき、不足分のうち体の貯水から補う割合
const float STORE_SURPLUS   = FRAC(PMF(0.5f,PLANT_STRIDE),PLANT_STRIDE);   // 余った土の水のうち、1 回の更新で貯水に回せる割合
const float FIRE_SHOOT_LOSS = LIN(PM(0.60f,PLANT_STRIDE),PLANT_STRIDE);  // 炎が幹と葉を焼く速さ(高さ/火勢)
const float N_UPTAKE        = LIN(PM(0.010f));    // 根が土から窒素を吸う速さの上限(同じマスの株で分ける)
const float CROWD_COST      = LIN(PM(0.020f));    // 自己間引き。混み合った株 1 本あたりの消耗
const float SAT_FLUX        = LIN(PM(1.0f));      // 光合成の飽和(saturation の値 1 が 1 tick に受けられる光)
const float CLIMB_RATE0     = LIN(PM(0.08f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 登る速さ(高さ/更新)
const float CLIMB_RATE1     = LIN(PM(0.25f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 登攀の巧さで足される速さ
const float DESCEND_RATE    = LIN(PM(0.15f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 降りる速さ
const float BURROW_RATE     = FRAC(PMF(0.15f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 潜り具合が遺伝子の値へ近づく割合
const float DIGEST_COST     = LIN(PM(0.018f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 長い腸の維持費
const float FAT_DRAW        = LIN(PM(0.05f,ANIMAL_STRIDE),ANIMAL_STRIDE);  // 飢えたとき脂肪を取り崩す速さ
const int   BANK_MIN        = (int)TICKS(40*DT0);     // 土中で待つ最短
const int   SEED_RETRY      = (int)TICKS(48*DT0);     // 乾いていて発芽できなかった種子が、次に土を見るまで
const int   RELEASE_TICK    = (int)TICKS(150*DT0);    // 植物が増えなくても、この時刻には動物を放つ

// --- 記録 ---
const int   HIST_MAX        = 200;   // 300tick刻みで6万tick。直近の変化が見える
const float SENESCE_START   = 0.55f;   // 寿命のこの割合を過ぎると衰えが始まる
const float SENESCE_RATE    = 2.4f;    // 衰えの速さ
const float LIFESPAN_COST   = LIN(PM(0.030f));  // 長命な体の維持費 [1分あたり 0.000292]
const int   LIFE_MIN        = (int)TICKS(900*DT0);  // 最短の寿命 [64.3 日]
const int   LIFE_MAX        = (int)TICKS(40000*DT0);  // 最長の寿命 [2.86e+03 日]
// 動物(地面を歩く昆虫程度)の寿命の幅。植物の幅(惑星版のまま)とは分ける。生理の条件
// 短い側: 卵から死ぬまで約 30 日(Drosophila melanogaster、25℃で発生 10 日+成虫の中央値 30〜40 日の下端)
// 長い側: 約 3 年(オサムシ科の成虫は 2〜3 年、まれに 4 年生きる。Lövei & Sunderland 1996)
const int   LIFE_MIN_A      = (int)TICKS(30*DAY_MIN);     // [30 日]
const int   LIFE_MAX_A      = (int)TICKS(3*YEAR_MIN);     // [1095 日]
const int   SAMPLE_INTERVAL = 300;
const int   AUTOSAVE_EVERY  = 7200;
const int   PRINT_EVERY     = 60;
