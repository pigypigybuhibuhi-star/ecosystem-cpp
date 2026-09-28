#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <cstdio>
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

// ==========================================================
//  パラメータ
// ==========================================================
const int   GRID_W          = 640;
const int   GRID_H          = 360;

// --- 地形 ---
const float SEA_LEVEL       = 0.42f;
const float MOUNTAIN_LEVEL  = 0.72f;
// プレートテクトニクス。山脈・海溝・海嶺・平坦な深海底は、
// すべてプレート境界での相対運動から出る
const int   N_PLATES        = 24;      // プレートの枚数
const float OCEANIC_FRAC    = 0.72f;   // 海洋プレートの割合
const float RIDGE_HEIGHT    = 0.045f;  // 海嶺。現実では海面下2,500mに留まる
const float ISLAND_ARC      = 0.30f;   // 島弧。日本列島やインドネシア
const int   N_HOTSPOT       = 5;       // ホットスポットの数
const int   HOTSPOT_CHAIN   = 9;       // 一列に並ぶ島の数
const float HOTSPOT_GAP     = 11.f;    // 島の間隔
const float HOTSPOT_R       = 6.5f;    // 島の大きさ
const float HOTSPOT_H       = 0.42f;   // 火山の高さ
const float TRENCH_DEPTH    = 0.15f;   // 海溝(沈み込み帯)
const float RIFT_DEPTH      = 0.10f;   // 大陸の裂け目
const float BOUNDARY_W      = 5.0f;    // 境界が効く幅(マス)
// 細部は標準偏差で正規化する。オクターブを増やしても振幅が崩れない
const float DETAIL_OCEAN    = 0.004f;  // 海底の細部。平坦さを保つため極小
const float DETAIL_LAND     = 0.075f;  // 陸の細部
const int   BASE_BLUR       = 16;      // 地殻の境界を均す
const int   FINAL_BLUR      = 1;       // 仕上げ。強くかけると海岸線が直線化する
const float WARP_AMP        = 32.f;    // 座標の歪み。Voronoiの直線を蛇行させる
const int   WARP_R          = 20;      // 歪みの波長
// --- 地形発達 ---
// 地殻の厚さから標高を決め(アイソスタシー)、低海面で川に刻ませ、
// そのあと海面を上げて谷を沈ませる。複雑な海岸線は溺れた谷から出る
const float CRUST_CONT      = 38.0f;   // 大陸地殻の厚さ
const float CRUST_OCEAN     = 7.0f;    // 海洋地殻の厚さ
const float CRUST_VARY      = 4.0f;    // 大陸ごとのばらつき
// 大陸内部の構造。一様な塊ではなく、古い核と若い帯と沈んだ盆地でできている
const float CRUST_INTERNAL  = 9.0f;    // 大陸内部の厚さの起伏
const int   CRUST_SCALE     = 26;      // その波長
const int   N_RIFT          = 7;       // 失敗した地溝の数
const float RIFT_THIN       = 7.5f;    // 地殻がどれだけ薄くなるか
const float RIFT_LEN        = 70.f;    // 地溝の長さ
const float RIFT_WIDTH      = 7.0f;    // 幅
const float CRATON_THICK    = 5.0f;    // 楯状地の地殻の厚み増し
const float CRATON_RESIST   = 0.25f;   // 侵食への抵抗(古い岩盤は硬い)
const int   CRATON_SCALE    = 34;
const float OROGEN_THICK    = 22.0f;   // 衝突でどれだけ厚くなるか
const float RHO_RATIO       = 0.1515f; // 1 - ρ地殻/ρマントル
const float SUBSIDE_K       = 0.9f;    // 熱的沈降。√年齢に比例
const float SUBSIDE_MAX     = 60.f;    // 沈降が飽和する年齢
// Cordonnier et al. 2016 のストリームパワー式。n=1, m=0.5
const float EROSION_M       = 0.5f;    // 流域面積の指数
const float EROSION_DT      = 2.5f;    // 時間刻み。陰的解法なので大きく取れる
const float TARGET_RELIEF   = 0.42f;   // 目指す最高峰の高さ(海面から)
const float MAX_SLOPE       = 0.022f;  // 熱的侵食。これ以上の勾配は崩れる
const float DIFFUSE_K       = 0.020f;  // 斜面拡散
const int   EROSION_STEPS   = 250;      // 侵食を回す回数
const float UPLIFT_CONT     = 0.0045f; // 大陸地殻全体の隆起(アイソスタシー)
const float UPLIFT_OROGEN   = 0.0085f; // 衝突境界での造山
const float GLACIAL_LEVEL   = 0.52f;   // 刻むときの低海面(分位)

// --- 海流 ---
const int   CURRENT_ITER    = 90;      // 熱が大洋を横断するまで回す
const float CURRENT_ADV     = 0.55f;
const float WBC_RANGE       = 14.f;    // 西岸強化が効く距離
const float CURRENT_HEAT    = 9.0f;    // 海流が動かす気温の幅(度)
const float OCEAN_FRACTION  = 0.64f;
// --- 気候 ---
const float TEMP_EQUATOR    = 30.0f;
const float TEMP_POLE       = -15.0f;
const float LAPSE_RATE      = 28.0f;
const int   SEASON_PERIOD   = 500;
const float SEASON_AMP      = 24.0f;   // 年較差。大きいほど亜寒帯が広がる
const float SEA_SEASON      = 0.30f;   // 海の年較差は陸の3割ほど
// 海面水温は気温と別の曲線を描く。凍結点より下がれないため、極でも−数度で止まる
const float SST_EQUATOR     = 29.0f;
const float SST_POLE        = -3.5f;
const float ICE_LAG         = 0.055f;  // 氷の成長と融解は季節に遅れる
const int   COAST_RANGE     = 7;       // 陸の冷気が海へ及ぶ距離
const float COAST_CHILL     = 0.55f;   // 沿岸の海が陸の気温に引かれる強さ

// --- 水循環 ---
const float EVAP_COEF       = 0.050f;  // 海面温度あたりの蒸発量
const float MOIST_CAP       = 8.00f;   // 大気が運べる水蒸気の上限
const float LAND_RECYCLE    = 0.040f;  // 陸からの蒸発散(水分の再循環)
const float RAIN_BASE       = 0.024f;  // 対流による降水率
const float OROGRAPHIC      = 10.0f;   // 地形性降雨(斜面を上る空気は雨を落とす)
const int   ORO_SPAN        = 6;       // 斜面を測る幅。広いほど雨陰がなだらかになる
const int   PRECIP_BLUR     = 6;       // 降水場のぼかし半径
const int   PRECIP_INTERVAL = 300;     // 降水場を計算し直す間隔
// ケッペンの乾燥限界(20T+140)にならい、限界は気温に比例させる
const float ARID_COEF       = 0.0010f; // 乾燥限界の係数
const float ARID_SOFT       = 1.20f;   // 大きいほど乾湿の差がなだらかになる
const float FIELD_CAP       = 1.00f;   // 土壌が保持できる上限
const float WATER_RELAX     = 0.020f;  // 土壌水分が目標へ近づく速さ
const float RIVER_THRESHOLD = 260.f;   // 侵食後の地形に合わせて下げる
const float RIVER_SUPPLY    = 0.60f;   // 河川が供給する水
const float DROUGHT_COST    = 0.055f;  // 耐乾装置の維持費
const float STORE_COST      = 0.030f;  // 貯水組織の維持費
const float STORE_FILL      = 0.020f;  // 湿っているとき貯める速さ
const float STORE_MAX       = 3.00f;   // 遺伝子が届く貯水量の上限

// --- 大気 ---
const float CO2_BASE        = 8.00f;   // 大気は巨大な貯蔵庫
const float O2_BASE         = 8.00f;
const float CO2_PER_ENERGY  = 0.020f;
const float O2_PER_METAB    = 0.020f;
const float N2_BASE         = 24.0f;   // O2分圧を地球並み(約21%)にする
const float SOIL_N_BASE     = 0.35f;   // 土壌の窒素
const float N_PER_GROWTH    = 0.055f;  // 体を作るのに要する窒素
const float N_FIX_COST      = 0.006f;  // 純益0.048に対して現実的な負担に
const float N_FIX_RATE      = 0.0180f; // 共生菌は効率がよい
const float N_MINERALIZE    = 0.035f;
const float N_SINK          = 0.022f;  // 沈降。外洋を貧栄養にするが、枯渇させない
const float N_RUNOFF        = 0.045f;  // 河川が沿岸へ運ぶ栄養
const float N_UPWELL        = 0.040f;  // 湧昇。深層水が栄養を持ち上げる
const int   SHELF_RANGE     = 6;

// --- 海氷 ---
// 海氷の形成には、継続した冷却が要る。瞬間の水温より低い閾値で表す
const float FREEZE_TEMP     = -5.0f;
const float THAW_TEMP       = -1.5f;
const float ICE_LIGHT       = 0.08f;   // 氷が通す光

const float PLANT_RESP_CO2  = 0.35f;   // 植物の呼吸(維持コストのうちCO2を出す割合)
const float RESP_EFFICIENCY = 0.60f;   // 呼吸効率の遺伝子が持つ効き
const float RESP_COST       = 0.032f;  // 高効率な呼吸器の維持費
const float DISSOLVED_O2    = 0.45f;   // 鰓は水から効率よく酸素を取り出す
const float ALT_O2_DECAY    = 2.20f;   // 標高による空気の薄まり
const float GAS_DIFFUSE     = 0.55f;   // 4tickに一度呼ぶぶん強める
const float GAS_RELAX       = 0.0060f;
const float WIND_STRENGTH   = 0.35f;   // 拡散に負けない強さにする
const float O2_PER_CO2      = 1.00f;

// --- 火 ---
// 酸素が17%未満では燃え広がらず、35%を超えると消火が困難になる(火の窓)
// 地球の21%を1.0とした比。17%未満では広がらず、35%超で消火困難
const float FIRE_O2_MIN     = 0.81f;   // 17% / 21%
const float FIRE_O2_WET     = 1.67f;   // 35% / 21%
const float FIRE_IGNITE     = 8.0e-8f; // 落雷による着火率/マス/tick
const float FIRE_SPREAD     = 0.20f;   // 隣へ移る強さ
const float FIRE_BURN       = 0.45f;   // 一tickで燃える燃料の割合
const float FIRE_DECAY      = 0.38f;   // 火勢の減衰。燃え尽きたら消える
const float FIRE_MOIST      = 0.55f;   // これ以上湿っていると燃えにくい
const float FIRE_O2_USE     = 0.025f;  // 燃焼が消費する酸素
const float FIRE_TOL_COST   = 0.075f;  // 厚い樹皮は高くつく
const float FIRE_DAMAGE     = 9.0f;    // 炎そのものの破壊力
const float FIRE_FUEL_MIN   = 1.5f;

// --- 撹乱 ---
const float SLIDE_SLOPE     = 0.045f;  // これ以上の傾斜で崩れうる
const float SLIDE_RATE      = 1.2e-5f; // 崩落率/マス/tick

// --- 変異 ---
const float MUT_STRENGTH    = 0.02f;
const float MUT_RATE_MAX    = 3.0f;    // 変異率の上限(基準の何倍まで)
const float MUT_RATE_COST   = 0.004f;  // 修復機構を緩めることの代償

// --- 植物 ---
const float GROW_RATE       = 0.02f;
const float GROW_COST       = 3.00f;
const float HEIGHT_COST     = 0.060f;
const float RESP_COEF       = 0.120f;  // 真の昼夜で取り込みが半減したぶんを補正
const float TOUGH_COST      = 0.070f;  // 防御を「今」張るための費用
const float INDUCE_RATE     = 0.22f;   // 食害の匂いに応じて立ち上がる速さ
const float INDUCE_DECAY    = 0.020f;  // 攻撃がなければ緩む
const float CLONAL_COST     = 0.010f;  // 地下茎の維持費
const float CLIMB_COST      = 0.022f;  // つるの維持費(巻きひげと気根)
const float CLIMB_GAIN      = 0.85f;   // 支えの高さをどれだけ利用できるか
const float CLIMB_PENALTY   = 0.45f;
const float COLD_COST_P     = 0.060f;
const float DORM_UPKEEP     = 0.12f;   // 休眠中の維持費の割合
const float DORM_HARDEN     = 30.0f;   // 休眠がもたらす耐凍性
const float DORM_HYST       = 3.0f;    // 覚醒に要する余裕(ばたつき防止)
const float GENET_SHARE     = 0.20f;
const float SEED_FAR        = 12.0f;
const float PLANT_DIVIDE_TH = 8.0f;
const float PLANT_INIT_E    = 6.0f;    // 立ち上がりの緩衝
const int   PLANT_INIT_N    = 400;
const int   PLANT_MAX       = 260000;  // これを超えると繁殖を止める。観察のための上限
const int   ANIMAL_MAX      = 60000;
const float SEX_COST        = 0.015f;  // 生殖器官の維持費
const float SEX_SEARCH      = 3.0f;    // 相手を探す範囲(マス)

// --- 病原体 ---
// 病原体はその土地で最も多い免疫型に適応する。
// 稀な型が有利になり、組み換えで新しい型を作れる有性生殖に価値が生まれる
const int   N_IMMUNE        = 3;       // 免疫の座位
const float PATH_ADAPT      = 0.030f;  // 病原体が宿主の型へ寄る速さ
const float PATH_GROW       = 0.020f;  // 宿主密度に応じた増殖
const float PATH_DECAY      = 0.015f;  // 宿主がいなければ衰える
const float PATH_LOAD_MAX   = 2.50f;
const float PATH_VIRULENCE  = 0.055f;  // 型が一致したときの被害
const float PATH_WIDTH      = 0.42f;   // 一致とみなす幅(狭いほど型が効く)
const int   PATH_INTERVAL   = 3;       // 何tickごとに更新するか

// --- 樹上生活 ---
const float ARBOREAL_COST   = 0.028f;  // 登攀に適した体の維持費
const float CLIMB_UP_COST   = 0.55f;   // 高さ1を登るのに要するエネルギー
const float PERCH_GAP       = 0.45f;   // 枝から枝へ渡れる高さの差
const float BRANCH_SUPPORT  = 1.20f;   // 枝が支えられる体重の目安
const float FALL_DAMAGE     = 0.9f;    // 落下で失うエネルギー(高さあたり)
const float TEMP_WIDTH_P    = 18.0f;

// --- 動物 ---
const float A_MAX_SPEED     = 0.20f;
const float A_METAB         = 0.040f;
const float A_MOVE_COST     = 0.060f;
const float A_REACH_COST    = 0.018f;
const float A_SIZE_COST     = 0.020f;
const float JAW_COST        = 0.012f;
const float COLD_COST_A     = 0.050f;
const float A_BITE          = 0.45f;
const float GRAZE_DEFOLIATE = 0.18f;   // 採食が葉を奪う量。局所的な枯渇を生む
const float A_DIVIDE_TH     = 6.0f;
const float A_EAT_LOSS      = 0.50f;   // 同化効率50%。若い葉や種子は消化しやすい
const float PRED_RATIO      = 1.15f;
const float PRED_EFF        = 0.85f;
const float SCAV_EFF        = 0.85f;   // 死骸は抵抗しないので取り分が大きい
const float MEAT_LOSS       = 0.30f;
const float DETRITUS_LOSS   = 0.80f;   // 分解物は栄養価が中間
const float DETRITUS_COST   = 0.008f;  // 分解酵素の維持費
// 分解者。動かず、死骸だけを食い、炭素と窒素を土に返す
const float MICROBE_COST    = 0.0004f; // 微生物の代謝は極小。休眠に近い状態で待てる
const float MICROBE_EFF     = 0.60f;
const float MICROBE_DIV     = 0.70f;   // 分裂の閾値を下げ、立ち上がりを速くする
const float MICROBE_RETURN  = 0.55f;   // 分解して大気と土に返す割合
const float MICROBE_DISP    = 2.0f;    // 胞子が飛ぶ距離
const int   MICROBE_INIT    = 3000;
const float BODY_GROW_RATE  = 0.010f;
const float BODY_GROW_COST  = 4.00f;
const float CORPSE_RETURN   = 2.50f;
const float MATURE_FRAC     = 0.80f;
const int   MAX_BROOD       = 12;      // 一度に産める上限
const float BROOD_TAX       = 0.22f;   // 多産による取り分の目減り
const float ROT_RATE        = 0.004f;  // 分解者が主役。これは物理的な風化だけ
const int   N_SENSOR        = 5;       // 感覚器の枠。使うかどうかは進化が決める
const int   RAY_STEPS       = 3;       // レイの段数。多いほど重い
const float SENSOR_MAX      = 14.0f;   // 射程の上限
const float SENSOR_MIN      = 0.8f;
const int   N_TUNE          = 13;      // 視覚7 + 嗅覚4 + 同種 + 水
const float SENSOR_COST     = 0.0010f; // 射程あたりの維持費
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
const float NODE_COST       = 0.0022f; // 軍拡で代謝が膨らむと、脳が真っ先に切られる
const float SYN_COST        = 0.00010f;
const int   ANIMAL_RELEASE  = 20000;
const int   ANIMAL_INIT_N   = 500;

// --- 水域と地形 ---
const float WATER_ATTEN     = 3.0f;    // 水深による光の減衰
const float BUOY_COST       = 0.030f;  // 浮くための組織の維持費
const float BUOY_DRIFT      = 0.55f;   // 海流に流される強さ
const float HABITAT_COST_P  = 0.35f;   // 生息環境の不適合(植物)
const float HABITAT_COST_A  = 0.10f;   // 生息環境の不適合(動物)
const float SLOPE_COEF      = 25.0f;    // 傾斜による減速
const float WATER_SPEED     = 1.30f;   // 水生が水中で得る速度
const float SHALLOW_EASE    = 0.30f;   // 浅瀬の緩さ(小さいほど岸が入りやすい)

// --- 匂い ---
// 4本のチャンネル。視覚と違い、時間に残り、遮蔽を回り込む
const int   N_ODOR          = 4;   // 0=生体 1=腐敗 2=食害 3=災害
const char* ODOR_NAME[4]    = {"生体","腐敗","食害","災害"};
const float ODOR_DIFFUSE    = 0.42f;  // 拡散
const float ODOR_DECAY      = 0.055f; // 散逸
const float ODOR_WIND       = 0.30f;  // 風による移流
const int   ODOR_INTERVAL   = 4;      // 何tickごとに更新するか
const float ODOR_SCALE      = 2.5f;   // 感覚器が飽和する濃さ
// 放出量
const float ODOR_PLANT      = 0.010f; // 植物の体臭(高さと毒に比例)
const float ODOR_ANIMAL     = 0.035f; // 動物の体臭(体格に比例)
const float ODOR_ESTRUS     = 0.055f; // 発情
const float ODOR_SICK       = 0.045f; // 病んだ個体
const float ODOR_ROT        = 0.008f; // 死骸
const float ODOR_WOUND      = 0.90f;  // 食害を受けた植物
const float ODOR_SMOKE      = 0.45f;  // 煙

// --- 昼夜 ---
// 太陽高度 sin(h) = sinφ·sinδ + cosφ·cosδ·cosH を各マスで解く。
// 昼夜の境界、極夜、白夜がすべてこの一式から出る
const int   DAY_PERIOD      = 14;      // 自転の周期(tick)。夜が貯蔵を超えないように
const float AXIAL_TILT      = 23.4f;   // 地軸の傾き(度)。極夜・白夜の範囲を決める
const float SOLAR_GAIN      = 2.30f;   // 正午に飽和で捨てるぶんを補う
const float NIGHT_RESP      = 0.35f;   // 暗所では光化学系が止まり、呼吸が落ちる
const float DIURNAL_LAND    = 9.0f;    // 陸の日較差(度)
const float DIURNAL_SEA     = 1.2f;    // 海は熱容量が大きく、ほとんど振れない

// --- 植物の新しい戦略 ---
const float PARASITE_COST   = 0.030f;  // 吸器の維持費
const float PARASITE_DRAIN  = 0.055f;  // 宿主から奪う割合
const float CARNIVORY_COST  = 0.045f;  // 捕虫器の維持費
const float CARNIVORY_CATCH = 0.14f;   // 捕獲率
const float CARNIVORY_N     = 0.09f;   // 獲物から得る窒素
const float FIRE_SEED_COST  = 0.012f;  // 休眠種子の維持費
const int   BURN_WINDOW     = 900;
const float BANK_COST       = 0.0035f;  // 土中でも呼吸し、食われ、朽ちる
const float BANK_PREDATION  = 0.0012f;  // 種子食者に見つかる率/tick
const int   BANK_MAX        = 1400;    // 待てる上限(tick)

// --- 動物の新しい戦略 ---
const float TOXIN_COST      = 0.035f;  // 毒の合成費
const float TOXIN_HARM      = 0.55f;   // 毒が捕食者に与える損害
const float TOXRES_COST     = 0.022f;  // 解毒の維持費
const float BURROW_COST     = 0.020f;  // 掘る体の維持費
const float BURROW_BUFFER   = 0.75f;   // 地中が和らげる気温の割合
const float FAT_COST        = 0.012f;  // 脂肪を担ぐ維持費
const float FAT_MAX         = 8.0f;
const float FAT_RATE        = 0.12f;   // 余剰を脂肪に回す速さ
const float NOCTURNAL_COST  = 0.014f;  // 夜目の維持費

// --- 系統 ---
const float SPECIES_THRESHOLD = 0.85f;
const int   REP_N             = 14;   // Lineage の rep / mean の長さ
const int   LINEAGE_INTERVAL  = 300;
const int   LINEAGE_SAMPLE    = 5000;
const int   LINEAGE_GRACE     = 3;      // 連続でこの回数見失って初めて絶滅とする
// --- 記録 ---
const int   HIST_MAX        = 4000;
const int   PLANT_STRIDE    = 2;       // 植物の更新を何tickに分けるか
const int   SEED_STRIDE     = 8;       // 土中の種子はさらに間引く
const float SENESCE_START   = 0.55f;   // 寿命のこの割合を過ぎると衰えが始まる
const float SENESCE_RATE    = 2.4f;    // 衰えの速さ
const float LIFESPAN_COST   = 0.030f;  // 長命な体の維持費
const int   LIFE_MIN        = 900;     // 最短の寿命(tick)
const int   LIFE_MAX        = 40000;   // 最長の寿命
const int   SAMPLE_INTERVAL = 300;
const int   AUTOSAVE_EVERY  = 7200;
const int   PRINT_EVERY     = 60;
// ==========================================================

float SCREEN_W = 1280.f, SCREEN_H = 720.f;

std::vector<float> elev(GRID_W*GRID_H, 0.f);
std::vector<char>  sea (GRID_W*GRID_H, 0);
std::vector<float> solar(GRID_W*GRID_H, 0.f);
std::vector<float> temper(GRID_W*GRID_H, 0.f);
std::vector<float> co2(GRID_W*GRID_H, CO2_BASE);
std::vector<float> o2 (GRID_W*GRID_H, O2_BASE);
std::vector<float> gbuf(GRID_W*GRID_H, 0.f);
std::vector<float> wind(GRID_H, 0.f);
std::vector<float> slope(GRID_W*GRID_H, 0.f);
std::vector<float> precip(GRID_W*GRID_H, 0.f);
std::vector<float> soil_water(GRID_W*GRID_H, 0.f);
std::vector<float> soil_n(GRID_W*GRID_H, SOIL_N_BASE);
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
std::vector<float> odor[N_ODOR], odor_buf;
double odor_total[N_ODOR]={0,0,0,0};
std::vector<float> flow_acc(GRID_W*GRID_H, 0.f);
std::vector<char>  river(GRID_W*GRID_H, 0);
std::vector<float> temp_annual(GRID_W*GRID_H, 0.f);
std::vector<float> aridity(GRID_W*GRID_H, 0.f);
std::vector<unsigned char> climate(GRID_W*GRID_H, 0);
std::vector<float> elev_smooth(GRID_W*GRID_H, 0.f);
std::vector<int>   plate_id(GRID_W*GRID_H, 0);     // どのプレートに属するか
std::vector<char>  plate_oce(N_PLATES, 0);         // 海洋プレートか
std::vector<float> bnd_type(GRID_W*GRID_H, 0.f);
std::vector<float> tect_f(GRID_W*GRID_H, 0.f);
std::vector<float> hotspot(GRID_W*GRID_H, 0.f);
std::vector<float> craton(GRID_W*GRID_H, 0.f);   // 古い安定した大陸核。侵食に抗う
std::vector<float> crust(GRID_W*GRID_H, 0.f);     // 地殻の厚さ
std::vector<float> sst_anom(GRID_W*GRID_H, 0.f);  // 海流による水温偏差
std::vector<float> cur_u(GRID_W*GRID_H, 0.f);     // 海流の東西成分
std::vector<float> cur_v(GRID_W*GRID_H, 0.f);     // 海流の南北成分
// 0=海 1=EF氷雪 2=ETツンドラ 3=BW砂漠 4=BSステップ 5=D亜寒帯 6=C温帯 7=Awサバナ 8=Af熱帯雨林
inline int idx(int x,int y){ return y*GRID_W + x; }

std::mt19937 rng(12344);
std::uniform_real_distribution<float> dist01(0.f,1.f);
inline float frand(float a,float b){ return a + (b-a)*dist01(rng); }
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
    float arboreal;   // 遺伝子: 登攀に適した体
    float brood;      // 遺伝子: 一度に産む数。多産少投資か、少産多投資か
    float lifespan;   // 遺伝子: 寿命の長さ
    float toxin;      // 遺伝子: 毒を持つ
    float tox_resist; // 遺伝子: 毒に耐える
    float burrow;     // 遺伝子: 穴を掘る
    float nocturnal;  // 遺伝子: 夜に活動する
    float store_fat;  // 遺伝子: 余剰を脂肪に回す傾向
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
    float rec[N_REC];
    // 本当の系譜。形質による分類(lin)とは別に、実際の親子関係を持つ
    int   mother, father;   // 親の id。無性なら father=-1
    int   gen;              // 世代数
    int   id,age,eaten,lin;
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
}
inline void sg_vec_p(const Plant& p,double* v){
    v[0]=p.max_height; v[1]=p.shade_tol; v[2]=p.tough;      v[3]=p.clonal;
    v[4]=p.disperse;   v[5]=p.cold_tol;  v[6]=p.aquatic;    v[7]=p.drought_tol;
    v[8]=p.climb;      v[9]=p.parasite;  v[10]=p.carnivory; v[11]=p.brood;
}
std::vector<std::vector<int>> animal_grid(GRID_W*GRID_H);
int next_animal_id = 0;

struct Corpse { float x,y,energy; int origin; };   // 0=植物由来 1=動物由来
std::vector<Corpse> corpses;
std::vector<std::vector<int>> corpse_grid(GRID_W*GRID_H);
std::vector<float> corpse_density(GRID_W*GRID_H,0.f);
std::vector<float> fuel(GRID_W*GRID_H,0.f);
std::vector<float> fire(GRID_W*GRID_H,0.f);
std::vector<float> fire_buf(GRID_W*GRID_H,0.f);
std::vector<float> canopy_h(GRID_W*GRID_H,0.f);   // そのマスで最も高い植物の高さ
// 病原体: 型(免疫座位ごと)と密度。植物と動物で別の病がつく
std::vector<float> pstrain[N_IMMUNE], astrain[N_IMMUNE];
std::vector<float> pload(GRID_W*GRID_H,0.f), aload(GRID_W*GRID_H,0.f);
double infect_p=0, infect_a=0;
long fire_events=0; double fire_area=0;
long slide_events=0;

// 選択勾配。繁殖した個体と集団全体で遺伝子の平均を比べる。
// 差がそのまま「その形質が今どちらへ押されているか」になる。
// 収支を手計算せずに、効いているかどうかを測定で言えるようにする
const int N_SG_A = 12;
const char* SG_A_NAME[N_SG_A]={
    "diet","detr","reach","size","speed","jaw",
    "cold","aqua","resp","arbor","sexual","brood"};
const int N_SG_P = 12;
const char* SG_P_NAME[N_SG_P]={
    "maxH","shade","tough","clonal","disp","cold",
    "aqua","dry","climb","paras","carni","brood"};
double sg_a_pop[N_SG_A]={0}, sg_a_rep[N_SG_A]={0};
double sg_p_pop[N_SG_P]={0}, sg_p_rep[N_SG_P]={0};
long   sg_a_np=0, sg_a_nr=0, sg_p_np=0, sg_p_nr=0;
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
Series s_plant, s_animal, s_maxh, s_reach, s_temp, s_co2, s_o2, s_soil, s_nut,
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
long d_starve=0, d_cold=0, d_preyed=0;
// 植物の収支を測る。推測ではなく実測で判断するため
double dbg_gain=0, dbg_cost=0; long dbg_n=0;
double dbg_light=0, dbg_tf=0, dbg_wf=0, dbg_co2=0;

// 最後に生きていた動物の遺伝子を保存しておく(再投入用)
struct AnimalGene {
    float speed,turn,reach,size,diet,jaw,cold_tol,aquatic,resp,detritus,sexual,arboreal;
    float brood,lifespan;
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
        g.sexual=a.sexual;   g.arboreal=a.arboreal; g.lin=a.lin;
        g.brood=a.brood;     g.lifespan=a.lifespan;
        for(int k=0;k<N_IMMUNE;k++) g.im[k]=a.im[k];
        for(int k=0;k<N_W;k++) g.w[k]=a.w[k];
        for(int k=0;k<N_NODE;k++) g.ngain[k]=a.ngain[k];
        for(int k=0;k<N_SENSOR;k++) g.sen[k]=a.sen[k];
        last_animals.push_back(g);
    }
}

// ===== 系統 =====
struct Lineage {
    int   id, kind, num;      // kind 0=植物 1=動物 / num は種別ごとの通し番号
    float rep[REP_N];
    int   nrep, birth_tick, death_tick, parent, missing;
    bool  alive;
    int   pop, max_pop;
    float mean[REP_N];        // 表示用の平均形質
    float final_mean[REP_N];  // 絶滅した瞬間の平均。個体が消えても残す
    int   final_pop;          // 絶滅直前の個体数
};
std::vector<Lineage> lineages;
int next_lineage_id = 0;
int next_plant_num = 1, next_animal_num = 1;

inline int plant_vec(const Plant& p, float* v){
    v[0]=p.absorb         *1.0f;
    v[1]=p.max_height/5.f *2.5f;
    v[2]=p.shade_tol      *2.0f;
    v[3]=p.tough          *2.0f;
    v[4]=p.clonal         *1.5f;
    v[5]=p.disperse       *1.0f;
    v[6]=p.cold_tol       *1.5f;
    v[7]=p.aquatic        *3.0f;   // 生息環境は最も重い軸
    v[8]=p.drought_tol    *2.0f;
    return 9;
}
inline int animal_vec(const Animal& a, float* v){
    v[0]=a.speed          *1.0f;
    v[1]=a.turn           *0.5f;
    v[2]=a.reach/5.f      *2.0f;
    v[3]=a.size/4.f       *2.5f;
    v[4]=a.diet           *3.0f;
    v[5]=a.jaw/2.f        *1.5f;
    v[6]=a.cold_tol       *1.5f;
    // 神経網の規模を、結線の強さの平均で代表させる
    { float s=0.f; for(int k=0;k<N_W;k++) s+=std::fabs(a.w[k]);
      v[7]=clamp01(s/(N_W*1.5f))*0.8f; }
    // 神経網の形。深さと幅が、生態的な位置を分ける
    { int nn=0; for(int k=0;k<N_NODE;k++) if(a.ngain[k]>NODE_ON) nn++;
      v[8]=(float)nn/N_NODE *1.5f; }
    v[9]=a.aquatic        *3.0f;
    // 感覚器の構成。何をどこまで見るかは、生態的な位置を強く表す
    float rs=0.f, tn=0.f; int nv=0;
    for(int k=0;k<N_SENSOR;k++)
        if(a.sen[k].range>=SENSOR_MIN){ rs+=a.sen[k].range; tn+=a.sen[k].tune; nv++; }
    v[10]=rs/(SENSOR_MAX*N_SENSOR) *2.0f;
    v[11]=(nv? tn/nv : 0.f)        *2.5f;
    v[12]=a.arboreal               *2.0f;   // 地上性か樹上性かは生態的な位置を分ける
    return 13;
}
inline float vdist(const float* a,const float* b,int n){
    float s=0.f; for(int i=0;i<n;i++){ float d=a[i]-b[i]; s+=d*d; } return std::sqrt(s);
}
inline float depthn(int i){ return std::max(0.f,(SEA_LEVEL-elev[i])/SEA_LEVEL); }
inline float habitat_p(const Plant& p,int ci){
    if(sea[ci]) return (1.f-p.aquatic)*(SHALLOW_EASE+(1.f-SHALLOW_EASE)*depthn(ci));
    return p.aquatic;
}
inline float habitat_a(const Animal& a,int ci){
    if(sea[ci]) return (1.f-a.aquatic)*(SHALLOW_EASE+(1.f-SHALLOW_EASE)*depthn(ci));
    return a.aquatic;
}
sf::Color lineage_color(int i){
    float h=std::fmod(0.13f+(float)i*0.61803399f,1.0f)*360.f;
    float s=0.70f,v=1.0f;
    float c=v*s,x=c*(1.f-std::fabs(std::fmod(h/60.f,2.f)-1.f)),m=v-c;
    float r,g,b;
    if(h<60){r=c;g=x;b=0;} else if(h<120){r=x;g=c;b=0;}
    else if(h<180){r=0;g=c;b=x;} else if(h<240){r=0;g=x;b=c;}
    else if(h<300){r=x;g=0;b=c;} else {r=c;g=0;b=x;}
    return sf::Color((int)((r+m)*255),(int)((g+m)*255),(int)((b+m)*255));
}

void update_lineages(int kind,int tick){
    int total = (kind==0)?(int)plants.size():(int)animals.size();
    if(total<=0){
        for(auto& L:lineages)
            if(L.alive && L.kind==kind){ L.alive=false; L.death_tick=tick; }
        return;
    }
    const int D = (kind==0)?9:13;
    std::vector<float> reps;
    float v[14];
    int stride=std::max(1,total/LINEAGE_SAMPLE), seen=0;
    for(int i=0;i<total;i++){
        if((seen++ % stride)!=0) continue;
        if(kind==0) plant_vec(plants[i],v); else animal_vec(animals[i],v);
        int best=-1; float bd=1e9f; int nc=(int)reps.size()/D;
        for(int c=0;c<nc;c++){ float d=vdist(v,&reps[c*D],D); if(d<bd){bd=d;best=c;} }
        if(best<0 || bd>=SPECIES_THRESHOLD)
            for(int k=0;k<D;k++) reps.push_back(v[k]);
    }
    int nc=(int)reps.size()/D;
    int oldN=(int)lineages.size();
    std::vector<int>  cl(nc,-1);
    std::vector<char> claimed(oldN,0);
    for(int c=0;c<nc;c++){
        int best=-1; float bd=1e9f;
        for(int L=0;L<oldN;L++){
            if(!lineages[L].alive || lineages[L].kind!=kind) continue;
            float d=vdist(&reps[c*D],lineages[L].rep,D);
            if(d<bd){ bd=d; best=L; }
        }
        if(best>=0 && bd<SPECIES_THRESHOLD && !claimed[best]){
            cl[c]=best; claimed[best]=1; lineages[best].missing=0;
            for(int k=0;k<D;k++) lineages[best].rep[k]=reps[c*D+k];
        } else {
            Lineage L; L.id=next_lineage_id++; L.kind=kind; L.nrep=D;
            L.num = (kind==0)? next_plant_num++ : next_animal_num++;
            for(int k=0;k<REP_N;k++){ L.rep[k]=0.f; L.mean[k]=0.f; L.final_mean[k]=0.f; }
            L.final_pop=0;
            for(int k=0;k<D;k++) L.rep[k]=reps[c*D+k];
            L.birth_tick=tick; L.death_tick=-1; L.missing=0;
            L.parent=(best>=0)?lineages[best].id:-1;
            L.alive=true; L.pop=0; L.max_pop=0;
            cl[c]=(int)lineages.size(); lineages.push_back(L);
        }
    }
    // 一度見失っただけでは絶滅としない(標本の揺らぎを吸収する)
    for(int L=0;L<oldN;L++)
        if(lineages[L].alive && lineages[L].kind==kind && !claimed[L]){
            lineages[L].missing++;
            if(lineages[L].missing >= LINEAGE_GRACE){
                lineages[L].alive=false; lineages[L].death_tick=tick;
                // 最後に見えていた姿を残す。個体がいなくなっても参照できる
                for(int k=0;k<REP_N;k++) lineages[L].final_mean[k]=lineages[L].mean[k];
                lineages[L].final_pop=lineages[L].pop;
            }
        }

    for(auto& L:lineages) if(L.kind==kind){ L.pop=0; for(int k=0;k<REP_N;k++) L.mean[k]=0.f; }
    for(int i=0;i<total;i++){
        if(kind==0) plant_vec(plants[i],v); else animal_vec(animals[i],v);
        int best=-1; float bd=1e9f;
        for(int c=0;c<nc;c++){ float d=vdist(v,&reps[c*D],D); if(d<bd){bd=d;best=c;} }
        if(best<0) continue;
        Lineage& L = lineages[cl[best]];
        if(kind==0){
            const Plant& p=plants[i];
            plants[i].lin=L.id;
            L.mean[0]+=p.absorb;   L.mean[1]+=p.max_height; L.mean[2]+=p.shade_tol;
            L.mean[3]+=p.tough;    L.mean[4]+=p.clonal;     L.mean[5]+=p.disperse;
            L.mean[6]+=p.cold_tol; L.mean[7]+=p.height; L.mean[8]+=p.aquatic;
            L.mean[9]+=p.drought_tol;
        } else {
            const Animal& a=animals[i];
            animals[i].lin=L.id;
            L.mean[0]+=a.speed; L.mean[1]+=a.turn;  L.mean[2]+=a.reach;
            L.mean[3]+=a.size;  L.mean[4]+=a.diet;  L.mean[5]+=a.jaw;
            L.mean[6]+=a.cold_tol; L.mean[7]+=a.body;
            L.mean[8]+=a.w[0];
            { int nn=0; for(int k=0;k<N_NODE;k++) if(a.ngain[k]>NODE_ON) nn++;
              L.mean[9]+=(float)nn; } L.mean[10]+=a.aquatic;
            { float rs=0.f; int nv=0;
              for(int k=0;k<N_SENSOR;k++)
                  if(a.sen[k].range>=SENSOR_MIN){ rs+=a.sen[k].range; nv++; }
              L.mean[11]+=rs; }
            L.mean[12]+=a.arboreal;
            L.mean[13]+=a.perch;
        }
        L.pop++;
    }
    for(auto& L:lineages){
        if(L.kind!=kind) continue;
        if(L.pop>0) for(int k=0;k<REP_N;k++) L.mean[k]/=(float)L.pop;
        if(L.pop>L.max_pop) L.max_pop=L.pop;
    }
}

// 生きていれば現在の平均、絶滅していれば最後に見えていた姿を返す
inline const float* view_mean(const Lineage& L){
    return L.alive ? L.mean : L.final_mean;
}

// フィルタ
const int FILTER_N = 10;
const char* FILTER_NAME[FILTER_N] = {
    "全て","現存","絶滅","耐寒型","肉食","草食","高木","硬い",
    "水生","陸生"
};
bool filter_on[FILTER_N] = {true,false,false,false,false,false,false,false,false,false};
bool filter_single(const Lineage& L,int f);
// 有効なフィルタ全てを満たすものだけ通す(ANDで合成)
bool lineage_match_multi(const Lineage& L){
    bool any=false;
    for(int f=0;f<FILTER_N;f++){
        if(!filter_on[f]) continue;
        any=true;
        if(!filter_single(L,f)) return false;
    }
    return any ? true : true;
}
bool filter_single(const Lineage& L,int f){
    switch(f){
        case 0: return true;
        case 1: return L.alive;
        case 2: return !L.alive;
        case 3: return view_mean(L)[6] > 0.60f;
        case 4: return L.kind==1 && view_mean(L)[4] > 0.60f;
        case 5: return L.kind==1 && view_mean(L)[4] < 0.30f;
        case 6: return L.kind==0 && view_mean(L)[1] > 1.50f;
        case 7: return L.kind==0 && view_mean(L)[3] > 0.40f;
        case 8: return (L.kind==0? view_mean(L)[8] : view_mean(L)[10]) > 0.60f;
        case 9: return (L.kind==0? view_mean(L)[8] : view_mean(L)[10]) < 0.30f;
    }
    return true;
}

// 分離型ボックスブラー(累積和を使うので半径に関係なく速い)
void box_blur(std::vector<float>& f,int R,int passes){
    static std::vector<float> tmp;
    tmp.assign(GRID_W*GRID_H,0.f);
    float inv=1.f/(2*R+1);
    for(int p=0;p<passes;p++){
        for(int y=0;y<GRID_H;y++){                    // 東西はループ
            float sum=0.f;
            for(int dx=-R;dx<=R;dx++){
                int nx=((dx%GRID_W)+GRID_W)%GRID_W;
                sum+=f[idx(nx,y)];
            }
            for(int x=0;x<GRID_W;x++){
                tmp[idx(x,y)]=sum*inv;
                int o=((x-R)%GRID_W+GRID_W)%GRID_W;
                int i=((x+R+1)%GRID_W+GRID_W)%GRID_W;
                sum+=f[idx(i,y)]-f[idx(o,y)];
            }
        }
        for(int x=0;x<GRID_W;x++){                    // 南北は端で止める
            float sum=0.f;
            for(int dy=-R;dy<=R;dy++){
                int ny=std::max(0,std::min(GRID_H-1,dy));
                sum+=tmp[idx(x,ny)];
            }
            for(int y=0;y<GRID_H;y++){
                f[idx(x,y)]=sum*inv;
                int oy=std::max(0,std::min(GRID_H-1,y-R));
                int iy=std::max(0,std::min(GRID_H-1,y+R+1));
                sum+=tmp[idx(x,iy)]-tmp[idx(x,oy)];
            }
        }
    }
}

// 窪地を埋める。どの点からも必ず海まで水が流れるようにする(Priority-Flood法)
void fill_depressions(std::vector<float>& z, float sea_lv){
    struct Node { float v; int i; };
    struct Cmp { bool operator()(const Node&a,const Node&b) const { return a.v>b.v; } };
    std::priority_queue<Node,std::vector<Node>,Cmp> pq;
    std::vector<char> done(GRID_W*GRID_H,0);
    for(int i=0;i<GRID_W*GRID_H;i++)
        if(z[i]<sea_lv){ done[i]=1; pq.push({z[i],i}); }
    for(int x=0;x<GRID_W;x++)
        for(int y : {0,GRID_H-1}){
            int i=idx(x,y);
            if(!done[i]){ done[i]=1; pq.push({z[i],i}); }
        }
    while(!pq.empty()){
        Node n=pq.top(); pq.pop();
        int x=n.i%GRID_W, y=n.i/GRID_W;
        for(int d=0;d<4;d++){
            const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
            int nx=x+DX[d], ny=y+DY[d];
            if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
            if(ny<0||ny>=GRID_H) continue;
            int j=idx(nx,ny);
            if(done[j]) continue;
            if(z[j]<n.v+1e-5f) z[j]=n.v+1e-5f;
            done[j]=1; pq.push({z[j],j});
        }
    }
}

// 最急降下で下流を決め、上流から流量を積む
void route_flow(const std::vector<float>& zf, const std::vector<float>& rain,
                std::vector<int>& rec, std::vector<float>& acc,
                std::vector<int>& order){
    int N=GRID_W*GRID_H;
    order.resize(N);
    for(int i=0;i<N;i++) order[i]=i;
    std::sort(order.begin(),order.end(),
              [&](int a,int b){ return zf[a]>zf[b]; });   // 高い順
    rec.assign(N,-1);
    acc=rain;
    const int DX[8]={1,-1,0,0,1,1,-1,-1}, DY[8]={0,0,1,-1,1,-1,1,-1};
    for(int i=0;i<N;i++){
        int x=i%GRID_W, y=i/GRID_W;
        int best=-1; float bs=0.f;
        for(int d=0;d<8;d++){
            int nx=x+DX[d], ny=y+DY[d];
            if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
            if(ny<0||ny>=GRID_H) continue;
            int j=idx(nx,ny);
            float len=(d<4)?1.f:1.41421f;
            float s=(zf[i]-zf[j])/len;
            if(s>bs){ bs=s; best=j; }
        }
        rec[i]=best;
    }
    for(int k=0;k<N;k++){                 // 高い順に流量を下流へ渡す
        int i=order[k];
        if(rec[i]>=0) acc[rec[i]]+=acc[i];
    }
}

// h_max = 2.244 u/k。望む起伏から侵食係数を逆算する
// h_max = 2.244 u/k。最大の隆起から、望む最高峰になる侵食係数を逆算する
inline float erosion_k(){
    return 2.244f*(UPLIFT_CONT+UPLIFT_OROGEN)/std::max(1e-6f,TARGET_RELIEF);
}

void erode_step(std::vector<float>& z, const std::vector<float>& rain,
                const std::vector<float>& uplift, float sea_lv){
    int N=GRID_W*GRID_H;
    for(int i=0;i<N;i++) z[i]+=uplift[i]*EROSION_DT;   // 隆起は海面下でも続く
    std::vector<float> zf=z;
    fill_depressions(zf,sea_lv);
    static std::vector<int> rec, order;
    static std::vector<float> acc;
    route_flow(zf,rain,rec,acc,order);
    const float K=erosion_k();

    // 陰的スキーム。下流の新しい高さを先に求め、それを使って上流を解く。
    // 大きな時間刻みでも安定するため、収束が速い
    for(int k=N-1;k>=0;k--){          // order は高い順。逆に辿れば下流から
        int i=order[k];
        int r=rec[i];
        if(r<0) continue;             // 河口。海面に固定され、削られない
        if(z[i]<sea_lv) continue;
        int x1=i%GRID_W, y1=i/GRID_W, x2=r%GRID_W, y2=r/GRID_W;
        int dx=x1-x2; if(dx>GRID_W/2) dx-=GRID_W; if(dx<-GRID_W/2) dx+=GRID_W;
        int dy=y1-y2;
        float dist=std::sqrt((float)(dx*dx+dy*dy));
        // 楯状地は古い硬い岩盤。侵食が遅く、平坦な地形が長く残る
        float resist=1.f-CRATON_RESIST*craton[i];
        float c=K*resist*std::pow(acc[i],EROSION_M)/std::max(0.5f,dist);
        z[i]=(z[i]+EROSION_DT*c*z[r])/(1.f+c*EROSION_DT);
    }

    // 熱的侵食。流域面積が小さい場所に立つ、非現実的な尖峰を崩す
    for(int pass=0;pass<2;pass++){
        for(int k=0;k<N;k++){
            int i=order[k];
            if(z[i]<sea_lv) continue;
            int x=i%GRID_W, y=i/GRID_W;
            const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
            for(int d=0;d<4;d++){
                int nx=x+DX[d], ny=y+DY[d];
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                int j=idx(nx,ny);
                float diff=z[i]-z[j];
                if(diff<=MAX_SLOPE) continue;
                float move=(diff-MAX_SLOPE)*0.5f;   // 斜面上の物質が落ちる
                z[i]-=move; z[j]+=move;
            }
        }
    }

    // 地殻均衡による戻り。削られた場所は浮き上がる
    // (陰的解法に含まれないので、別に扱う)
    static std::vector<float> sm;
    sm=z;
    box_blur(sm,1,1);
    for(int i=0;i<N;i++) z[i]+=(sm[i]-z[i])*DIFFUSE_K;
}

// 海流。風に押された表層水が海岸で曲がり、西岸で強化される。
// 大陸東岸に暖流、西岸に寒流と湧昇という実際の配置が出る
void compute_currents(){
    std::fill(sst_anom.begin(),sst_anom.end(),0.f);
    std::fill(cur_u.begin(),cur_u.end(),0.f);
    std::fill(cur_v.begin(),cur_v.end(),0.f);
    // 流速場。風に従い、岸にぶつかれば岸沿いに向きを変える
    for(int y=0;y<GRID_H;y++){
        float latn=((GRID_H-1)/2.f-(float)y)/((GRID_H-1)/2.f);
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(!sea[i]) continue;
            float u=wind[y]/WIND_STRENGTH, v=0.f;
            // 東に陸があれば、ここは海盆の西岸。極向きに強化される
            int look=(int)WBC_RANGE, e_land=0, w_land=0;
            for(int k=1;k<=look;k++){
                if(!e_land && !sea[idx((x+k)%GRID_W,y)]) e_land=k;
                if(!w_land && !sea[idx((x-k+GRID_W)%GRID_W,y)]) w_land=k;
            }
            if(e_land){ float s=1.f-(float)e_land/look; v+= latn*1.8f*s; u*=0.5f; }
            if(w_land){ float s=1.f-(float)w_land/look; v-= latn*1.1f*s; }
            cur_u[i]=u; cur_v[i]=v;
        }
    }
    // 熱を流れに乗せて運ぶ。低緯度の暖かさが高緯度へ届く
    std::vector<float> heat(GRID_W*GRID_H,0.f), buf(GRID_W*GRID_H,0.f);
    for(int y=0;y<GRID_H;y++){
        float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
        for(int x=0;x<GRID_W;x++)
            if(sea[idx(x,y)]) heat[idx(x,y)]=1.f-lat;    // 低緯度が熱源
    }
    for(int it=0;it<CURRENT_ITER;it++){
        buf=heat;
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(!sea[i]) continue;
            int ux=(cur_u[i]>0)?(x-1+GRID_W)%GRID_W:(x+1)%GRID_W;
            int uy=y-((cur_v[i]>0)?1:-1);
            if(uy<0) uy=0; if(uy>=GRID_H) uy=GRID_H-1;
            float fu=std::min(0.7f,std::fabs(cur_u[i]))*CURRENT_ADV;
            float fv=std::min(0.7f,std::fabs(cur_v[i]))*CURRENT_ADV;
            int iu=idx(ux,y), iv=idx(x,uy);
            float src=0.f, wsum=0.f;
            if(sea[iu]){ src+=heat[iu]*fu; wsum+=fu; }
            if(sea[iv]){ src+=heat[iv]*fv; wsum+=fv; }
            buf[i]=heat[i]*(1.f-wsum)+src;
        }
        heat.swap(buf);
        box_blur(heat,1,1);
        for(int y=0;y<GRID_H;y++){        // 熱源と放熱で平衡へ引き戻す
            float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
            for(int x=0;x<GRID_W;x++){
                int i=idx(x,y);
                if(sea[i]) heat[i]+=((1.f-lat)-heat[i])*0.015f;
            }
        }
    }
    // 緯度ごとの平均からのずれが、海流による気温の偏差になる
    for(int y=0;y<GRID_H;y++){
        double m=0; int c=0;
        for(int x=0;x<GRID_W;x++) if(sea[idx(x,y)]){ m+=heat[idx(x,y)]; c++; }
        if(!c) continue;
        m/=c;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]) sst_anom[i]=(float)(heat[i]-m)*CURRENT_HEAT;
        }
    }
    box_blur(sst_anom,3,2);              // 沿岸の陸にも及ぶ
}

void init_terrain(){
    const float PI=3.14159265f;
    std::vector<float> psx(N_PLATES),psy(N_PLATES),pvx(N_PLATES),pvy(N_PLATES),
                       pthick(N_PLATES);
    // プレートの種を撒き、それぞれに地殻の厚さと運動を与える
    for(int k=0;k<N_PLATES;k++){
        psx[k]=frand(0.f,(float)GRID_W);
        psy[k]=frand(0.f,(float)GRID_H);
        plate_oce[k]=(dist01(rng)<OCEANIC_FRAC)?1:0;
        pthick[k]= plate_oce[k] ? CRUST_OCEAN
                                : CRUST_CONT+(dist01(rng)-0.5f)*CRUST_VARY;
        float a=frand(0.f,2.f*PI), sp=0.4f+frand(0.f,0.6f);
        pvx[k]=std::cos(a)*sp; pvy[k]=std::sin(a)*sp;
    }

    // 座標そのものをノイズでずらす。プレート境界が直線でなくなる
    std::vector<float> wpx(GRID_W*GRID_H), wpy(GRID_W*GRID_H);
    {
        std::vector<float> a1(GRID_W*GRID_H),a2(GRID_W*GRID_H),
                           b1f(GRID_W*GRID_H),b2f(GRID_W*GRID_H);
        for(int i=0;i<GRID_W*GRID_H;i++){
            a1[i]=frand(0.f,1.f); a2[i]=frand(0.f,1.f);
            b1f[i]=frand(0.f,1.f); b2f[i]=frand(0.f,1.f);
        }
        box_blur(a1,WARP_R,3);  box_blur(a2,WARP_R,3);
        box_blur(b1f,std::max(2,WARP_R/3),3); box_blur(b2f,std::max(2,WARP_R/3),3);
        auto nrm=[&](std::vector<float>& v){
            float lo=1e9f,hi=-1e9f;
            for(float t:v){ if(t<lo)lo=t; if(t>hi)hi=t; }
            for(auto& t:v) t=(t-lo)/std::max(1e-6f,hi-lo);
        };
        nrm(a1); nrm(a2); nrm(b1f); nrm(b2f);
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            wpx[i]=(float)x+((a1[i]-0.5f)*2.f+(b1f[i]-0.5f)*0.9f)*WARP_AMP;
            wpy[i]=(float)y+((a2[i]-0.5f)*2.f+(b2f[i]-0.5f)*0.9f)*WARP_AMP;
            if(wpy[i]<0.f) wpy[i]=0.f;
            if(wpy[i]>GRID_H-1.f) wpy[i]=GRID_H-1.f;
        }
    }

    // プレートの割り当て: ランダム塗りつぶし。
    // キューから無作為に選んで広げると、境界が直線でなくフラクタルになる
    {
        int N=GRID_W*GRID_H;
        std::fill(plate_id.begin(),plate_id.end(),-1);
        std::vector<int> queue;
        for(int k=0;k<N_PLATES;k++){
            int gx=(int)psx[k], gy=(int)psy[k];
            gx=std::max(0,std::min(GRID_W-1,gx));
            gy=std::max(0,std::min(GRID_H-1,gy));
            int i=idx(gx,gy);
            if(plate_id[i]>=0) continue;
            plate_id[i]=k; queue.push_back(i);
        }
        const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
        for(size_t out=0; out<queue.size(); out++){
            // 残りの中から無作為に一つ選び、先頭と入れ替える
            size_t pos=out+(size_t)(frand(0.f,(float)(queue.size()-out)-0.001f));
            int cur=queue[pos];
            queue[pos]=queue[out]; queue[out]=cur;
            int x=cur%GRID_W, y=cur/GRID_W;
            for(int d=0;d<4;d++){
                int nx=x+DX[d], ny=y+DY[d];
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                int j=idx(nx,ny);
                if(plate_id[j]>=0) continue;
                plate_id[j]=plate_id[cur];
                queue.push_back(j);
            }
        }
        for(int i=0;i<N;i++) if(plate_id[i]<0) plate_id[i]=0;

        // 境界を見つけ、そこでの相対運動を求める
        std::vector<float> conv_f(N,0.f);
        std::vector<int> bq;
        std::vector<float> bdist(N,1e9f);
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y); int a=plate_id[i];
            float nx=0.f, ny=0.f, cs=0.f; int cnt=0;
            for(int d=0;d<4;d++){
                int mx=x+DX[d], my=y+DY[d];
                if(mx<0) mx+=GRID_W; if(mx>=GRID_W) mx-=GRID_W;
                if(my<0||my>=GRID_H) continue;
                int j=idx(mx,my); int b=plate_id[j];
                if(b==a) continue;
                // 相手のプレートへ向かう向きが法線
                float ux=(float)DX[d], uy=(float)DY[d];
                float c=(pvx[a]-pvx[b])*ux+(pvy[a]-pvy[b])*uy;  // 正なら衝突
                nx+=ux; ny+=uy; cs+=c; cnt++;
            }
            if(cnt){ conv_f[i]=cs/cnt; bdist[i]=0.f; bq.push_back(i); }
        }
        // 境界から内陸へ、影響を減衰させながら広げる
        for(size_t h=0;h<bq.size();h++){
            int i=bq[h];
            if(bdist[i]>=BOUNDARY_W) continue;
            int x=i%GRID_W, y=i/GRID_W;
            for(int d=0;d<4;d++){
                int nx=x+DX[d], ny=y+DY[d];
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                int j=idx(nx,ny);
                if(bdist[j]<=bdist[i]+1.f) continue;
                bdist[j]=bdist[i]+1.f;
                conv_f[j]=conv_f[i];
                bq.push_back(j);
            }
        }
        for(int i=0;i<N;i++){
            float w=clamp01(1.f-bdist[i]/BOUNDARY_W); w*=w;
            bnd_type[i]=conv_f[i]*w;
        }
        // 境界ごとの地形。海溝・海嶺・島弧・地溝は、ここで作られる
        for(int i=0;i<N;i++){
            float bt=bnd_type[i];
            if(std::fabs(bt)<1e-4f){ tect_f[i]=0.f; continue; }
            float ac=std::fabs(bt);
            bool oce=plate_oce[plate_id[i]]!=0;
            // 隣接する相手のプレートを探す
            int x=i%GRID_W, y=i/GRID_W, other=-1;
            for(int d=0;d<4&&other<0;d++){
                const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
                for(int step=1;step<=(int)BOUNDARY_W+1&&other<0;step++){
                    int nx=x+DX[d]*step, ny=y+DY[d]*step;
                    if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                    if(ny<0||ny>=GRID_H) continue;
                    int p=plate_id[idx(nx,ny)];
                    if(p!=plate_id[i]) other=p;
                }
            }
            bool oce2 = (other>=0) ? (plate_oce[other]!=0) : oce;
            float t=0.f;
            if(bt>0.f){
                if(!oce && !oce2)      t=0.f;                   // 大陸衝突は地殻の厚さで表す
                else if(oce && !oce2)  t=-TRENCH_DEPTH*ac;      // 沈み込む側 → 海溝
                else if(!oce && oce2)  t=0.f;                   // 乗り上げる側も地殻で表す
                else {
                    // 海洋同士。速度の大きい方が沈み込み、他方に島弧が立つ
                    bool sub=(pvx[plate_id[i]]*pvx[plate_id[i]]
                             +pvy[plate_id[i]]*pvy[plate_id[i]])
                           > (other>=0? pvx[other]*pvx[other]+pvy[other]*pvy[other] : 0.f);
                    t = sub ? -TRENCH_DEPTH*ac : ISLAND_ARC*ac;
                }
            } else {
                if(oce) t= RIDGE_HEIGHT*ac;    // 海洋の発散 → 海嶺
                else    t=-RIFT_DEPTH*ac;      // 大陸の発散 → 地溝
            }
            tect_f[i]=t;
        }
    }

    // ホットスポット。マントル深部の熱の柱は動かないが、
    // その上をプレートが動くので、火山島が一列に並ぶ
    for(int h=0;h<N_HOTSPOT;h++){
        float hx=frand(0.f,(float)GRID_W), hy=frand(GRID_H*0.1f,GRID_H*0.9f);
        int seed=idx(std::max(0,std::min(GRID_W-1,(int)hx)),
                     std::max(0,std::min(GRID_H-1,(int)hy)));
        int pl=plate_id[seed];
        float dx=-pvx[pl], dy=-pvy[pl];          // プレートが去った方向へ島が残る
        float nl=std::sqrt(dx*dx+dy*dy)+1e-9f; dx/=nl; dy/=nl;
        for(int k=0;k<HOTSPOT_CHAIN;k++){
            float px=hx+dx*k*HOTSPOT_GAP, py=hy+dy*k*HOTSPOT_GAP;
            if(px<0.f) px+=GRID_W; if(px>=GRID_W) px-=GRID_W;
            if(py<0.f||py>=GRID_H) break;
            // 古い島ほど沈む。海山になり、やがて見えなくなる
            float age=(float)k/HOTSPOT_CHAIN;
            float amp=HOTSPOT_H*(1.f-age*0.75f);
            int r=(int)HOTSPOT_R;
            for(int ddy=-r;ddy<=r;ddy++) for(int ddx=-r;ddx<=r;ddx++){
                float d=std::sqrt((float)(ddx*ddx+ddy*ddy));
                if(d>HOTSPOT_R) continue;
                int gx=(int)px+ddx, gy=(int)py+ddy;
                if(gx<0) gx+=GRID_W; if(gx>=GRID_W) gx-=GRID_W;
                if(gy<0||gy>=GRID_H) continue;
                float f=1.f-d/HOTSPOT_R;
                hotspot[idx(gx,gy)]+=amp*f*f;
            }
        }
    }

    // --- 地殻の厚さ。標高はアイソスタシーで決まる ---
    for(int i=0;i<GRID_W*GRID_H;i++) crust[i]=pthick[plate_id[i]];

    // --- 大陸内部の構造 ---
    {
        int N=GRID_W*GRID_H;
        // 厚さの起伏。高原と内陸盆地を生む
        std::vector<float> f(N);
        for(int oc=0;oc<3;oc++){
            int R=CRUST_SCALE>>oc;
            if(R<2) R=2;
            for(int i=0;i<N;i++) f[i]=frand(0.f,1.f);
            box_blur(f,R,2);
            float lo=1e9f,hi=-1e9f;
            for(float v:f){ if(v<lo)lo=v; if(v>hi)hi=v; }
            float amp=CRUST_INTERNAL/(1<<oc);
            for(int i=0;i<N;i++)
                if(!plate_oce[plate_id[i]])
                    crust[i]+=((f[i]-lo)/std::max(1e-6f,hi-lo)-0.5f)*amp;
        }
        // 楯状地。古く厚く、そして侵食に抗う
        for(int i=0;i<N;i++) f[i]=frand(0.f,1.f);
        box_blur(f,CRATON_SCALE,3);
        {
            float lo=1e9f,hi=-1e9f;
            for(float v:f){ if(v<lo)lo=v; if(v>hi)hi=v; }
            for(int i=0;i<N;i++){
                if(plate_oce[plate_id[i]]){ craton[i]=0.f; continue; }
                float t=(f[i]-lo)/std::max(1e-6f,hi-lo);
                craton[i]=clamp01((t-0.55f)/0.35f);     // 上位のみが楯状地
                crust[i]+=CRATON_THICK*craton[i];
            }
        }
        // 失敗した地溝。大陸が割れかけて止まった跡が、内陸の細長い低地になる
        for(int k=0;k<N_RIFT;k++){
            int sx=-1,sy=-1;
            for(int t=0;t<200;t++){
                int gx=(int)frand(0.f,(float)GRID_W), gy=(int)frand(0.f,(float)GRID_H);
                if(plate_oce[plate_id[idx(gx,gy)]]) continue;
                sx=gx; sy=gy; break;
            }
            if(sx<0) continue;
            float ang=frand(0.f,6.2832f);
            float px=(float)sx, py=(float)sy;
            for(int s=0;s<(int)RIFT_LEN;s++){
                ang+=frand(-0.14f,0.14f);         // 蛇行させる
                px+=std::cos(ang); py+=std::sin(ang);
                if(px<0.f) px+=GRID_W; if(px>=GRID_W) px-=GRID_W;
                if(py<1.f||py>=GRID_H-1.f) break;
                float fade=1.f-(float)s/RIFT_LEN;  // 先へ行くほど浅くなる
                int r=(int)RIFT_WIDTH;
                for(int dy=-r;dy<=r;dy++) for(int dx=-r;dx<=r;dx++){
                    float d=std::sqrt((float)(dx*dx+dy*dy));
                    if(d>RIFT_WIDTH) continue;
                    int gx=(int)px+dx, gy=(int)py+dy;
                    if(gx<0) gx+=GRID_W; if(gx>=GRID_W) gx-=GRID_W;
                    if(gy<0||gy>=GRID_H) continue;
                    int j=idx(gx,gy);
                    if(plate_oce[plate_id[j]]) continue;
                    float w=(1.f-d/RIFT_WIDTH);
                    crust[j]-=RIFT_THIN*w*w*fade;
                    craton[j]*=(1.f-w*0.8f);      // 割れた場所は古い核ではない
                }
            }
        }
    }
    box_blur(crust,BASE_BLUR,2);
    // 収束境界で地殻が厚くなる。これが山脈になる
    for(int i=0;i<GRID_W*GRID_H;i++)
        if(bnd_type[i]>0.f && !plate_oce[plate_id[i]])
            crust[i]+=OROGEN_THICK*bnd_type[i];
    // アイソスタシー
    for(int i=0;i<GRID_W*GRID_H;i++)
        elev[i]=crust[i]*RHO_RATIO+hotspot[i]+tect_f[i];
    // 熱的沈降。海洋底は海嶺から離れるほど深い
    {
        std::vector<float> dr(GRID_W*GRID_H,1e9f);
        std::vector<int> q;
        for(int i=0;i<GRID_W*GRID_H;i++)
            if(plate_oce[plate_id[i]] && bnd_type[i]<-0.05f){ dr[i]=0.f; q.push_back(i); }
        for(size_t h=0;h<q.size();h++){
            int i=q[h]; int x=i%GRID_W, y=i/GRID_W;
            for(int d=0;d<4;d++){
                const int DX[4]={1,-1,0,0}, DY[4]={0,0,1,-1};
                int nx=x+DX[d], ny=y+DY[d];
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                int j=idx(nx,ny);
                if(dr[j]>dr[i]+1.f){ dr[j]=dr[i]+1.f; q.push_back(j); }
            }
        }
        for(int i=0;i<GRID_W*GRID_H;i++)
            if(plate_oce[plate_id[i]])
                elev[i]-=SUBSIDE_K*std::sqrt(std::min(dr[i],SUBSIDE_MAX));
        // 海嶺は海面を超えない。海洋地殻が陸として顔を出すのは例外的な場合だけ
        float ocean_top=CRUST_OCEAN*RHO_RATIO+RIDGE_HEIGHT*0.35f;
        for(int i=0;i<GRID_W*GRID_H;i++)
            if(plate_oce[plate_id[i]] && elev[i]>ocean_top+hotspot[i])
                elev[i]=ocean_top+hotspot[i];
    }
    // 地殻の厚さは km 単位で扱ったので、表示と侵食の尺度 0〜1 に写す
    {
        float lo=1e9f,hi=-1e9f;
        for(float v:elev){ if(v<lo)lo=v; if(v>hi)hi=v; }
        float sp=std::max(1e-6f,hi-lo);
        for(auto& v:elev) v=(v-lo)/sp;
    }

    // 細部。海底は平坦に保つため、海洋プレートでは極小にする
    {
        std::vector<float> det(GRID_W*GRID_H,0.f), layer(GRID_W*GRID_H);
        const int   RR[6]={40,20,10,5,2,1};
        const float RA[6]={0.42f,0.40f,0.38f,0.36f,0.34f,0.30f};
        for(int oc=0;oc<6;oc++){
            for(int i=0;i<GRID_W*GRID_H;i++) layer[i]=frand(0.f,1.f);
            box_blur(layer,RR[oc],3);
            float a=1e9f,b=-1e9f;
            for(float v:layer){ if(v<a)a=v; if(v>b)b=v; }
            for(int i=0;i<GRID_W*GRID_H;i++){
                float t=(layer[i]-a)/std::max(1e-6f,b-a);
                det[i]+=(t-0.5f)*RA[oc];
            }
        }
        // 標準偏差で揃える。オクターブ数を変えても振幅が崩れない
        double m=0; for(float v:det) m+=v; m/=GRID_W*GRID_H;
        double s2=0; for(float v:det){ double d=v-m; s2+=d*d; }
        float sd=(float)std::sqrt(s2/(GRID_W*GRID_H))+1e-6f;
        for(int i=0;i<GRID_W*GRID_H;i++)
            elev[i]+= det[i]/sd
                    * (plate_oce[plate_id[i]]?DETAIL_OCEAN:DETAIL_LAND);
    }
    if(FINAL_BLUR>0) box_blur(elev,FINAL_BLUR,1);

    // --- 川に刻ませる ---
    // 低い海面で谷を刻み、そのあと海面を上げて沈ませる。
    // 複雑な海岸線は、溺れた谷から生まれる
    {
        std::vector<float> s(elev);
        std::sort(s.begin(),s.end());
        float glacial=s[(size_t)(s.size()*GLACIAL_LEVEL)];
        // 降水の見込み。緯度と標高から粗く見積もる
        std::vector<float> rain(GRID_W*GRID_H,1.f);
        for(int y=0;y<GRID_H;y++){
            float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
            float band=0.35f+0.65f*(0.5f+0.5f*std::cos(lat*3.14159265f*3.f));
            for(int x=0;x<GRID_W;x++) rain[idx(x,y)]=band;
        }
        // 造山は続いている。衝突境界では隆起が侵食に抗う
        // 隆起は滑らかな勾配より、一定値の段を並べた方が自然な谷を生む
        std::vector<float> up(GRID_W*GRID_H,0.f);
        std::vector<float> band(GRID_W*GRID_H);
        for(int i=0;i<GRID_W*GRID_H;i++) band[i]=frand(0.f,1.f);
        box_blur(band,9,2);
        { float lo=1e9f,hi=-1e9f;
          for(float v:band){ if(v<lo)lo=v; if(v>hi)hi=v; }
          for(auto& v:band) v=std::floor((v-lo)/std::max(1e-6f,hi-lo)*4.f)/3.f; }
        // 大陸地殻は全体が隆起し続ける。衝突境界ではそこに造山が上乗せされる
        for(int i=0;i<GRID_W*GRID_H;i++){
            if(plate_oce[plate_id[i]]){
                // 海洋地殻は沈むが、火山島は供給が続くので持ちこたえる
                up[i]=-UPLIFT_CONT*0.15f+hotspot[i]*0.8f;
                continue;
            }
            float u=UPLIFT_CONT;
            if(bnd_type[i]>0.f) u+=UPLIFT_OROGEN*bnd_type[i];
            up[i]=u*(0.55f+0.45f*band[i]);
        }
        for(int it=0;it<EROSION_STEPS;it++)
            erode_step(elev,rain,up,glacial);
    }

    // 極を冷たい海にする
    for(int y=0;y<GRID_H;y++){
        float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
        float edge=(lat>0.90f)?(lat-0.90f)/0.10f:0.f;
        for(int x=0;x<GRID_W;x++) elev[idx(x,y)]-=edge*0.45f;
    }

    // 海陸比を強制し、陸の起伏を引き伸ばす
    {
        std::vector<float> s(elev);
        std::sort(s.begin(),s.end());
        float cut=s[(size_t)(s.size()*OCEAN_FRACTION)];
        float shift=SEA_LEVEL-cut;
        for(auto& v:elev){ v+=shift; if(v<0.02f)v=0.02f; if(v>1.f)v=1.f; }
        // 起伏はアイソスタシーと侵食が作る。人為的な引き伸ばしは不要
        float mx=SEA_LEVEL;
        for(float v:elev) if(v>mx) mx=v;
        float span=std::max(1e-6f,mx-SEA_LEVEL);
        for(auto& v:elev){
            if(v<=SEA_LEVEL) continue;
            v=SEA_LEVEL+(v-SEA_LEVEL)/span*(1.f-SEA_LEVEL);
        }
    }

    for(int i=0;i<GRID_W*GRID_H;i++) sea[i]=(elev[i]<SEA_LEVEL)?1:0;
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
        int xr=(x+1)%GRID_W, xl=(x-1+GRID_W)%GRID_W;
        int yu=std::max(0,y-1), yd=std::min(GRID_H-1,y+1);
        float gx=(elev[idx(xr,y)]-elev[idx(xl,y)])*0.5f;
        float gy=(elev[idx(x,yd)]-elev[idx(x,yu)])*0.5f;
        slope[idx(x,y)]=std::sqrt(gx*gx+gy*gy);
    }
    // 地形性降雨は数十キロの規模で起きる。細かい凹凸を均した標高を使う
    elev_smooth=elev;
    box_blur(elev_smooth,4,2);
    for(int y=0;y<GRID_H;y++){
        float latn=((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
        wind[y]=-std::cos(latn*3.14159265f*3.f)*WIND_STRENGTH;
    }
    // 陸からの距離。沿岸は河川の栄養が届き、外洋は届かない
    std::fill(shelf.begin(),shelf.end(),0.f);
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
        int i=idx(x,y);
        if(!sea[i]){ shelf[i]=1.f; continue; }
        float best=0.f;
        for(int dy=-SHELF_RANGE;dy<=SHELF_RANGE;dy++)
        for(int dx=-SHELF_RANGE;dx<=SHELF_RANGE;dx++){
            int nx=x+dx, ny=y+dy;
            if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
            if(ny<0||ny>=GRID_H) continue;
            if(sea[idx(nx,ny)]) continue;
            float d=std::sqrt((float)(dx*dx+dy*dy));
            float v=1.f-d/(SHELF_RANGE+1.f);
            if(v>best) best=v;
        }
        shelf[i]=best;
    }
    // 湧昇: 風が岸から沖へ吹き、表層水を運び去る場所に深層水が上がる
    for(int y=0;y<GRID_H;y++){
        int d=(wind[y]>0)?1:-1;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(!sea[i]) continue;
            int back=((x-d)%GRID_W+GRID_W)%GRID_W;   // 風上側
            if(!sea[idx(back,y)] && shelf[i]>0.4f) upwell[i]=1;
        }
    }
    // 陸からの近さ。沿岸の海は、大陸の気温に引かれる
    std::fill(coast_w.begin(),coast_w.end(),0.f);
    for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
        int i=idx(x,y);
        if(!sea[i]) continue;
        float best=0.f;
        for(int dy=-COAST_RANGE;dy<=COAST_RANGE;dy++)
        for(int dx=-COAST_RANGE;dx<=COAST_RANGE;dx++){
            int nx=x+dx, ny=y+dy;
            if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
            if(ny<0||ny>=GRID_H) continue;
            if(sea[idx(nx,ny)]) continue;
            float d=std::sqrt((float)(dx*dx+dy*dy));
            float v=1.f-d/(COAST_RANGE+1.f);
            if(v>best) best=v;
        }
        coast_w[i]=best*best;          // 岸に近いほど強く引かれる
    }
}

inline float altitude(int i){ return std::max(0.f, elev[i]-SEA_LEVEL); }

// 各セルから最も低い隣へ水を流し(D8)、上流面積を累積して河川を決める
void compute_rivers(){
    std::vector<int> order(GRID_W*GRID_H);
    for(int i=0;i<GRID_W*GRID_H;i++) order[i]=i;
    std::sort(order.begin(),order.end(),
              [](int a,int b){ return elev[a]>elev[b]; });
    std::fill(flow_acc.begin(),flow_acc.end(),1.f);
    for(int i : order){
        if(sea[i]) continue;
        int x=i%GRID_W, y=i/GRID_W;
        int best=-1; float lowest=elev[i];
        for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
            if(dx==0&&dy==0) continue;
            int nx=x+dx, ny=y+dy;
            if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
            if(ny<0||ny>=GRID_H) continue;
            int ni=idx(nx,ny);
            if(elev[ni]<lowest){ lowest=elev[ni]; best=ni; }
        }
        if(best>=0) flow_acc[best]+=flow_acc[i];   // 下流へ積み上げる
    }
    for(int i=0;i<GRID_W*GRID_H;i++)
        river[i] = (!sea[i] && flow_acc[i]>RIVER_THRESHOLD) ? 1 : 0;
}

// 降水: 緯度帯(赤道で湿り、亜熱帯で乾き、中緯度で再び湿る)に地形性降雨を重ねる
// 水蒸気を風で運ぶ。海で蒸発し、斜面を上るか空気が冷えると降る。
// 降水量は緯度ではなく、海からの距離と地形で決まる。
void compute_precip(){
    for(int y=0;y<GRID_H;y++){
        float u=wind[y];
        int d=(u>0)?1:-1;
        int x0=(d>0)?0:GRID_W-1;
        float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
        // ハドレー循環: 赤道で上昇し、亜熱帯で下降する
        float lift=0.5f+0.5f*std::cos(lat*3.14159265f*3.f);
        float m=0.f;
        // 東西はループしているので、数周させて定常状態に落ち着かせる
        for(int lap=0;lap<3;lap++){
            for(int st=0;st<GRID_W;st++){
                int x=((x0+d*st)%GRID_W+GRID_W)%GRID_W;
                int i=idx(x,y);
                if(sea[i]){
                    float sst=std::max(0.f,temper[i]);
                    m += EVAP_COEF*sst*(1.f-m/MOIST_CAP);   // 暖かい海ほど蒸発する
                } else {
                    m += LAND_RECYCLE*soil_water[i];        // 陸からの再循環
                }
                if(m>MOIST_CAP) m=MOIST_CAP;
                // 空気が保持できる量は気温で決まる。冷えれば溢れて降る
                float cap=MOIST_CAP*clamp01((temper[i]+20.f)/50.f);
                // 斜面は広い幅で測る。1セルの凹凸で雨を落とし切らせない
                int up=((x-d*ORO_SPAN)%GRID_W+GRID_W)%GRID_W;
                float rise=std::max(0.f,(elev_smooth[i]-elev_smooth[idx(up,y)])
                                        /(float)ORO_SPAN);
                float rate=RAIN_BASE*lift+OROGRAPHIC*rise;
                float rain=m*rate+std::max(0.f,m-cap)*0.5f;
                if(rain>m) rain=m;
                m-=rain;
                precip[i]=rain;
            }
        }
    }
    // 降水は数十キロの範囲に広がる。行ごとの独立計算による筋も、これで消える
    box_blur(precip,PRECIP_BLUR,2);
}

// 水収支: AW = [(1-RO)*AW + P - ET]+
// ケッペンの乾燥指数にもとづいて土壌水分の目標を決め、そこへ緩やかに近づける
void update_water(int frame){
    float ph=std::sin(2.f*3.14159265f*frame/SEASON_PERIOD);
    for(int y=0;y<GRID_H;y++){
        float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
        float sign=((float)y<(GRID_H-1)/2.f)?1.f:-1.f;
        float wet=1.f+0.45f*lat*ph*sign;          // 雨季と乾季
        if(wet<0.f) wet=0.f;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]){
                soil_water[i]=FIELD_CAP;
                // 有機物が沈降して栄養が抜ける。沿岸と湧昇域だけが補われる
                soil_n[i]-=soil_n[i]*N_SINK;
                soil_n[i]+=N_RUNOFF*shelf[i]*shelf[i];
                if(upwell[i]) soil_n[i]+=N_UPWELL;
                if(soil_n[i]<0.f) soil_n[i]=0.f;
                if(soil_n[i]>3.f) soil_n[i]=3.f;
                continue;
            }
            // 乾燥限界は気温に比例する。寒い土地は少ない雨でも乾かない
            float thr=ARID_COEF*std::max(0.5f,temper[i]+7.f);
            float A=precip[i]*wet/std::max(1e-9f,thr);
            aridity[i]=A;
            float target=A/(A+ARID_SOFT);          // 飽和させて差をなだらかにする
            if(river[i]) target=std::max(target,RIVER_SUPPLY);
            soil_water[i]+=(target-soil_water[i])*WATER_RELAX;
            soil_water[i]=clamp01(soil_water[i]);
        }
    }
}

// 年平均気温と降水から気候区分を決める(ケッペンの判定順に従う)
void classify_climate(){
    for(int y=0;y<GRID_H;y++){
        float lat=std::fabs((float)y-(GRID_H-1)/2.f)/((GRID_H-1)/2.f);
        float amp=SEASON_AMP*lat;
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]){ climate[i]=0; continue; }
            float Tann=TEMP_EQUATOR+(TEMP_POLE-TEMP_EQUATOR)*lat
                      -LAPSE_RATE*altitude(i);
            temp_annual[i]=Tann;
            float Tw=Tann+amp, Tc=Tann-amp;
            float thr=ARID_COEF*std::max(0.5f,Tann+7.f);
            float A=precip[i]/std::max(1e-9f,thr);
            if(Tw<0.f)       climate[i]=1;   // EF 氷雪
            else if(Tw<10.f) climate[i]=2;   // ET ツンドラ
            else if(A<0.5f)  climate[i]=3;   // BW 砂漠
            else if(A<1.0f)  climate[i]=4;   // BS ステップ
            else if(Tc>=18.f) climate[i]=(A<2.0f)?7:8;  // Aw サバナ / Af 熱帯雨林
            else if(Tc>-3.f)  climate[i]=6;  // C 温帯
            else              climate[i]=5;  // D 亜寒帯
        }
    }
}

void update_climate(int frame){
    const float PI=3.14159265f;
    float ph=std::sin(2.f*PI*frame/SEASON_PERIOD);
    // 太陽赤緯。夏の半球へ傾く
    sun_dec = AXIAL_TILT*PI/180.f*ph;
    float sd=std::sin(sun_dec), cd=std::cos(sun_dec);
    // 自転の位相。時間とともに昼の帯が西へ移る
    float rot=2.f*PI*(float)(frame%DAY_PERIOD)/DAY_PERIOD;
    double dsum=0;

    for(int y=0;y<GRID_H;y++){
        // 符号つき緯度。+1が北極、-1が南極
        float latn=((GRID_H-1)/2.f-(float)y)/((GRID_H-1)/2.f);
        float lat=std::fabs(latn);
        float phi=latn*PI/2.f;
        float sp=std::sin(phi), cp=std::cos(phi);
        float base=TEMP_EQUATOR+(TEMP_POLE-TEMP_EQUATOR)*lat;
        float base_sea=SST_EQUATOR+(SST_POLE-SST_EQUATOR)*lat;
        float seas_land=SEASON_AMP*latn*ph;
        float seas_sea =SEASON_AMP*SEA_SEASON*latn*ph;   // 海は熱容量で振れにくい
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            // 時角。経度と自転から決まる
            float H=2.f*PI*(float)x/GRID_W+rot;
            float sinh=sp*sd+cp*cd*std::cos(H);
            float sun=(sinh>0.f)?sinh:0.f;      // 地平線の下なら日射はゼロ
            solar[i]=sun*SOLAR_GAIN;
            float dl=clamp01(sinh*5.f);          // 薄明を挟んだ昼夜
            daylit[i]=dl;
            dsum+=dl;
            // 日較差。海は熱容量が大きくほとんど振れない
            float amp=sea[i]?DIURNAL_SEA:DIURNAL_LAND;
            if(sea[i]){
                float t_sea=base_sea+seas_sea+amp*(dl-0.5f)+sst_anom[i];
                // 沿岸の海は、隣接する大陸の気温に引かれる。
                // 冬に大陸から吹き出す寒気が、外洋より低緯度で海を凍らせる
                if(coast_w[i]>0.01f){
                    float t_land=base+seas_land+amp*(dl-0.5f)+sst_anom[i];
                    float pull=COAST_CHILL*coast_w[i];
                    if(t_land<t_sea) t_sea+=(t_land-t_sea)*pull;        // 冷気は強く効く
                    else             t_sea+=(t_land-t_sea)*pull*0.35f;  // 暖気は弱く
                }
                temper[i]=t_sea;
            } else {
                temper[i]=base+seas_land+amp*(dl-0.5f)+sst_anom[i]
                         -LAPSE_RATE*altitude(i);
            }
        }
    }
    // 海氷。凍れば光を遮り、陸の生物には橋になる
    // 氷は一日では張らないし、一日では融けない。厚さとして蓄える
    static std::vector<float> ice_amt(GRID_W*GRID_H,0.f);
    for(int i=0;i<GRID_W*GRID_H;i++){
        if(!sea[i]){ ice[i]=0; ice_amt[i]=0.f; continue; }
        float target = (temper[i]<FREEZE_TEMP) ? 1.f
                     : (temper[i]>THAW_TEMP)   ? 0.f
                     : (THAW_TEMP-temper[i])/(THAW_TEMP-FREEZE_TEMP);
        ice_amt[i]+=(target-ice_amt[i])*ICE_LAG;
        ice[i]=(ice_amt[i]>0.5f)?1:0;
    }
    day_light=(float)(dsum/(GRID_W*GRID_H));
}

void update_gases(){
    for(int y=0;y<GRID_H;y++){
        float u=wind[y];
        if(std::fabs(u)<1e-4f) continue;
        for(int pass=0;pass<2;pass++){
            std::vector<float>& g=(pass==0)?co2:o2;
            for(int x=0;x<GRID_W;x++){
                int up=(u>0)?(x-1+GRID_W)%GRID_W:(x+1)%GRID_W;
                gbuf[idx(x,y)]=g[idx(x,y)]+(g[idx(up,y)]-g[idx(x,y)])*std::fabs(u);
            }
            for(int x=0;x<GRID_W;x++) g[idx(x,y)]=gbuf[idx(x,y)];
        }
    }
    for(int pass=0;pass<2;pass++){
        std::vector<float>& g=(pass==0)?co2:o2;
        float base=(pass==0)?CO2_BASE:O2_BASE;
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            float s=0; int c=0;
            for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
                int nx=x+dx,ny=y+dy;
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                s+=g[idx(nx,ny)]; c++;
            }
            float mean=s/c; int i=idx(x,y);
            gbuf[i]=g[i]+(mean-g[i])*GAS_DIFFUSE+(base-g[i])*GAS_RELAX;
            if(gbuf[i]<0.f) gbuf[i]=0.f;
        }
        g.swap(gbuf);
    }
}

// 老化。寿命の後半に入ると、代謝が上がり、繁殖の力が落ちる
inline float senesce(int age,float lifespan){
    float L=LIFE_MIN+(LIFE_MAX-LIFE_MIN)*lifespan;
    float u=(float)age/L;
    if(u<SENESCE_START) return 0.f;
    float t=(u-SENESCE_START)/(1.f-SENESCE_START);
    return t*t*SENESCE_RATE;
}

inline float saturation(const Plant& p){ return 0.25f+0.75f*(1.f-p.shade_tol); }
inline float opt_temp_p(const Plant& p){ return 26.f-26.f*p.cold_tol; }
inline float kill_temp_p(const Plant& p){
    return 8.f-62.f*p.cold_tol-(p.dormant?DORM_HARDEN:0.f);
}
inline float opt_temp_a(const Animal& a){ return 26.f-26.f*a.cold_tol; }
inline float kill_temp_a(const Animal& a){ return  8.f-62.f*a.cold_tol; }

// 病原体は、その土地で最も多い免疫型を追いかける。
// 型が一致した宿主は病み、稀な型は逃れる(負の頻度依存選択)
inline float immune_match(const float* im,int ci,std::vector<float>* strain){
    float d=0.f;
    for(int k=0;k<N_IMMUNE;k++){ float t=im[k]-strain[k][ci]; d+=t*t; }
    d=std::sqrt(d/N_IMMUNE);
    float z=d/PATH_WIDTH;
    return std::exp(-z*z);          // 近いほど強く感染する
}

void update_pathogens(){
    static std::vector<float> sum_im[N_IMMUNE];
    static std::vector<float> host_n;
    for(int k=0;k<N_IMMUNE;k++) sum_im[k].assign(GRID_W*GRID_H,0.f);
    host_n.assign(GRID_W*GRID_H,0.f);

    // --- 植物 ---
    for(const auto& p:plants){
        if(!p.alive) continue;
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int ci=idx(gx,gy);
        for(int k=0;k<N_IMMUNE;k++) sum_im[k][ci]+=p.im[k];
        host_n[ci]+=1.f;
    }
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        float n=host_n[ci];
        if(n>0.f){
            for(int k=0;k<N_IMMUNE;k++){
                float mean=sum_im[k][ci]/n;
                pstrain[k][ci]+=(mean-pstrain[k][ci])*PATH_ADAPT;   // 多数派へ適応
            }
            pload[ci]+=PATH_GROW*std::min(6.f,n)-PATH_DECAY*pload[ci];
        } else pload[ci]-=PATH_DECAY*3.f*pload[ci];
        if(pload[ci]<0.f) pload[ci]=0.f;
        if(pload[ci]>PATH_LOAD_MAX) pload[ci]=PATH_LOAD_MAX;
    }

    // --- 動物 ---
    for(int k=0;k<N_IMMUNE;k++) sum_im[k].assign(GRID_W*GRID_H,0.f);
    host_n.assign(GRID_W*GRID_H,0.f);
    for(const auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int ci=idx(gx,gy);
        for(int k=0;k<N_IMMUNE;k++) sum_im[k][ci]+=a.im[k];
        host_n[ci]+=1.f;
    }
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        float n=host_n[ci];
        if(n>0.f){
            for(int k=0;k<N_IMMUNE;k++){
                float mean=sum_im[k][ci]/n;
                astrain[k][ci]+=(mean-astrain[k][ci])*PATH_ADAPT;
            }
            aload[ci]+=PATH_GROW*std::min(6.f,n)-PATH_DECAY*aload[ci];
        } else aload[ci]-=PATH_DECAY*3.f*aload[ci];
        if(aload[ci]<0.f) aload[ci]=0.f;
        if(aload[ci]>PATH_LOAD_MAX) aload[ci]=PATH_LOAD_MAX;
    }
}

// 匂いを放ち、拡散させ、風で流す。
// 視覚と違って時間に残り、遮蔽を回り込み、発信者の意図と無関係に漏れる
void update_odor(){
    // --- 放出 ---
    for(const auto& p:plants){
        if(!p.alive||p.seed_wait>0) continue;   // 種子は匂わない
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int i=idx(gx,gy);
        odor[0][i]+=ODOR_PLANT*(p.height+p.tough*1.5f+p.carnivory);
        if(pload[i]>0.5f) odor[0][i]+=ODOR_SICK*0.4f;
    }
    for(const auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int i=idx(gx,gy);
        float hide=(a.depth>0.5f)?0.25f:1.f;      // 潜れば漏れにくい
        odor[0][i]+=ODOR_ANIMAL*(a.body+a.toxin*1.2f)*hide;
        // 発情。繁殖できる状態の個体が出す
        if(a.sexual>0.15f && a.body>=a.size*MATURE_FRAC
           && a.energy>A_DIVIDE_TH*0.7f)
            odor[0][i]+=ODOR_ESTRUS*a.sexual;
        if(aload[i]>0.3f) odor[0][i]+=ODOR_SICK;
    }
    for(const auto& cp:corpses){
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        odor[1][idx(gx,gy)]+=ODOR_ROT*cp.energy;
    }
    for(int i=0;i<GRID_W*GRID_H;i++){
        if(fire[i]>0.05f)            odor[3][i]+=ODOR_SMOKE*fire[i];
        else if(burn_age[i]<BURN_WINDOW) odor[3][i]+=ODOR_SMOKE*0.05f;
    }

    // --- 移流・拡散・散逸 ---
    for(int k=0;k<N_ODOR;k++){
        std::vector<float>& g=odor[k];
        // 風で流れる
        for(int y=0;y<GRID_H;y++){
            float u=wind[y]*ODOR_WIND/WIND_STRENGTH;
            if(std::fabs(u)<1e-4f) continue;
            float f=std::min(0.85f,std::fabs(u));
            for(int x=0;x<GRID_W;x++){
                int up=(u>0)?(x-1+GRID_W)%GRID_W:(x+1)%GRID_W;
                odor_buf[idx(x,y)]=g[idx(x,y)]*(1.f-f)+g[idx(up,y)]*f;
            }
            for(int x=0;x<GRID_W;x++) g[idx(x,y)]=odor_buf[idx(x,y)];
        }
        // 拡散して散逸する
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            float s=0.f; int c=0;
            for(int dy=-1;dy<=1;dy++) for(int dx=-1;dx<=1;dx++){
                int nx=x+dx, ny=y+dy;
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                s+=g[idx(nx,ny)]; c++;
            }
            int i=idx(x,y);
            float v=g[i]+(s/c-g[i])*ODOR_DIFFUSE;
            v*=(1.f-ODOR_DECAY);
            odor_buf[i]=(v<1e-5f)?0.f:v;
        }
        g.swap(odor_buf);
        double t=0; for(int i=0;i<GRID_W*GRID_H;i++) t+=g[i];
        odor_total[k]=t/(GRID_W*GRID_H);
    }
}

// 崖崩れ。急斜面が滑り、その場の生物を押し流す
void update_slides(){
    // 崖崩れ。急斜面が滑り、その場の生物を押し流す
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        if(sea[ci]||slope[ci]<SLIDE_SLOPE) continue;
        if(dist01(rng)>=SLIDE_RATE*(slope[ci]/SLIDE_SLOPE)) continue;
        slide_events++;
        for(int t=pgrid.b(ci);t<pgrid.e(ci);t++){
            Plant& p=plants[pgrid.items[t]];
            if(!p.alive) continue;
            p.alive=false;
            corpses.push_back({p.x,p.y,p.height*0.8f,0});
        }
        soil_n[ci]*=0.35f;                 // 表土が流される
    }
}

// 分解者。死骸を食い、炭素と窒素を還す
void update_microbes(){
    std::vector<Microbe> babies;
    for(auto& m:microbes){
        if(!m.alive) continue;
        m.age++;
        int gx=(int)m.x, gy=(int)m.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H){ m.alive=false; continue; }
        int ci=idx(gx,gy);
        // 気温。低温では活動が落ちる
        float kill=-8.f-30.f*m.cold_tol;   // 胞子は凍結にも耐える
        if(temper[ci]<kill){ m.alive=false; d_micro++; continue; }
        float dT=(temper[ci]-(24.f-20.f*m.cold_tol))/22.f;
        float act=std::exp(-dT*dT);

        // 死骸を分解する。速いほど多く取れるが、酵素の維持費も上がる
        int cb=cgrid.b(ci), ce=cgrid.e(ci);
        if(cb<ce){
            int pk=cgrid.items[cb+(int)(frand(0.f,(float)(ce-cb)-0.001f))];
            float dec=corpses[pk].energy*m.rate*0.10f*act;
            corpses[pk].energy-=dec;
            m.energy+=dec*MICROBE_EFF;
            // 残りは大気と土へ
            co2[ci]   += dec*(1.f-MICROBE_EFF)*MICROBE_RETURN*CO2_PER_ENERGY;
            soil_n[ci]+= dec*(1.f-MICROBE_EFF)*MICROBE_RETURN*N_MINERALIZE*3.f;
            if(soil_n[ci]>3.f) soil_n[ci]=3.f;
        }
        m.energy -= MICROBE_COST*(0.5f+m.rate) + m.cold_tol*0.0015f;
        if(m.energy<=0.f){ m.alive=false; d_micro++; continue; }
        if(m.energy>=MICROBE_DIV){
            m.energy*=0.5f;
            Microbe c=m;
            c.energy=m.energy;
            float a=frand(0.f,6.2832f), d=frand(0.f,MICROBE_DISP);
            c.x=m.x+std::cos(a)*d; c.y=m.y+std::sin(a)*d;
            if(c.x<0.f) c.x+=GRID_W; if(c.x>=GRID_W) c.x-=GRID_W;
            if(c.y<0.f) c.y=0.f;     if(c.y>=GRID_H) c.y=GRID_H-1.f;
            c.rate    =mutate(m.rate,    0.05f,1.0f);
            c.cold_tol=mutate(m.cold_tol,0.0f,1.0f);
            c.id=next_microbe_id++; c.age=0;
            babies.push_back(c);
        }
    }
    for(auto& b:babies) microbes.push_back(b);
    microbes.erase(std::remove_if(microbes.begin(),microbes.end(),
        [](const Microbe& m){ return !m.alive; }),microbes.end());
    // 分解者が絶えた場合、死骸の山から再び湧く
    if(microbes.size()<50 && corpses.size()>1000){
        for(int t=0;t<200;t++){
            const Corpse& cp=corpses[(size_t)frand(0.f,(float)corpses.size()-0.001f)];
            Microbe m; m.x=cp.x; m.y=cp.y;
            m.energy=0.4f; m.alive=true; m.age=0;
            m.rate=frand(0.1f,0.8f); m.cold_tol=frand(0.f,0.6f);
            m.id=next_microbe_id++;
            microbes.push_back(m);
        }
    }
}

// 山火事。落雷で着火し、乾いた燃料の上を風下へ広がる。
// 埋もれた炭素を大気に戻し、焼け跡という新しい場所を開く。
void update_fire(int frame){
    // 酸素の分圧で判定する。窒素のない世界なので、他のガスとの比で見る
    float o2m=0.f, tot=0.f;
    for(int i=0;i<GRID_W*GRID_H;i++){ o2m+=o2[i]; tot+=o2[i]+co2[i]+N2_BASE; }
    float frac = (tot>1e-6f) ? o2m/tot : 0.f;    // 0..1
    // 地球の21%を基準に、火の窓(17%〜35%)へ写す
    float o2r = frac/0.21f;
    if(o2r<FIRE_O2_MIN){                      // 酸素が足りず火が広がらない
        std::fill(fire.begin(),fire.end(),0.f);
        return;
    }
    // 酸素が多いほど、湿った燃料まで燃えるようになる
    float wet_ok=FIRE_MOIST*std::min(2.0f,o2r/FIRE_O2_WET);

    for(int y=0;y<GRID_H;y++){
        for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(sea[i]||fuel[i]<FIRE_FUEL_MIN){ fire_buf[i]=0.f; continue; }
            if(soil_water[i]>wet_ok){ fire_buf[i]=fire[i]*(1.f-FIRE_DECAY*2.f); continue; }
            float f=fire[i];
            // 風下へ移る。隣の火勢を集める
            int up=(wind[y]>0)?(x-1+GRID_W)%GRID_W:(x+1)%GRID_W;
            float nb=fire[idx(up,y)]*1.6f;
            if(y>0)          nb+=fire[idx(x,y-1)];
            if(y<GRID_H-1)   nb+=fire[idx(x,y+1)];
            nb+=fire[idx((x+1)%GRID_W,y)]+fire[idx((x-1+GRID_W)%GRID_W,y)];
            float dryness=1.f-clamp01(soil_water[i]/std::max(0.05f,wet_ok));
            f+=nb*FIRE_SPREAD*dryness;
            if(dist01(rng)<FIRE_IGNITE*dryness*o2r) { f+=1.f; fire_events++; }
            f*=(1.f-FIRE_DECAY);
            if(f>3.f) f=3.f;
            fire_buf[i]=f;
        }
    }
    fire.swap(fire_buf);

    // 燃焼: 死骸を灰にし、炭素を大気へ戻す
    double area=0;
    for(auto& cp:corpses){
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        int i=idx(gx,gy);
        if(fire[i]<0.05f) continue;
        float burn=cp.energy*FIRE_BURN*std::min(1.f,fire[i]);
        cp.energy-=burn;
        co2[i]+=burn*CO2_PER_ENERGY;
        o2[i]=std::max(0.f,o2[i]-burn*FIRE_O2_USE);
    }
    for(int i=0;i<GRID_W*GRID_H;i++) if(fire[i]>=0.05f) area+=1;
    for(int i=0;i<GRID_W*GRID_H;i++){
        if(fire[i]>=0.05f) burn_age[i]=0;
        else if(burn_age[i]<999999) burn_age[i]++;
    }
    fire_area=area;
}

void update_plants(){
    pgrid.reset();
    for(const auto& p:plants){
        if(!p.alive||p.seed_wait>0) continue;   // 種子は光も遮らず、食われもしない
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        pgrid.count(idx(gx,gy));
    }
    pgrid.finalize();
    for(int i=0;i<(int)plants.size();i++){
        if(!plants[i].alive||plants[i].seed_wait>0) continue;
        int gx=(int)plants[i].x, gy=(int)plants[i].y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        pgrid.put(idx(gx,gy),i);
        fuel[idx(gx,gy)]+=plants[i].height*1.2f;
    }
    std::fill(canopy_h.begin(),canopy_h.end(),0.f);
    static std::vector<int> order;
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        int b0=pgrid.start[ci], e0=pgrid.start[ci+1];
        if(b0==e0) continue;
        const int* use=&pgrid.items[b0];
        int cnt=e0-b0;
        // このマスで最も高い株。つるはこれを支えにして登る
        float tallest=0.f;
        for(int t=b0;t<e0;t++)
            if(plants[pgrid.items[t]].height>tallest)
                tallest=plants[pgrid.items[t]].height;
        for(int t=b0;t<e0;t++){
            Plant& q=plants[pgrid.items[t]];
            float support=std::max(0.f,tallest-q.height);
            // つるは自立する幹を作らない。支えがなければ地面を這うしかない
            float self=q.height*(1.f-CLIMB_PENALTY*q.climb);
            q.eff_height=self+support*q.climb*CLIMB_GAIN;
            if(q.eff_height>tallest) tallest=q.eff_height;   // つるも林冠を作る
        }
        canopy_h[ci]=tallest;
        if(cnt>1){                               // 1株なら背比べは要らない
            order.assign(use,use+cnt);
            std::sort(order.begin(),order.end(),[](int a,int b){
                return plants[a].eff_height>plants[b].eff_height; });
            use=order.data();
        }
        float remaining=solar[ci];
        if(ice[ci]) remaining *= ICE_LIGHT;     // 氷が光を遮る
        for(int t=0;t<cnt;t++){
            if(remaining<=0.f) break;
            int pi=use[t];
            Plant& p=plants[pi];
            if(p.dormant) continue;
            // 浮く個体は水面近くに留まり、深さに関わらず光を受ける
            float att=1.f;
            if(sea[ci]) att=std::exp(-WATER_ATTEN*depthn(ci)*(1.f-p.buoyancy));
            float got=std::min(remaining*att,saturation(p));
            remaining -= got/std::max(0.05f,att);
            float dT=(temper[ci]-opt_temp_p(p))/TEMP_WIDTH_P;
            float tf=std::exp(-dT*dT);
            float gain=got*p.absorb*tf;
            // 水: 耐乾性が高いほど乏しい水で足りるが、装置に維持費がかかる
            float wneed = 0.55f*(1.f-0.75f*p.drought_tol);
            // 土壌が乏しければ、体に貯めた水で補う
            float avail_w = soil_water[ci]
                          + p.water_store/std::max(0.05f,wneed)*0.25f;
            float wf = clamp01(avail_w/std::max(0.02f,wneed));
            gain *= wf;
            if(soil_water[ci]<wneed && p.water_store>0.f){
                float draw=std::min(p.water_store,(wneed-soil_water[ci])*0.02f);
                p.water_store-=draw;
            }
            soil_water[ci] = std::max(0.f, soil_water[ci]-gain*0.010f);  // 蒸散
            // 炭素(リービッヒの最小律)
            float need=gain*CO2_PER_ENERGY;
            float avail=std::min(need,co2[ci]);
            if(need>1e-9f) gain*=avail/need;
            co2[ci]-=avail; o2[ci]+=avail*O2_PER_CO2;
            p.energy+=gain;
            dbg_gain+=gain; dbg_n++;
            dbg_light+=got; dbg_tf+=tf; dbg_wf+=wf;
            dbg_co2+=(need>1e-9f? avail/need : 1.f);
        }
    }
    // 寄生: 吸器を同じマスの最も大きな株に刺し、エネルギーを直接奪う。
    // 宿主の tough が防御として働く
    static int para_tick=0;
    if(++para_tick%4==0)
    for(int ci=0;ci<GRID_W*GRID_H;ci++){
        int b0=pgrid.start[ci], e0=pgrid.start[ci+1];
        if(e0-b0<2) continue;
        int host=-1; float best=-1.f;
        for(int t=b0;t<e0;t++){
            const Plant& q=plants[pgrid.items[t]];
            if(q.parasite>0.30f) continue;          // 寄生者は宿主にしない
            if(q.energy>best){ best=q.energy; host=pgrid.items[t]; }
        }
        if(host<0) continue;
        for(int t=b0;t<e0;t++){
            int pi=pgrid.items[t];
            if(pi==host) continue;
            Plant& q=plants[pi];
            if(q.parasite<0.10f) continue;
            float pierce=q.parasite/(q.parasite+plants[host].defense);
            float drain=plants[host].energy*PARASITE_DRAIN*pierce*4.f;
            plants[host].energy-=drain;
            q.energy+=drain*0.85f;
            q.n_store=std::min(1.5f,q.n_store+drain*0.02f);
        }
    }
    static int genet_tick=0;
    if(++genet_tick%10==0){          // 地下茎の融通は10tickに一度で足りる
        static std::unordered_map<int,std::pair<double,int>> genet;
        genet.clear();
        for(const auto& p:plants){
            if(!p.alive||p.clonal<0.05f) continue;
            auto& e=genet[p.genet]; e.first+=p.energy; e.second++;
        }
        for(auto& p:plants){
            if(!p.alive||p.clonal<0.05f) continue;
            auto it=genet.find(p.genet);
            if(it==genet.end()||it->second.second<2) continue;
            float mean=(float)(it->second.first/it->second.second);
            p.energy+=(mean-p.energy)*std::min(1.f,GENET_SHARE*10.f)*p.clonal;
        }
    }
    static int pphase=0;
    pphase=(pphase+1)%PLANT_STRIDE;
    std::vector<Plant> babies;
    int pidx=-1;
    for(auto& p:plants){
        pidx++;
        if(!p.alive) continue;
        // 土中で待つ種子。光合成もせず、食われもせず、ほぼ無費用で耐える
        if(p.seed_wait>0){
            // 待機中の種子は数が多い。8tickに一度だけ、まとめて処理する
            if((pidx+pphase)%SEED_STRIDE!=0) continue;
            p.seed_wait-=SEED_STRIDE;
            p.energy-=BANK_COST*SEED_STRIDE;
            if(dist01(rng)<BANK_PREDATION*SEED_STRIDE){ p.alive=false; continue; }
            if(p.energy<=0.f){ p.alive=false; continue; }
            if(p.seed_wait<=0){
                // 発芽の時だけ環境を見る
                int gx0=(int)p.x, gy0=(int)p.y;
                int ci0=(gx0>=0&&gx0<GRID_W&&gy0>=0&&gy0<GRID_H)?idx(gx0,gy0):0;
                if(soil_water[ci0]<0.25f){ p.seed_wait=SEED_STRIDE*6; continue; }
                p.seed_wait=0; p.height=0.05f;
            }
            continue;
        }
        if((pidx+pphase)%PLANT_STRIDE!=0) continue;
        p.age++;
        int gx=(int)p.x, gy=(int)p.y;
        int ci=(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H)?idx(gx,gy):0;
        // 休眠の切り替え。覚醒には余裕をもたせ、境界でばたつかせない
        if(p.dormant){ if(temper[ci]>p.dorm_temp+DORM_HYST) p.dormant=false; }
        else         { if(temper[ci]<p.dorm_temp)            p.dormant=true;  }

        if(temper[ci]<kill_temp_p(p)){
            p.alive=false;
            corpses.push_back({p.x,p.y,p.energy+p.height*2.0f,0});
            continue;
        }
        // 炎。厚い樹皮を持つ個体は耐える
        if(fire[ci]>0.05f){
            // 厚い樹皮でも完全には守れない。炎は幹を焼き、葉を奪う
            float dmg=fire[ci]*(1.f-p.fire_tol*0.65f);
            p.energy-=dmg*FIRE_DAMAGE;
            p.height=std::max(0.03f,p.height-dmg*0.60f);
            if(p.energy<=0.f){
                p.alive=false;
                corpses.push_back({p.x,p.y,p.height*0.25f,0});   // 大半は灰になる
                continue;
            }
        }
        // 窒素を取り込む。固定できる個体は大気から、できない個体は土壌から
        if(!p.dormant){
            p.n_store += N_FIX_RATE*p.n_fix*PLANT_STRIDE;
            // 同じマスの株数で分け合う。密集すれば一株あたりの取り分が減る
            int nn2=std::max(1,pgrid.n_at(ci));
            float take=std::min(soil_n[ci],0.010f*PLANT_STRIDE/(float)nn2);
            soil_n[ci]-=take; p.n_store+=take;
            if(p.n_store>1.5f) p.n_store=1.5f;
        }
        if(!p.dormant && daylit[ci]>0.3f && p.height<p.max_height && p.energy>0.5f){
            float g=std::min(GROW_RATE*PLANT_STRIDE,p.max_height-p.height);
            float nn=g*N_PER_GROWTH;
            if(p.n_store<nn) g*=p.n_store/std::max(1e-6f,nn);   // 窒素が律速する
            p.n_store=std::max(0.f,p.n_store-g*N_PER_GROWTH);
            p.height+=g; p.energy-=g*GROW_COST;
            // 建設は呼吸を伴う
            float need=g*GROW_COST*CO2_PER_ENERGY*0.5f;
            float got=std::min(need,o2[ci]);
            o2[ci]-=got; co2[ci]+=got;
        }
        float hm=habitat_p(p,ci);
        float upkeep = RESP_COEF*saturation(p)+HEIGHT_COST*p.height
                 +TOUGH_COST*p.defense+CLONAL_COST*p.clonal+COLD_COST_P*p.cold_tol
                     +HABITAT_COST_P*hm*hm+DROUGHT_COST*p.drought_tol;
        // 植物も呼吸する: 維持のために酸素を使い、二酸化炭素を出す
        {
            float need = upkeep*PLANT_RESP_CO2*CO2_PER_ENERGY
                       * (p.dormant?DORM_UPKEEP:1.f);
            float got  = std::min(need,o2[ci]);
            o2[ci]-=got; co2[ci]+=got;
            if(need>1e-9f && got<need*0.5f) upkeep *= 1.3f;   // 酸欠は割高になる
        }
        // 余裕があるとき、組織に水を貯める
        if(!p.dormant && p.water_store<p.store_cap && soil_water[ci]>0.5f){
            float f=std::min(STORE_FILL,p.store_cap-p.water_store);
            f=std::min(f,(soil_water[ci]-0.5f)*0.5f);
            if(f>0.f){ p.water_store+=f; soil_water[ci]-=f*0.5f; }
        }
        // 誘導防御。食害の匂いを感じて初めて武装する。常時ではコストが高い
        float threat=clamp01(odor[2][ci]/1.5f);
        p.defense += (p.tough*threat - p.defense)*INDUCE_RATE;
        p.defense -= p.defense*INDUCE_DECAY;
        if(p.defense<0.f) p.defense=0.f;
        if(p.defense>p.tough) p.defense=p.tough;
        // 浮く個体は海流に流される。根を張れないので固定できない
        if(sea[ci] && p.buoyancy>0.4f){
            p.x+=cur_u[ci]*BUOY_DRIFT*p.buoyancy;
            p.y+=cur_v[ci]*BUOY_DRIFT*p.buoyancy;
            if(p.x<0.f) p.x+=GRID_W; if(p.x>=GRID_W) p.x-=GRID_W;
            if(p.y<0.f) p.y=0.f;     if(p.y>=GRID_H) p.y=GRID_H-1.f;
        }
        upkeep += BUOY_COST*p.buoyancy
                + STORE_COST*p.store_cap + FIRE_TOL_COST*p.fire_tol
                + CLIMB_COST*p.climb + N_FIX_COST*p.n_fix
                + PARASITE_COST*p.parasite + CARNIVORY_COST*p.carnivory
                + FIRE_SEED_COST*p.fire_seed;
        // 暗所では維持費が落ちる。葉の光化学系を回す必要がないため
        upkeep *= (NIGHT_RESP+(1.f-NIGHT_RESP)*daylit[ci]);
        float sen=senesce(p.age,p.lifespan);
        upkeep *= (1.f+sen);
        upkeep += LIFESPAN_COST*p.lifespan;
        // 自己間引き。混み合うほど、光以外の要因でも消耗する
        {
            int nn=pgrid.n_at(ci);
            if(nn>4) upkeep += (nn-4)*0.020f;
        }
        if(p.dormant) upkeep*=DORM_UPKEEP;   // 葉を落として代謝を絞る
        // 発病。この土地の病原体と型が合ってしまった個体が病む
        if(pload[ci]>0.02f){
            float inf=pload[ci]*immune_match(p.im,ci,pstrain);
            upkeep+=inf*PATH_VIRULENCE;
            infect_p+=inf;
        }
        p.energy-=upkeep*PLANT_STRIDE;
        dbg_cost+=upkeep;      // 更新を飛ばしたぶんをまとめて払う
        if(p.energy<=0.f){
            p.alive=false;
            corpses.push_back({p.x,p.y,p.height*2.0f});
            continue;
        }
        if(!p.dormant && p.energy>=PLANT_DIVIDE_TH
           && (int)plants.size()+(int)babies.size() < PLANT_MAX){
            // 有性生殖を選ぶなら、まず相手を見つける必要がある
            const Plant* mate=nullptr;
            if(dist01(rng)<p.sexual){
                int r=(int)SEX_SEARCH;
                for(int dy=-r;dy<=r&&!mate;dy++) for(int dx=-r;dx<=r&&!mate;dx++){
                    int nx=gx+dx, ny=gy+dy;
                    if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                    if(ny<0||ny>=GRID_H) continue;
                    int ni=idx(nx,ny);
                    for(int t=pgrid.b(ni);t<pgrid.e(ni);t++){
                        const Plant& q=plants[pgrid.items[t]];
                        if(!q.alive||&q==&p||q.lin!=p.lin) continue;
                        mate=&q; break;
                    }
                }
                if(!mate) continue;      // 相手がいなければ、この機会を逃す
            }
            { double v[N_SG_P]; sg_vec_p(p,v);
              for(int k=0;k<N_SG_P;k++) sg_p_rep[k]+=v[k];
              sg_p_nr++; }
            int nb=1+(int)(p.brood*(MAX_BROOD-1)+0.5f);
            float share=1.f/(1.f+ (nb-1)*BROOD_TAX*nb );
            float invest=p.energy*(0.20f+0.55f*p.repro_alloc);
            p.energy-=invest;
            for(int bi=0;bi<nb;bi++){
            Plant c=p;
            c.energy=invest*share;
            // 組み換え: 各形質を、どちらかの親から受け取る
            if(mate){
                #define PICK(F) c.F = (dist01(rng)<0.5f)? p.F : mate->F
                PICK(absorb); PICK(max_height); PICK(shade_tol); PICK(tough);
                PICK(clonal); PICK(disperse);   PICK(cold_tol);  PICK(aquatic);
                PICK(drought_tol); PICK(dorm_temp); PICK(store_cap);
                PICK(fire_tol); PICK(climb); PICK(sexual); PICK(brood);
                PICK(n_fix); PICK(lifespan); PICK(parasite);
                PICK(carnivory); PICK(fire_seed); PICK(repro_alloc); PICK(seed_bank);
                PICK(buoyancy);
                #undef PICK
                // 免疫は座位ごとに独立に受け継ぐ。組み合わせで新しい型が生まれる
                for(int k=0;k<N_IMMUNE;k++)
                    c.im[k]=(dist01(rng)<0.5f)? p.im[k] : mate->im[k];
            }
            bool stay=(dist01(rng)<p.clonal);
            float d=stay?1.5f:(1.5f+SEED_FAR*p.disperse);
            float ang=frand(0.f,6.2832f);
            c.x=p.x+std::cos(ang)*d; c.y=p.y+std::sin(ang)*d;
            if(c.x<0.f) c.x+=GRID_W; if(c.x>=GRID_W) c.x-=GRID_W;
            if(c.y<0.f) c.y=0.f;     if(c.y>=GRID_H) c.y=GRID_H-1.f;
            int cgx=(int)c.x, cgy=(int)c.y;
            if(cgx<0||cgx>=GRID_W||cgy<0||cgy>=GRID_H) continue;   // この子は失われる
            // 火後発芽: 焼け跡にしか芽吹かない代わり、競争相手がいない
            if(p.fire_seed>0.35f && burn_age[idx(cgx,cgy)]>BURN_WINDOW) continue;
            c.height=0.05f; c.age=0; c.lin=p.lin;
            c.defense=0.f;
            // 土中で待つか、すぐ芽吹くか
            c.seed_wait = (dist01(rng)<p.seed_bank)
                        ? (int)(frand(40.f,(float)BANK_MAX)) : 0;
            c.id=next_plant_id++;
            c.genet=stay?p.genet:c.id;
            c.absorb    =mutate(p.absorb,    0.05f,1.0f);
            c.max_height=mutate(p.max_height,0.05f,5.0f);
            c.shade_tol =mutate(p.shade_tol, 0.0f,1.0f);
            c.tough     =mutate(p.tough,     0.0f,1.0f);
            c.clonal    =mutate(p.clonal,    0.0f,1.0f);
            c.disperse  =mutate(p.disperse,  0.0f,1.0f);
            c.cold_tol  =mutate(p.cold_tol,  0.0f,1.0f);
            c.aquatic     =mutate(p.aquatic,     0.0f,1.0f);
            c.drought_tol =mutate(p.drought_tol, 0.0f,1.0f);
            c.dorm_temp   =mutate(p.dorm_temp,  -20.0f,20.0f);
            c.store_cap   =mutate(p.store_cap,   0.0f,STORE_MAX);
            c.fire_tol    =mutate(p.fire_tol,    0.0f,1.0f);
            c.climb       =mutate(p.climb,       0.0f,1.0f);
            c.sexual      =mutate(p.sexual,      0.0f,1.0f);
            for(int k=0;k<N_IMMUNE;k++) c.im[k]=mutate(p.im[k],0.0f,1.0f);
            c.brood   =mutate(p.brood,   0.0f,1.0f);
            c.n_fix   =mutate(p.n_fix,   0.0f,1.0f);
            c.lifespan   =mutate(p.lifespan,   0.0f,1.0f);
            c.parasite   =mutate(p.parasite,   0.0f,1.0f);
            c.carnivory  =mutate(p.carnivory,  0.0f,1.0f);
            c.fire_seed  =mutate(p.fire_seed,  0.0f,1.0f);
            c.repro_alloc=mutate(p.repro_alloc,0.0f,1.0f);
            c.seed_bank  =mutate(p.seed_bank,  0.0f,1.0f);
            c.buoyancy   =mutate(p.buoyancy,   0.0f,1.0f);
            c.mut_rate   =mutate(p.mut_rate,   0.0f,1.0f);
            c.mother=p.id; c.father=mate?mate->id:-1; c.gen=p.gen+1;
            c.n_store =0.f;
            c.dormant     =p.dormant;
            c.water_store =0.f;
            babies.push_back(c);
            }
        }
    }
    for(auto& b:babies) plants.push_back(b);
    plants.erase(std::remove_if(plants.begin(),plants.end(),
        [](const Plant& p){ return !p.alive; }),plants.end());
}



void update_animals(){
    agrid.reset();
    for(const auto& a:animals){
        if(!a.alive) continue;
        int gx=(int)a.x, gy=(int)a.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        agrid.count(idx(gx,gy));
    }
    agrid.finalize();
    for(int i=0;i<(int)animals.size();i++){
        if(!animals[i].alive) continue;
        int gx=(int)animals[i].x, gy=(int)animals[i].y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        agrid.put(idx(gx,gy),i);
    }
    cgrid.reset();
    for(const auto& cp:corpses){
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        cgrid.count(idx(gx,gy));
    }
    cgrid.finalize();
    for(int i=0;i<(int)corpses.size();i++){
        int gx=(int)corpses[i].x, gy=(int)corpses[i].y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        cgrid.put(idx(gx,gy),i);
    }
    std::fill(fuel.begin(),fuel.end(),0.f);
    for(int ci=0;ci<GRID_W*GRID_H;ci++)
        corpse_density[ci]=(float)cgrid.n_at(ci);
    for(const auto& cp:corpses){
        if(cp.origin!=0) continue;            // 燃えるのは枯れた植物
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        fuel[idx(gx,gy)]+=cp.energy;
    }
    std::vector<Animal> babies;
    int n0=(int)animals.size();
    for(int i=0;i<n0;i++){
        Animal& a=animals[i];
        if(!a.alive) continue;
        a.age++;
        // それぞれの感覚器が、自分の向いた方角に、自分の射程だけ探る
        float in[N_IN];
        float sensor_span=0.f;
        for(int k=0;k<N_SENSOR;k++){
            in[k]=0.f;
            const Animal::Sensor& S=a.sen[k];
            if(S.range<SENSOR_MIN) continue;      // 器官として成立していない
            sensor_span+=S.range;
            float ang=a.dir+S.angle;
            float cs=std::cos(ang), sn=std::sin(ang);
            for(int s=1;s<=RAY_STEPS;s++){
                float d=S.range*s/(float)RAY_STEPS;
                int rx=(int)(a.x+cs*d), ry=(int)(a.y+sn*d);
                if(rx<0) rx+=GRID_W; if(rx>=GRID_W) rx-=GRID_W;
                if(ry<0||ry>=GRID_H) continue;
                int ri=idx(rx,ry);
                int band=(int)(S.tune*N_TUNE); if(band>N_TUNE-1) band=N_TUNE-1;
                float sig=0.f;
                switch(band){
                case 0: if(!pgrid.empty_at(ri)) sig=1.f; break;      // 植物
                case 1:                                              // 自分より小さい動物
                    { int e=std::min(agrid.e(ri),agrid.b(ri)+4);     // 先頭数体だけ見る
                      for(int t=agrid.b(ri);t<e;t++){
                        int j=agrid.items[t];
                        if(j!=i && animals[j].body < a.body*0.9f){ sig=1.f; break; }
                      } } break;
                case 2:                                              // 自分より大きい動物
                    { int e=std::min(agrid.e(ri),agrid.b(ri)+4);
                      for(int t=agrid.b(ri);t<e;t++){
                        int j=agrid.items[t];
                        if(j!=i && animals[j].body > a.body*1.1f){ sig=1.f; break; }
                      } } break;
                case 3: if(!cgrid.empty_at(ri)) sig=1.f; break;      // 死骸
                case 4: sig=std::min(1.f,fire[ri]); break;           // 炎
                case 5: sig=clamp01(canopy_h[ri]/3.f); break;        // 林冠の高さ
                case 6: sig=clamp01(aload[ri]/PATH_LOAD_MAX); break; // 病原体
                case 11:                                             // 同種。群れの前提
                    { int c=0, e=std::min(agrid.e(ri),agrid.b(ri)+6);
                      for(int t=agrid.b(ri);t<e;t++){
                          int j=agrid.items[t];
                          if(j!=i && animals[j].lin==a.lin) c++;
                      }
                      sig=std::min(1.f,c*0.45f); }
                    break;
                case 12: sig=sea[ri]?1.f:0.f; break;                 // 水際
                // 嗅覚。濃さそのものを読む。遮蔽を越え、時間に残った跡も追える
                default: sig=clamp01(odor[band-7][ri]/ODOR_SCALE); break;
                }
                if(sig>0.02f){ in[k]=sig*(1.f-d/S.range); break; }
            }
        }
        int sci=idx(std::max(0,std::min(GRID_W-1,(int)a.x)),
                    std::max(0,std::min(GRID_H-1,(int)a.y)));
        in[N_SENSOR+0]=std::min(1.f,a.energy/A_DIVIDE_TH);
        in[N_SENSOR+1]=clamp01((temper[sci]+30.f)/70.f);
        in[N_SENSOR+2]=soil_water[sci];
        in[N_SENSOR+3]=std::min(1.f,a.body/3.f);   // 自分の体格。戦うか逃げるかを分ける
        in[N_SENSOR+4]=1.f;                        // バイアス
        // 再帰入力。tanh を通しているので本来は有界だが、念のため縛る
        for(int k=0;k<N_REC;k++){
            float r=a.rec[k];
            if(!(r>-1.5f&&r<1.5f)) r=0.f;   // NaN もここで潰れる
            in[N_SENSOR+5+k]=r;
        }
        // 神経網。働いているニューロンだけを通す。
        // すべて無効なら入力から出力への直結だけが残り、単層として振る舞う
        float bus[MAX_NODE]={0}, nxt[MAX_NODE];
        int   n_node=0, n_syn=0;
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
        float out[N_OUT]={0};
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
        float o0=out[0], o1=out[1];
        float want_eat   =out[2];   // 採食するか
        float want_attack=out[3];   // 襲うか
        for(int k=0;k<N_REC;k++) a.rec[k]=out[4+k];   // 次のtickへ持ち越す
        // 神経組織は高くつく。ニューロンと結線の数に比例する
        float brain_use=n_node*(NODE_COST/SYN_COST)+n_syn;
        float turn_mult=(o0+1.f)*0.5f, thrust_mult=(o1+1.f)*0.5f;

        if(a.body<a.size && a.energy>1.0f){
            float g=std::min(BODY_GROW_RATE,a.size-a.body);
            a.body+=g; a.energy-=g*BODY_GROW_COST;
        }
        int hci=idx(std::max(0,std::min(GRID_W-1,(int)a.x)),
                    std::max(0,std::min(GRID_H-1,(int)a.y)));

        // 枝が支えられる高さ。重い体は高くまで登れない
        float bearable = canopy_h[hci]*BRANCH_SUPPORT/std::max(0.15f,a.body);
        float reachable = std::min(canopy_h[hci],bearable);

        // 登り降り。登攀が巧いほど安く登れる
        float want = reachable*a.arboreal;
        if(want>a.perch+0.02f){
            float up=std::min(want-a.perch,0.08f+0.25f*a.arboreal);
            a.perch+=up;
            a.energy-=up*CLIMB_UP_COST/(0.25f+a.arboreal);
        } else if(a.perch>reachable){
            // 支えが失われた。滑り落ちる
            float drop=a.perch-reachable;
            a.perch=reachable;
            if(drop>0.5f) a.energy-=drop*FALL_DAMAGE*(1.f-a.arboreal*0.7f);
        } else if(a.perch>want){
            a.perch=std::max(want,a.perch-0.15f);    // 降りるのは安い
        }
        if(a.perch<0.f) a.perch=0.f;

        float terr;
        if(a.perch>0.25f){
            // 樹上。地面の起伏に縛られない。登攀が巧いほど速く渡れる
            terr = 0.35f+0.75f*a.arboreal;
        } else if(sea[hci]) terr = (0.25f + WATER_SPEED*a.aquatic);   // 泳ぐ
        else         terr = (1.f - 0.75f*a.aquatic)                    // 陸を這う
                          / (1.f + SLOPE_COEF*slope[hci]);             // 急坂は遅い
        if(terr<0.05f) terr=0.05f;

        // 酸素の分圧。窒素がないので、他のガスとの比で見る
        float o2f = o2[hci]/std::max(1e-6f,o2[hci]+co2[hci]+N2_BASE);
        if(sea[hci]) o2f *= DISSOLVED_O2;              // 水に溶ける酸素は僅か
        o2f *= std::exp(-altitude(hci)*ALT_O2_DECAY);  // 高所は空気が薄い
        float o2rel = o2f/0.21f;                       // 地球の21%を1.0とする
        // 呼吸効率が高いほど、薄い酸素でも動ける
        float thr_half = 0.55f*(1.f-0.80f*a.resp);
        float o2fit = std::min(1.f, o2rel*(1.f+thr_half)/(o2rel+thr_half));
        a.aerobic = o2fit;

        a.dir+=frand(-a.turn*turn_mult,a.turn*turn_mult);
        float sp=a.speed*A_MAX_SPEED*thrust_mult*terr*o2fit;
        float nx=a.x+std::cos(a.dir)*sp, ny=a.y+std::sin(a.dir)*sp;
        if(nx<0.f) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
        if(ny<0.f){ ny=0.f; a.dir=-a.dir; }
        if(ny>=GRID_H){ ny=GRID_H-1.f; a.dir=-a.dir; }
        int ngx=(int)nx, ngy=(int)ny;
        bool blocked=false;
        if(ngx>=0&&ngx<GRID_W&&ngy>=0&&ngy<GRID_H){
            int ni=idx(ngx,ngy);
            // 氷の上は陸として歩ける。冬に大陸が繋がる
            bool passable_land = !sea[ni] || ice[ni];
            if(!passable_land && a.aquatic<0.20f) blocked=true;   // 泳げない
            if( passable_land && a.aquatic>0.85f && !sea[ni]) blocked=true;
            if(sea[ni] && ice[ni] && a.aquatic>0.85f) blocked=true;  // 氷に蓋をされた
            // 樹上を渡る: 隣の林冠が近い高さでないと渡れず、降りるしかない
            if(a.perch>0.25f && canopy_h[ni] < a.perch-PERCH_GAP){
                a.perch=std::max(0.f,canopy_h[ni]);         // 木が途切れた。降りる
                a.energy-=0.04f*(1.f-a.arboreal*0.6f);
            }
        }
        if(blocked) a.dir+=3.14159f;
        else { a.x=nx; a.y=ny; }
        int gx=(int)a.x, gy=(int)a.y;
        int ci=(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H)?idx(gx,gy):0;
        // 潜り具合。地中は気温が和らぐ
        a.depth += (a.burrow - a.depth)*0.15f;
        float felt = temper[ci];
        if(a.depth>0.05f){
            float year_mean=TEMP_EQUATOR
                +(TEMP_POLE-TEMP_EQUATOR)*std::fabs((float)((int)a.y)-(GRID_H-1)/2.f)
                 /((GRID_H-1)/2.f) - LAPSE_RATE*altitude(ci);
            felt = temper[ci]+(year_mean-temper[ci])*BURROW_BUFFER*a.depth;
        }
        if(felt<kill_temp_a(a)){
            a.alive=false; d_cold++;
            corpses.push_back({a.x,a.y,a.body*CORPSE_RETURN+a.energy,1});
            continue;
        }
        float dT=(felt-opt_temp_a(a))/25.f;
        float act=std::exp(-dT*dT)*o2fit;   // 酸素が薄ければ動きも鈍る
        // 昼夜。夜行性なら夜に本領を発揮する
        float night=1.f-daylit[ci];      // その土地が今、昼か夜か
        act *= (0.35f + 0.65f*(a.nocturnal*night + (1.f-a.nocturnal)*daylit[ci]));
        // 潜っていれば地表のものを食えない
        if(a.depth>0.5f) act*=0.15f;

        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H){
            // 三つの消化能力。合計に上限があるため、全部は伸ばせない
            float cap_sum = a.diet + a.detritus*0.5f;   // 分解は草食と干渉しにくい
            if(cap_sum > 1.f) cap_sum = 1.f;
            float plant_dig = 1.f - cap_sum*cap_sum;
            float detr_dig  = a.detritus;
            float meat_dig  = a.diet;
            if(plant_dig>0.01f && want_eat>ACT_THRESH){
                // 摂食量は体格に比例する。大きい体は多く処理できる
                int max_bites=1+(int)(a.body*2.0f);
                if(max_bites>6) max_bites=6;
                int bites=0;
                for(int t=pgrid.b(ci);t<pgrid.e(ci)&&bites<max_bites;t++){
                    int pi=pgrid.items[t];
                    Plant& p=plants[pi];
                    if(!p.alive||p.eff_height>a.perch+a.reach) continue;
                    float chew=a.jaw/(a.jaw+p.defense*2.0f);
                    float bite=p.energy*A_BITE*chew*act;
                    p.energy-=bite;
                    // 齧られれば葉を失う。光合成能力が落ち、回復に時間がかかる。
                    // これがあって初めて、その場が枯れ、探す必要が生まれる
                    p.height=std::max(0.03f,p.height-bite*GRAZE_DEFOLIATE);
                    a.energy+=bite*(1.f-A_EAT_LOSS)*plant_dig;
                    odor[2][ci]+=ODOR_WOUND*bite;   // 傷ついた植物が匂いを放つ
                    a.eaten++; bites++;
                    if(p.energy<=0.f) p.alive=false;
                }
            }
            // 死骸を漁る: 植物由来は分解能力で、動物由来は肉食能力で処理する
            // 死骸を漁る。分解者は栄養価の低い落葉を大量に処理する
            {
                int cb=cgrid.b(ci), ce=cgrid.e(ci);
                if(cb<ce && a.perch<0.5f && want_eat>ACT_THRESH){
                    int nd=ce-cb;
                    int tries=1+(int)(3.f*detr_dig);        // 能力が高いほど多く処理する
                    for(int t=0;t<tries;t++){
                        int pk=cgrid.items[cb+(int)(frand(0.f,(float)nd-0.001f))];
                        bool is_plant=(corpses[pk].origin==0);
                        float cap =is_plant?detr_dig:meat_dig;
                        float loss=is_plant?DETRITUS_LOSS:MEAT_LOSS;
                        if(cap<=0.005f) continue;
                        float bite=corpses[pk].energy*SCAV_EFF*cap;
                        corpses[pk].energy-=bite;
                        a.energy+=bite*(1.f-loss); a.eaten++;
                    }
                }
            }
            // 生きた動物を襲う
            if(meat_dig>0.02f && want_attack>ACT_THRESH){
                int bx=(int)a.x, by=(int)a.y, pk=-1;
                for(int dy=-1;dy<=1&&pk<0;dy++) for(int dx=-1;dx<=1&&pk<0;dx++){
                    int nx=bx+dx, ny=by+dy;
                    if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                    if(ny<0||ny>=GRID_H) continue;
                    int ni=idx(nx,ny);
                    for(int t=agrid.b(ni);t<agrid.e(ni);t++){
                        int j=agrid.items[t];
                        if(j==i||!animals[j].alive) continue;
                        if(a.body<=animals[j].body*PRED_RATIO) continue;
                        if(std::fabs(a.perch-animals[j].perch)>=0.6f) continue;
                        if(animals[j].depth>0.5f) continue;      // 穴に逃げ込まれた
                        pk=j; break;
                    }
                }
                if(pk>=0){
                    Animal& v=animals[pk];
                    // わずかな肉食性でも噛めるようにする。線形だと谷を渡れない
                    float bite=v.energy*PRED_EFF*std::sqrt(meat_dig)*act;
                    v.energy-=bite;
                    // 毒。解毒できなければ、食っても損をする
                    float tox=v.toxin*(1.f-a.tox_resist);
                    a.energy+=bite*(1.f-MEAT_LOSS)*(1.f-tox);
                    life_note(a.id,g_frame,"捕食",bite);
                    a.energy-=v.toxin*TOXIN_HARM*(1.f-a.tox_resist);
                    a.eaten++;
                    if(v.energy<=0.f){
                        v.alive=false; d_preyed++;
                        corpses.push_back({v.x,v.y,v.body*CORPSE_RETURN+v.energy,1});
                    }
                }
            }
        }
        float hma=habitat_a(a,ci);
        // 食虫植物に捕らえられる。土壌の窒素が乏しい場所ほど罠が育つ
        if(a.depth<0.2f){
            for(int t=pgrid.b(ci);t<pgrid.e(ci);t++){
                Plant& tp=plants[pgrid.items[t]];
                if(!tp.alive||tp.carnivory<0.12f) continue;
                float catch_p=tp.carnivory*CARNIVORY_CATCH/(1.f+a.body*1.5f);
                if(dist01(rng)>=catch_p) continue;
                float taken=a.energy*0.55f+a.fat*0.4f;
                a.energy-=taken; a.fat*=0.6f;
                tp.energy  += taken*0.70f;                 // 直接エネルギーになる
                tp.n_store  = std::min(1.5f,tp.n_store+taken*CARNIVORY_N);
                a.eaten--;
                break;
            }
            if(a.energy<=0.f){
                a.alive=false; d_preyed++;
                corpses.push_back({a.x,a.y,a.body*CORPSE_RETURN,1});
                continue;
            }
        }
        float metab=A_METAB+sp*A_MOVE_COST+a.reach*A_REACH_COST
                   +a.body*A_SIZE_COST+a.jaw*JAW_COST+COLD_COST_A*a.cold_tol
                   +HABITAT_COST_A*hma*hma
                   +a.detritus*DETRITUS_COST+sensor_span*SENSOR_COST
                   +brain_use*SYN_COST
                   +a.sexual*SEX_COST+a.arboreal*ARBOREAL_COST
                   +(1.f-a.mut_rate)*MUT_RATE_COST   // 精密な修復は高くつく
                   +a.toxin*TOXIN_COST+a.tox_resist*TOXRES_COST
                   +a.burrow*BURROW_COST+a.nocturnal*NOCTURNAL_COST
                   +a.fat*FAT_COST;
        float sen=senesce(a.age,a.lifespan);
        metab *= (1.f+sen);
        metab += LIFESPAN_COST*a.lifespan;
        // 発病
        if(aload[ci]>0.02f){
            float inf=aload[ci]*immune_match(a.im,ci,astrain);
            metab+=inf*PATH_VIRULENCE;
            infect_a+=inf;
        }
        metab+=0.f
                   +HABITAT_COST_A*hma*hma;
        // 呼吸: 酸素の分圧で効率が決まる。水中と高所は薄い
        float o2_need = metab*O2_PER_METAB*(1.f-RESP_EFFICIENCY*a.resp);
        float avail   = std::min(o2_need,o2[ci]);
        o2[ci]-=avail; co2[ci]+=avail;
        float eff=(o2_need>1e-9f)?avail/o2_need:1.f;
        if(eff<0.15f) eff=0.15f;
        a.energy -= metab/eff + a.resp*RESP_COST;
        // 余剰を脂肪に回す。飢えれば取り崩す
        if(a.energy>A_DIVIDE_TH*0.6f && a.fat<a.store_fat*FAT_MAX){
            float put=std::min((a.energy-A_DIVIDE_TH*0.6f)*FAT_RATE,
                               a.store_fat*FAT_MAX-a.fat);
            a.energy-=put; a.fat+=put*0.85f;
        } else if(a.energy<A_DIVIDE_TH*0.25f && a.fat>0.f){
            float take=std::min(a.fat,0.05f);
            a.fat-=take; a.energy+=take;
        }
        if(a.energy<=0.f){
            a.alive=false; d_starve++;
            life_note(a.id,g_frame,"餓死",(float)a.age);
            corpses.push_back({a.x,a.y,a.body*CORPSE_RETURN+a.energy,1});
            continue;
        }
        if((int)animals.size()+(int)babies.size() < ANIMAL_MAX
           && a.energy>=A_DIVIDE_TH*(1.4f-0.6f*a.repro_alloc)
           && a.body>=a.size*MATURE_FRAC){
            const Animal* mate=nullptr;
            if(dist01(rng)<a.sexual){
                int r=(int)SEX_SEARCH;
                int bx=(int)a.x, by=(int)a.y;
                for(int dy=-r;dy<=r&&!mate;dy++) for(int dx=-r;dx<=r&&!mate;dx++){
                    int nx=bx+dx, ny=by+dy;
                    if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                    if(ny<0||ny>=GRID_H) continue;
                    int ni=idx(nx,ny);
                    for(int t=agrid.b(ni);t<agrid.e(ni);t++){
                        const Animal& q=animals[agrid.items[t]];
                        if(!q.alive||&q==&a||q.lin!=a.lin) continue;
                        if(q.body<q.size*MATURE_FRAC) continue;
                        mate=&q; break;
                    }
                }
                if(!mate) continue;
            }
            // 何匹産むか。多く産めば一匹あたりの持ち分が減る
            { double v[N_SG_A]; sg_vec_a(a,v);
              for(int k=0;k<N_SG_A;k++) sg_a_rep[k]+=v[k];
              sg_a_nr++; }
            int nb=1+(int)(a.brood*(MAX_BROOD-1)+0.5f);
            float share=1.f/(1.f+ (nb-1)*BROOD_TAX*nb );
            float invest=a.energy*(0.20f+0.55f*a.repro_alloc);
            a.energy-=invest;
            for(int bi=0;bi<nb;bi++){
            Animal c=a;
            c.energy=invest*share;
            if(mate){
                #define PICK(F) c.F = (dist01(rng)<0.5f)? a.F : mate->F
                PICK(speed); PICK(turn); PICK(reach); PICK(size); PICK(diet);
                PICK(jaw); PICK(cold_tol); PICK(aquatic); PICK(resp);
                PICK(detritus); PICK(sexual); PICK(arboreal); PICK(brood);
                PICK(lifespan); PICK(toxin); PICK(tox_resist);
                PICK(burrow); PICK(nocturnal); PICK(store_fat); PICK(repro_alloc);
                #undef PICK
                for(int k=0;k<N_W;k++)
                    c.w[k]=(dist01(rng)<0.5f)? a.w[k] : mate->w[k];
                for(int k=0;k<N_NODE;k++)
                    c.ngain[k]=(dist01(rng)<0.5f)? a.ngain[k] : mate->ngain[k];
                for(int k=0;k<N_SENSOR;k++)
                    c.sen[k]=(dist01(rng)<0.5f)? a.sen[k] : mate->sen[k];
                for(int k=0;k<N_IMMUNE;k++)
                    c.im[k]=(dist01(rng)<0.5f)? a.im[k] : mate->im[k];
            }
            c.x=a.x+frand(-2.f,2.f); c.y=a.y+frand(-2.f,2.f);
            if(c.x<0.f) c.x+=GRID_W; if(c.x>=GRID_W) c.x-=GRID_W;
            if(c.y<0.f) c.y=0.f;     if(c.y>=GRID_H) c.y=GRID_H-1.f;
            c.dir=frand(0.f,6.2832f); c.age=0; c.eaten=0;
            c.body=0.10f; c.lin=a.lin;
            c.id=next_animal_id++;
            c.speed   =mutate(a.speed,   0.0f,1.0f);
            c.turn    =mutate(a.turn,    0.0f,1.0f);
            c.reach   =mutate(a.reach,   0.05f,5.0f);
            c.size    =mutate(a.size,    0.10f,4.0f);
            c.diet    =mutate(a.diet,    0.0f,1.0f);
            c.jaw     =mutate(a.jaw,     0.0f,2.0f);
            c.cold_tol=mutate(a.cold_tol,0.0f,1.0f);
            c.aquatic =mutate(a.aquatic, 0.0f,1.0f);
            c.resp     =mutate(a.resp,     0.0f,1.0f);
            c.detritus =mutate(a.detritus, 0.0f,1.0f);
            c.sexual   =mutate(a.sexual,   0.0f,1.0f);
            c.arboreal =mutate(a.arboreal, 0.0f,1.0f);
            c.perch    =0.f;                       // 子は地面から始める
            c.aerobic  =1.f;
            for(int k=0;k<N_IMMUNE;k++) c.im[k]=mutate(a.im[k],0.0f,1.0f);
            // 変異率そのものが遺伝する。環境が荒れた後に上がりうる
            float mr=0.4f+a.mut_rate*MUT_RATE_MAX;
            c.mut_rate=mutate_r(a.mut_rate,0.f,1.f,mr);
            for(int k=0;k<N_REC;k++) c.rec[k]=0.f;   // 記憶は引き継がない
            c.mother=a.id; c.father=mate?mate->id:-1; c.gen=a.gen+1;
            life_note(a.id,g_frame,mate?"有性繁殖":"無性繁殖",(float)nb);
            c.brood   =mutate(a.brood,   0.0f,1.0f);
            c.lifespan   =mutate(a.lifespan,   0.0f,1.0f);
            c.toxin      =mutate(a.toxin,      0.0f,1.0f);
            c.tox_resist =mutate(a.tox_resist, 0.0f,1.0f);
            c.burrow     =mutate(a.burrow,     0.0f,1.0f);
            c.nocturnal  =mutate(a.nocturnal,  0.0f,1.0f);
            c.store_fat  =mutate(a.store_fat,  0.0f,1.0f);
            c.repro_alloc=mutate(a.repro_alloc,0.0f,1.0f);
            c.fat=0.f; c.depth=0.f;
            for(int k=0;k<N_W;k++) c.w[k]=mutate(a.w[k],-3.0f,3.0f);
            // ニューロンの有効・無効。稀に大きく飛んで、層の追加や削除が起きる
            for(int k=0;k<N_NODE;k++)
                c.ngain[k]=(dist01(rng)<0.015f) ? frand(0.f,1.6f)
                                                : mutate(a.ngain[k],0.f,1.6f);
            for(int k=0;k<N_SENSOR;k++){
                c.sen[k].angle=mutate(a.sen[k].angle,-3.1416f,3.1416f);
                c.sen[k].range=mutate(a.sen[k].range, 0.0f,SENSOR_MAX);
                // 受容体の特異性は、稀に大きく変わる
                c.sen[k].tune = (dist01(rng)<0.02f) ? frand(0.f,1.f)
                                                    : mutate(a.sen[k].tune,0.0f,1.0f);
            }
            babies.push_back(c);
            }
        }
    }
    for(auto& b:babies) animals.push_back(b);
    animals.erase(std::remove_if(animals.begin(),animals.end(),
        [](const Animal& a){ return !a.alive; }),animals.end());
    // 分解: 死骸が朽ちると、炭素は大気へ戻る
    for(auto& cp:corpses){
        float rot = cp.energy*ROT_RATE;
        cp.energy -= rot;
        int gx=(int)cp.x, gy=(int)cp.y;
        if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H){
            int ri=idx(gx,gy);
            co2[ri]   += rot * CO2_PER_ENERGY;
            // 海では分解された栄養の一部が有光層に残る
            soil_n[ri]+= rot * N_MINERALIZE * (sea[ri]?0.45f:1.f);
            if(soil_n[ri]>3.f) soil_n[ri]=3.f;
        }
    }
    corpses.erase(std::remove_if(corpses.begin(),corpses.end(),
        [](const Corpse& c){ return c.energy<0.02f; }),corpses.end());
}

// ===== セーブ形式 v2 =====
// 各レコードの並びを、ファイル冒頭に名前で宣言する。
// 読み込みは名前で対応づけるので、遺伝子を足しても消しても古いセーブが読める。
//   足した遺伝子 → 既定値で埋める
//   消した遺伝子 → 読み飛ばす
//   並べ替え    → 影響なし

template<class T>
struct Fld {
    std::string name;
    std::function<double(const T&)> get;
    std::function<void(T&,double)> set;
    double def;
};
// 普通のメンバ(float / int / bool)
template<class T,class M>
Fld<T> F(const char* n, M T::* p, double d){
    return { n, [p](const T& o){ return (double)(o.*p); },
                [p](T& o,double v){ o.*p=(M)v; }, d };
}
// 配列メンバの一要素
template<class T,class M,size_t N>
Fld<T> FA(const std::string& n, M (T::*p)[N], size_t k, double d){
    return { n, [p,k](const T& o){ return (double)((o.*p)[k]); },
                [p,k](T& o,double v){ (o.*p)[k]=(M)v; }, d };
}

// --- 各レコードの定義。遺伝子を足すときは、ここに一行足すだけ ---
std::vector<Fld<Plant>> plant_schema(){
    std::vector<Fld<Plant>> s={
        F("x",&Plant::x,0), F("y",&Plant::y,0),
        F("energy",&Plant::energy,1), F("height",&Plant::height,0.05),
        F("absorb",&Plant::absorb,0.9), F("max_height",&Plant::max_height,0.5),
        F("shade_tol",&Plant::shade_tol,0.1), F("tough",&Plant::tough,0),
        F("clonal",&Plant::clonal,0), F("disperse",&Plant::disperse,0.2),
        F("cold_tol",&Plant::cold_tol,0.05), F("aquatic",&Plant::aquatic,0),
        F("drought_tol",&Plant::drought_tol,0), F("dorm_temp",&Plant::dorm_temp,-15),
        F("dormant",&Plant::dormant,0), F("store_cap",&Plant::store_cap,0),
        F("water_store",&Plant::water_store,0), F("eff_height",&Plant::eff_height,0.05),
        F("fire_tol",&Plant::fire_tol,0), F("climb",&Plant::climb,0),
        F("sexual",&Plant::sexual,0.2), F("brood",&Plant::brood,0.1),
        F("n_fix",&Plant::n_fix,0), F("n_store",&Plant::n_store,0.5),
        F("lifespan",&Plant::lifespan,0.05), F("parasite",&Plant::parasite,0),
        F("carnivory",&Plant::carnivory,0), F("fire_seed",&Plant::fire_seed,0),
        F("repro_alloc",&Plant::repro_alloc,0.5), F("seed_bank",&Plant::seed_bank,0),
        F("defense",&Plant::defense,0), F("seed_wait",&Plant::seed_wait,0),
        F("id",&Plant::id,0), F("genet",&Plant::genet,0),
        F("age",&Plant::age,0), F("lin",&Plant::lin,-1),
    };
    // 後から足した遺伝子。抜けていると読み込み時に未初期化になる
    s.push_back(F("seed_bank",&Plant::seed_bank,0));
    s.push_back(F("buoyancy",&Plant::buoyancy,0.6));
    s.push_back(F("mut_rate",&Plant::mut_rate,0.3));
    s.push_back(F("mother",&Plant::mother,-1));
    s.push_back(F("father",&Plant::father,-1));
    s.push_back(F("gen",&Plant::gen,0));
    for(int k=0;k<N_IMMUNE;k++) s.push_back(FA("im"+std::to_string(k),&Plant::im,k,0.5));
    return s;
}
std::vector<Fld<Animal>> animal_schema(){
    std::vector<Fld<Animal>> s={
        F("x",&Animal::x,0), F("y",&Animal::y,0), F("dir",&Animal::dir,0),
        F("energy",&Animal::energy,1), F("speed",&Animal::speed,0.5),
        F("turn",&Animal::turn,0.2), F("reach",&Animal::reach,0.7),
        F("size",&Animal::size,0.5), F("diet",&Animal::diet,0),
        F("jaw",&Animal::jaw,1), F("cold_tol",&Animal::cold_tol,0.3),
        F("body",&Animal::body,0.1), F("aquatic",&Animal::aquatic,0),
        F("resp",&Animal::resp,0.3), F("detritus",&Animal::detritus,0),
        F("sexual",&Animal::sexual,0.2), F("arboreal",&Animal::arboreal,0),
        F("perch",&Animal::perch,0), F("brood",&Animal::brood,0.1),
        F("lifespan",&Animal::lifespan,0.2), F("toxin",&Animal::toxin,0),
        F("tox_resist",&Animal::tox_resist,0), F("burrow",&Animal::burrow,0),
        F("nocturnal",&Animal::nocturnal,0), F("store_fat",&Animal::store_fat,0.2),
        F("repro_alloc",&Animal::repro_alloc,0.5), F("fat",&Animal::fat,0),
        F("depth",&Animal::depth,0),
        F("id",&Animal::id,0), F("age",&Animal::age,0),
        F("eaten",&Animal::eaten,0), F("lin",&Animal::lin,-1),
    };
    for(int k=0;k<N_W;k++) s.push_back(FA("w"+std::to_string(k),&Animal::w,k,0));
    for(int k=0;k<N_NODE;k++)
        s.push_back(FA("ng"+std::to_string(k),&Animal::ngain,k,k<3?1.0:0.05));
    s.push_back(F("mut_rate",&Animal::mut_rate,0.3));
    s.push_back(F("mother",&Animal::mother,-1));
    s.push_back(F("father",&Animal::father,-1));
    s.push_back(F("gen",&Animal::gen,0));
    for(int k=0;k<N_REC;k++) s.push_back(FA("rec"+std::to_string(k),&Animal::rec,k,0));
    for(int k=0;k<N_SENSOR;k++){
        std::string b="sen"+std::to_string(k)+"_";
        s.push_back({b+"angle",[k](const Animal& a){return (double)a.sen[k].angle;},
                               [k](Animal& a,double v){a.sen[k].angle=(float)v;},0});
        s.push_back({b+"range",[k](const Animal& a){return (double)a.sen[k].range;},
                               [k](Animal& a,double v){a.sen[k].range=(float)v;},0});
        s.push_back({b+"tune", [k](const Animal& a){return (double)a.sen[k].tune;},
                               [k](Animal& a,double v){a.sen[k].tune=(float)v;},0});
    }
    for(int k=0;k<N_IMMUNE;k++) s.push_back(FA("im"+std::to_string(k),&Animal::im,k,0.5));
    return s;
}
std::vector<Fld<Microbe>> microbe_schema(){
    return { F("x",&Microbe::x,0), F("y",&Microbe::y,0),
             F("energy",&Microbe::energy,0.4), F("rate",&Microbe::rate,0.3),
             F("cold_tol",&Microbe::cold_tol,0.2),
             F("id",&Microbe::id,0), F("age",&Microbe::age,0) };
}
std::vector<Fld<Corpse>> corpse_schema(){
    return { F("x",&Corpse::x,0), F("y",&Corpse::y,0),
             F("energy",&Corpse::energy,0), F("origin",&Corpse::origin,1) };
}
std::vector<Fld<Lineage>> lineage_schema(){
    std::vector<Fld<Lineage>> s={
        F("id",&Lineage::id,0), F("kind",&Lineage::kind,0), F("num",&Lineage::num,0),
        F("nrep",&Lineage::nrep,0), F("birth",&Lineage::birth_tick,0),
        F("death",&Lineage::death_tick,-1), F("parent",&Lineage::parent,-1),
        F("alive",&Lineage::alive,1), F("pop",&Lineage::pop,0),
        F("max_pop",&Lineage::max_pop,0), F("missing",&Lineage::missing,0),
    };
    for(int k=0;k<REP_N;k++) s.push_back(FA("rep"+std::to_string(k),&Lineage::rep,k,0));
    for(int k=0;k<REP_N;k++) s.push_back(FA("mean"+std::to_string(k),&Lineage::mean,k,0));
    for(int k=0;k<REP_N;k++)
        s.push_back(FA("fmean"+std::to_string(k),&Lineage::final_mean,k,0));
    s.push_back(F("final_pop",&Lineage::final_pop,0));
    return s;
}
// マスはインデックスで引く
std::vector<Fld<int>> cell_schema(){
    std::vector<Fld<int>> s;
    auto fv=[&](const std::string& n,std::vector<float>& v,double d){
        s.push_back({n,[&v](const int& i){return (double)v[i];},
                       [&v](int& i,double x){v[i]=(float)x;},d}); };
    auto iv=[&](const std::string& n,std::vector<int>& v,double d){
        s.push_back({n,[&v](const int& i){return (double)v[i];},
                       [&v](int& i,double x){v[i]=(int)x;},d}); };
    fv("elev",elev,0.4); fv("co2",co2,CO2_BASE); fv("o2",o2,O2_BASE);
    fv("soil_n",soil_n,SOIL_N_BASE); fv("soil_water",soil_water,0.4);
    fv("pload",pload,0); fv("aload",aload,0);
    iv("burn_age",burn_age,999999);
    for(int k=0;k<N_ODOR;k++) fv("odor"+std::to_string(k),odor[k],0);
    for(int k=0;k<N_IMMUNE;k++){
        fv("pstrain"+std::to_string(k),pstrain[k],0.5);
        fv("astrain"+std::to_string(k),astrain[k],0.5);
    }
    // 地形の構造。再生成すると別の世界になるので、そのまま残す
    iv("plate_id",plate_id,0);
    fv("craton",craton,0); fv("crust",crust,0);
    fv("bnd_type",bnd_type,0); fv("hotspot",hotspot,0);
    return s;
}

// 書き出しと読み込みを同じ関数で記述するための器。
// フィールドの並びが一箇所にしかないので、片方だけ直し忘れることがない
struct SaveW {
    std::ofstream& f;
    template<class T> void operator()(T& v){ f << v << ' '; }
    void operator()(bool& v){ f << (v?1:0) << ' '; }
    void nl(){ f << '\n'; }
};
struct SaveR {
    std::ifstream& f;
    template<class T> void operator()(T& v){ f >> v; }
    void operator()(bool& v){ int t=0; f >> t; v=(t!=0); }
    void nl(){}
};

// --- 書き出しと読み込みの汎用部 ---
template<class T>
void write_schema(std::ostream& f,const char* tag,const std::vector<Fld<T>>& sc){
    f<<"SCHEMA "<<tag<<" "<<sc.size();
    for(auto& d:sc) f<<' '<<d.name;
    f<<'\n';
}
template<class T>
void write_rec(std::ostream& f,const std::vector<Fld<T>>& sc,const T& o){
    for(size_t k=0;k<sc.size();k++){ if(k) f<<' '; f<<sc[k].get(o); }
    f<<'\n';
}
// ファイルの列を、今のフィールドに対応づける
struct ColMap { std::vector<int> col; std::vector<char> seen; };
template<class T>
ColMap build_map(const std::vector<std::string>& names,
                 const std::vector<Fld<T>>& sc,const char* tag){
    ColMap m; m.seen.assign(sc.size(),0);
    std::unordered_map<std::string,int> by;
    for(size_t k=0;k<sc.size();k++) by[sc[k].name]=(int)k;
    std::string extra;
    for(auto& nm:names){
        auto it=by.find(nm);
        int fi=(it==by.end())?-1:it->second;
        m.col.push_back(fi);
        if(fi>=0) m.seen[fi]=1; else extra+=" "+nm;
    }
    std::string miss;
    for(size_t k=0;k<sc.size();k++) if(!m.seen[k]) miss+=" "+sc[k].name;
    if(!miss.empty())  printf("  %s: filled with defaults:%s\n",tag,miss.c_str());
    if(!extra.empty()) printf("  %s: ignored (no longer used):%s\n",tag,extra.c_str());
    return m;
}
template<class T>
void read_rec(std::istream& f,const std::vector<Fld<T>>& sc,const ColMap& m,T& o){
    for(auto& d:sc) d.set(o,d.def);            // まず既定値で埋める
    for(int fi:m.col){ double v; f>>v; if(fi>=0) sc[fi].set(o,v); }
}

// 乱数の内部状態。これが無いと、ロード後の乱数列が元の続きにならず、
// 同じ地点から再開しても別の歴史になる
template<class S> void io_rng(S& s, std::ostream* w, std::istream* r){
    if(w){ *w << rng << '\n'; }
    if(r){ *r >> rng; }
}

// グラフの履歴。世界の状態ではないが、変遷を見せるために残す
template<class S> void io_hist(S& s){
    auto fv=[&](std::vector<float>& v){
        int n=(int)v.size(); s(n);
        v.resize(n);
        for(int i=0;i<n;i++) s(v[i]);
        s.nl();
    };
    auto iv=[&](std::vector<int>& v){
        int n=(int)v.size(); s(n);
        v.resize(n);
        for(int i=0;i<n;i++) s(v[i]);
        s.nl();
    };
    fv(h_maxH); fv(h_reach); fv(h_diet); fv(h_shade); fv(h_tough); fv(h_jaw);
    fv(h_coldP); fv(h_coldA); fv(h_co2); fv(h_o2); fv(h_temp);
    fv(h_soil); fv(h_nut); fv(h_ice); fv(h_path); fv(h_odor); fv(h_sst);
    iv(h_plant); iv(h_animal); iv(h_corpse);
    Series* ss[13]={&s_plant,&s_animal,&s_maxh,&s_reach,&s_temp,
                    &s_co2,&s_o2,&s_soil,&s_nut,
                    &s_shade,&s_diet,&s_coldP,&s_coldA};
    for(int k=0;k<13;k++)
        for(int sc=0;sc<N_SCALE;sc++){
            fv(ss[k]->v[sc]);
            s(ss[k]->acc[sc]); s(ss[k]->cnt[sc]); s.nl();
        }
}

void save_state(const std::string& fn,int frame){
    auto t0=std::chrono::steady_clock::now();
    std::string tmp=fn+".tmp";                 // 書き終えてから差し替える
    {
        std::ofstream f(tmp);
        if(!f){ printf("save failed: %s\n",fn.c_str()); return; }
        f.precision(8);
        auto cs=cell_schema(); auto ps=plant_schema(); auto as=animal_schema();
        auto ks=corpse_schema(); auto ls=lineage_schema();
        auto ms=microbe_schema();
        SaveW sw{f};
        f<<"P3SAVE2\n";
        f<<"HEADER "<<frame<<' '<<next_plant_id<<' '<<next_animal_id<<' '
         <<next_lineage_id<<' '<<next_plant_num<<' '<<next_animal_num<<' '
         <<GRID_W<<' '<<GRID_H<<'\n';
        write_schema(f,"cell",cs);   write_schema(f,"plant",ps);
        write_schema(f,"animal",as); write_schema(f,"microbe",ms);
        write_schema(f,"corpse",ks);
        write_schema(f,"lineage",ls);
        f<<"CELLS "<<GRID_W*GRID_H<<'\n';
        for(int i=0;i<GRID_W*GRID_H;i++) write_rec(f,cs,i);
        f<<"PLANTS "<<plants.size()<<'\n';
        for(auto& p:plants) write_rec(f,ps,p);
        f<<"ANIMALS "<<animals.size()<<'\n';
        for(auto& a:animals) write_rec(f,as,a);
        f<<"MICROBES "<<microbes.size()<<'\n';
        for(auto& m:microbes) write_rec(f,ms,m);
        f<<"CORPSES "<<corpses.size()<<'\n';
        for(auto& c:corpses) write_rec(f,ks,c);
        f<<"LINEAGES "<<lineages.size()<<'\n';
        for(auto& L:lineages) write_rec(f,ls,L);
        f<<"RNG\n"<<rng<<'\n';
        f<<"HISTORY\n";
        io_hist(sw);
        f<<"END\n";
        if(!f){ printf("save failed while writing\n"); return; }
    }
    std::remove(fn.c_str());
    if(std::rename(tmp.c_str(),fn.c_str())!=0){
        printf("save failed: could not replace %s\n",fn.c_str()); return;
    }
    double ms=std::chrono::duration<double,std::milli>(
        std::chrono::steady_clock::now()-t0).count();
    printf("saved: %s (plants %d animals %d lineages %d, %.0f ms)\n",fn.c_str(),
           (int)plants.size(),(int)animals.size(),(int)lineages.size(),ms);
}

bool load_state(const std::string& fn,int& frame){
    std::ifstream f(fn);
    if(!f){ printf("load failed (no file): %s\n",fn.c_str()); return false; }
    std::string tag; f>>tag;
    if(tag!="P3SAVE2"){
        printf("load failed: '%s' is the old format and cannot be read\n",fn.c_str());
        return false;
    }
    try {
        SaveR sr{f};
        std::string sec;
        int gw=0,gh=0;
        f>>sec;
        if(sec!="HEADER") throw std::runtime_error("missing HEADER");
        f>>frame>>next_plant_id>>next_animal_id>>next_lineage_id
         >>next_plant_num>>next_animal_num>>gw>>gh;
        if(gw!=GRID_W||gh!=GRID_H) throw std::runtime_error("grid size differs");

        auto cs=cell_schema(); auto ps=plant_schema(); auto as=animal_schema();
        auto ks=corpse_schema(); auto ls=lineage_schema();
        auto ms=microbe_schema();
        ColMap mc,mp,ma,mk,ml,mm;
        auto read_schema=[&](){
            std::string s,kind; int n; f>>s>>kind>>n;
            if(s!="SCHEMA") throw std::runtime_error("expected SCHEMA, got "+s);
            std::vector<std::string> names(n);
            for(auto& x:names) f>>x;
            if(kind=="cell")    mc=build_map(names,cs,"cell");
            else if(kind=="plant")  mp=build_map(names,ps,"plant");
            else if(kind=="animal") ma=build_map(names,as,"animal");
            else if(kind=="microbe")mm=build_map(names,ms,"microbe");
            else if(kind=="corpse") mk=build_map(names,ks,"corpse");
            else if(kind=="lineage")ml=build_map(names,ls,"lineage");
        };
        for(int k=0;k<6;k++) read_schema();

        auto want=[&](const char* name)->long{
            std::string s; long n=-1; f>>s>>n;
            if(!f||s!=name) throw std::runtime_error(std::string("expected ")+name+", got "+s);
            return n;
        };
        long n=want("CELLS");
        for(int i=0;i<n;i++){ int ii=i; read_rec(f,cs,mc,ii); }
        for(int i=0;i<GRID_W*GRID_H;i++) sea[i]=(elev[i]<SEA_LEVEL)?1:0;

        n=want("PLANTS"); plants.clear(); plants.reserve(n);
        for(long i=0;i<n;i++){
            Plant p; read_rec(f,ps,mp,p);
            if(!f) throw std::runtime_error("plant read failed at "+std::to_string(i));
            p.alive=true; plants.push_back(p);
        }
        n=want("ANIMALS"); animals.clear(); animals.reserve(n);
        for(long i=0;i<n;i++){
            Animal a; read_rec(f,as,ma,a);
            if(!f) throw std::runtime_error("animal read failed at "+std::to_string(i));
            a.alive=true; a.aerobic=1.f; animals.push_back(a);
        }
        n=want("MICROBES"); microbes.clear(); microbes.reserve(n);
        for(long i=0;i<n;i++){
            Microbe m; read_rec(f,ms,mm,m);
            m.alive=true; microbes.push_back(m);
        }
        n=want("CORPSES"); corpses.clear(); corpses.reserve(n);
        for(long i=0;i<n;i++){ Corpse c; read_rec(f,ks,mk,c); corpses.push_back(c); }
        n=want("LINEAGES"); lineages.clear(); lineages.reserve(n);
        for(long i=0;i<n;i++){ Lineage L; read_rec(f,ls,ml,L); lineages.push_back(L); }
        f>>sec;
        if(sec=="RNG"){ f>>rng; f>>sec; }
        if(sec=="HISTORY"){
            io_hist(sr);
            f>>sec;
        }
        if(sec!="END") throw std::runtime_error("END marker missing — file truncated?");

        // 地形から導かれるものを作り直す
        elev_smooth=elev; box_blur(elev_smooth,4,2);
        // 沿岸の近さを作り直す
        std::fill(coast_w.begin(),coast_w.end(),0.f);
        for(int y=0;y<GRID_H;y++) for(int x=0;x<GRID_W;x++){
            int i=idx(x,y);
            if(!sea[i]) continue;
            float best=0.f;
            for(int dy=-COAST_RANGE;dy<=COAST_RANGE;dy++)
            for(int dx=-COAST_RANGE;dx<=COAST_RANGE;dx++){
                int nx=x+dx, ny=y+dy;
                if(nx<0) nx+=GRID_W; if(nx>=GRID_W) nx-=GRID_W;
                if(ny<0||ny>=GRID_H) continue;
                if(sea[idx(nx,ny)]) continue;
                float d=std::sqrt((float)(dx*dx+dy*dy));
                float v=1.f-d/(COAST_RANGE+1.f);
                if(v>best) best=v;
            }
            coast_w[i]=best*best;
        }
        compute_currents();
        compute_rivers();
        update_climate(frame);
        compute_precip();
        classify_climate();
        pgrid.reset();
        for(const auto& p:plants){
            int gx=(int)p.x, gy=(int)p.y;
            if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) pgrid.count(idx(gx,gy));
        }
        pgrid.finalize();
        for(int i=0;i<(int)plants.size();i++){
            int gx=(int)plants[i].x, gy=(int)plants[i].y;
            if(gx>=0&&gx<GRID_W&&gy>=0&&gy<GRID_H) pgrid.put(idx(gx,gy),i);
        }
        printf("loaded: %s (plants %d animals %d lineages %d)\n",fn.c_str(),
               (int)plants.size(),(int)animals.size(),(int)lineages.size());
        return true;
    } catch(const std::exception& e){
        printf("load failed: %s\n",e.what());
        printf("  -> world cleared.\n");
        plants.clear(); animals.clear(); corpses.clear(); lineages.clear();
        return false;
    }
}

// 遺伝子の適性に合う場所を探して動物を放つ
void release_animals(int per_species,int frame){
    if(last_animals.empty()){ printf("no stored genes\n"); return; }
    int total=0;
    for(const auto& g : last_animals){
        // この遺伝子にとって最も条件の良いマスを、無作為に探して選ぶ
        for(int n=0;n<per_species;n++){
            int bx=-1,by=-1; float bscore=-1e9f;
            for(int t=0;t<400;t++){
                int gx=(int)frand(0.f,(float)GRID_W);
                int gy=(int)frand(0.f,(float)GRID_H);
                int ci=idx(gx,gy);
                bool s=sea[ci];
                // 生息環境の不適合
                float hm = s ? (1.f-g.aquatic)*(SHALLOW_EASE+(1.f-SHALLOW_EASE)*depthn(ci))
                             : g.aquatic;
                if(hm>0.55f) continue;
                // 移動できない場所は除く
                if( s && g.aquatic<0.20f) continue;
                if(!s && g.aquatic>0.85f) continue;
                // 気温
                float kill=8.f-40.f*g.cold_tol;
                if(temper[ci]<kill+3.f) continue;
                float opt=26.f-26.f*g.cold_tol;
                float dT=(temper[ci]-opt)/25.f;
                float tfit=std::exp(-dT*dT);
                // 餌: 届く高さの植物と、食える死骸
                float food=0.f;
                for(int t=pgrid.b(ci);t<pgrid.e(ci);t++)
                    if(plants[pgrid.items[t]].height<=g.reach) food+=1.f;
                float detr=(g.detritus>0.05f)? corpse_density[ci]*0.5f : 0.f;
                float score = tfit*2.f + food + detr - hm*3.f;
                if(score>bscore){ bscore=score; bx=gx; by=gy; }
            }
            if(bx<0) continue;
            Animal a;
            a.x=bx+frand(0.f,1.f); a.y=by+frand(0.f,1.f);
            a.dir=frand(0.f,6.2832f);
            a.energy=A_DIVIDE_TH*0.7f; a.alive=true; a.age=0; a.eaten=0;
            a.speed=g.speed; a.turn=g.turn; a.reach=g.reach; a.size=g.size;
            a.diet=g.diet;   a.jaw=g.jaw;   a.cold_tol=g.cold_tol;
            a.aquatic=g.aquatic; a.resp=g.resp; a.detritus=g.detritus;
            a.sexual=g.sexual;   a.arboreal=g.arboreal; a.perch=0.f; a.aerobic=1.f;
            a.mut_rate=0.3f;                    a.mother=-1; a.father=-1; a.gen=0;
            for(int k=0;k<N_REC;k++) a.rec[k]=0.f;
            a.brood=g.brood;     a.lifespan=g.lifespan;
            for(int k=0;k<N_IMMUNE;k++) a.im[k]=g.im[k];
            a.body=g.size*MATURE_FRAC;      // 成熟した状態で放つ
            a.id=next_animal_id++; a.lin=g.lin;
            for(int k=0;k<N_W;k++) a.w[k]=g.w[k];
            for(int k=0;k<N_NODE;k++) a.ngain[k]=g.ngain[k];
            for(int k=0;k<N_SENSOR;k++) a.sen[k]=g.sen[k];
            animals.push_back(a); total++;
        }
    }
    printf("=== released %d animals from %d species at t=%d ===\n",
           total,(int)last_animals.size(),frame);
}

// 長時間の放置中に画面や本体が眠らないようにする
#ifdef __APPLE__
static IOPMAssertionID g_sleep_assert = 0;
void prevent_sleep_begin(){
    CFStringRef reason = CFSTR("evosim long run");
    IOReturn r = IOPMAssertionCreateWithName(
        kIOPMAssertionTypeNoDisplaySleep, kIOPMAssertionLevelOn,
        reason, &g_sleep_assert);
    printf(r==kIOReturnSuccess ? "sleep prevented\n" : "sleep prevention failed\n");
}
void prevent_sleep_end(){
    if(g_sleep_assert) IOPMAssertionRelease(g_sleep_assert);
}
#elif defined(_WIN32)
void prevent_sleep_begin(){
    if(SetThreadExecutionState(ES_CONTINUOUS|ES_SYSTEM_REQUIRED|ES_DISPLAY_REQUIRED))
        printf("sleep prevented\n");
    else
        printf("sleep prevention failed\n");
}
void prevent_sleep_end(){
    SetThreadExecutionState(ES_CONTINUOUS);
}
#else
void prevent_sleep_begin(){}
void prevent_sleep_end(){}
#endif

// C++17 では file_time_type から time_t への標準の変換がないため、
// 現在時刻との差分を取って system_clock に移す
std::time_t file_time_to_t(std::filesystem::file_time_type ft){
    auto sys = std::chrono::time_point_cast<std::chrono::system_clock::duration>(
        ft - std::filesystem::file_time_type::clock::now()
           + std::chrono::system_clock::now());
    return std::chrono::system_clock::to_time_t(sys);
}

// 植物の姿。遺伝子の組み合わせから、どの型に属するかを決める
const int N_FORM = 8;
const char* FORM_NAME[N_FORM]={
    "プランクトン","寄生植物","食虫植物","着生植物","つる","竹型","樹木","草本"};
inline int plant_form(const Plant& p){
    if(p.aquatic > 0.50f) return 0;
    if(p.parasite  > 0.35f) return 1;
    if(p.carnivory > 0.30f) return 2;
    if(p.climb>0.45f && p.max_height<0.8f) return 3;   // 登るが幹を持たない
    if(p.climb     > 0.45f) return 4;
    if(p.clonal    > 0.50f) return 5;
    if(p.max_height> 1.50f) return 6;
    return 7;
}
// 型ごとの、世界に現存する最大と最小の背丈
float form_min[N_FORM], form_max[N_FORM];
int   form_cnt[N_FORM];
void scan_forms(){
    for(int k=0;k<N_FORM;k++){ form_min[k]=1e9f; form_max[k]=-1e9f; form_cnt[k]=0; }
    for(const auto& p:plants){
        if(!p.alive) continue;
        int f=plant_form(p);
        if(p.max_height<form_min[f]) form_min[f]=p.max_height;
        if(p.max_height>form_max[f]) form_max[f]=p.max_height;
        form_cnt[f]++;
    }
    for(int k=0;k<N_FORM;k++) if(form_cnt[k]==0){ form_min[k]=0.f; form_max[k]=1.f; }
}

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
                const float* raw,int n,sf::Color col)
{
    const float TAU=6.2831853f, TOP=-1.5707963f;
    auto vtx=[&](float a,float r){ return sf::Vector2f{cx+std::cos(a)*r,cy+std::sin(a)*r}; };

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
        sf::Text t(font,jp(lb),11);
        t.setFillColor(sf::Color(206,206,214));
        sf::FloatRect bb=t.getLocalBounds();
        sf::Vector2f lp=vtx(a,R+17.f);
        t.setPosition({lp.x-bb.size.x*0.5f,lp.y-7.f});
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

int main(){
    prevent_sleep_begin();
    for(int k=0;k<N_ODOR;k++) odor[k].assign(GRID_W*GRID_H,0.f);
    odor_buf.assign(GRID_W*GRID_H,0.f);
    for(int k=0;k<N_IMMUNE;k++){
        pstrain[k].assign(GRID_W*GRID_H,0.5f);
        astrain[k].assign(GRID_W*GRID_H,0.5f);
    }
    init_terrain();
    compute_currents();      // 海流を求める。気温の計算がこれを使う
    compute_rivers();
    update_climate(0);
    for(int i=0;i<GRID_W*GRID_H;i++) soil_water[i]=sea[i]?FIELD_CAP:0.35f;
    for(int k=0;k<40;k++){                   // 降水と土壌水分を交互に収束させる
        compute_precip();
        update_water(0);
    }
    classify_climate();
    int placed=0,guard=0;
    while(placed<PLANT_INIT_N && guard++<200000){
        Plant p;
        p.x=frand(0.f,(float)GRID_W); p.y=GRID_H/2.f+frand(-40.f,40.f);
        int gx=(int)p.x, gy=(int)p.y;
        if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
        if(!sea[idx(gx,gy)]) continue;              // 生命は海から始まる
        if(depthn(idx(gx,gy))>0.35f) continue;      // ただし浅い海に
        p.energy=PLANT_INIT_E; p.height=0.05f; p.alive=true; p.age=0;
        p.absorb=0.92f; p.max_height=0.25f; p.shade_tol=0.10f;
        p.tough=0.0f; p.clonal=0.0f; p.disperse=0.2f; p.cold_tol=0.05f;
        p.aquatic=1.0f; p.drought_tol=0.0f; p.dorm_temp=-15.0f; p.dormant=false;
        p.store_cap=0.0f; p.water_store=0.0f; p.eff_height=0.05f;
        p.fire_tol=0.0f; p.climb=0.0f; p.sexual=0.2f;
        for(int k=0;k<N_IMMUNE;k++) p.im[k]=frand(0.f,1.f);
        p.brood=0.1f; p.n_fix=0.0f; p.n_store=1.0f; p.lifespan=0.05f;
        p.parasite=0.f; p.carnivory=0.f; p.fire_seed=0.f; p.repro_alloc=0.5f;
        p.seed_bank=0.f; p.defense=0.f; p.seed_wait=0; p.buoyancy=0.6f;
        p.mut_rate=0.3f;        p.mother=-1; p.father=-1; p.gen=0;
        p.id=next_plant_id++; p.genet=p.id; p.lin=-1;
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
    sf::RenderWindow window(desktop,"Phase 3",sf::State::Fullscreen);
    window.setFramerateLimit(60);
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
    if(!texture.loadFromImage(image)){ printf("texture failed\n"); return 1; }
    sf::Sprite sprite(texture);

    float zoom=std::min(SCREEN_W/(float)GRID_W,SCREEN_H/(float)GRID_H);
    float cam_x=GRID_W/2.f, cam_y=GRID_H/2.f;
    // 東西は無限にループする。一番近い複製の位置を返す
    auto w2sx=[&](float wx){
        float dx=wx-cam_x;
        if(dx> GRID_W*0.5f) dx-=GRID_W;
        if(dx<-GRID_W*0.5f) dx+=GRID_W;
        return dx*zoom+SCREEN_W/2.f;
    };
    auto w2sy=[&](float wy){ return (wy-cam_y)*zoom+SCREEN_H/2.f; };
    auto s2wx=[&](float sx){
        float wx=(sx-SCREEN_W/2.f)/zoom+cam_x;
        wx=std::fmod(wx,(float)GRID_W);
        if(wx<0.f) wx+=GRID_W;
        return wx;
    };
    auto s2wy=[&](float sy){ return (sy-SCREEN_H/2.f)/zoom+cam_y; };
    auto clamp_cam=[&](){
        while(cam_x< 0.f)      cam_x+=GRID_W;   // 東西は巻く
        while(cam_x>=GRID_W)   cam_x-=GRID_W;
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

    while(window.isOpen()){
        while(const std::optional event=window.pollEvent()){
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
                    if(k->code==sf::Keyboard::Key::Z) tree_zoom=std::min(12.f,tree_zoom*1.2f);
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
                            p.cold_tol=0.10f; p.aquatic=1.0f; p.drought_tol=0.f;
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
                        }
                        printf("=== released %d plankton at t=%d ===\n",n,frame);
                    }
                    if(k->code==sf::Keyboard::Key::N){ show_night=!show_night; g_last_view=-1; }
                    if(k->code==sf::Keyboard::Key::T){ show_tree=true; sel_lineage=-1; }
                    if(k->code==sf::Keyboard::Key::Z) zoom_at(1.25f,SCREEN_W/2.f,SCREEN_H/2.f);
                    if(k->code==sf::Keyboard::Key::X) zoom_at(0.80f,SCREEN_W/2.f,SCREEN_H/2.f);
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
            if(key_w) tree_oy+=tp;
            if(key_s) tree_oy-=tp;
            if(key_a) tree_ox+=tp;
            if(key_d) tree_ox-=tp;
            float lx=SCREEN_W*0.5f*(tree_zoom-1.f)+40.f;
            float ly=SCREEN_H*0.5f*(tree_zoom-1.f)+40.f;
            if(tree_ox> lx) tree_ox= lx;  if(tree_ox<-lx) tree_ox=-lx;
            if(tree_oy> ly) tree_oy= ly;  if(tree_oy<-ly) tree_oy=-ly;
        } else if(ui_mode==0){
            float pan=12.f/zoom*6.f;
            if(ui_mode==0){
                if(key_w) cam_y-=pan;
                if(key_s) cam_y+=pan;
                if(key_a) cam_x-=pan;
                if(key_d) cam_x+=pan;
            }
            clamp_cam();
        }

        if(!paused){
            g_frame=frame;
            prof_begin(); update_climate(frame);                     prof_end(0);
            prof_begin();
            if(frame%PRECIP_INTERVAL==0){ compute_precip(); classify_climate(); }
                                                                      prof_end(1);
            prof_begin(); update_water(frame);                        prof_end(2);
            prof_begin(); update_animals();                           prof_end(3);
            prof_begin();
            if(frame%PATH_INTERVAL==0) update_pathogens();
            update_plants(); update_microbes(); update_fire(frame);
            if(frame%ODOR_INTERVAL==0) update_odor();
            update_slides();
            prof_end(4);
            prof_begin();
            if(frame%4==0) update_gases();      // 拡散は4tickに一度で足りる
            prof_end(5);
            frame++;

            if((int)plants.size() >peak_plant)  peak_plant =(int)plants.size();
            if((int)animals.size()>peak_animal) peak_animal=(int)animals.size();
            if((int)corpses.size()>peak_corpse) peak_corpse=(int)corpses.size();

            // 植物の数だけを待つと、その間に抑制なしで増え続ける。
            // 時間でも切って、暴走する前に消費者を入れる
            if(!released && ((int)plants.size()>ANIMAL_RELEASE || frame>=150)){
                int put=0,g2=0;
                while(put<ANIMAL_INIT_N && g2++<200000){
                    Animal a;
                    a.x=frand(0.f,(float)GRID_W); a.y=GRID_H/2.f+frand(-50.f,50.f);
                    int gx=(int)a.x, gy=(int)a.y;
                    if(gx<0||gx>=GRID_W||gy<0||gy>=GRID_H) continue;
                    if(!sea[idx(gx,gy)]) continue;
        if(depthn(idx(gx,gy))>0.12f) continue;
        if(shelf[idx(gx,gy)]<0.5f) continue;   // 岸のそばに置く
                    a.dir=frand(0.f,6.2832f); a.energy=3.0f; a.alive=true;
                    a.age=0; a.eaten=0;
                    a.speed=0.5f; a.turn=0.2f; a.reach=0.7f;
                    a.size=0.5f; a.diet=0.0f; a.jaw=1.0f; a.cold_tol=0.3f;
                    a.body=0.10f; a.aquatic=1.0f; a.resp=0.3f; a.detritus=0.0f;
                    a.sexual=0.2f; a.arboreal=0.0f; a.perch=0.0f; a.aerobic=1.0f;
                    a.mut_rate=0.3f;                    a.mother=-1; a.father=-1; a.gen=0;
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
                }
                released=true;
                printf("=== animals released at t=%d ===\n",frame);
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
                { double sfi=0,scl2=0,sx=0;
                  for(const auto& p:plants){ sfi+=p.fire_tol; scl2+=p.climb; sx+=p.sexual; }
                  double asx=0; for(const auto& a:animals) asx+=a.sexual;
                  printf("      fire:%.0f cells (%ld ignitions) | fireTol:%.3f climb:%.3f | sexual P:%.3f A:%.3f\n",
                         fire_area,fire_events,n?sfi/n:0.0,n?scl2/n:0.0,
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
                  for(int i=0;i<GRID_W*GRID_H;i++){ o2m2+=o2[i]; tot2+=o2[i]+co2[i]+N2_BASE; }
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
                    printf("      defense %.3f (tough %.3f) | seedBank %.3f waiting %d | ice %d cells | slides %ld\n",
                           n?sdf/n:0.0, n?stg/n:0.0, n?sbk/n:0.0, nseed, nice,
                           slide_events); }
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
                printf("      A:%d speed:%.3f reach:%.3f body:%.3f diet:%.3f detr:%.3f jaw:%.3f cold:%.3f aq:%.3f | CO2:%.3f O2:%.3f corpse:%d(plant %d) | sp P%d A%d\n",
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
                  printf("      microbe: %d  rate %.3f  cold %.3f  deaths %ld\n",
                         nm, nm?mr2/nm:0.0, nm?mc2/nm:0.0, d_micro); }
                printf("      deaths: starve %ld  cold %ld  preyed %ld\n",
                       d_starve,d_cold,d_preyed);
            }
        }

        // ===== 背景 =====
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
                float ca=clamp01(co2[i]/(CO2_BASE*1.5f));
                float cb=clamp01(o2[i]/(O2_BASE*1.5f));
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
                int pl=plate_id[i];
                float hue=std::fmod(0.11f+pl*0.61803399f,1.0f)*6.f;
                int hi=(int)hue; float fr=hue-hi;
                float q=plate_oce[pl]?0.55f:1.0f;
                float r,g,b;
                switch(hi){
                    case 0: r=1;g=fr;b=0; break;   case 1: r=1-fr;g=1;b=0; break;
                    case 2: r=0;g=1;b=fr; break;   case 3: r=0;g=1-fr;b=1; break;
                    case 4: r=fr;g=0;b=1; break;   default:r=1;g=0;b=1-fr; break;
                }
                float bt=bnd_type[i];
                if(craton[i]>0.35f && std::fabs(bt)<0.15f){
                    float cr=craton[i];
                    c=sf::Color((std::uint8_t)(90+80*cr),(std::uint8_t)(70+50*cr),
                                (std::uint8_t)(40+30*cr));
                }
                else if(bt>0.15f) c=sf::Color(255,90,60);
                else if(bt<-0.15f)c=sf::Color(90,200,255);     // 発散境界
                else c=sf::Color((std::uint8_t)(60+150*r*q),
                                 (std::uint8_t)(60+150*g*q),
                                 (std::uint8_t)(60+150*b*q));
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
                    // 高所は岩肌が露出し、頂は雪を頂く
                    if(elev[i]>MOUNTAIN_LEVEL){
                        float m=clamp01((elev[i]-MOUNTAIN_LEVEL)/(1.f-MOUNTAIN_LEVEL));
                        r=r*(1.f-m)+(96+150*m)*m;
                        g=g*(1.f-m)+(94+150*m)*m;
                        b=b*(1.f-m)+(90+155*m)*m;
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
        nrep  = (int)std::ceil(SCREEN_W/std::max(1.f,wrapW))+1;
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
            int klo=(int)std::ceil((-20.f-sx0)/wrapW);
            int khi=(int)std::floor((SCREEN_W+20.f-sx0)/wrapW);
            for(int k=klo;k<=khi;k++) quad(sx0+k*wrapW,sy,r,col);
        };
        // 引きの画面では点で描く。頂点数が6分の1になる
        auto point_wrapped=[&](float sx0,float sy,sf::Color col){
            if(sy<-20.f||sy>SCREEN_H+20.f) return;
            int klo=(int)std::ceil((-20.f-sx0)/wrapW);
            int khi=(int)std::floor((SCREEN_W+20.f-sx0)/wrapW);
            for(int k=klo;k<=khi;k++){
                sf::Vertex v; v.position={sx0+k*wrapW,sy}; v.color=col;
                vbuf.push_back(v);
            }
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
        for(int i=0;i<(int)plants.size();i++){
            if(plants[i].seed_wait>0) continue;       // 土中の種子は見えない
            if(plants[i].y<vy0||plants[i].y>vy1) continue;
            dorder.push_back(i);
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
                float PX=18.f,PY=18.f,PW=330.f,PH=sa?520.f:560.f;
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
                        draw_radar(window,font,PX+PW*0.5f,by+100.f,82.f,
                                   AN,av,ar,13,sf::Color(255,190,110));
                        by+=208.f;
                        sf::Text nt(font,jp("神経網"),12);
                        nt.setFillColor(sf::Color(200,210,230));
                        nt.setPosition({PX+10.f,by}); window.draw(nt);
                        by+=16.f;
                        draw_nn(window,font,*sa,PX+6.f,by,PW-12.f,116.f);
                        by+=124.f;
                    }
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
                        draw_radar(window,font,PX+PW*0.5f,by+100.f,82.f,
                                   PN,pv,pr,13,sf::Color(140,235,140));
                        by+=208.f;
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
                for(int i=0;i<N&&i<(int)v.size();i++){
                    sf::Vertex vt;
                    vt.position={px0+pw*i/(float)(N-1),
                                 py0+ph*(1.f-(float)v[i]/std::max(1,pk))};
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
                lineF(py,h_co2,0.f,CO2_BASE*1.5f,sf::Color(230,90,200));
                lineF(py,h_o2 ,0.f,O2_BASE*1.5f ,sf::Color(70,180,240));
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
                                    "プレート","海流"};
                int span=N*SAMPLE_INTERVAL;
                snprintf(b,240,"span -%.1fk .. now   t %d   zoom %.1fx   [V] %s%s",
                         span/1000.f,frame,zoom,vm[view_mode],paused?"   [PAUSED]":"");
                // 今どの地図を見ているかを、画面上部にはっきり出す
                if(font_ok){
                    char mb[64];
                    snprintf(mb,64,"[V]  %s",vm[view_mode]);
                    sf::Text mt(font,jp(mb),17);
                    sf::FloatRect mr=mt.getLocalBounds();
                    sf::RectangleShape mbg({mr.size.x+24.f,28.f});
                    mbg.setPosition({SCREEN_W*0.5f-(mr.size.x+24.f)*0.5f,8.f});
                    mbg.setFillColor(sf::Color(0,0,0,200));
                    mbg.setOutlineColor(sf::Color(150,170,190));
                    mbg.setOutlineThickness(1.f);
                    window.draw(mbg);
                    mt.setFillColor(sf::Color(235,240,250));
                    mt.setPosition({SCREEN_W*0.5f-mr.size.x*0.5f,11.f});
                    window.draw(mt);
                }
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
                float PX=18.f,PY=110.f,PW=330.f,PH=330.f;
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
                    draw_radar(window,font,PX+PW*0.5f,by+94.f,78.f,
                               PN,pv,pr,9,lineage_color(L.id));
                    by+=196.f;
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
                    draw_radar(window,font,PX+PW*0.5f,by+94.f,78.f,
                               AN,av,ar,9,lineage_color(L.id));
                    by+=196.f;
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
            float lx=18.f, ly=SCREEN_H-8.f-8*20.f;
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
            float PW2=SCREEN_W*0.80f, PH2=SCREEN_H*0.66f;
            float PX2=(SCREEN_W-PW2)*0.5f, PY2=(SCREEN_H-PH2)*0.5f;
            sf::RectangleShape bg2({PW2,PH2});
            bg2.setPosition({PX2,PY2});
            bg2.setFillColor(sf::Color(8,14,12,245));
            bg2.setOutlineColor(sf::Color(130,150,160)); bg2.setOutlineThickness(1.f);
            window.draw(bg2);

            struct GDef { const char* name; Series* a; Series* b;
                          const char* la; const char* lb; float lo,hi; };
            GDef gd[7]={
                {"個体数",     &s_plant,&s_animal,"植物","動物",     0.f,-1.f},
                {"植物の形質", &s_maxh, &s_shade, "最大樹高","耐陰性", 0.f,5.f},
                {"動物の形質", &s_reach,&s_diet,  "届く高さ","肉食性", 0.f,5.f},
                {"陸の平均気温",&s_temp, nullptr,  "気温",nullptr,   -30.f,40.f},
                {"耐寒性",     &s_coldP,&s_coldA, "植物","動物",     0.f,1.f},
                {"大気",       &s_co2,  &s_o2,    "CO2","O2",       0.f,CO2_BASE*1.5f},
                {"生産者と消費者",&s_plant,&s_animal,"植物","動物",  0.f,-1.f},
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
            float ax=PX2+70.f, ay=PY2+60.f, aw2=PW2-100.f, ah=PH2-140.f;
            sf::RectangleShape fr2({aw2,ah});
            fr2.setPosition({ax,ay});
            fr2.setFillColor(sf::Color(0,0,0,0));
            fr2.setOutlineColor(sf::Color(70,80,90)); fr2.setOutlineThickness(1.f);
            window.draw(fr2);
            const std::vector<float>& va=g.a->v[graph_scale];
            const std::vector<float>* vb=g.b? &g.b->v[graph_scale] : nullptr;
            int N2=(int)va.size();
            int win=std::min(N2,300);
            int maxoff=std::max(0,N2-win);
            if(graph_off>maxoff) graph_off=maxoff;
            if(graph_off<0) graph_off=0;
            int st3=N2-win-graph_off; if(st3<0) st3=0;
            float lo=g.lo, hi=g.hi;
            if(hi<0.f){                        // 自動目盛
                hi=1.f;
                for(int i=st3;i<st3+win&&i<N2;i++){
                    if(va[i]>hi) hi=va[i];
                    if(vb&&i<(int)vb->size()&&(*vb)[i]>hi) hi=(*vb)[i];
                }
                hi*=1.1f;
            }
            auto plot=[&](const std::vector<float>& v,sf::Color col){
                if((int)v.size()<2) return;
                sf::VertexArray ln(sf::PrimitiveType::LineStrip);
                for(int i=st3;i<st3+win&&i<(int)v.size();i++){
                    float t=(v[i]-lo)/std::max(1e-6f,hi-lo);
                    sf::Vertex vt;
                    vt.position={ax+aw2*(i-st3)/(float)std::max(1,win-1),
                                 ay+ah*(1.f-clamp01(t))};
                    vt.color=col; ln.append(vt);
                }
                window.draw(ln);
            };
            plot(va,sf::Color(120,235,140));
            if(vb) plot(*vb,sf::Color(255,170,90));
            // 目盛
            char lb2[64];
            for(int k=0;k<=4;k++){
                float t=k/4.f, yy=ay+ah*(1.f-t);
                sf::RectangleShape gl2({aw2,1.f});
                gl2.setPosition({ax,yy}); gl2.setFillColor(sf::Color(50,58,66));
                window.draw(gl2);
                snprintf(lb2,64,"%.1f",lo+(hi-lo)*t);
                sf::Text t2(font,jp(lb2),11);
                t2.setFillColor(sf::Color(150,160,170));
                t2.setPosition({PX2+18.f,yy-8.f}); window.draw(t2);
            }
            // 凡例と操作の案内
            snprintf(lb2,64,"%s",g.la);
            { sf::Text t2(font,jp(lb2),14); t2.setFillColor(sf::Color(120,235,140));
              t2.setPosition({ax,ay+ah+10.f}); window.draw(t2); }
            if(g.lb){
                snprintf(lb2,64,"%s",g.lb);
                sf::Text t2(font,jp(lb2),14); t2.setFillColor(sf::Color(255,170,90));
                t2.setPosition({ax+110.f,ay+ah+10.f}); window.draw(t2);
            }
            snprintf(lb2,64,"スクロールで時間移動   %d / %d",graph_off,maxoff);
            { sf::Text t2(font,jp(lb2),12); t2.setFillColor(sf::Color(140,150,160));
              t2.setPosition({ax,ay+ah+32.f}); window.draw(t2); }
            ui_hits.push_back({ax,ay,aw2,ah,220});   // スクロール領域
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
                                     "病原体","匂い","昼夜","栄養","プレート","海流"};
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
    }
    prevent_sleep_end();
    return 0;
}
