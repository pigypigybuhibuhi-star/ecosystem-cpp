# evosim phase3 の仕様書(ODD プロトコル)

ODD プロトコル(Overview, Design concepts, Details。Grimm et al. 2020)の形で、シミュレータの中身を書く。
コードを読まずに分かることを目標にする。式と定数はコード(phase3.cpp と src/)の値をそのまま写した。

- 版: 最初の版(2026-10-28 夜)。main 3e4e230 の phase3.cpp に基づく。7.14 病原体を 2026-10-28 22 時に書き足した。7.10 の植物の特殊な形質と文献の一覧を 2026-10-28 23 時に書き足した
- 書き終えていない所は「(次の版で書く)」と記した
- 印: **[法則]** は既知の物理・生物の法則から来る規則。**[条件]** は外から与える条件(惑星規模の状態・始まりの状態・区画の外)。**[未換算]** は、移行の段階 3 より前の「1 回の更新あたりの値」を時間の単位に換算しただけで、まだ実際の 1 分あたりの物理量に置き換えていない定数(移行の段階 5 の対象)

### 時間の換算の書き方
速さの定数はすべて「1 分あたりの物理量 × DT」で持つ(DT = 1 分)。次の換算が全体で使われる。

| 書き方 | 意味 |
|---|---|
| PM(v, n) | 段階 3 より前に n tick(1 tick = 約 103 分)ごとの 1 回の更新で v だけ増減していた量を、1 分あたり v/(n·DT0) に直したもの(DT0 = 1440/14 分) |
| PMF(f, n) | 同じく、1 回の更新で割合 f だけ減っていた量を、1 分あたりの率 −ln(1−f)/(n·DT0) に直したもの |
| LIN(r, n) | 1 分あたり r の量を n tick 分にしたもの(r·n·DT) |
| FRAC(k, n) | 率 k で n tick の間に起きる(減る)割合 1 − exp(−k·n·DT) |

PM と PMF を通った値が **[未換算]** に当たる。1 回の更新に何回起きるか(齧る回数、繁殖の判定など)は「1 分あたりの率 × DT」の期待値で決め、端数は確率で決める(端数が 0 なら乱数を引かない)。向きのゆらぎは角度 × √(DT/DT0)。

---

## 1. 目的と型(Purpose and patterns)

(作者が書く)

---

## 2. 世界にいるもの・状態の変数・尺度(Entities, state variables, and scales)

### 2.1 尺度
- 空間: 900 × 500 マス、1 マス = 10 m(9 km × 5 km の一区画の土地)。**[条件]**(作者の決定)
- 時間: 1 tick = 1 分。1 日 = 1440 tick、1 年 = 365 日 = 525600 tick。昼夜と季節はこの二つから導く
- 区画の端: 開いた境界。輪につながない。風上の端から入る空気の湿度・温度・組成は条件として与える
- 生き物の位置はマス単位の連続座標(実数)で持つ。マスは環境の場に使う
- 体の尺度: 動物の「体 1」= 生重 100 mg **[条件]**。組織 1 あたりの炭素(エネルギー)は CORPSE_RETURN = 2.5 単位。炭素の単位は「1 単位 = 6 mg C」に揃える方針(D-2026-10-24-2 の作者の答え)

### 2.2 環境の場(マスごと)
| 場 | 意味 |
|---|---|
| elev, slope | 標高(尺度 0〜1、海面 0.42、尺度 1 ≈ 690 m)、勾配(m/m) |
| sea, ice, depth | 海かどうか、海氷、水深 |
| river | 川のマスかどうか(上流の流量から決める) |
| soil_water | 根の層の土の水(1.0 = 100 mm) |
| 地下水の飽和した厚さ | 帯水層(厚さ 3 m)の中の水の厚さ |
| salt_soil, salt_gw | 根の層の土の水と帯水層の水に溶けた塩の量(陸のマス。水の量で割ると塩分 g/kg。海は 35 g/kg で一定) |
| soil_n | 土(海では有光層)の無機の窒素 |
| temper | 気温(°C) |
| solar, daylit | 日射、昼かどうか |
| co2, o2 | 大気の二酸化炭素と酸素の、マスの上の空気の柱の量からのずれ(気体の単位 = 0.3 g C = 0.025 mol。7.8) |
| odor[4] | 匂い: 0 生体・1 腐敗・2 食害・3 災害(煙) |
| fire, fuel, 燃料の湿り | 火勢、燃料(枯れた植物の死骸と株の高さ)、燃料の含水率 |
| pstrain, astrain, pload, aload | 植物と動物の病原体の型(3 座位)と量 |
| canopy_h | そのマスの最も高い株の高さ(林冠) |
| corpse_density | 死骸の数 |

### 2.3 植物(Plant)
状態: 位置、エネルギー(体の外の糖の貯え)、高さ(体)、体内の窒素、貯水、実効的な高さ(つるの分を含む)、今張っている防御、土の中で待つ残り時間(種子)、休眠しているか、年齢、系譜(母・父・世代)、ジェネット。

遺伝子(すべて 0〜1 などの連続値。変異する): 光の吸収(absorb)、最大の高さ、耐陰性、防御の上限(tough)、常に張る防御の割合(def_const)、地下茎(clonal)、散布(disperse)、耐寒性、水に浸かることへの強さ(aquatic。通気組織)、塩への強さ(salt_tol。適合溶質で浸透圧を合わせられる塩分、海水を 1 とする比)、耐乾性、休眠の気温、貯水の容量、耐火性(厚い樹皮)、つる(climb)、有性生殖の度合い、一度に作る種子の数(brood)、窒素固定、寿命、寄生、食虫、火後発芽、繁殖に回す割合、土中で待つ傾向、浮く強さ(buoyancy)、変異率、免疫の型(3 座位)。

### 2.4 動物(Animal)
地面を歩く昆虫程度を想定する(歩く速さ 0.5〜10 cm/s、一世代およそ 2 週間。CLAUDE.md「世界の形」)。

状態: 位置、向き、エネルギー、体(body、成長して size に近づく)、脂肪、潜り具合、今いる高さ(perch)、酸素の充足度、再帰の値 2 つ、年齢、食べた回数、系譜。

遺伝子: 速さ、旋回、届く高さ(reach)、成体の大きさ(size)、肉食(diet)、顎、耐寒性、水中での呼吸と泳ぎ(aquatic)、塩への強さ(salt_tol。体液の浸透圧を保てる外の塩分、海水を 1 とする比)、呼吸の効率、分解物食(detritus)、有性生殖の度合い、登攀(arboreal)、一度に産む数(brood)、寿命、毒、解毒、穴掘り(burrow)、夜目(nocturnal)、脂肪を貯める日数(store_fat × 60 日)、繁殖に回す割合、性フェロモンを出す量(estrus)、変異率、一口の大きさ、消化効率、免疫の型(3 座位)、感覚器 5 つ(向き・射程・何に反応するか)、神経網の重み(N_W 個)とニューロンの有効・無効の値(24 個)。

### 2.5 微生物(Microbe)
動かない分解者。状態: 位置、体の炭素、年齢。遺伝子: 分解の速さ(rate)、耐寒性。

### 2.6 死骸(Corpse)
位置、炭素(エネルギー)、窒素、由来(植物か動物か)。糞と落とした葉も植物由来の死骸に足す。

### 2.7 雨雲(Storm)
位置(区画の外の風上を含む)、年齢、中心の降水強度。一生は約 60 分。

### 2.7.1 本震(Mainshock)
起きた tick、規模、震央の位置(区画の外の断層の上)。余震の率の計算のために残す(7.7.1)。

### 2.8 帳簿(Budget)
水・炭素・窒素・塩の出入りを過程ごとに数える観察のための量。生き物と環境の規則の入力には使わない。

### 2.9 系統(lin)
観察のための分類。形質の距離 SPECIES_THRESHOLD = 0.85 で分ける。**注意**: 今のコードでは、交配の相手を選ぶとき(動物・植物)に lin を使っている。これは「観察の量を生き物の規則に使わない」に反し、「上の階層の書き込みを外す」の項目として作者が洗い出し済み(未修正)。仲間の感覚(帯 11)は lin を使わず、免疫の型の似ている度合いで感じ取る(下の「感覚」)。

---

## 3. 処理の順番(Process overview and scheduling)

1 tick(1 分)ごとに、次の順で処理する。

1. 気候(update_climate): 太陽の位置・日射・気温・海氷
2. 雨雲(update_storms): 発生・移動・雨・冷気・雷。雨を土と地下水へ入れる
3. 300 tick ごとに: 降水の場の計算(compute_precip)と気候区分(classify_climate)
4. 土の水(update_water)、地下水(update_groundwater、60 tick ごとに流れを解く)
5. 動物(update_animals): 2 tick に分けて半数ずつ。感覚 → 神経網 → 成長・登り降り → 移動 → 気温による死 → 採食・漁る・襲う → 食虫植物 → 維持費と呼吸 → 脂肪 → 繁殖。最後に死んだ個体を除き、糞を死骸に、死骸の風化
6. 病原体(3 tick ごと)、植物(update_plants): 光合成は毎 tick、体の更新は 3 tick に分けて 1/3 ずつ、土の中の種子は 8 tick に分けて
7. 微生物(4 tick に分けて 1/4 ずつ)
8. 火(update_fire)
9. 匂い(4 tick ごと)
10. 地震(update_quakes)と崖崩れ(update_slides)。再生のときは記録の出来事を置く
11. 大気の移流と拡散(4 tick ごと)

索引(植物・動物・死骸のマスごとの一覧)は使う直前に作り直す。動物と植物の配列は生まれた順に並び、子は tick の終わりに後ろへ足す。

---

## 4. 設計の考え方(Design concepts)

- **基本の原則**: 物理法則に忠実であること。地球との一致は目標ではなく結果。惑星規模で決まるもの(緯度・日射・卓越風・背景の気温と湿度)は条件として外から与え、区画の中で起きること(雨・川・地下水・火・煙・崖崩れ・生き物)は法則で動かす
- **出てくるもの(emergence)**: 砂漠、動物の行動、種は直接書かず、下の階層の規則から出す。気候区分は観察のために判定するだけ
- **適応**: 動物の行動はすべて神経網の出力で決める。遺伝子が決めるのは能力・費用・感覚の範囲だけ
- **目的(objectives)**: 個体は何も最適化しない。残るかどうかは繁殖とエネルギーの収支から決まる
- **感じ取ること**: 動物は感覚器(視覚・嗅覚・水際・似た体)と体の内の状態(エネルギー・気温・土の水・体格)を入力にする。感覚器が何に反応するかは遺伝子で、進化する。感じ取れない前触れは前触れにならない
- **相互作用**: 食う・食われる(植物、動物、死骸)、寄生、食虫、光をめぐる背比べ、窒素の取り合い、病原体(最も多い免疫の型に適応する。負の頻度依存)、常に張る防御と匂いによる誘導防御
- **確率**: 起きる確率は 1 − exp(−率 × 時間)。乱数は一本の列で、同じ種と同じ引数なら同じ結果になる(決定性)
- **集団(collectives)**: 地下茎でつながった株(ジェネット)はエネルギーを融通する
- **観察**: 帳簿、系統、選択勾配、全体の記録(--census)、神経網の書き出し(--dump-brains)と試験器(--assay)、野外の記録(--track・--trigger-log)。観察の量は生き物の規則の入力に使わない(2.9 の注意を除く)
- **保存則**: 水・炭素・窒素は明示した出入り(区画の端・海・雨・流出など)以外で増減しない。帳簿で確かめる

---

## 5. 始まりの状態(Initialization)

1. 地形を作る(7.1)。300 万年分の隆起と侵食を回す
2. 吹送流、川を決める
3. 気候を t=0 で評価する。土の水は仮の値 0.35(陸)で降水を決め、その降水の釣り合いの値を置き直して降水を決め直す **[条件]**
4. 植物 400 株を置く **[条件]**(区画のまわりの陸と海から種が入ってくる。D-2026-10-24-1 の答え (b)・(c))。遺伝子の初期値は sea_plant_genes(光の吸収 0.92、最大の高さ 0.25、水生 1.0 など)
   - 最初の 200 株(PLANT_INIT_LAND)は陸の、土の水が発芽の下限 0.25 以上のマスに置き、陸に合う遺伝子にする: 水生 0、根を張る(浮く強さ 0.3)、休眠の気温 5 °C(PLANT_LAND_DORM。温帯の草本の生育の基準温度、McMaster & Wilhelm 1997)
   - 残りの 200 株は浅い海(深さの尺度 0.35 未満)に置く。その半数は海底に根を張る(浮く強さ 0.3)。風上の端から入る浮く植物の密度の条件(6 の表)は、この海に置いた株の、海全体の平均の密度とする
5. 微生物 3000 個体
6. 動物は、植物が 20000 株に達したとき、または RELEASE_TICK(約 10.7 日)に放つ。遺伝子の適性に合うマスを探して置く。500 匹 **[条件]**。放す個体は育ち切った体で放す

---

## 6. 外から与える条件(Input data)

| 条件 | 値 | 根拠 |
|---|---|---|
| 緯度 | 北緯 40° | 作者の決定 |
| 背景の気温 | 赤道 30 °C・極 −15 °C の間を緯度で補う。気温減率 6.5 °C/km。年較差 24 °C × 緯度、日較差 陸 9 °C・海 1.2 °C | 惑星規模の条件 |
| 地軸の傾き | 23.4° | 地球 |
| 卓越風 | 地上 10 m で 3 m/s、向きは緯度の風の帯 | Vautard et al. 2010 |
| 中層の風(雨雲を流す) | 10 m/s | 中緯度 500 hPa |
| 風上から入る空気 | 相対湿度 0.75、CO2 と O2 は柱の量(ずれ 0)、匂いは持たない | Dai 2006 |
| 陸の日平均の相対湿度(燃料の湿り用) | 0.70 | Dai 2006 |
| 隆起 | 北の端 0.3 mm/年、場所ごとのばらつき 40%、300 万年 | 中程度に活発な海岸の山地 |
| 区画の南の端の海の深さ | 60 m、南から 25% の帯が沈降 | 作者の決定 |
| 帯水層 | 厚さ 3 m、比産出率 0.2 | Heimsath et al. 1997、Johnson 1967 |
| 海の水の塩分 | 35 g/kg(陸の帯水層へ入る海水も同じ) | Millero et al. 2008 |
| 風上の端から入る浮く植物 | 最初に海に置いた株(200 株)の海の平均の密度 × 漂流の速さ | D-2026-10-26-1 の答え (b) |
| 区画の外の断層 | 北の端から 10 km 北を東西に走る長さ 60 km、地震発生層の幅 15 km、すべりの速さ = 隆起の速さ 0.3 mm/年、剛性率 3e10 Pa | 傾動地塊を持ち上げる断層とみなす(作者の CLAUDE.md「震源の位置と規模は条件」) |
| 強制の本震(--quake T M) | tick T に規模 M の本震を断層の中央(区画の中央の北)に置く | 研究の比較・試験のため |
| 記録した環境の出来事(--replay) | 火・崖崩れ・雨雲・雷・地震を、記録の tick・位置・強さで置く | 研究の比較のため |

---

## 7. 個々の仕組み(Submodels)

### 7.1 地形の生成 [法則]+[条件]
- 初期の面: 南の端で海の深さ 60 m の緩い傾きに、起伏(標準偏差 2 m)を重ねる **[条件]**
- 隆起 U: 北ほど大きい傾動する地塊。北の端 0.3 mm/年 **[条件]**
- 川の侵食: ストリームパワー則 dz/dt = U − K·A^m·S(m = 0.5、n = 1、K = 1e-5 /年、A は m²)。Cordonnier et al. 2016 の陰的解法、刻み 1 万年 × 300 回 **[法則]**(K は Stock & Montgomery 1999、Harel et al. 2016 の範囲の中ほど)
- 川が削り始める流域: 100 マス(0.01 km²)**[法則]**(Montgomery & Dietrich 1988)
- 斜面拡散: dz/dt = D·∇²z、D = 0.0032 m²/年 **[法則]**(Roering et al. 1999)
- 安息角: 勾配 0.70(35°)を超える斜面は崩れて下へ移す。質量は保つ **[法則]**
- 窪地は Priority-Flood 法で埋め、南の端(海)だけを出口にする。北と東西の端からは川が出入りしない **[条件]**
- 出てきた高さと傾斜は写し直さずに使う

### 7.2 気候 [条件]+[法則]
- 太陽高度 sin h = sin φ·sin δ + cos φ·cos δ·cos H を解く。赤緯 δ は季節、時角 H は自転から出す **[法則]**。区画の中で緯度と太陽時は一様
- 気温 = 背景の気温(緯度・季節・標高)+ 日較差。沿岸の海は陸の気温に引かれる(7 マス、強さ 0.55)**[未換算に近い経験式]**
- 雨雲の冷気で気温が下がる(7.5)
- 海氷: −5 °C で張り、−1.5 °C で融ける。季節に遅れて厚さとして蓄える。氷は光を 8% だけ通し、陸の動物には橋になる

### 7.3 降水の場と土の水 [条件]+経験式
- 降水の場(300 tick ごと): 風上の端から相対湿度 0.75 の空気が入り、海で蒸発(海面の温度に比例)、陸で再循環(0.04)、空気が冷えるか斜面を上ると保持できる量を超えた分が降る(地形性降雨)。6 マスでぼかす
- 土の水: ケッペンの乾燥限界(20T + 140 にならい気温に比例)から目標 A/(A + 1.2) を決め、その目標へ率 1.96e-4 /分 で近づける(緩和)**[未換算]**。雨季と乾季は wet = 1 + 0.45·緯度·sin(季節)
- **注意**: この緩和は「降水 − 蒸発散 − 排水」の流れではなく、気候区分の指数を使う上の階層の書き方で、雨からの涵養が 0 になる原因になっている。D-2026-10-28-2 で作者の判断を待っている
- 川のマスは河川の供給 0.6 を受ける。川は D8 の最急降下で流量を上流から積んで決める(閾値 260 マス)。**注意**: 今の川は「上流に降った雨を積み上げた流量」ではなく流域のマスの数で決まっていて、CLAUDE.md「水の循環」の項目(未実装)

### 7.4 地下水 [法則]
- ダルシーの式とブシネスクの近似。水頭 = 岩盤の高さ + 飽和した厚さ。隣のマスとの流れ = K × 上流側の飽和した厚さ × 水頭の差、K = 1e-5 m/s(Freeze & Cherry 1979)。60 tick ごとに解く(陽解法の安定の限界は約 1100 分)
- 海のマスは海面の水頭で一定で、帯水層は海につながり満ちているとする。海との間は両向きに流れ、海から陸へ入る水は塩分 35 g/kg の海水(海水の侵入。D-2026-10-28-1 の作者の答え (b))**[条件]**。区画の端は、端の傾きがそのまま外へ続き、端の外の水の塩分は端のマスと同じとする **[条件]**
- 地下水面が地表に届いた分は湧き出し、まず根の層の土を満たし、残りは地表を流れ去る
- **塩** [法則]: 塩は水に溶けて運ばれるだけで、生まれも消えもしない。帯水層の中では水と一緒に動く(移流。辺を通る塩 = 水の流量 × 上流側の塩分)。湧いた水は帯水層の塩分の塩を持って出て、根の層に入った分は土に残り、流れ去った分は帳簿の out に付ける。雨は塩を持たない。根の層の水が蒸発・蒸散・土の水の緩和で減っても塩は残り、濃くなる。雨雲の雨が根の層を満たして下へ染みるとき、染みる水は根の層の水と混ざり、その塩分の塩を帯水層へ運ぶ(洗い出し)。塩分が NaCl の飽和 264 g/kg を超えた分は結晶として土に残る(在庫には数え、濃さは 264 で頭打ち)。帯水層と根の層の間の拡散は、分子拡散(1e-9 m²/s)で数 m を渡るのに数十年かかるので入れていない(Freeze & Cherry 1979)。分散も入れていない。崖崩れは塩を運ばない(崩れたマスに残る)。海の塩(35 g/kg、Millero et al. 2008)は帳簿の外の溜め **[条件]**
- 涵養: 雨雲の雨が根の層を満たした残りだけが下へ染みる。地表の浸透の上限 13 mm/h(ローム、Rawls et al. 1982)を超える雨は地表を流れ去る(ホートン型)

### 7.5 天気: 対流の雨雲・突風前線・雷 [法則]+[条件]
- 発生: 1 km² あたり 1 分あたり 1.3e-6 × 日射の正弦 × 飽和水蒸気圧の 20 °C との比 **[条件]**(一地点が年に約 30 回雨雲の下に入る値。Changnon 1988)。区画と、一生のうちに流れ込みうる風上の範囲で起こす
- 移動: 中層の風 10 m/s で流れる。一生 60 分(Byers & Braham 1949)。強さは一生の中で sin の形に盛り上がって衰える
- 雨: 半径 2.5 km のガウス形。成熟期の中心の強さは雲ごとに指数分布、平均 30 mm/h × 飽和水蒸気圧の比(クラウジウス・クラペイロン、Lenderink & van Meijgaard 2008)、上限 150 mm/h
- 冷気: 半径 5 km、中心で −6 °C(Engerer et al. 2008)。雨より広く、流れてくる雲の前触れになる
- 突風前線: 冷気の密度流 c = 0.7·√(g'·h)、g' = g·ΔT/T、冷気の厚さ h = 1000 m(Benjamin 1968、Charba 1974)。冷気 6 °C で約 10 m/s。風は冷気の気圧の高まりの勾配に沿って外へ吹く。匂いを運ぶ
- 雷: 対地放電の回数の期待値 = その tick の雨の質量 ÷ 1e8 kg(rain yield、Petersen & Rutledge 1998)。落ちる場所は雨の芯と同じ広がり。長い連続電流を伴うもの 20%(Fuquay et al. 1972)が燃料に着火しうる
- 雲と雷の乱数は別の列にし、雷を足しても雲の列は変わらない

### 7.6 燃料の湿りと火 [法則]+[未換算]
- 燃料の湿り: Simard (1968) の平衡含水率に、時定数 60 分(1 時間燃料、Fosberg & Deeming 1971)で追従する。水蒸気圧は日平均の気温と相対湿度 0.70 で決まり、一日の中では変わらない。昼に気温が上がれば相対湿度が下がり、燃料が乾く。雨雲の冷気で下がった分は水蒸気圧を変えないので、冷えた分だけ相対湿度が上がる
- 消炎含水率 0.25(Anderson 1982 の中央値)。湿りによる燃焼の弱まりは Rothermel (1972) の η_M
- 酸素: 地球の 21% を 1.0 とし、17% 未満では広がらない。酸素が多いほど湿った燃料まで燃える
- 着火: 落雷(7.5)。燃え広がりは風下の隣から火勢を集める(率 2.17e-3 /分 **[未換算]**)、燃える割合 5.81e-3 /分 **[未換算]**、減衰 4.65e-3 /分 **[未換算]**
- 燃焼で死骸を灰にし、炭素は CO2 になる(O2 を 1:1 で使う)。焼けた組織の窒素は揮発する。炭として残る分(高さあたり 0.25)は死骸に戻す
- 燃料の量: 枯れた植物の死骸と、株の高さ × 1.2。**注意**: 燃料を kg/m² の乾燥質量で数えて燃え広がりの下限と比べる仕組みは未実装(D-2026-10-26-3 で作者の判断を経た「植物の尺度」の作業に含まれる)
- 再生のときは落雷と燃え広がりを止め、記録の各件だけを置く。置いた火は隣へ移らない

### 7.7 崖崩れ [法則]
- 起こり方: 無限斜面の安全率 FS = c'/(γ z sinα cosα) + (1 − m γw/γ) tanφ'/tanα が 1 を下回ったマスが、その tick に崩れる(Skempton & DeLory 1957、Montgomery & Dietrich 1994)。z = 崩れる厚さ 1 m、γ = 17 kN/m³、γw = 9.81 kN/m³、φ' = 安息角 35°、c' = 5 kPa(土 2 kPa と根の補強 3 kPa。Schmidt et al. 2001)。m は崩れる層のうち地下水で飽和した割合で、帯水層(地表から 3 m、7.4)の飽和した厚さが 2 m を超えた分。乾いた斜面は最も急な所(勾配約 1)でも FS > 1 で崩れず、飽和すると勾配約 0.76(37°)より急な斜面が崩れる。雨雲の雨が帯水層を満たすこと(7.4 の涵養)が引き金になる。乱数は引かない
- 以前は勾配 0.58 以上の斜面に勾配に比例する一定の率(隆起と釣り合う削られる速さ 5.71e-10 /マス/分)で起こしていた。雨や地下水と無関係だったので、D-2026-10-26-2 の答え (a) のとおり安全率に置き換えた
- 厚さ 1 m(Larsen et al. 2010)。崩土は最も急な下りを伝い、出発点からの見通し角が H/L を下回る所で止まる。H/L は log10(H/L) = −0.085·log10(V) + 0.047(Corominas 1996)で、1 マスの崩壊で約 0.75
- 崩れた跡は 1 m 低くなり傾斜が緩む。表土の窒素と土砂は到達点へ。植物は死骸として運ばれ、巻き込まれた動物は埋まって死ぬ
- **注意**: 今は雨からの涵養がほぼ 0 で(7.3 の注意、D-2026-10-28-2)、陸の帯水層はほとんど空なので、崩れはほぼ起きない。土の厚さは区画のどこでも 1 m で、根の補強は植物の量によらない
- 地震の揺れによる崩れは 7.7.1 の仕組みで起き、崩れそのもの(到達・運ばれる物・地形の変化)はこの節と同じ

### 7.7.1 地震の揺れ [法則]+[条件]
- **断層** [条件]: 区画の外、北の端から 10 km 北を東西に走る長さ 60 km の断層(区画の中央の北に中心)。北の端を頂にする傾動地塊を持ち上げる断層とみなし、すべりの速さは隆起の速さ 0.3 mm/年と同じ。地震発生層の幅 15 km、剛性率 3e10 Pa
- **規模の分布** [法則]: グーテンベルク・リヒター則 log10 N(≥M) = a − bM、b = 1(Gutenberg & Richter 1944)。最小 M4(これより小さい揺れは区画の斜面の限界加速度に届かない)、最大は断層の長さから M = 5.08 + 1.16 log10(60) = 7.14(Wells & Coppersmith 1994)。a は断層に蓄えられるモーメントの速さ μLWs = 8.1e15 N·m/年(Brune 1968)と、打ち切った G-R 則が解放するモーメント b/(1.5−b)·10^(a+9.1+(1.5−b)Mmax)(M0 = 10^(1.5M+9.1) N·m、Hanks & Kanamori 1979、Youngs & Coppersmith 1985)を等しく置いて決める。a ≈ 2.94、M4 以上の本震は年 0.086 回、M6 以上は年 0.00086 回
- **本震** [法則]: 1 分あたりの率 × DT の期待値で起こす(端数は確率)。震央は断層に沿って一様。規模は打ち切った G-R 則から引く
- **余震** [法則]: 本震の規模 Mm の後、時刻 t(日)に規模 M 以上が起きる率 10^(a'+b'(Mm−M))(t+c)^−p(改良大森則、Reasenberg & Jones 1989 の一般的な値 a' = −1.67、b' = 0.91、p = 1.08、c = 0.05 日)。各 tick の期待値は、その tick の区間でこの式を積分した値。震央は本震の破壊の範囲に一様、規模は b' の G-R 則(M4〜7.14)から引く。余震はさらに余震を生まない(式が二次の余震を含めた列全体への当てはめのため)。M6.8 の後の 1 日に M4 以上はおよそ 26 回
- **破壊の長さ** [経験式]: log10 L(km) = −2.44 + 0.59 M(Wells & Coppersmith 1994、地下の破壊の長さ、全種類)。震央を中心に断層に沿って置き、断層の端で切る
- **揺れ** [経験式]: 最大加速度(g)log10 A = −1.02 + 0.249 M − log10 r − 0.00255 r、r = √(d² + 7.3²) km、d は破壊の地表への投影までの水平距離(Joyner & Boore 1981)。ばらつきは入れず中央の値。区画の最も近いマスで d = 10 km なので、M6.8 で 0.35 g、M4 で 0.07 g
- **斜面の限界加速度** [法則]: 7.7 の無限斜面の安全率 FS(地下水の間隙水圧を含む)から、限界加速度 ac = (FS − 1) g sinα(Newmark 1965)。z = 崩れる厚さ 1 m、γ = 17 kN/m³、φ' = 地形を作るときの安息角 35°(粘着の無い粒の安息角は内部摩擦角に等しい)、c' = 5 kPa(土 2 kPa と根の補強 3 kPa。Montgomery & Dietrich 1994、Schmidt et al. 2001)。この値で、乾いた斜面は区画の最も急な所(勾配の大きさ約 1)でも FS > 1(ac ≈ 0.2 g)で、揺れが無ければ崩れない
- **崩れ** [経験式]: A > ac のマスで、ニューマークの変位 log10 Dn(cm) = 0.215 + log10[(1 − ac/A)^2.341 (ac/A)^−1.438](Jibson 2007)を求め、崩れる確率 P = 0.335[1 − exp(−0.048 Dn^1.565)](Jibson et al. 2000)で崩す。A ≤ ac のマスでは乱数を引かない。崩れは 7.7 の slide_at と同じ
- 乱数は地震だけの列を使う。地震が起きなければ、生き物と環境の結果は地震の無い版とバイト単位で同じ
- 揺れは今は崩れを起こすだけで、生き物の感覚の入力にはなっていない(神経網の入力の変更は、明るさ・登る・潜る・前触れの感覚とまとめて一度に行う)
- **注意**: 根の補強は植物の量によらない定数。地下水で湿った斜面ほど限界加速度が下がり、小さな揺れで崩れる

### 7.8 大気の気体 [法則]+[未換算]
- 量 **[条件]**: 各マスの上の空気の柱の全量 = 地上の気圧 101325 Pa ÷ 重力 9.81 m/s² ÷ 空気のモル質量 0.028965 kg/mol = 3.57e5 mol/m²。体積の割合は CO2 420 ppm、O2 20.95%、N2 と Ar 79.01%(地球の今の値。惑星規模の条件)。気体の単位はエネルギー 1(= 6 mg C)を酸化して出る CO2 の量 × 1/0.02 = 0.3 g C = 0.025 mol。1 マス(100 m²)の柱は CO2 6.0e5、O2 3.0e8、N2 と Ar 1.13e9 単位
- 場が持つのは柱の量からのずれだけ。マスにある量 = 柱の量 + ずれ。生き物の出入り(1 分に 1e-4 単位ほど)は柱の量の 1e-10 ほどなので、柱の量をそのまま float に持つと丸めで消えるため。ずれへの端数は持ち越して足す
- 酸素の分圧 = O2 ÷(CO2 + O2 + N2 と Ar)。動物の移動と活動の酸素の充足(7.11 の「移動」)と、火の広がりの下限(17%)に使う。以前は CO2 と O2 を 8、N2 を 24 と置いた比(20%)だった
- CO2 と O2 は 4 tick ごとに、風 3 m/s(1 分に 18 マス)のセミラグランジュ法で移流し、拡散する。風上の端からは柱の量(ずれ 0)の空気が入り、風下の端から出る **[条件]**。端の出入りは帳簿に付ける
- 光合成で CO2 を取り O2 を出す(1:1)。呼吸・分解・燃焼は O2 を使い CO2 を出す(呼吸商 1)
- 標高による空気の薄まり: スケールハイト 8.4 km

### 7.9 匂い [法則]
- 4 本: 生体(植物の体臭は高さと毒に比例、動物は体格に比例。繁殖できる状態の動物の性フェロモンも 0 番に入る)、腐敗(死骸)、食害(齧られた植物)、災害(煙)
- 生体・腐敗・食害の匂いは OH ラジカルとオゾンに酸化されて消える。率 = (k_OH·[OH] + k_O3·[O3])·60 /分。代表の物質 α-ピネン・二硫化ジメチル・(Z)-3-ヘキセノールの 298 K の速度定数(Atkinson & Arey 2003)。[OH] は太陽高度の正弦に比例し真上で 4e6 個/cm³(Prinn et al. 2005)、[O3] は 30 ppb で一定 **[条件]**
- 煙は粒子で、沈着で消える(率 1.37e-4 /分 **[未換算]**)
- 運ぶ風: 地表から 0.1 m までの層の平均の速さ。対数分布 u(z) = u*/κ·ln(z/z0)、z0 = 0.01 m(Wieringa 1993)から約 0.68 m/s。突風前線の風を重ねる
- 横への広がり: 渦拡散 K = κ·u*·z(Prandtl)
- 広がる流れ(∇·u > 0)では上から匂いの無い空気が降りてきて、濃さは exp(−∇·u·Δt) に薄まる

### 7.10 植物 [法則]+[未換算]
- **光**: マスの株を実効的な高さの順に並べ、上から日射を受け取る。一株が受けられる量は飽和(0.25 + 0.75·(1 − 耐陰性))まで。水中では深さで減衰(浮く株は減衰が小さい)。氷の下は 8%
- **光合成**: 取り込み = 受けた光 × 吸収 × 気温の適合 exp(−((T − T_opt)/18)²) × 水の適合 × CO2 の足り具合(リービッヒの最小律)。T_opt = 26 − 26·耐寒性
- **水**: 必要な水 0.55·(1 − 0.75·耐乾性)。土の水が足りなければ体の貯水で補う。蒸散で土の水が減る
- **成長**: 高さを最大の高さへ伸ばす(率 1.94e-4 /分 **[未換算]**)。建設に高さあたり 3.0 の炭素と、その半分の建設の呼吸(Penning de Vries 1974)。窒素が足りなければ成長が律速される(体を作るのに窒素 0.055/高さ)
- **維持の呼吸**: 葉の暗呼吸(1.17e-3 /分 × 飽和 **[未換算]**)+ 茎と根(体の炭素あたり 1%/日、McCree 1974・Ryan 1991・Amthor 2000)+ 各形質の維持費 **[未換算]**。暗所では 0.35 倍、休眠中は 0.12 倍、酸欠では 1.3 倍。払った炭素は CO2 になる
- **窒素**: 固定できる株は大気から(費用を払う)、できない株は土から吸う。同じマスの株で分け合う
- **塩と水に浸かること** [法則]+[未換算]: 生息環境の不適合を二つに分ける。(1) 水に浸かることへの不適合 hm: 海では (1 − aquatic) × (0.30 + 0.70 × 水深の比)、陸では aquatic(通気組織と浮く体の代償)。(2) 塩への不適合 hs = 根の層の塩分 / 35 − salt_tol。まわりの水が塩への強さより塩辛ければ水を吸えず(生理的な乾燥)、薄ければ溜めた適合溶質と輸送の仕組みが無駄になる(塩生植物の多くは淡水で育ちが落ちる。Flowers & Colmer 2008)。維持費に 0.0034 /分 × (hm² + hs²) **[未換算]** を足す。最初の海の株は salt_tol 1、陸の株は 0。動物の aquatic はまだ一つのまま(水中での呼吸と泳ぎ・塩への強さに分けるのは次の段階)
- **防御** [法則]+[未換算]: 防御の目標 = tough × (def_const + (1 − def_const) × 脅威)。脅威 = min(1, 食害の匂い / 1.5)。今張っている防御は目標へ率 8.05e-4 /分で近づき、6.55e-5 /分で緩む **[未換算]**。維持費 6.81e-4 /分 × 今張っている防御 **[未換算]**。遺伝子 def_const(0〜1)は、常に張る分(構成的防御)の割合。常に張れば食われ始めから守れるが費用が続き、匂いで張れば費用は食害のあるときだけで済むが、立ち上がるまで食われる(最適防御の理論。Karban & Baldwin 1997)。どちらの生き方が残るかは決め打ちせず進化に任せる。最初の株と区画の外から入る株は 0.5 **[条件]**。古いセーブは 0(匂いで張るだけ。この変更の前の規則と同じ)で読む
- **休眠**: 気温が休眠の気温を下回ると休眠し、3 °C の余裕を超えて上がると覚める
- **つる** [法則]+[未換算]: 遺伝子 climb(0〜1)。マスの中で最も高い株を支えにする。実効的な高さ = 自分の高さ × (1 − 0.45·climb) + 支えの余り × climb × 0.85(支えの余り = max(0, そのマスで最も高い株 − 自分の高さ))。自立する幹を作らないので、支えの無いマスでは自分の高さより低くなる。つるの実効的な高さも林冠になり、後に見る株の支えになる(マスの株を見る順に依存する)。維持費 2.14e-4 /分 × climb **[未換算]**(巻きひげと気根)
- **寄生** [法則]+[未換算]: 遺伝子 parasite(0〜1)。4 分に一度、2 株以上いるマスで、最もエネルギーの多い株を宿主に選ぶ(寄生者も宿主になりうる)。parasite > 0 の他の株はすべて吸器を刺し、宿主のエネルギー × (1 − exp(−5.5e-4 /分 × parasite × 刺さり具合 × 4 分)) を奪う。吸器が宿主の維管束とつながる面の広さは吸器への投資に比例するので、奪う率に parasite が掛かる(Press & Phoenix 2005)。刺さり具合は宿主の防御との競り合いで、比なので投資の大きさを表さない。前の規則の「parasite ≥ 0.10 の株だけが刺す」「0.30 を超える株は宿主にしない」の閾値は外した。刺さり具合 = parasite / (parasite + 宿主の防御)。奪った量の 85% を自分のエネルギーにし、15% は呼吸で CO2 になる(O2 を 1:1 で使う)。窒素も、宿主の貯え・奪った量の 2%・自分の空き の最も小さい分を移す。同じ宿主に複数の寄生者がいれば、順に奪う。維持費 2.92e-4 /分 × parasite **[未換算]**
- **地下茎(クローン)** [法則]+[未換算]: 遺伝子 clonal(0〜1)。子を産むとき、確率 clonal で親の近く(1.5 マス)に置き、親と同じ株(genet)の一員にする。それ以外は 1.5 + 12·disperse マス先に散らし、新しい株にする。10 分に一度、clonal > 0 の株は地下茎を通してエネルギーを融通する(閾値は置かない。clonal が小さければ通りやすさ a も小さく、融通は連続に小さくなる。地下茎の師管の断面が投資に比例するため。Pitelka & Ashmun 1985、Alpert 1996)。各株は通りやすさ a = 1 − exp(−2.17e-3 /分 × clonal × 10 分) で共通の節につながり(clonal が率に掛かる)、節の値は a で重みを付けた平均 Σa·E / Σa になる(節は蓄えを持たない。キルヒホッフの電流則)。各株は E += a·(節 − E) だけ節へ寄る。Σ a·(節 − E) = 0 なので、同じ genet の中で clonal が違っても(子は変異する)融通の前後で genet のエネルギーの合計が保たれる。維持費 9.72e-5 /分 × clonal **[未換算]**。土中で待つ種子も同じ genet に数えられる
- **火後発芽** [法則]: 遺伝子 fire_seed(0〜1)。種子ごとに確率 fire_seed で、焼け跡でしか芽吹かない種子になる(球果や果実の閉じ方の度合い、serotiny は連続で、部分的に閉じる株が多い。Lamont et al. 1991)。そういう種子は、落ちたマスが焼けてから 64.3 日(BURN_WINDOW)より後なら芽吹かずに朽ちる(炭素は死骸へ)。焼けて間もないマスにだけ芽生える。残りの種子は普通に芽吹く。どちらに賭けるかの割合は決め打ちせず進化に任せる(前の規則の「0.35 を超える株はすべての種子が焼け跡を待ち、以下はすべて普通に芽吹く」は外した)。維持費 1.17e-4 /分 × fire_seed **[未換算]**。燃えたことは煙と焼け跡の匂い(7.9)で感じ取れる
- **耐火(厚い樹皮)** [未換算]: 遺伝子 fire_tol。炎の強さ × (1 − 0.65·fire_tol) の害を受け、エネルギーを失い、高さが縮む(下限 0.03)。焼けた組織の炭素は CO2 に、窒素は揮発する(帳簿の fire)。焼け死んだ株は高さ × 0.25 を炭として死骸に残し、窒素も同じ割合で残す。維持費 7.29e-4 /分 × fire_tol **[未換算]**
- **土中の種子** [法則]+[未換算]: 遺伝子 seed_bank(0〜1)。種子は確率 seed_bank で土中に入り、2.86〜100 日(一様)待つ。待つ間は光合成も成長もせず、寒さや火でも死なない。8 分に一度まとめて、呼吸 3.4e-5 /分(CO2 になる)を払い、種子食者に見つかる率 1.17e-5 /分 **[未換算]** で死ぬ(1 回の判定の確率 1 − exp(−1.17e-5 /分 × 8 分))(体は死骸へ。種子食者は動物の個体ではなく、率として与えている)。待ち終えたとき土の水が 0.25 未満なら、3.43 日後にもう一度土を見る。土中の種子も上限 PLANT_MAX の数に入る
- **食虫** [法則]+[未換算]: 遺伝子 carnivory(0〜1)。0.12 以上の株のいるマスで、潜っていない(深さ 0.2 未満の)動物は、動物の更新ごとに確率 carnivory × 7.33e-4 /分(を更新の間隔に換算) ÷ (1 + 1.5·体) で捕らえられる。最初に捕らえた株が、動物のエネルギーの 55% と脂肪の 40% を奪う。奪った量の 70% が株のエネルギーに、30% が消化の呼吸で CO2 になる。株は奪った量 × 0.09 の窒素を得る(貯えの上限まで)。**注意**: 動物は窒素の在庫を持たないので、この窒素は帳簿の外から入る量として帳簿の carn に付けている。捕らえられてエネルギーが尽きた動物は死ぬ。維持費 4.37e-4 /分 × carnivory **[未換算]**
- 時間の換算: 寄生の奪う割合・地下茎の融通の割合・土中の種子の死ぬ確率は、1 分あたりの率 r と判定の間隔 n(4・10・8 分)から指数の形 1 − exp(−r·n·DT) で決める(「守り続けること」の換算の規則)
- 寄生・食虫・つる・地下茎の遺伝子の閾値(0.35・0.30・0.45・0.50)は、画面の分類と要約にだけ使い、規則の入力には使わない。規則の中の閾値は上に書いたもの(食虫 0.12)。地下茎は閾値を持たず、clonal に比例する通りやすさで融通する。寄生は閾値を持たず、投資 parasite に比例して奪う。火後発芽は閾値を持たず、種子ごとの確率 fire_seed で働く
- **繁殖**: エネルギーが 8.0 を超え、休眠していなければ判定する。有性なら同じ lin の相手を探す(2.9 の注意)。種子の数 = 1 + brood·11。投資 = エネルギー × (0.20 + 0.55·繁殖に回す割合)。多産ほど一つあたりの構造に余分(0.22)がかかり、余分は呼吸で CO2 になる。芽生えの体は種子のエネルギーと親の窒素で作り、足りなければ朽ちる
- **遺伝**: 有性なら各形質と各免疫の座位をどちらかの親から半々で受け取り、そのあと変異(7.13)を選んだ値に掛ける(変異率は母から受け取る)。種子の散らばり方(clonal・disperse・seed_bank・fire_seed)は母の体の性質として母の値で決める
- **死**: エネルギーが尽きる、寿命、寒さ、火。体の炭素と窒素はすべて死骸に
- **漂流**: 浮く株は水面の漂流(風の 3%、Wu 1983)× 浮く強さで流される。東西の端を越えた株は区画の外へ出る(帳簿に付ける)
- 上限 PLANT_MAX = 260000(計算の安全のため)。土の中の種子も数に入る

### 7.11 動物 [法則]+[未換算]
- **感覚**: 感覚器 5 つ。それぞれ自分の向きに射程(最大 14 マス)だけ 3 段で探り、反応する帯(13: 植物・小さい動物・大きい動物・死骸・炎・林冠・病原体・嗅覚 4・似た体・水際)で最初に見つけた信号 × (1 − 距離/射程)を入力にする。射程が 0.8 マス未満なら器官として成立しない。帯 11(似た体)[法則]: そのマスにいる他の動物(先頭 6 体まで)それぞれについて、免疫の型 im(3 座位、各 0〜1)の似ている度合い 1 − Σ|im_自分 − im_相手|/3 を足し、× 0.45 を 1 で頭打ちにした値を信号とする。相手の体の表面の化学的な手がかりを自分の型と照らし合わせて仲間を見分ける表現型の照合(Lacy & Sherman 1983)。観察の分類 lin は使わない
- **入力**: 感覚 5 + エネルギー(繁殖の閾値で割る)・気温・土の水・体格・バイアス + 再帰 2 = 12
- **神経網**: brain_forward。隠れ層は最大 3 層 × 8。各ニューロンの有効・無効の値が 0.35 を超えたものだけが働く。tanh。入力から出力への直結もある。出力 6 = 旋回・前進・採食・攻撃・再帰 2。出力が 0.15 を超えた行為だけを行う。試験器(--assay)も同じ計算を使う
- **成長**: 体を size へ伸ばす(率 4.86e-5 /分 **[未換算]**)。体 1 を作るのに 4.0、組織に残るのは 2.5、残りは呼吸で CO2
- **登り降り**: 林冠の高さと、枝が支えられる体重から届く高さを決め、arboreal の分だけ登る。登る費用 0.55/高さ。支えを失って滑り落ちるときも、樹上を渡る先の林冠が 0.45 より低く途切れて下の林冠(無ければ地面)まで落ちるときも、同じ落下の規則に従う: 高さ 0.5 を超える落下は、落ちた高さ × 0.9 × (1 − 0.7·arboreal) を傷の修復のエネルギーとして使う(呼吸で CO2 に)。以前は、林冠が途切れて落ちるときだけ高さによらず 0.04 × (1 − 0.6·arboreal) だった。**注意**: 落下の怪我を重力と終端速度から決める(追加の例「落下」)には、植物の高さの単位をメートルに結びつける必要がある(D-2026-10-30-1 の植物の尺度の作業)。**注意**: 登るかどうかを神経網に通す作業が「上の階層の書き込みを外す」に残っている(今は遺伝子の分だけ自動で登る)
- **移動**: 速さ = 遺伝子 × A_MAX_SPEED × 前進の出力 × 地形(陸・水・樹上、勾配で 1/(1 + 0.5·勾配))× 酸素の充足。水と陸の境に通れない閾値は置かない(以前は aquatic < 0.20 の個体は水に入れず、aquatic > 0.85 の個体は陸と海氷のある海に入れなかった)。水中の速さ 0.25 + 1.30·aquatic、陸の速さ 1 − 0.75·aquatic と、水に浸かることへの不適合の費用(維持費の hm²)がどちらも aquatic に連続に従う。水に落ちた陸の昆虫も水面で遅く漕いで岸に戻ることがあり、水生の幼虫も陸を這い出せる(泳ぐ能力と陸を歩く能力は連続に変わる)。入るかどうかは神経網の向きと、この速さと費用の結果で決まる。**注意**: A_MAX_SPEED は 1 分あたり 0.000972 マス **[未換算]** で、速さの遺伝子 1 でも約 0.016 cm/s。条件の 0.5〜10 cm/s から大きくずれている(優先順位 1 の作業)
- **気温による死**: 感じる気温(潜っていれば年平均へ 0.75 × 潜り具合だけ寄る)が 8 − 62·耐寒性 を下回ると死ぬ
- **活動**: 気温の適合 × 酸素 × 昼夜(0.35 + 0.65·(夜目 × 夜 + (1 − 夜目)× 昼))。深く潜っていれば 0.15 倍。**注意**: 昼に動くか夜に動くかを神経網に決めさせる作業が残っている
- **採食**: 植物・分解物・肉の三つの消化能力の合計に上限がある(植物の消化 = 1 − (diet + 0.5·detritus)²)。一回の更新に齧れる株の数は体格で決まる(最大 6)。一口 = 株のエネルギー × 0.45 × 一口の大きさ × 噛み砕き(顎/(顎 + 2·防御))× 活動。同化効率 50% を基準に消化効率と一口の大きさで変わる。同化しなかった分は糞、落とした葉は落葉(死骸)
- **漁る**: 死骸を一回に 0.85 × 能力だけ取る。植物由来は 20%、動物由来は 70% を同化
- **襲う**: 隣の 3×3 で、自分の体が相手の 1.15 倍を超え、高さが近く、穴に潜っていない相手。一口 = 相手のエネルギー × 率 9.22e-3 /分 × √肉食 × 活動。毒は解毒できなければ同化を減らし、解毒の費用がかかる
- **食虫植物**: 地表にいれば、食虫の株に捕らえられうる
- **維持費**: 基礎代謝 1.8e-5 /分 × 組織の炭素(Chown & Nicolson 2004、Addo-Bediako et al. 2002)×(1 + 形質の維持費の和 / 基準)+ 歩く費用。老化で上がる。病原体に型が合えば上がる
- **塩と水に浸かること** [法則]+[未換算]: 生息環境の不適合を二つに分ける。(1) 水中での呼吸と泳ぎ aquatic: 海では (1 − aquatic) × (0.30 + 0.70 × 水深の比)、陸では aquatic を不適合とする(今までの式のまま)。(2) 塩への不適合 hs = 今いるマスの水の塩分(陸は根の層の土の水、海は 35 g/kg)/ 35 − salt_tol。外の水が体液より塩辛ければ塩を汲み出す輸送(直腸の塩の分泌など)、薄ければ溶質を取り戻す輸送に費用がかかる(Bradley 1987)。維持費の形質の和に 0.000486 /分 × hs² **[未換算]**(水に浸かることへの不適合と同じ係数)を足す。最初に放す動物は salt_tol 1。外から入れ直す動物の置き場所は、塩への不適合 0.55 を超えるマスを除く **[条件]**
- **呼吸**: 必要な O2 = 代謝 × 0.02 × (1 − 0.6·呼吸の効率)。マスの O2 で足りない分は発酵で、炭素は有機物のまま糞に出る
- **脂肪** [法則]+[条件の幅]: 体の維持費(基礎代謝と形質の維持費。歩く費用は含まない)の 1 日分を day_cost とする。エネルギーが糖の貯え GLYCO_DAYS = 1 日分(絶食した昆虫は先に糖を約 1 日で使い切る。Arrese & Soulages 2010)を超えれば、超えた分 × 6.21e-4 /分(を更新の間隔に換算)を、脂肪が目標 store_fat × 60 日分(FAT_DAYS_MAX。昆虫の絶食に耐える日数は数日〜数か月。Rion & Kawecki 2007、Hahn & Denlinger 2011)に届くまで脂肪へ回す(合成で 15% を呼吸)。エネルギーが 1 日分を下回れば、脂肪を 1 回の更新で最大 FAT_DRAW 取り崩す。繁殖の閾値は使わない(D-2026-10-24-3 の作者の答え)。何日分貯めるかは遺伝子 store_fat で進化する
- **繁殖**: エネルギーが A_DIVIDE_TH × (1.4 − 0.6·繁殖に回す割合)(A_DIVIDE_TH = 6.0)を超え、体が size の 80% 以上なら判定する。有性なら半径 3 マスで同じ lin の成熟した相手を探す(2.9 の注意)。産む数 = 1 + brood·11、投資 = エネルギー × (0.20 + 0.55·繁殖に回す割合)、一匹の取り分 = 1/(数 × (1 + (数 − 1)·0.22))。子の体は 0.10(作る費用 0.4)で、取り分がそれに届かなければ死骸になる。**注意**: 繁殖の閾値を組織の炭素に比例させ、子の体を「投資 ÷ 産む数」から決める作業(移行の段階 5 の作者の判断)は、優先順位 1 で試験中
- **性フェロモン** [法則]+[条件に近い推定]: 体が size の 80% 以上で、エネルギーが A_DIVIDE_TH × 0.7 を超える個体は、匂いの場の 0 番に 1.34e-4 /分 × estrus **[未換算]** を出す。遺伝子 estrus(0〜1)が出す量を決め、出している間は合成の費用 3.37e-7 /分 × estrus × 体 を代謝に足す(費用の炭素は CO2 として出す)。費用は、estrus = 1 で体 1(生重 100 mg)が 100 ng/時を出すとして(蛾の雌の放出は 1〜100 ng/時の桁、Cardé & Baker 1984。その上端)、脂質 1 g の合成にブドウ糖 3.03 g(Penning de Vries et al. 1974)から決めた。基礎代謝の約 0.75%。匂いは捕食者にも届くので、多く出すほど相手にも捕食者にも見つかりやすい。出す量は決め打ちせず進化に任せる。最初に放す動物は 0.2 **[条件]**(前の規則で sexual 0.2 の個体が出していた量)。古いセーブも 0.2 で読む。前の規則の「sexual が 0.15 を超える個体だけが sexual に比例して出す」は外した。相手の探し方は今のまま(半径 3 マスの探索)で、匂いは神経網の入力を通してだけ働く
- **遺伝**: 有性なら各形質・各重み・各感覚器・各免疫の座位をどちらかの親から受け取る。そのあと変異(7.13)を、選んだ値に掛ける(一口・消化・変異率は母から受け取る)
- 上限 ANIMAL_MAX = 60000(計算の安全のため)。上限に達したときに配列の前(古く生まれた個体)から枠が埋まる偏りは D-2026-10-28-3 で作者の判断を待っている

### 7.12 微生物と死骸 [法則]
- 分解: 一次反応。マスの死骸のすべてに同じ割合で取りつく。率 = 1.02e-3 /分 **[未換算]** × 分解の速さ × 活動(気温)。取り込みの上限 4.17e-3 /分(0.25 /時)× 体の炭素(最大比増殖速度 μmax。基質を加えた呼吸の解析で 0.1〜0.4 /時、Blagodatsky et al. 2000、Panikov 1995)
- 炭素利用効率 0.40(Manzoni et al. 2012、Sinsabaugh et al. 2013): 分解した炭素の 40% が体に、残りは CO2。分解した割合だけ窒素が土へ無機化する
- 維持の呼吸 1.7e-5 /分 × 体の炭素 × 酵素の係数(Anderson & Domsch 1985、1993)
- 入れ替わり: 基質のあるマスにいる個体は 0.05 /日 で死ぬ(Rousk & Bååth 2011)。死んだ細胞の炭素は CO2 へ
- 胞子の定着: 微生物のいない死骸のマスに、平均 173 時間で定着する(胞子 7.5e-8 mg C から群落 2.4 mg C まで、最大比増殖速度 0.1 /時で育つ時間)
- 胞子の散布: 平均 40 m の指数分布(高さ 5 cm × 風 0.8 m/s ÷ 終端速度 1 mm/s。Gregory 1973)
- 分解者が 50 個体未満になり死骸が 1000 を超えると、外から飛んできた胞子として 200 個体を死骸の上に置く(持ち込む炭素は帳簿の移入に付ける)**[条件]**
- 死骸の風化(光分解): 0.3 /年(Austin & Vivanco 2006)。炭素は CO2、窒素は土へ(海では 45% だけ有光層に残り、残りは沈む)
- 炭素 0.02 を下回った死骸は消し、残りを CO2 と土の窒素に出す

### 7.13 変異 [条件(加速)]
- 子の各遺伝子に、範囲の幅 × 一様乱数 ±0.02(MUT_STRENGTH)を足す。変異率の遺伝子を持つ形質(変異率・一口・消化)は、さらに 0.4 + 変異率 × 3 倍
- ニューロンの有効・無効は 1.5% の確率で大きく飛ぶ(層の追加・削除)。感覚器の反応先は 2% の確率で大きく変わる
- 進化の速さは実際より大幅に速い。変異の強さによる意図した加速で、時間の換算で遅くしない

### 7.14 病原体 [法則(負の頻度依存)]+[未換算]
- **状態**: 植物の病原体と動物の病原体を別々の場として持つ。マスごとに型(3 座位 pstrain/astrain、各 0〜1、始まりは 0.5)と量(pload/aload、0〜2.5)。量は個体の数でも炭素でもない抽象的な濃さで、帳簿(炭素・窒素・水)には入らない
- **宿主の型**: 植物と動物はそれぞれ免疫の型 im(3 座位、各 0〜1)を持つ。始まりは一様乱数。子の型は、有性生殖では座位ごとに父母のどちらかから半々の確率で受け継ぎ、無性生殖では親から受け継ぐ。受け継いだ値に変異(7.13)を足す
- **更新**(3 tick ごと): マスの中の生きた宿主の数 n と、型の平均 m_k を数える
  - 宿主がいるマス: 型が平均へ寄る s_k ← s_k + (m_k − s_k)·(1 − exp(−9.87e-5·3))。量は L ← L + 6.48e-5·3·min(6, n) − L·(1 − exp(−4.9e-5·3))
  - 宿主のいないマス: L ← L − 3·L·(1 − exp(−4.9e-5·3))(3 倍の速さで衰える)
  - 量は 0〜2.5 に切る
  - 定数はいずれも段階 3 より前の 1 回の更新あたりの値(寄る 3%・増える 0.02・衰える 1.5%)を換算したもの **[未換算]**
- **発病**: 量が 0.02 を超えるマスにいる個体は、型の近さ match = exp(−(d/0.42)²)(d = 座位の差の二乗平均の平方根)を掛けた感染 inf = L·match を受け、維持費に inf × 5.35e-4 /分 **[未換算]** を足す。この分のエネルギーは他の維持費と同じく呼吸で CO2 になるので、保存則は宿主の側で保たれる
- **何が起きるか**: 病原体はそのマスで多い型に寄るので、多数派の型の個体ほど病みやすく、少数派の型が有利になる(負の頻度依存選択。Haldane 1949、Hamilton 1980 の赤の女王)。病原体の場はマスに固定で、移動や拡散はしない(宿主の移動と子の分散だけが型の空間の模様を動かす)
- **匂いと感覚**: 病原体の量が多いマスでは「病んだ個体」の匂いが出る(植物は量 0.5 超で 1.09e-4 × 0.4 /分、動物は量 0.3 超で 1.09e-4 /分。7.9 の匂いの場の 0 番)。動物の感覚器の帯 6 は、そのマスの動物の病原体の量 / 2.5 を信号として受ける
- **注意**: 量の増え方は宿主の数で決まり、宿主の体の大きさや接触の頻度は使っていない。濃さに単位が無いので、実際の感染率や致死率と比べられない(段階 5 の対象)

### 7.15 記録と再生
- --record: 生き物と無関係に起きる環境の出来事(火の着火・崖崩れ・雨雲・雷・地震)を tick・種類・位置・強さで書き出す。地震の位置は区画の外の震央、強さは規模
- --replay(--shuffle-seed で混ぜる): 記録の各件を、混ぜた時刻と場所で起こす。再生した種類の出来事は物理の仕組みからは起こさず、再生した出来事から連鎖を起こさない。生き物の仕組みと感覚の物理は通常通り動く
- 記録を付けても付けなくても、シミュレーションの結果は変わらない

### 7.16 測定の道具 [観察。世界の規則ではない]
どの道具も乱数を使わず、世界の状態に触れず、標準出力にも書かない。だから道具を付けても付けなくても、シミュレーションの結果とログは変わらない。道具が書き出す量(lin・選択勾配など)は、生き物や環境の規則の入力には使わない。

- **神経網の書き出し**(--dump-brains 間隔、--dump-max 上限 1000): 間隔ごとに、生きている動物を生まれた順に等間隔で上限まで選び、id・世代・lin・体・大きさ・感覚器(向き・射程・同調)・節の利得・重みを brains_t########.txt に書く。浮動小数は 16 進で書き、正確に読み戻せる
- **試験器**(--assay ファイル、または --assay-save セーブ。--assay-max・--assay-out): 世界を作らず、神経網だけを、世界の動物と同じ brain_forward で動かす。1 段 = 動物の 1 回の更新(ANIMAL_STRIDE = 2 tick = 2 分)
  - 刺激: 捕食者(感覚の帯 2)・火(4)・煙(10)を脅威、植物(0)・死骸(3)を餌、獲物(1)・仲間(11)をその他とする。刺激は向きによらず、帯の合う感覚器すべてに 強さ × 2/3 で入る(世界で最初のレイの段で見つけたときの入力に合わせる)
  - 他の入力は固定: エネルギー 0.5(繁殖の閾値の半分)、気温 15 °C、土の水 0.4、体の大きさ、バイアス 1
  - 手順: 刺激なしで 60 段落ち着かせる → 刺激を 20 段 → 止めて 120 段。刺激なしの試験(偽の試験)でも再帰回路が 0.02 より大きく動く神経網は「落ち着かない」として除く
  - **価**: 反応した脅威の再帰回路の変化の和と、餌の変化の和の余弦が −0.5 より小さい(逆向き)
  - **持続**: 刺激を止めた後、変化が刺激中の 10% より大きく残る段数が 3 段以上
  - **強度**: 最も強く反応した脅威を強さ 0.25・0.5・0.75・1 で与え、変化の大きさとの相関が 0.9 より大きい
  - **一般化**: 違う脅威どうしの変化の向きの余弦の平均が 0.5 より大きく(入口)、かつ刺激を止めた後の 10 段で、再帰回路を通して行動の出力(旋回・前進・採食・攻撃)のうち 2 つ以上が平均 0.05 より大きくずれる(出口)
  - 根拠: Anderson & Adolphs (2014) の、感情に似た状態の四つの性質
- **陽性対照**: 手で組んだ神経網。感覚器を帯 2・4・10・0・3 に置き、脅威で再帰回路 0 を +3、餌で −3、自分への結合 0.9 でゆっくり戻り、その値が旋回 +1・前進 +2・採食 −2 を変える。試験器で四つすべてが検出されることを確かめる
- **陰性対照**: 同じ神経網の再帰の入力を 0 に切ったもの(--assay では各個体について切ったものも測る)。四つが検出されないことを確かめる
- **野外の対照**(--positive-control N・--no-recurrent N): 最初に放す動物の先頭から、乱数を使わずに N 匹を陽性対照の神経網と感覚器に置き換え、続く N 匹の再帰の入力を切る。切ったことは遺伝子ではなく、子にそのまま写る
- **個体の追跡**(--track N、--track-from T): t ≥ T に生まれた(放された)最初の N 匹を、生まれてから死ぬまで毎回の更新で track.csv に書く(位置・体力・感覚器の入力・出力・再帰回路)
- **きっかけの前後の記録**(--trigger-log W、--trigger-max M): 脅威(帯 2・4・10)か餌(帯 0・3)の感覚器の入力が 0 から立ち上がったとき、その前 W 回と後 W 回の更新を trigger.csv に書く。記録中に次のきっかけが来ても延ばさない
- **全体の記録**(--census N): N tick ごとに census.csv へ 1 行 1 値(t,group,name,stat,value)で書く。個体数(植物・種子・動物・微生物・死骸・系統の数)、繁殖の回数、死因ごとの累計(飢え・寒さ・捕食・微生物・崖崩れ)、遺伝子ごとの分布(平均・標準偏差・10/50/90 百分位)、選択勾配(前の記録から繁殖した親の遺伝子の平均 − 集団の平均、標準偏差で割った値も)、環境の出来事の件数と強さの和(--record と同じ種類)、陸の平均気温・日照・太陽の赤緯
- **進化の時間変化**(--assay-series ファイル、--dump-brains と一緒に): 書き出すたびに、同じ動物をそのままと再帰を切ったものの両方で試験器にかけ、四つの性質それぞれ・四つすべてが検出された割合と持続の平均を 1 行ずつ足す

---

## 文献
本文に出てくる文献。書誌の細部を確かめきれていないものは「(書誌は次の版で確かめる)」と記した。

- Addo-Bediako, A., Chown, S. L., Gaston, K. J. (2002) Metabolic cold adaptation in insects: a large-scale perspective. Functional Ecology 16: 332–338
- Alpert, P. (1996) Nutrient sharing in natural clonal fragments of Fragaria chiloensis. Journal of Ecology 84: 395–406
- Amthor, J. S. (2000) The McCree–de Wit–Penning de Vries–Thornley respiration paradigms: 30 years later. Annals of Botany 86: 1–20
- Anderson, D. J., Adolphs, R. (2014) A framework for studying emotions across species. Cell 157: 187–200
- Anderson, H. E. (1982) Aids to determining fuel models for estimating fire behavior. USDA Forest Service General Technical Report INT-122
- Anderson, T.-H., Domsch, K. H. (1985) Maintenance carbon requirements of actively-metabolizing microbial populations under in situ conditions. Soil Biology and Biochemistry 17: 197–203
- Anderson, T.-H., Domsch, K. H. (1993) The metabolic quotient for CO2 (qCO2) as a specific activity parameter to assess the effects of environmental conditions on the microbial biomass of forest soils. Soil Biology and Biochemistry 25: 393–395
- Atkinson, R., Arey, J. (2003) Atmospheric degradation of volatile organic compounds. Chemical Reviews 103: 4605–4638
- Austin, A. T., Vivanco, L. (2006) Plant litter decomposition in a semi-arid ecosystem controlled by photodegradation. Nature 442: 555–558
- Benjamin, T. B. (1968) Gravity currents and related phenomena. Journal of Fluid Mechanics 31: 209–248
- Blagodatsky, S. A., Heinemeyer, O., Richter, J. (2000) Estimating the active and total soil microbial biomass by kinetic respiration analysis. Biology and Fertility of Soils 32: 73–81. doi:10.1007/s003740000219
- Brune, J. N. (1968) Seismic moment, seismicity, and rate of slip along major fault zones. Journal of Geophysical Research 73: 777–784
- Bradley, T. J. (1987) Physiology of osmoregulation in mosquitoes. Annual Review of Entomology 32: 439–462
- Byers, H. R., Braham, R. R. (1949) The Thunderstorm. U.S. Weather Bureau
- Cardé, R. T., Baker, T. C. (1984) Sexual communication with pheromones. In: Bell, W. J., Cardé, R. T. (eds.) Chemical Ecology of Insects. Chapman and Hall, pp. 355–383
- Changnon, S. A. (1988) Climatography of thunder events in the conterminous United States. Part I: Temporal aspects. Journal of Climate 1: 389–398
- Charba, J. (1974) Application of gravity current model to analysis of squall-line gust front. Monthly Weather Review 102: 140–156
- Chown, S. L., Nicolson, S. W. (2004) Insect Physiological Ecology. Oxford University Press
- Cordonnier, G. et al. (2016) Large scale terrain generation from tectonic uplift and fluvial erosion. Computer Graphics Forum 35(2)
- Corominas, J. (1996) The angle of reach as a mobility index for small and large landslides. Canadian Geotechnical Journal 33: 260–271
- Dai, A. (2006) Recent climatology, variability, and trends in global surface humidity. Journal of Climate 19: 3589–3606
- Engerer, N. A., Stensrud, D. J., Coniglio, M. C. (2008) Surface characteristics of observed cold pools. Monthly Weather Review 136: 4839–4849
- Flowers, T. J., Colmer, T. D. (2008) Salinity tolerance in halophytes. New Phytologist 179: 945–963
- Fosberg, M. A., Deeming, J. E. (1971) Derivation of the 1- and 10-hour timelag fuel moisture calculations for fire-danger rating. USDA Forest Service Research Note RM-207
- Freeze, R. A., Cherry, J. A. (1979) Groundwater. Prentice-Hall
- Fuquay, D. M., Taylor, A. R., Hawe, R. G., Schmid, C. W. Jr. (1972) Lightning discharges that caused forest fires. Journal of Geophysical Research 77(12): 2156–2158
- Gregory, P. H. (1973) The Microbiology of the Atmosphere, 2nd ed.
- Grimm, V. et al. (2020) The ODD protocol for describing agent-based and other simulation models: a second update to improve clarity, replication, and structural realism. Journal of Artificial Societies and Social Simulation 23(2): 7
- Gutenberg, B., Richter, C. F. (1944) Frequency of earthquakes in California. Bulletin of the Seismological Society of America 34: 185–188
- Haldane, J. B. S. (1949) Disease and evolution. La Ricerca Scientifica 19 (Suppl.): 68–76
- Hamilton, W. D. (1980) Sex versus non-sex versus parasite. Oikos 35: 282–290
- Hanks, T. C., Kanamori, H. (1979) A moment magnitude scale. Journal of Geophysical Research 84: 2348–2350
- Harel, M.-A., Mudd, S. M., Attal, M. (2016) Global analysis of the stream power law parameters based on worldwide 10Be denudation rates. Geomorphology 268: 184–196
- Heimsath, A. M. et al. (1997) The soil production function and landscape equilibrium. Nature 388: 358–361
- Jibson, R. W. (2007) Regression models for estimating coseismic landslide displacement. Engineering Geology 91: 209–218
- Jibson, R. W., Harp, E. L., Michael, J. A. (2000) A method for producing digital probabilistic seismic landslide hazard maps. Engineering Geology 58: 271–289
- Johnson, A. I. (1967) Specific yield — compilation of specific yields for various materials. U.S. Geological Survey Water-Supply Paper 1662-D
- Joyner, W. B., Boore, D. M. (1981) Peak horizontal acceleration and velocity from strong-motion records including records from the 1979 Imperial Valley, California, earthquake. Bulletin of the Seismological Society of America 71: 2011–2038
- Karban, R., Baldwin, I. T. (1997) Induced Responses to Herbivory. University of Chicago Press
- Lamont, B. B., Le Maitre, D. C., Cowling, R. M., Enright, N. J. (1991) Canopy seed storage in woody plants. The Botanical Review 57: 277–317
- Larsen, I. J., Montgomery, D. R., Korup, O. (2010) Landslide erosion controlled by hillslope material. Nature Geoscience 3: 247–251
- Lenderink, G., van Meijgaard, E. (2008) Increase in hourly precipitation extremes beyond expectations from temperature changes. Nature Geoscience 1: 511–514
- Manzoni, S. et al. (2012) Environmental and stoichiometric controls on microbial carbon-use efficiency in soils. New Phytologist 196: 79–91
- McCree, K. J. (1974) Equations for the rate of dark respiration of white clover and grain sorghum, as functions of dry weight, photosynthetic rate, and temperature. Crop Science 14: 509–514
- Millero, F. J., Feistel, R., Wright, D. G., McDougall, T. J. (2008) The composition of Standard Seawater and the definition of the Reference-Composition Salinity Scale. Deep-Sea Research I 55: 50–72
- Montgomery, D. R., Brandon, M. T. (2002) Topographic controls on erosion rates in tectonically active mountain ranges. Earth and Planetary Science Letters 201: 481–489
- Montgomery, D. R., Dietrich, W. E. (1988) Where do channels begin? Nature 336: 232–234
- Montgomery, D. R., Dietrich, W. E. (1994) A physically based model for the topographic control on shallow landsliding. Water Resources Research 30: 1153–1171
- Newmark, N. M. (1965) Effects of earthquakes on dams and embankments. Géotechnique 15: 139–160
- Panikov, N. S. (1995) Microbial Growth Kinetics. Chapman & Hall, London
- Penning de Vries, F. W. T., Brunsting, A. H. M., van Laar, H. H. (1974) Products, requirements and efficiency of biosynthesis: a quantitative approach. Journal of Theoretical Biology 45: 339–377
- Penning de Vries, F. W. T. (1974) Substrate utilization and respiration in relation to growth and maintenance in higher plants. Netherlands Journal of Agricultural Science 22: 40–44
- Petersen, W. A., Rutledge, S. A. (1998) On the relationship between cloud-to-ground lightning and convective rainfall. Journal of Geophysical Research 103: 14025–14040
- Pitelka, L. F., Ashmun, J. W. (1985) Physiology and integration of ramets in clonal plants. In: Jackson, Buss, Cook (eds) Population Biology and Evolution of Clonal Organisms. Yale University Press, 399–435
- Press, M. C. & Phoenix, G. K. (2005) Impacts of parasitic plants on natural communities. New Phytologist 166: 737–751
- Prinn, R. G. et al. (2005) Evidence for variability of atmospheric hydroxyl radicals over the past quarter century. Geophysical Research Letters 32: L07809
- Rawls, W. J., Brakensiek, D. L., Saxton, K. E. (1982) Estimation of soil water properties. Transactions of the ASAE 25: 1316–1320
- Roering, J. J., Kirchner, J. W., Dietrich, W. E. (1999) Evidence for nonlinear, diffusive sediment transport on hillslopes and implications for landscape morphology. Water Resources Research 35: 853–870
- Reasenberg, P. A., Jones, L. M. (1989) Earthquake hazard after a mainshock in California. Science 243: 1173–1176
- Rothermel, R. C. (1972) A mathematical model for predicting fire spread in wildland fuels. USDA Forest Service Research Paper INT-115
- Lacy, R.C., Sherman, P.W. (1983) Kin recognition by phenotype matching. American Naturalist 121: 489–512
- Rousk, J., Bååth, E. (2011) Growth of saprotrophic fungi and bacteria in soil. FEMS Microbiology Ecology 78: 17–30
- Ryan, M. G. (1991) Effects of climate change on plant respiration. Ecological Applications 1: 157–167
- Schmidt, K. M. et al. (2001) The variability of root cohesion as an influence on shallow landslide susceptibility in the Oregon Coast Range. Canadian Geotechnical Journal 38: 995–1024
- Skempton, A. W., DeLory, F. A. (1957) Stability of natural slopes in London Clay. Proc. 4th Int. Conf. Soil Mechanics and Foundation Engineering 2: 378–381
- Simard, A. J. (1968) The moisture content of forest fuels — I. A review of the basic concepts. Information Report FF-X-14, Forest Fire Research Institute, Department of Forestry and Rural Development, Ottawa
- Sinsabaugh, R. L. et al. (2013) Carbon use efficiency of microbial communities: stoichiometry, methodology and modelling. Ecology Letters 16: 930–939
- Stock, J. D., Montgomery, D. R. (1999) Geologic constraints on bedrock river incision using the stream power law. Journal of Geophysical Research 104: 4983–4993
- Vautard, R. et al. (2010) Northern Hemisphere atmospheric stilling partly attributed to an increase in surface roughness. Nature Geoscience 3: 756–761
- Wells, D. L., Coppersmith, K. J. (1994) New empirical relationships among magnitude, rupture length, rupture width, rupture area, and surface displacement. Bulletin of the Seismological Society of America 84: 974–1002
- Wieringa, J. (1993) Representative roughness parameters for homogeneous terrain. Boundary-Layer Meteorology 63: 323–363
- Wu, J. (1983) Sea-surface drift currents induced by wind and waves. Journal of Physical Oceanography 13: 1441–1451
- Youngs, R. R., Coppersmith, K. J. (1985) Implications of fault slip rates and earthquake recurrence models to probabilistic seismic hazard estimates. Bulletin of the Seismological Society of America 75: 939–964
