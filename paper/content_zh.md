::title
面向5G-A的短距毫米波链路降雨感知误差预算：<br>湿天线偏移的识别与校正

::abstract
5G-A毫米波链路路径短，降雨衰减常只有零点几分贝，与路径长度无关的衰减误差会按路径长度的倒数放大为雨强误差。本文建立短距链路降雨反演的误差预算，把误差分为湿期长度无关偏移、晴空参考漂移、功率报告量化和导频估计噪声，给出量化报告的共享参考偏差和导频共享参考误差的解析式，并用不同长度链路的跨链路回归识别长度无关偏移。基于OpenMRG数据集中{{n_use}}条{{f_min}}–{{f_max}} GHz子链路、天气雷达和雨量计的分析表明：长度无关偏移随雨强增大，雨停后约半小时内消退到0.1 dB以下，且消退过程与链路长度无关，具有湿天线衰减的特征，是1 km以下链路的主导误差；参考漂移和1 dB报告量化处于0.1–0.3 dB量级；毫秒级晴空参考下导频估计噪声可以忽略。用随雨强变化的湿天线模型校正后，1 km以下链路相对雨量计的小时雨强均方根误差由{{eg_radar_none_le1}} mm/h降至{{eg_radar_kralt_le1}} mm/h，低于天气雷达的{{eg_none_radar_le1}} mm/h。在意大利艾米利亚-罗马涅的第二个地区（{{r2_n_use}}条{{r2_f_min}}–{{r2_f_max}} GHz子链路，电平分辨率1 dB），跨链路回归给出相同的长度无关偏移（分别为{{W}}和{{r2_W}} dB），湿天线校正同样使2 km以下链路的误差低于雷达。5G-A短链路降雨感知的首要任务是校正湿天线偏移，导频测量精度不是限制因素。

::keywords
5G-A；通感一体化；毫米波链路；降雨反演；湿天线衰减；误差预算。

# 引言

地面微波链路的雨致衰减提供路径平均的降雨信息。利用商用微波链路（CML）的接收电平测雨，已在城市和国家尺度得到验证（Messer et al., 2006；Leijnse et al., 2007；Overeem et al., 2011, 2013），并形成了干湿判别、基线估计、湿天线修正和雨强转换的处理流程（Schleiss & Berne, 2010；Overeem et al., 2016；Chwala & Kunstmann, 2019；Graf et al., 2020）。

5G-A为降雨感知提供了更密集的毫米波链路。FR2频段的小站回传和固定无线接入链路多在1 km以内，频率高，对降雨敏感（Han et al., 2020）。通感一体化使基站可以利用解调参考信号（DM-RS）等导频的信道估计（Liu et al., 2022；3GPP, 2024），也有研究把基站作为机会式天气雷达估计降雨（Tornielli Bellini et al., 2026）；而标准测量报告中的参考信号接收功率（RSRP）以1 dB为步长上报（3GPP, 2026）。这使人自然地把测量分辨率和导频估计精度视为5G-A链路降雨感知的关键。

路径变短改变了误差的构成。雨致衰减与路径长度成正比，而湿天线衰减、晴空基线误差和功率量化误差与路径长度无关（Kharadly & Ross, 2001；Leijnse et al., 2008；Zinevich et al., 2010）。换算为雨强后，这些误差随路径缩短而放大。已有湿天线研究（Schleiss et al., 2013；Fencl et al., 2019；Pastorek et al., 2022）和误差来源试验（van Leth et al., 2018）主要针对数公里的回传链路或少数实验链路。对于0.1–1 km的链路，各误差项的相对大小、湿天线校正能否使其可用，以及导频级测量能带来多少改进，还缺少基于大量真实链路的定量结果。

本文研究一个问题：在短距毫米波链路上，哪些与路径长度无关的误差决定降雨反演精度，主导项能否被识别和校正。主要贡献有三点：（1）建立短距链路降雨反演的误差预算，给出量化报告的共享参考偏差解析式和导频共享参考误差的信息下界，并分别用逐事件实测和NR链路仿真检验；（2）用OpenMRG数据集中{{n_use}}条{{f_min}}–{{f_max}} GHz、{{L_min}}–{{L_max}} km的子链路，通过跨链路回归识别长度无关偏移，给出它随雨强和频率的变化以及雨停后的消退，证实其为湿天线衰减并主导1 km以下链路；（3）以雷达和雨量计为参考比较偏移校正方法，表明随雨强变化的湿天线模型使短链路的小时雨强优于雷达，而1 dB报告量化和导频噪声的影响很小，并在设备、频率、分辨率和气候都不同的第二个地区（OpenRainER）复现了主要结论。

# 误差模型

## 2.1 短路径对衰减误差的放大

降雨事件中的衰减相对于事件前的晴空参考计算。设$\hat{P}_t$为第$t$分钟的接收电平（dB），$\hat{P}_{\mathrm{ref}}$为事件前干期的平均电平，则测得衰减为

$$\hat{A}_t = \hat{P}_{\mathrm{ref}} - \hat{P}_t = kLR_t^{\alpha} + W_t + d_t + \eta_{\mathrm{ref}} - \eta_t$$ {#eq:decomp}

^ 式中，$kLR_t^{\alpha}$为长度$L$（km）的路径上雨强$R_t$（mm/h）引起的衰减，$k$和$\alpha$由频率和极化按ITU-R P.838-3确定（ITU-R, 2005）；$W_t$为与路径长度无关的湿期偏移，主要来自天线和天线罩上的水膜；$d_t$为晴空电平相对参考时段的漂移；$\eta_{\mathrm{ref}}$和$\eta_t$分别为参考和当前观测的测量误差，即量化误差或导频估计噪声。

同一事件中，$\eta_{\mathrm{ref}}$对所有分钟相同，$W_t$和$d_t$变化缓慢，三者在事件内近似共享，不随时间平均而减小。一个衰减误差$\delta$（dB）相当于雨强

$$R_{\mathrm{eq}}(\delta) = \left(\frac{\delta}{kL}\right)^{1/\alpha}$$ {#eq:req}

^ 它与路径长度近似成反比。28 GHz垂直极化时$k$ = {{k28}}，$\alpha$ = {{a28}}，0.5 km链路上0.1 dB的误差即相当于{{req_05_01}} mm/h的雨强。

## 2.2 量化报告的共享参考偏差

功率报告把电平取整到步长为$\Delta$的网格：

$$Q(x) = u + \Delta\,\operatorname{round}\!\left(\frac{x-u}{\Delta}\right)$$ {#eq:quant}

^ 式中$u$为网格偏置。若窗口内电平围绕均值$c$以标准差$\sigma$起伏，平均后残留的取整误差为（Widrow & Kollár, 2008）

$$b(c;\sigma) = \sum_{m=1}^{\infty}\frac{(-1)^{m}\Delta}{\pi m}\,\mathrm{e}^{-2\pi^{2}m^{2}\sigma^{2}/\Delta^{2}}\sin\frac{2\pi m(c-u)}{\Delta}$$ {#eq:bias}

^ $b$随$c$在网格中的位置周期变化，起伏越小，衰减越慢。对随机的网格偏置，其均方根为

$$\sigma_{b}(\sigma) = \left[\sum_{m=1}^{\infty}\frac{\Delta^{2}}{2\pi^{2}m^{2}}\,\mathrm{e}^{-4\pi^{2}m^{2}\sigma^{2}/\Delta^{2}}\right]^{1/2}$$ {#eq:brms}

^ $\sigma$ = 0时$\sigma_b = \Delta/\sqrt{12}$，$\sigma > 0.42\Delta$时$\sigma_b < 0.01\Delta$。晴空参考窗口内电平平稳（本文数据中$\sigma$的中位数为{{sig_ref_med}} dB），取整误差几乎不被平均，$\eta_{\mathrm{ref}} = b(c_{\mathrm{ref}};\sigma_{\mathrm{ref}})$成为整个事件共享的偏差；降雨期电平随雨强变化（$\sigma$的中位数为{{sig_ev_med}} dB），取整误差被部分平均。

## 2.3 导频测量的共享参考误差

NR下行DM-RS配置类型1位于OFDM符号2和11，每隔一个子载波放置一个导频（3GPP, 2024），66个资源块时每个时隙有$J$ = {{J}}个导频。用$M$个晴空时隙平均得到参考$r_j = h_j + e_{0j}$，降雨时隙的最小二乘导频估计为$z_{tj} = a_t\exp(\mathrm{i}\phi_t)h_j + e_{tj}$，其中幅度$a_t$满足$A_t = -20\log_{10}a_t$，$\phi_t$为未知公共相位，$e_{0j}$和$e_{tj}$为方差$v/M$和$v$的复高斯噪声。消去未知信道和相位后，$n$次观测平均衰减的Cramér–Rao下界为（Kay, 1993）

$$\operatorname{Var}(\bar{A}_n) \ge \frac{K^{2}}{2J}\left(\frac{1}{\gamma_{\mathrm{r}}} + \frac{1}{n\gamma_{\mathrm{t}}}\right),\quad K = \frac{20}{\ln 10}$$ {#eq:crb}

^ 式中$\gamma_{\mathrm{r}} = M\,\mathrm{SNR}$和$\gamma_{\mathrm{t}} = a^{2}\mathrm{SNR}$分别为单导频的参考信噪比和当前信噪比（信道平均功率归一化为1）。第一项为共享参考误差，不随$n$减小。子载波间隔120 kHz时每个时隙0.125 ms，参考时长为$T_{\mathrm{ref}}$时，共享误差的标准差为

$$\sigma_{\mathrm{ref}} = \frac{K}{\sqrt{2J\gamma_{\mathrm{r}}}},\quad \gamma_{\mathrm{r}} = \frac{T_{\mathrm{ref}}}{0.125\ \mathrm{ms}}\,\mathrm{SNR}$$ {#eq:pilotref}

^ 对剖面幅度估计器作一阶传播，其共享分量为$K^2/(2J\gamma_{\mathrm{r}}) + K^2/[J(\gamma_{\mathrm{r}} + \gamma_{\mathrm{t}})^2]$，在$\gamma_{\mathrm{r}} \gg 1$时与下界一致。导频估计在物理层以浮点数表示，不经过报告量化。

## 2.4 长度无关偏移的识别与校正

同一城市的链路经历相近的降雨，雨致衰减与$kL$成正比，而$W$不随路径变化。按参考雨强（本文用雷达路径平均雨强）把湿期分钟分档，在第$b$档内对全部链路回归

$$\bar{a}_{i,b} = W_b + \beta_b\,\bar{e}_{i,b}$$ {#eq:regress}

^ 式中$\bar{a}_{i,b}$为链路$i$在该档的平均测得衰减，$\bar{e}_{i,b}$为按链路$i$的$k$、$\alpha$和$L$由雷达雨强计算的平均雨致衰减$kLR^{\alpha}$。截距$W_b$即长度无关偏移；参考降雨的乘性偏差只改变斜率$\beta_b$，不改变截距。雷达雨强与路径真实雨强之间存在采样和$Z$–$R$误差，同一档内的真实雨强较分散，$W_b$随雨强的变化因此被压平（回归稀释）。分档曲线$W(R)$适合诊断偏移的存在、量级和随雨强的变化。

校正时比较三种偏移：常数偏移、分档曲线$W(R)$，以及随雨强变化的湿天线衰减模型KR-alt（Pastorek et al., 2022），后者源自Kharadly和Ross（2001）：

$$W(R) = A_{\max}\left[1 - \mathrm{e}^{-dR^{\zeta}}\right]$$ {#eq:kralt}

^ 本文采用文献参数$A_{\max}$ = {{kr_Amax}} dB、$\zeta$ = {{kr_zeta}}、$d$ = {{kr_d}}，不作重新标定。对随雨强变化的偏移，雨强由下式求解：

$$W(\hat{R}) + kL\hat{R}^{\alpha} = \hat{A}$$ {#eq:invert}

^ 左端随$\hat{R}$单调增加，解唯一；$\hat{A} < W(0)$时取$\hat{R}$ = 0。

# 数据与方法

## 3.1 OpenMRG数据

OpenMRG数据集（Andersson et al., 2022）包含瑞典哥德堡728条子链路2015年6–8月的发射和接收电平（TSL、RSL），采样间隔10 s，分辨率约1/3 dB，并配有SMHI天气雷达和10个城市雨量计的1 min资料。本文选取24–40 GHz的{{n_sel}}条子链路（实际频率{{f_min}}–{{f_max}} GHz），其中{{n_use}}条（{{n_phys}}条物理链路）的湿期分钟数不少于600，用于分析。链路长度为{{L_min}}–{{L_max}} km，中位数{{L_med}} km，其中1 km以下{{n_short}}条，与5G-A毫米波链路的长度相当。主要设置见表@tab:data。

第二个地区采用OpenRainER数据集（Covi & Roversi, 2026；Covi et al., 2026），包含意大利艾米利亚-罗马涅Lepida网络的151条链路在2021—2022年的1 min收发电平，以及C波段雷达15 min降雨量和气象站15 min雨量（Roversi et al., 2020）。数据提供方在记录前已将电平按1 dB取整，与5G的RSRP上报步长相同；本文按数据说明再取整到1 dB。链路频率为{{r2_f_min}}–{{r2_f_max}} GHz（落在5G的n258频段内），长度{{r2_L_min}}–{{r2_L_max}} km（中位数{{r2_L_med}} km），湿期分钟数满足要求的有{{r2_n_use}}条子链路（{{r2_n_phys}}条物理链路），其中2 km以下{{r2_n_le2}}条。为避开降雪，只用5–9月的数据。

::table tab:data | 数据与处理设置 | 3.0,11.4
项目 | 设置
链路 | OpenMRG，2015年6–8月；24–40 GHz子链路{{n_sel}}条，分析{{n_use}}条（{{n_phys}}条物理链路）
电平 | TSL、RSL，10 s采样，约1/3 dB分辨率，按1 min平均；衰减取TSL − RSL
参考降雨 | SMHI雷达，2 km、5 min，$Z = 200R^{1.5}$，沿路径20点平均；10个城市雨量计，1 min
湿期 | 雷达路径雨强在前后30 min内超过0.1 mm/h；连续湿期构成事件（不短于10 min）
晴空参考 | 事件开始前6 h内干期分钟的平均衰减，不少于60 min，在事件内保持不变
雨强分档 | 0、0.1、0.5、1、2、5、10、200 mm/h
报告量化 | 10 s RSL按1 dB和2 dB取整，各子链路的网格偏置随机
偏移校正 | 常数偏移、分档$W(R)$（6–7月标定）；KR-alt（文献参数）
检验 | 8月；雷达参考用全部子链路，雨量计参考用距雨量计{{gauge_km}} km以内的子链路
区间估计 | 以物理链路为单元的集群自助法，2000次
::end

## 3.2 处理流程

信号电平按1 min平均，衰减取TSL − RSL。雷达反射率按$Z = 200R^{1.5}$换算为雨强，沿路径等间隔取20个点平均。湿期由雷达判别，使各试验的干湿分类完全相同，从而只比较测量误差和偏移。每个事件以开始前6 h内干期分钟的平均衰减作为晴空参考，事件内各分钟减去同一参考得到降雨期衰减。共得到{{n_events}}个事件，持续时间的中位数为{{ev_med_min}} min。

为检验不依赖雷达的应用情形，另用链路自身信号判别干湿（Schleiss & Berne, 2010）：衰减的{{lm_window}} min滑动标准差超过阈值时记为湿期，阈值按6–7月的湿期比例为每条子链路确定，8月不再使用雷达。

## 3.3 试验设计

（1）长度无关偏移。按式(@eq:regress)对全部湿期分钟和7个雨强档分别回归，以物理链路为单元做集群自助，得到95%区间。另以$k_iL_i$代替$\bar{e}_{i,b}$作不依赖雷达雨强数值的回归，以Theil–Sen回归检验离群链路的影响，并检验湿期边界（15、30、60 min）和参考窗口（3、6、12 h）的影响。对每个事件，从雷达路径雨强最后一次超过0.1 mm/h起，保持事件参考，跟踪其后90 min内路径无雨时的衰减，得到雨停后的消退曲线。

（2）参考漂移。在干期中重现一次降雨事件：以6 h参考预测其后一个完整干期小时的衰减（真值为0），参考窗口结束至检验小时结束的时间取1、2、4、6 h，共{{n_dry_hours}}个检验小时。

（3）报告量化。把10 s RSL按$\Delta$ = 1 dB（RSRP的上报步长）和2 dB（差分L1-RSRP的步长）取整（3GPP, 2026），各子链路的网格偏置$u$随机抽取；量化后的报告用相同的干期分钟重建参考。逐事件比较参考和降雨期的量化误差与式(@eq:bias)的预测值。

（4）导频仿真。在NR DM-RS时域链路中（子载波间隔120 kHz，1024点FFT，66个资源块，固定多径），取信噪比5、15、25 dB和参考长度1、8、80个时隙；每种条件生成400个参考，每个参考复用于128个降雨时隙（衰减0.5 dB），用同一参考下估计值的协方差估计共享分量。

（5）偏移校正的效果。常数偏移和分档$W(R)$用6–7月数据标定，KR-alt采用文献参数；在8月比较各校正的小时雨强。以雷达路径平均雨强为参考时使用全部子链路；以雨量计为参考时使用距最近雨量计不超过{{gauge_km}} km的{{g_n_le1}}条1 km以下和{{g_n_gt1}}条1 km以上子链路，并同时评价雷达本身。另对链路自身判别干湿和1 dB、2 dB报告重复KR-alt校正。各长度组以子链路的中位数报告，降幅的区间按物理链路集群自助计算。

（6）第二个地区。对OpenRainER采用相同的处理和参数：2021年5–9月标定常数偏移和分档$W(R)$，2022年5–9月检验。雷达和气象站时间戳表示累计时段的终点还是起点，由其与1 min链路衰减的相关性确定（均为终点）。雷达为15 min分辨率，雨停后的消退曲线按15 min分档。对两个地区，均在每个时滞上对各链路的残余衰减作关于$k_iL_i$的回归，以截距表示与长度无关的部分。

# 结果

## 4.1 长度无关偏移主导短链路误差

图@fig:example为一条{{ex_L}} km、{{ex_f}} GHz子链路在一次{{ex_hours}} h降雨中的衰减。该子链路按固定规则选取：其8月小时雨强均方根误差在KR-alt校正后的降幅最接近0.5 km以下链路的中位数（{{ex_red_med}}%），事件为其8月雷达累计雨量最大的一次。测得衰减始终明显高于雷达对应的雨致衰减，差值随雨强增大；加上KR-alt湿天线衰减后，两者基本吻合。该事件雷达路径累计雨量为{{ex_radar_mm}} mm；链路在不扣除偏移、扣除常数偏移、扣除分档$W(R)$和扣除KR-alt湿天线衰减时分别为{{ex_link_none_mm}}、{{ex_link_const_mm}}、{{ex_link_rain_mm}}和{{ex_link_kralt_mm}} mm。

跨链路回归给出湿期平均的长度无关偏移$W$ = {{W}} dB（95%区间{{W_ci}} dB），斜率$\beta$ = {{beta}}（{{beta_ci}}）（图@fig:offset(a)）。不依赖雷达雨强数值的回归给出{{W_rf}} dB（{{W_rf_ci}} dB），Theil–Sen回归给出{{W_ts}} dB。28–29、32–33和37–39 GHz三个频段分别为{{W_26_30}}、{{W_32_34}}和{{W_37_40}} dB，37–40 GHz最大。改变湿期边界和参考窗口时，$W$在{{W_sens}} dB之间变化。湿期平均衰减中，这一偏移在0.5 km以下链路占{{share_05}}%，0.5–1 km占{{share_1}}%，4 km以上仍占{{share_4}}%。

偏移随雨强显著增大（图@fig:offset(b)）：雷达雨强低于{{bin_low}} mm/h时为{{W_low}} dB，超过{{bin_high}} mm/h时达{{W_high}} dB（95%区间{{W_high_ci}} dB）。1 mm/h时$W$ = {{W_at_1}} dB，是0.5 km以下链路雨致衰减的{{W_over_rain_1_05}}倍。只有雨强超过一定值，雨致衰减才大于偏移；这一临界雨强的中位数在0.5 km以下链路为{{Rc_05}} mm/h，0.5–1 km为{{Rc_1}} mm/h，1–2 km为{{Rc_2}} mm/h，4 km以上仅{{Rc_20}} mm/h。在除最高档以外的各档，按$k_iL_i$回归和Theil–Sen回归的截距与雷达回归相差不超过{{W_alt_maxdiff_low}} dB。在0.1–2 mm/h各档，KR-alt模型比分档估计高{{kr_minus_bin_lo}}–{{kr_minus_bin_hi}} dB，最高档为分档估计的{{kr_over_bin_top}}倍，与回归稀释的方向一致。

![一条短链路的降雨事件（子链路{{ex_sublink}}，{{ex_date}}）。(a) 衰减；(b) 雨强。灰底为湿期。](fig1_example.png){#fig:example width=13.2}

雨停后的消退（图@fig:offset(c)）提供了偏移来源的直接证据。雷达路径雨强降到0.1 mm/h以下后的前5 min，各子链路残余衰减的中位数为{{dry_first}} dB，在雨停后{{dry_half_lo}}–{{dry_half_hi}} min的分档中减半，在{{dry_below01_lo}}–{{dry_below01_hi}} min的分档中降到0.1 dB以下，60 min时为{{dry_60}} dB。1 km以下、1–4 km和4 km以上三组链路的消退曲线在雨停10 min后相差不超过{{dry_spread_after10}} dB。按跨链路回归的截距计，雨停后7.5 min和22.5 min时与长度无关的残余衰减分别为{{di_mrg_7}}和{{di_mrg_22}} dB，15 min后关于$kL$的斜率不超过{{di_mrg_slope_after15}}。残余路径降雨引起的衰减会随链路长度增加，而这一消退与长度无关，符合天线和天线罩上水膜的蒸发和流失。

![长度无关偏移。(a) 各子链路湿期平均衰减与雷达对应雨致衰减，直线为跨链路回归，截距即$W$；(b) 各雨强档的截距$W(R)$（阴影为95%区间）、KR-alt模型和28 GHz垂直极化下0.3 km、1 km路径的雨致衰减；(c) 雷达雨停后残余衰减的中位数。](fig2_offset.png){#fig:offset width=14.6}

## 4.2 参考漂移与报告量化

干期检验中，参考误差的小时平均标准差（各子链路的中位数）在参考结束1 h后为{{drift_1h}} dB，6 h后增至{{drift_6h}} dB（图@fig:drift(a)），平均偏差仅{{drift_bias}} dB。逐分钟均方根误差为{{drift_min}} dB，与小时平均值接近，说明这一误差主要来自缓慢漂移：它在事件内近似共享，不能靠时间平均消除。

1 dB报告下，参考的量化误差均方根为{{q1_ref_rms}} dB，降雨期为{{q1_ev_rms}} dB；小雨时两者落在同一量化间隔内，相互部分抵消，事件级衰减误差的均方根为{{q1_tot_rms}} dB。用式(@eq:bias)逐事件预测参考误差，回归斜率为{{q1_slope}}，$R^2$ = {{q1_r2}}（图@fig:drift(b)）。2 dB报告的对应值为{{q2_ref_rms}}、{{q2_ev_rms}}和{{q2_tot_rms}} dB，斜率{{q2_slope}}，$R^2$ = {{q2_r2}}。1 dB时预测的均方根（{{q1_ref_pred}} dB）略小于实测值，差异来自窗口内电平分布偏离高斯以及原始数据自身约1/3 dB的取整。

![参考漂移与报告量化。(a) 干期检验中参考误差随参考时效的变化，阴影为各子链路的四分位范围；(b) 1 dB报告下各事件参考量化误差的实测值与式(@eq:bias)预测值，红线为1:1线。](fig3_drift_quant.png){#fig:drift width=14.5}

## 4.3 导频测量的共享参考误差

图@fig:pilot比较了仿真得到的共享参考误差与闭式结果。9种条件下，仿真值与闭式值相差{{pilot_rel_range}}，闭式值均落在仿真的95%区间内。15 dB时，1 ms晴空导频的共享误差标准差为{{pilot_sd_15_1ms}} dB，10 ms时为{{pilot_sd_15_10ms}} dB；5 dB、10 ms时为{{pilot_sd_5_10ms}} dB。毫秒级参考下，导频测量的共享误差比参考漂移和1 dB量化误差小一到两个数量级，比长度无关偏移小两个数量级以上。

![导频测量的共享参考误差随晴空参考时长的变化。线为闭式结果，点为NR DM-RS链路仿真（误差线为95%区间）；水平线为OpenMRG中的长度无关偏移（1 mm/h）、1 h参考漂移和1 dB报告的事件级量化误差。](fig4_pilot.png){#fig:pilot width=10.5}

## 4.4 误差预算

表@tab:budget以1 mm/h雨强为例，把各误差项换算为等效雨强。0.5 km以下链路上，雨致衰减只有{{b_L05_signal_db}} dB，长度无关偏移相当于{{b_L05_W1_mmh}} mm/h；扣除平均偏移后，链路间残余偏移的离散仍相当于{{b_L05_resid_mmh}} mm/h；参考漂移和1 dB量化分别相当于{{b_L05_drift_mmh}}和{{b_L05_quant_mmh}} mm/h；导频噪声只相当于{{b_L05_pilot_mmh}} mm/h。路径越长，各项越小。

::table tab:budget | 1 mm/h时各误差项的等效雨强（mm/h，各子链路的中位数） | 1.4,1.45,1.8,2.3,1.8,1.8,1.8,1.8
长度 (km) | 链路数 | 雨致衰减 (dB) | 长度无关偏移 | 残余偏移 | 参考漂移 | 1 dB量化 | 导频噪声
≤0.5 | {{b_L05_n}} | {{b_L05_signal_db}} | {{b_L05_W1_mmh}} | {{b_L05_resid_mmh}} | {{b_L05_drift_mmh}} | {{b_L05_quant_mmh}} | {{b_L05_pilot_mmh}}
0.5–1 | {{b_L1_n}} | {{b_L1_signal_db}} | {{b_L1_W1_mmh}} | {{b_L1_resid_mmh}} | {{b_L1_drift_mmh}} | {{b_L1_quant_mmh}} | {{b_L1_pilot_mmh}}
1–2 | {{b_L2_n}} | {{b_L2_signal_db}} | {{b_L2_W1_mmh}} | {{b_L2_resid_mmh}} | {{b_L2_drift_mmh}} | {{b_L2_quant_mmh}} | {{b_L2_pilot_mmh}}
2–4 | {{b_L4_n}} | {{b_L4_signal_db}} | {{b_L4_W1_mmh}} | {{b_L4_resid_mmh}} | {{b_L4_drift_mmh}} | {{b_L4_quant_mmh}} | {{b_L4_pilot_mmh}}
>4 | {{b_L20_n}} | {{b_L20_signal_db}} | {{b_L20_W1_mmh}} | {{b_L20_resid_mmh}} | {{b_L20_drift_mmh}} | {{b_L20_quant_mmh}} | {{b_L20_pilot_mmh}}
::end

^ 注：长度无关偏移取分档$W$(1 mm/h) = {{W_at_1}} dB（KR-alt为{{kr_at_1}} dB）；残余偏移取1–2 mm/h档内链路间的离散{{resid_1}} dB；参考漂移为1 h干期检验的标准差；1 dB量化为事件级误差的均方根；导频噪声为15 dB、10 ms参考下的共享误差{{pilot_sd_15_10ms}} dB。

## 4.5 湿天线校正的效果

以雷达为参考时（表@tab:eval，图@fig:eval(a)），KR-alt校正使0.5 km以下链路的8月小时雨强均方根误差由{{er_radar_none_L05}} mm/h降至{{er_radar_kralt_L05}} mm/h，降幅中位数{{red_kr_L05}}%（95%区间{{red_kr_L05_ci}}%），{{impr_kr_L05}}%的子链路得到改善；0.5–1 km和1–2 km分别降低{{red_kr_L1}}%和{{red_kr_L2}}%，4 km以上降低{{red_kr_L20}}%。常数偏移只降低{{red_const_L20}}%–{{red_const_L05}}%。分档$W(R)$介于两者之间（0.5 km以下为{{er_radar_rain_L05}} mm/h），与第2.4节的回归稀释一致。0.5 km以下链路的累计雨量比（链路/雷达）由{{eq_radar_none_L05}}降至{{eq_radar_kralt_L05}}。改用链路自身判别干湿后（湿期命中率{{lm_hit}}%，误判率{{lm_false}}%），KR-alt的均方根误差几乎不变（0.5 km以下为{{er_link_kralt_L05}} mm/h）。

::table tab:eval | 8月小时雨强的均方根误差（mm/h，各子链路的中位数） | 3.7,1.3,1.3,1.3,1.3,1.5,1.5,1.75
方法 | ≤0.5 km | 0.5–1 km | 1–2 km | >4 km | ≤1 km† | >1 km† | 雨量比†
不扣除偏移 | {{er_radar_none_L05}} | {{er_radar_none_L1}} | {{er_radar_none_L2}} | {{er_radar_none_L20}} | {{eg_radar_none_le1}} | {{eg_radar_none_gt1}} | {{egq_radar_none_le1}}
常数偏移 | {{er_radar_constant_L05}} | {{er_radar_constant_L1}} | {{er_radar_constant_L2}} | {{er_radar_constant_L20}} | {{eg_radar_constant_le1}} | {{eg_radar_constant_gt1}} | {{egq_radar_constant_le1}}
分档$W(R)$ | {{er_radar_rain_L05}} | {{er_radar_rain_L1}} | {{er_radar_rain_L2}} | {{er_radar_rain_L20}} | {{eg_radar_rain_le1}} | {{eg_radar_rain_gt1}} | {{egq_radar_rain_le1}}
KR-alt | {{er_radar_kralt_L05}} | {{er_radar_kralt_L1}} | {{er_radar_kralt_L2}} | {{er_radar_kralt_L20}} | {{eg_radar_kralt_le1}} | {{eg_radar_kralt_gt1}} | {{egq_radar_kralt_le1}}
KR-alt，链路判别干湿 | {{er_link_kralt_L05}} | {{er_link_kralt_L1}} | {{er_link_kralt_L2}} | {{er_link_kralt_L20}} | — | — | —
KR-alt，1 dB报告 | {{er_radar_q1_kralt_L05}} | {{er_radar_q1_kralt_L1}} | {{er_radar_q1_kralt_L2}} | {{er_radar_q1_kralt_L20}} | {{eg_radar_q1_kralt_le1}} | {{eg_radar_q1_kralt_gt1}} | {{egq_radar_q1_kralt_le1}}
KR-alt，2 dB报告 | {{er_radar_q2_kralt_L05}} | {{er_radar_q2_kralt_L1}} | {{er_radar_q2_kralt_L2}} | {{er_radar_q2_kralt_L20}} | {{eg_radar_q2_kralt_le1}} | {{eg_radar_q2_kralt_gt1}} | {{egq_radar_q2_kralt_le1}}
天气雷达 | — | — | — | — | {{eg_none_radar_le1}} | {{eg_none_radar_gt1}} | {{egq_none_radar_le1}}
::end

^ 注：不带†的列以雷达路径平均雨强为参考；带†的列以距子链路{{gauge_km}} km以内的雨量计为参考，雨量比为1 km以下子链路的累计雨量之比（估计/雨量计）。

以雨量计为参考时（图@fig:eval(b)），1 km以下链路经KR-alt校正后的均方根误差为{{eg_radar_kralt_le1}} mm/h，低于雷达的{{eg_none_radar_le1}} mm/h；在{{g_pairs_le1}}条物理链路中有{{g_kr_better_le1}}条优于雷达（符号检验$p$ = ${{g_sign_p_le1}}$），累计雨量比为{{egq_radar_kralt_le1}}。1 km以上链路为{{eg_radar_kralt_gt1}} mm/h，{{g_pairs_gt1}}条物理链路中{{g_kr_better_gt1}}条优于雷达。不扣除偏移时，1 km以下链路的误差为{{eg_radar_none_le1}} mm/h，约为雷达的10倍。

湿天线校正后，报告分辨率的影响仍然很小。以雷达为参考，1 dB报告使0.5 km以下链路的均方根误差增加{{kr_q1_inc_05}}%，2 dB报告增加{{kr_q2_inc_05}}%；以雨量计为参考，1 km以下链路分别为{{eg_radar_q1_kralt_le1}}和{{eg_radar_q2_kralt_le1}} mm/h（原始分辨率为{{eg_radar_kralt_le1}} mm/h）。

![8月的小时雨强检验。(a) 以雷达为参考，点为各长度组子链路的中位数；(b) 以{{gauge_km}} km以内的雨量计为参考，括号内为子链路数。常数偏移和分档$W(R)$在6–7月标定，KR-alt采用文献参数。](fig5_evaluation.png){#fig:eval width=14.5}

## 4.6 第二个地区的检验

表@tab:regions和图@fig:region2对比了两个地区的结果。OpenRainER的跨链路回归给出长度无关偏移{{r2_W}} dB（95%区间{{r2_W_ci}} dB），与OpenMRG的{{W}} dB几乎相同；不依赖雷达雨强数值的回归给出{{r2_W_rf}} dB，2021年单独标定为{{r2_W_cal}} dB。这一偏移在2 km以下链路占湿期平均衰减的{{r2_share_le2}}%，在10 km以上链路仍占{{r2_share_gt10}}%。回归斜率$\beta$ = {{r2_beta}}，明显小于OpenMRG，原因是该地区未经雨量计校正的雷达降雨量偏高：在雨量计处，雷达累计雨量为雨量计的{{r2_radar_over_gauge}}倍。截距不受这一乘性偏差影响。

![第二个地区（OpenRainER）的检验。(a) 跨链路回归；(b) 两个地区各雨强档的截距$W(R)$（阴影为95%区间）与KR-alt模型；(c) 两个地区雨停后与长度无关的残余衰减（各时滞上关于$kL$回归的截距）；(d) 2022年5–9月以1 km以内雨量计为参考的小时雨强均方根误差，括号内为子链路数。](fig6_region2.png){#fig:region2 width=13.6}

偏移同样随雨强增大（图@fig:region2(b)），由{{r2_W_low}} dB增至{{r2_W_high}} dB；由于雷达偏高，同一雷达雨强对应的实际雨强较低，曲线相对OpenMRG右移。雨停后与长度无关的残余衰减在7.5 min和22.5 min时为{{di_er_7}}和{{di_er_22}} dB（图@fig:region2(c)），比OpenMRG略高，消退过程相似；残余衰减关于$kL$的斜率中位数为{{di_er_slope_med}}，长链路上还有少量路径残余降雨。1 h参考漂移为{{r2_drift_1h}} dB，约为OpenMRG的两倍，可能与1 dB取整、更长的路径和设备差异有关。

以雨量计为参考（2022年5–9月，{{r2_n_gauges}}个雨量计1 km以内的子链路），2 km以下链路的均方根误差在不扣除偏移时为{{r2_gauge_none_le2}} mm/h，KR-alt校正后为{{r2_gauge_kralt_le2}} mm/h（降幅中位数{{r2_g_red_le2}}%），累计雨量比为{{r2q_gauge_kralt_le2}}；同一位置雷达的误差为{{r2_gauge_radar_le2}} mm/h。{{r2_g_phys_le2}}条物理链路中有{{r2_g_better_le2}}条优于雷达（符号检验$p$ = {{r2_g_p_le2}}，样本较少）。2 km以上链路校正前后分别为{{r2_gauge_none_gt2}}和{{r2_gauge_kralt_gt2}} mm/h，{{r2_g_phys_gt2}}条物理链路中{{r2_g_better_gt2}}条优于雷达（$p$ = ${{r2_g_p_gt2}}$）。该地区的电平本身就是1 dB分辨率，因此这一结果直接说明：校正湿天线偏移后，1 dB报告足以支持优于雷达的小时雨强。以雷达为参考时，2 km以下链路的均方根误差由{{r2_radar_none_le2}} mm/h降至{{r2_radar_kralt_le2}} mm/h；2 km以上链路由{{r2_radar_none_gt2}} mm/h变为{{r2_radar_kralt_gt2}} mm/h，略有增加，原因是该地区的雷达本身偏高，以它为参考会把偏高的雷达降雨量当作真值。

::table tab:regions | 两个地区的结果对比 | 4.4,4.9,4.9
项目 | OpenMRG（瑞典哥德堡） | OpenRainER（意大利艾米利亚-罗马涅）
时段 | 2015年6–8月 | 2021、2022年5–9月
频率 | {{f_min}}–{{f_max}} GHz | {{r2_f_min}}–{{r2_f_max}} GHz
子链路（物理链路） | {{n_use}}（{{n_phys}}） | {{r2_n_use}}（{{r2_n_phys}}）
长度（中位数） | {{L_min}}–{{L_max}} km（{{L_med}} km） | {{r2_L_min}}–{{r2_L_max}} km（{{r2_L_med}} km）
电平记录 | 10 s，约1/3 dB | 1 min，1 dB
雷达 | 2 km，5 min | 约1 km，15 min
常数偏移$W$（95%区间） | {{W}}（{{W_ci}}）dB | {{r2_W}}（{{r2_W_ci}}）dB
雨停后7.5/22.5 min的长度无关残余 | {{di_mrg_7}}/{{di_mrg_22}} dB | {{di_er_7}}/{{di_er_22}} dB
1 h参考漂移 | {{drift_1h}} dB | {{r2_drift_1h}} dB
短链路相对雨量计的误差（不扣除→KR-alt） | {{eg_radar_none_le1}}→{{eg_radar_kralt_le1}} mm/h（≤1 km） | {{r2_gauge_none_le2}}→{{r2_gauge_kralt_le2}} mm/h（≤2 km）
同一位置雷达的误差 | {{eg_none_radar_le1}} mm/h | {{r2_gauge_radar_le2}} mm/h
雷达/雨量计累计雨量比 | {{egq_none_radar_le1}} | {{r2_radar_over_gauge}}
::end


# 讨论

## 5.1 对5G-A降雨感知的含义

在0.1–1 km的毫米波链路上，误差的主次依次为：随雨强变化的湿天线偏移，链路间的残余偏移和参考漂移，功率报告量化，导频估计噪声。两个地区的设备、频率（{{r2_f_min}}–{{f_max}} GHz）、电平分辨率和气候都不同，长度无关偏移却分别为{{W}}和{{r2_W}} dB，说明它主要取决于湿天线这一普遍机制。湿天线偏移在1 mm/h时相当于雨致衰减的数倍，必须随雨强校正；常数偏移只能去掉其中一小部分。KR-alt的参数来自其他城市的链路观测（Pastorek et al., 2022），在本文{{f_min}}–{{f_max}} GHz、{{L_min}}–{{L_max}} km的链路上不经重新标定即可把1 km以下链路的误差降到雷达以下。分档$W(R)$受回归稀释影响，宜用于诊断：跨链路回归的截距和雨停后的消退曲线只依赖网络自身的衰减数据和降雨的时间信息，可以用来检验湿天线模型是否适用于新的设备和天线罩。

校正湿天线偏移后，1 dB报告对小时雨强的影响只有几个百分点，2 dB报告在短链路上的影响更明显。毫秒级晴空参考下的导频噪声比其他误差小一到两个数量级。因此，就雨强精度而言，5G-A链路降雨感知的投入应首先放在湿天线校正和参考漂移控制上，提高导频估计精度的收益很小。导频级测量的价值在于避免报告量化，并有可能为湿天线消退等分钟级过程提供更密的观测。

由式(@eq:req)，在雨强$R$处要达到雨强误差$\varepsilon$，允许的共享衰减误差约为$\alpha kLR^{\alpha-1}\varepsilon$。28 GHz、1 mm/h、$\varepsilon$ = 1 mm/h时，0.3 km链路允许{{sigma_req_03}} dB，1 km链路允许{{sigma_req_10}} dB；而本文测得的1 h参考漂移为{{drift_1h}} dB。0.5 km以下链路要在弱降雨中达到1 mm/h的精度，还需要把参考漂移控制在0.1 dB以下，例如缩短参考与降雨之间的间隔。

## 5.2 局限与后续工作

（1）OpenMRG为回传链路设备，其天线和天线罩与5G-A有源阵列不同，湿天线偏移的量值和消退时间需要在5G-A链路上重新检验，但跨链路回归和消退曲线的方法可以直接使用。（2）雷达的2 km、5 min采样会进入分档$W(R)$和以雷达为参考的检验；雨量计是点观测，与路径平均存在代表性差异，本文只用距其1 km以内的子链路。（3）KR-alt采用文献参数，在8月以外和更强降雨下的表现有待更长数据检验。（4）干期检验中的漂移包含CML接收电平检测的误差，导频测量能否减小这一部分，需要5G-A实测数据检验。（5）两个数据集共覆盖三个暖季，冬季的雪和冰未作检验；OpenRainER中2 km以下且靠近雨量计的物理链路只有{{r2_g_phys_le2}}条，短链路的结论主要来自OpenMRG。

# 结论

（1）短距毫米波链路降雨反演的主导误差是湿天线偏移。OpenMRG的{{n_use}}条{{f_min}}–{{f_max}} GHz子链路中，跨链路回归给出的长度无关偏移随雨强由{{W_low}} dB增至{{W_high}} dB，在0.5 km以下链路占湿期平均衰减的{{share_05}}%；雨停后它在约半小时内消退到0.1 dB以下，且消退与链路长度无关。OpenRainER的{{r2_n_use}}条{{r2_f_min}}–{{r2_f_max}} GHz子链路给出相同的偏移（{{r2_W}} dB）。

（2）用随雨强变化的湿天线模型校正后，0.5 km以下链路以雷达为参考的小时雨强均方根误差降低{{red_kr_L05}}%；1 km以下链路以雨量计为参考的误差为{{eg_radar_kralt_le1}} mm/h，低于天气雷达的{{eg_none_radar_le1}} mm/h；在OpenRainER中，2 km以下链路由{{r2_gauge_none_le2}} mm/h降至{{r2_gauge_kralt_le2}} mm/h，同样低于雷达。常数偏移校正不能达到这一效果。

（3）参考漂移（1 h为{{drift_1h}} dB）和1 dB报告量化（事件级{{q1_tot_rms}} dB）是次要的共享误差，其中量化误差主要来自晴空参考，可由式(@eq:bias)逐事件预测；湿天线校正后，1 dB报告只使误差增加几个百分点。毫秒级参考下导频测量的共享误差约为{{pilot_sd_15_10ms}} dB，不构成限制。5G-A短链路降雨感知应首先校正湿天线偏移并控制参考漂移。

#* 数据与代码可用性

OpenMRG数据集（Andersson et al., 2022）以CC BY-SA 4.0许可发布于https://doi.org/10.5281/zenodo.7107689；OpenRainER数据集（Covi & Roversi, 2026）以CC BY 4.0许可发布于https://doi.org/10.5281/zenodo.10593848。全部分析代码、配置、汇总结果和排版脚本见https://github.com/KJTWangyw/Rain-5G-A（分支paper/short-link-error-budget），文中数字均由代码生成。

#* 参考文献

::refs
3GPP. (2024). NR; Physical channels and modulation (3GPP TS 38.211 version 18.3.0 Release 18). ETSI TS 138 211 V18.3.0.
3GPP. (2026). NR; Requirements for support of radio resource management (3GPP TS 38.133 version 18.15.0 Release 18). ETSI TS 138 133 V18.15.0.
Andersson, J. C. M., Olsson, J., van de Beek, R., & Hansryd, J. (2022). OpenMRG: Open data from Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden. Earth System Science Data, 14(12), 5411-5426.
Chwala, C., & Kunstmann, H. (2019). Commercial microwave link networks for rainfall observation: Assessment of the current status and future challenges. WIREs Water, 6(2), e1337.
Covi, E., & Roversi, G. (2026). OpenRainER (Version 2.0.3) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.10593848
Covi, E., Roversi, G., & Nebuloni, R. (2026). OpenRainER: An open-source dataset for studying the opportunistic sensing of rainfall in Emilia-Romagna, Italy. Earth System Science Data Discussions. https://doi.org/10.5194/essd-2026-160
Fencl, M., Valtr, P., Kvicera, M., & Bareš, V. (2019). Quantifying wet antenna attenuation in 38-GHz commercial microwave links of cellular backhaul. IEEE Geoscience and Remote Sensing Letters, 16(4), 514-518.
Graf, M., Chwala, C., Polz, J., & Kunstmann, H. (2020). Rainfall estimation from a German-wide commercial microwave link network: Optimized processing and validation for 1 year of data. Hydrology and Earth System Sciences, 24(6), 2931-2950.
Han, C., Huo, J., Gao, Q., Su, G., & Wang, H. (2020). Rainfall monitoring based on next-generation millimeter-wave backhaul technologies in a dense urban environment. Remote Sensing, 12(6), 1045.
ITU-R. (2005). Specific attenuation model for rain for use in prediction methods (Recommendation ITU-R P.838-3). International Telecommunication Union.
Kay, S. M. (1993). Fundamentals of statistical signal processing: Estimation theory. Prentice Hall.
Kharadly, M. M. Z., & Ross, R. (2001). Effect of wet antenna attenuation on propagation data statistics. IEEE Transactions on Antennas and Propagation, 49(8), 1183-1191.
Leijnse, H., Uijlenhoet, R., & Stricker, J. N. M. (2007). Rainfall measurement using radio links from cellular communication networks. Water Resources Research, 43(3), W03201.
Leijnse, H., Uijlenhoet, R., & Stricker, J. N. M. (2008). Microwave link rainfall estimation: Effects of link length and frequency, temporal sampling, power resolution, and wet antenna attenuation. Advances in Water Resources, 31(11), 1481-1493.
Liu, F., Cui, Y., Masouros, C., Xu, J., Han, T. X., Eldar, Y. C., & Buzzi, S. (2022). Integrated sensing and communications: Toward dual-functional wireless networks for 6G and beyond. IEEE Journal on Selected Areas in Communications, 40(6), 1728-1767.
Messer, H., Zinevich, A., & Alpert, P. (2006). Environmental monitoring by wireless communication networks. Science, 312(5774), 713.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2011). Measuring urban rainfall using microwave links from commercial cellular communication networks. Water Resources Research, 47(12), W12505.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2013). Country-wide rainfall maps from cellular communication networks. Proceedings of the National Academy of Sciences, 110(8), 2741-2745.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2016). Retrieval algorithm for rainfall mapping from microwave links in a cellular communication network. Atmospheric Measurement Techniques, 9(5), 2425-2444.
Pastorek, J., Fencl, M., Rieckermann, J., & Bareš, V. (2022). Precipitation estimates from commercial microwave links: Practical approaches to wet-antenna correction. IEEE Transactions on Geoscience and Remote Sensing, 60, 1-9.
Roversi, G., Alberoni, P. P., Fornasiero, A., & Porcù, F. (2020). Commercial microwave links as a tool for operational rainfall monitoring in Northern Italy. Atmospheric Measurement Techniques, 13(11), 5779-5797.
Schleiss, M., & Berne, A. (2010). Identification of dry and rainy periods using telecommunication microwave links. IEEE Geoscience and Remote Sensing Letters, 7(3), 611-615.
Schleiss, M., Rieckermann, J., & Berne, A. (2013). Quantification and modeling of wet-antenna attenuation for commercial microwave links. IEEE Geoscience and Remote Sensing Letters, 10(5), 1195-1199.
Tornielli Bellini, D., Montopoli, M., Tagliaferri, D., Baldini, L., Adirosi, E., Duque, S., Resteghini, L., & Spagnolini, U. (2026). Mobile radio networks and weather radars dualism: Rainfall measurement revolution in densely populated areas. arXiv:2603.19153.
van Leth, T. C., Overeem, A., Leijnse, H., & Uijlenhoet, R. (2018). A measurement campaign to assess sources of error in microwave link rainfall estimation. Atmospheric Measurement Techniques, 11(8), 4645-4669.
Widrow, B., & Kollár, I. (2008). Quantization noise: Roundoff error in digital computation, signal processing, control, and communications. Cambridge University Press.
Zinevich, A., Messer, H., & Alpert, P. (2010). Prediction of rainfall intensity measurement errors using commercial microwave communication links. Atmospheric Measurement Techniques, 3(5), 1385-1402.
::end
