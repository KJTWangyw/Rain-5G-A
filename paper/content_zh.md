::title
面向5G-A的短距毫米波链路降雨感知误差预算：<br>长度无关偏移的识别与标定

::abstract
5G-A毫米波链路路径短，降雨衰减常只有零点几分贝，与路径长度无关的衰减误差会按路径长度的倒数放大为雨强误差。本文建立短距链路降雨反演的误差预算，把误差分为湿期长度无关偏移、晴空参考漂移、功率报告量化和导频估计噪声四类，给出量化报告的共享参考偏差和导频共享参考误差的解析式，并利用不同长度链路的跨链路回归分离长度无关偏移。基于OpenMRG数据集中{{n_use}}条26–40 GHz子链路与天气雷达的分析表明，长度无关偏移随雨强由{{W_low}} dB增至{{W_high}} dB，是0.5 km以下链路的主导误差；参考漂移与1 dB报告量化处于0.1–0.3 dB量级；导频估计噪声在毫秒级晴空参考下比后两者小一到两个数量级。以6–7月标定、8月检验，按雨强标定的偏移使0.5 km以下链路的小时雨强均方根误差降低{{ev_red_05}}%。5G-A短链路降雨感知的首要问题是长度无关偏移的标定，导频测量精度不是限制因素。

::keywords
5G-A；通感一体化；毫米波链路；降雨反演；湿天线衰减；误差预算。

# 引言

地面微波链路的雨致衰减提供路径平均的降雨信息。利用商用微波链路（CML）的接收电平测雨，已在城市和国家尺度得到验证（Messer et al., 2006；Leijnse et al., 2007；Overeem et al., 2011, 2013），并形成了干湿判别、基线估计、湿天线修正和雨强转换的处理流程（Overeem et al., 2016；Chwala & Kunstmann, 2019；Graf et al., 2020）。

5G-A为降雨感知提供了更密集的毫米波链路。FR2频段的小站回传和固定无线接入链路多在1 km以内，频率高，对降雨敏感（Han et al., 2020）。通感一体化使基站可以直接利用解调参考信号（DM-RS）等导频的信道估计（Liu et al., 2022；3GPP, 2024），而标准测量报告中的参考信号接收功率（RSRP）以1 dB为步长上报（3GPP, 2026）。这使人自然地把测量分辨率和导频估计精度视为5G-A降雨感知的关键。

路径变短改变了误差的构成。雨致衰减与路径长度成正比，而湿天线衰减、晴空基线误差和功率量化误差与路径长度无关（Kharadly & Ross, 2001；Leijnse et al., 2008；Zinevich et al., 2010）。换算为雨强后，这些误差随路径缩短而放大。已有湿天线研究（Schleiss et al., 2013；Fencl et al., 2019；Pastorek et al., 2022）和误差来源试验（van Leth et al., 2018）主要针对数公里的回传链路或少数实验链路。对于0.1–1 km的链路，各误差项的相对大小以及导频级测量能带来多少改进，还缺少基于大量真实链路的定量结果。

本文研究一个问题：在短距毫米波链路上，哪些与路径长度无关的误差决定降雨反演精度，各占多大比例，主导项能否从网络数据中标定。主要贡献有三点：（1）建立短距链路降雨反演的误差预算，给出量化报告的共享参考偏差解析式和导频共享参考误差的信息下界；（2）用OpenMRG数据集中{{n_use}}条26–40 GHz、{{L_min}}–{{L_max}} km的子链路定量各误差项，表明长度无关偏移随雨强增大并主导短链路误差，参考漂移和报告量化次之，导频估计噪声可以忽略；（3）利用链路长度差异的跨链路回归标定随雨强变化的偏移，并在独立月份中检验其效果。

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

## 2.4 利用链路长度差异识别偏移

同一城市的链路经历相近的降雨，雨致衰减与$kL$成正比，而$W$不随路径变化。按参考雨强（本文用雷达路径平均雨强）把湿期分钟分档，在第$b$档内对全部链路回归

$$\bar{a}_{i,b} = W_b + \beta_b\,\bar{e}_{i,b}$$ {#eq:regress}

^ 式中$\bar{a}_{i,b}$为链路$i$在该档的平均测得衰减，$\bar{e}_{i,b}$为按链路$i$的$k$、$\alpha$和$L$由雷达雨强计算的平均雨致衰减$kLR^{\alpha}$。截距$W_b$即长度无关偏移；参考降雨的乘性偏差只改变斜率$\beta_b$，不改变截距。把各档的$W_b$按档内平均雨强插值为单调函数$W(R)$，雨强由下式求解：

$$W(\hat{R}) + kL\hat{R}^{\alpha} = \hat{A}$$ {#eq:invert}

^ 左端随$\hat{R}$单调增加，解唯一；$\hat{A} < W(0)$时取$\hat{R}$ = 0。

# 数据与方法

## 3.1 OpenMRG数据

OpenMRG数据集（Andersson et al., 2022）包含瑞典哥德堡728条子链路2015年6–8月的发射和接收电平（TSL、RSL），采样间隔10 s，分辨率约1/3 dB，并配有SMHI天气雷达和城市雨量计资料。本文选取24–40 GHz的{{n_sel}}条子链路（实际频率{{f_min}}–{{f_max}} GHz），其中{{n_use}}条（{{n_phys}}条物理链路）的湿期分钟数不少于600，用于分析。链路长度为{{L_min}}–{{L_max}} km，中位数{{L_med}} km，其中1 km以下{{n_short}}条，与5G-A毫米波链路的长度相当。主要设置见表@tab:data。

::table tab:data | 数据与处理设置 | 3.0,11.4
项目 | 设置
链路 | OpenMRG，2015年6–8月；24–40 GHz子链路{{n_sel}}条，分析{{n_use}}条（{{n_phys}}条物理链路）
电平 | TSL、RSL，10 s采样，约1/3 dB分辨率，按1 min平均；衰减取TSL − RSL
参考降雨 | SMHI雷达，2 km、5 min，$Z = 200R^{1.5}$，沿路径20点平均
湿期 | 雷达路径雨强在前后30 min内超过0.1 mm/h；连续湿期构成事件（不短于10 min）
晴空参考 | 事件开始前6 h内干期分钟的平均衰减，不少于60 min，在事件内保持不变
雨强分档 | 0、0.1、0.5、1、2、5、10、200 mm/h
报告量化 | 10 s RSL按1 dB和2 dB取整，各子链路的网格偏置随机
标定与检验 | 6–7月标定，8月检验
区间估计 | 以物理链路为单元的集群自助法，2000次
::end

## 3.2 处理流程

信号电平按1 min平均，衰减取TSL − RSL。雷达反射率按$Z = 200R^{1.5}$换算为雨强，沿路径等间隔取20个点平均。湿期由雷达判别，使各试验的干湿分类完全相同，从而只比较测量误差和偏移。每个事件以开始前6 h内干期分钟的平均衰减作为晴空参考，事件内各分钟减去同一参考得到降雨期衰减。共得到{{n_events}}个事件，持续时间的中位数为{{ev_med_min}} min。

## 3.3 试验设计

（1）长度无关偏移。按式(@eq:regress)对全部湿期分钟和7个雨强档分别回归，以物理链路为单元做集群自助，得到95%区间。另以$k_iL_i$代替$\bar{e}_{i,b}$作不依赖雷达雨强数值的回归，以Theil–Sen回归检验离群链路的影响，并检验湿期边界（15、30、60 min）和参考窗口（3、6、12 h）的影响。

（2）参考漂移。在干期中重现一次降雨事件：以6 h参考预测其后一个完整干期小时的衰减（真值为0），参考窗口结束至检验小时结束的时间取1、2、4、6 h，共{{n_dry_hours}}个检验小时。

（3）报告量化。把10 s RSL按$\Delta$ = 1 dB（RSRP的上报步长）和2 dB（差分L1-RSRP的步长）取整（3GPP, 2026），各子链路的网格偏置$u$随机抽取；量化后的报告用相同的干期分钟重建参考。逐事件比较参考和降雨期的量化误差与式(@eq:bias)的预测值。

（4）导频仿真。在NR DM-RS时域链路中（子载波间隔120 kHz，1024点FFT，66个资源块，固定多径），取信噪比5、15、25 dB和参考长度1、8、80个时隙；每种条件生成400个参考，每个参考复用于128个降雨时隙（衰减0.5 dB），用同一参考下估计值的协方差估计共享分量。

（5）偏移标定的效果。用6–7月数据标定常数偏移和$W(R)$，在8月比较三种处理的小时雨强与雷达路径平均雨强：不扣除偏移、扣除常数偏移、按式(@eq:invert)扣除$W(R)$。对1 dB和2 dB报告重复$W(R)$处理。

# 结果

## 4.1 长度无关偏移主导短链路误差

图@fig:example为一条{{ex_L}} km、{{ex_f}} GHz子链路在一次{{ex_hours}} h降雨中的衰减。该子链路按固定规则选取：其8月小时雨强均方根误差在扣除$W(R)$后的降幅最接近0.5 km以下链路的中位数（{{ex_red_med}}%），事件为其8月雷达累计雨量最大的一次。测得衰减始终明显高于雷达对应的雨致衰减，差值随雨强增大；加上$W(R)$后两者基本吻合。该事件雷达路径累计雨量为{{ex_radar_mm}} mm，链路在不扣除偏移、扣除常数偏移和扣除$W(R)$时分别为{{ex_link_none_mm}}、{{ex_link_const_mm}}和{{ex_link_rain_mm}} mm。

跨链路回归给出湿期平均的长度无关偏移$W$ = {{W}} dB（95%区间{{W_ci}} dB），斜率$\beta$ = {{beta}}（{{beta_ci}}）（图@fig:offset(a)）。不依赖雷达雨强数值的回归给出{{W_rf}} dB（{{W_rf_ci}} dB），Theil–Sen回归给出{{W_ts}} dB。26–30、32–34和37–40 GHz三个频段分别为{{W_26_30}}、{{W_32_34}}和{{W_37_40}} dB，37–40 GHz最大。改变湿期边界和参考窗口时，$W$在{{W_sens}} dB之间变化。湿期平均衰减中，这一偏移在0.5 km以下链路占{{share_05}}%，0.5–1 km占{{share_1}}%，4 km以上仍占{{share_4}}%。

![一条短链路的降雨事件（子链路{{ex_sublink}}，{{ex_date}}）。(a) 衰减；(b) 雨强。灰底为湿期。](fig1_example.png){#fig:example width=13.2}

偏移随雨强显著增大（图@fig:offset(b)）：雷达雨强低于{{bin_low}} mm/h时为{{W_low}} dB，超过{{bin_high}} mm/h时达{{W_high}} dB（95%区间{{W_high_ci}} dB）。1 mm/h时$W$ = {{W_at_1}} dB，是0.5 km以下链路雨致衰减的{{W_over_rain_1_05}}倍。只有雨强超过一定值，雨致衰减才大于偏移；这一临界雨强的中位数在0.5 km以下链路为{{Rc_05}} mm/h，0.5–1 km为{{Rc_1}} mm/h，1–2 km为{{Rc_2}} mm/h，4 km以上仅{{Rc_20}} mm/h。在除最高档以外的各档，按$k_iL_i$回归和Theil–Sen回归的截距与雷达回归相差不超过{{W_alt_maxdiff_low}} dB。

![长度无关偏移。(a) 各子链路湿期平均衰减与雷达对应雨致衰减，直线为跨链路回归，截距即$W$；(b) 各雨强档的截距$W(R)$（阴影为95%区间），细线为28 GHz垂直极化下0.3、1、3 km路径的雨致衰减$kLR^{\alpha}$。](fig2_offset.png){#fig:offset width=14.5}

## 4.2 参考漂移与报告量化

干期检验中，参考误差的小时平均标准差（各子链路的中位数）在参考结束1 h后为{{drift_1h}} dB，6 h后增至{{drift_6h}} dB（图@fig:drift(a)），平均偏差仅{{drift_bias}} dB。逐分钟均方根误差为{{drift_min}} dB，与小时平均值接近，说明这一误差主要来自缓慢漂移：它在事件内近似共享，不能靠时间平均消除。

1 dB报告下，参考的量化误差均方根为{{q1_ref_rms}} dB，降雨期为{{q1_ev_rms}} dB；小雨时两者落在同一量化间隔内，相互部分抵消，事件级衰减误差的均方根为{{q1_tot_rms}} dB。用式(@eq:bias)逐事件预测参考误差，回归斜率为{{q1_slope}}，$R^2$ = {{q1_r2}}（图@fig:drift(b)）。2 dB报告的对应值为{{q2_ref_rms}}、{{q2_ev_rms}}和{{q2_tot_rms}} dB，斜率{{q2_slope}}，$R^2$ = {{q2_r2}}。1 dB时预测的均方根（{{q1_ref_pred}} dB）略小于实测值，差异来自窗口内电平分布偏离高斯以及原始数据自身约1/3 dB的取整。

![参考漂移与报告量化。(a) 干期检验中参考误差随参考时效的变化，阴影为各子链路的四分位范围；(b) 1 dB报告下各事件参考量化误差的实测值与式(@eq:bias)预测值，红线为1:1线。](fig3_drift_quant.png){#fig:drift width=14.5}

## 4.3 导频测量的共享参考误差

图@fig:pilot比较了仿真得到的共享参考误差与闭式结果。9种条件下，仿真值与闭式值相差{{pilot_rel_range}}，闭式值均落在仿真的95%区间内。15 dB时，1 ms晴空导频的共享误差标准差为{{pilot_sd_15_1ms}} dB，10 ms时为{{pilot_sd_15_10ms}} dB；5 dB、10 ms时为{{pilot_sd_5_10ms}} dB。毫秒级参考下，导频测量的共享误差比参考漂移和1 dB量化误差小一到两个数量级，比长度无关偏移小两个数量级以上。

![导频测量的共享参考误差随晴空参考时长的变化。线为闭式结果，点为NR DM-RS链路仿真（误差线为95%区间）；水平线为OpenMRG中的长度无关偏移（1 mm/h）、1 h参考漂移和1 dB报告的事件级量化误差。](fig4_pilot.png){#fig:pilot width=10.5}

## 4.4 误差预算与偏移标定的效果

表@tab:budget以1 mm/h雨强为例，把各误差项换算为等效雨强。0.5 km以下链路上，雨致衰减只有{{b_L05_signal_db}} dB，长度无关偏移相当于{{b_L05_W1_mmh}} mm/h；扣除平均偏移后，链路间残余偏移的离散仍相当于{{b_L05_resid_mmh}} mm/h；参考漂移和1 dB量化分别相当于{{b_L05_drift_mmh}}和{{b_L05_quant_mmh}} mm/h；导频噪声只相当于{{b_L05_pilot_mmh}} mm/h。路径越长，各项越小。

::table tab:budget | 1 mm/h时各误差项的等效雨强（mm/h，各子链路的中位数） | 1.4,1.45,1.8,2.3,1.8,1.8,1.8,1.8
长度 (km) | 链路数 | 雨致衰减 (dB) | 长度无关偏移 | 残余偏移 | 参考漂移 | 1 dB量化 | 导频噪声
≤0.5 | {{b_L05_n}} | {{b_L05_signal_db}} | {{b_L05_W1_mmh}} | {{b_L05_resid_mmh}} | {{b_L05_drift_mmh}} | {{b_L05_quant_mmh}} | {{b_L05_pilot_mmh}}
0.5–1 | {{b_L1_n}} | {{b_L1_signal_db}} | {{b_L1_W1_mmh}} | {{b_L1_resid_mmh}} | {{b_L1_drift_mmh}} | {{b_L1_quant_mmh}} | {{b_L1_pilot_mmh}}
1–2 | {{b_L2_n}} | {{b_L2_signal_db}} | {{b_L2_W1_mmh}} | {{b_L2_resid_mmh}} | {{b_L2_drift_mmh}} | {{b_L2_quant_mmh}} | {{b_L2_pilot_mmh}}
2–4 | {{b_L4_n}} | {{b_L4_signal_db}} | {{b_L4_W1_mmh}} | {{b_L4_resid_mmh}} | {{b_L4_drift_mmh}} | {{b_L4_quant_mmh}} | {{b_L4_pilot_mmh}}
>4 | {{b_L20_n}} | {{b_L20_signal_db}} | {{b_L20_W1_mmh}} | {{b_L20_resid_mmh}} | {{b_L20_drift_mmh}} | {{b_L20_quant_mmh}} | {{b_L20_pilot_mmh}}
::end

^ 注：长度无关偏移取$W$(1 mm/h) = {{W_at_1}} dB；残余偏移取1–2 mm/h档内链路间的离散{{resid_1}} dB；参考漂移为1 h干期检验的标准差；1 dB量化为事件级误差的均方根；导频噪声为15 dB、10 ms参考下的共享误差{{pilot_sd_15_10ms}} dB。

8月的独立检验中（图@fig:eval），扣除$W(R)$使0.5 km以下链路的小时雨强均方根误差由{{ev_rmse_none_05}} mm/h降至{{ev_rmse_rain_05}} mm/h，降低{{ev_red_05}}%；0.5–1 km和1–2 km分别降低{{ev_red_1}}%和{{ev_red_2}}%。扣除常数偏移只能降至{{ev_rmse_constant_05}} mm/h。0.5 km以下链路的累计雨量比（链路/雷达）由{{ev_ratio_none_05}}降至{{ev_ratio_rain_05}}，4 km以上由{{ev_ratio_none_4}}降至{{ev_ratio_rain_4}}。小时相关系数基本不变（0.5 km以下为{{ev_r_none_05}}和{{ev_r_rain_05}}），偏移造成的主要是系统偏差。扣除$W(R)$后，1 dB报告使各长度组的均方根误差变化{{q1_inc_range}}，2 dB报告变化{{q2_inc_range}}。短链路的累计雨量仍然偏高（0.5 km以下为{{ev_ratio_rain_05}}倍），说明残余偏移和参考漂移（表@tab:budget）在短链路上仍不可忽略。

![8月的小时雨强检验。(a) 均方根误差；(b) 链路与雷达的累计雨量比。点为各长度组子链路的中位数；偏移在6–7月标定。](fig5_evaluation.png){#fig:eval width=14.5}

# 讨论

## 5.1 对5G-A降雨感知的含义

在0.1–1 km的毫米波链路上，误差的主次依次为：随雨强变化的长度无关偏移，链路间的残余偏移和参考漂移，功率报告量化，导频估计噪声。前三类都在事件内共享，提高测量分辨率不能消除。偏移标定后，1 dB报告对小时雨强的影响只有几个百分点，2 dB报告在短链路上的影响增至十几个百分点。导频级测量去除了报告量化，其噪声在毫秒级参考下可以忽略；进一步提高导频估计精度，对雨强精度几乎没有贡献。

由式(@eq:req)，在雨强$R$处要达到雨强误差$\varepsilon$，允许的共享衰减误差约为$\alpha kLR^{\alpha-1}\varepsilon$。28 GHz、1 mm/h、$\varepsilon$ = 1 mm/h时，0.3 km链路允许{{sigma_req_03}} dB，1 km链路允许{{sigma_req_10}} dB；而本文测得的1 h参考漂移为{{drift_1h}} dB，残余偏移离散为{{resid_1}} dB。因此，0.5 km以下链路的降雨感知需要把偏移和漂移控制在0.1 dB以下，这比提高导频精度更关键。

5G-A基站连接距离不同的终端和回传节点，网络中天然存在路径长度的差异。式(@eq:regress)的截距不受参考降雨乘性偏差的影响；对常数偏移，不依赖雷达雨强数值的回归给出相同结果（{{W_rf}} dB对{{W}} dB）。这为在5G-A网络中自标定长度无关偏移提供了条件。

## 5.2 局限与后续工作

（1）OpenMRG为回传链路设备，其天线和天线罩与5G-A有源阵列不同，偏移的量值需要在5G-A链路上重新标定，但式(@eq:regress)的方法可以直接使用。（2）雷达的2 km、5 min采样和$Z$–$R$关系误差进入$\beta_b$和8月检验；截距对乘性偏差不敏感，但会受加性偏差影响。（3）干期检验中的漂移包含CML接收电平检测的误差，导频测量能否减小这一部分，需要5G-A实测数据检验。（4）数据只覆盖一个城市的一个夏季；6–7月标定的偏移在最高雨强档为{{W_cal_top}} dB，低于全时段的{{W_high}} dB。

# 结论

（1）短距毫米波链路降雨反演的主导误差是随雨强增大的长度无关偏移。OpenMRG的{{n_use}}条26–40 GHz子链路中，偏移由{{W_low}} dB增至{{W_high}} dB，在0.5 km以下链路占湿期平均衰减的{{share_05}}%。

（2）参考漂移（1 h为{{drift_1h}} dB）和1 dB报告量化（事件级{{q1_tot_rms}} dB）是次要的共享误差，其中量化误差主要来自晴空参考，可由式(@eq:bias)逐事件预测。毫秒级参考下导频测量的共享误差约为{{pilot_sd_15_10ms}} dB，不构成限制。

（3）利用链路长度差异的跨链路回归可以标定随雨强变化的偏移。在独立月份中，它使0.5 km以下链路的小时雨强均方根误差降低{{ev_red_05}}%。5G-A短链路降雨感知应优先标定长度无关偏移并控制参考漂移。

#* 数据与代码可用性

OpenMRG数据集（Andersson et al., 2022）以CC BY-SA 4.0许可发布于https://doi.org/10.5281/zenodo.7107689。全部分析代码、配置、汇总结果和排版脚本见https://github.com/KJTWangyw/Rain-5G-A（分支paper/short-link-error-budget），文中数字均由代码生成。

#* 参考文献

::refs
3GPP. (2024). NR; Physical channels and modulation (3GPP TS 38.211 version 18.3.0 Release 18). ETSI TS 138 211 V18.3.0.
3GPP. (2026). NR; Requirements for support of radio resource management (3GPP TS 38.133 version 18.15.0 Release 18). ETSI TS 138 133 V18.15.0.
Andersson, J. C. M., Olsson, J., van de Beek, R., & Hansryd, J. (2022). OpenMRG: Open data from Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden. Earth System Science Data, 14(12), 5411-5426.
Chwala, C., & Kunstmann, H. (2019). Commercial microwave link networks for rainfall observation: Assessment of the current status and future challenges. WIREs Water, 6(2), e1337.
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
Schleiss, M., Rieckermann, J., & Berne, A. (2013). Quantification and modeling of wet-antenna attenuation for commercial microwave links. IEEE Geoscience and Remote Sensing Letters, 10(5), 1195-1199.
van Leth, T. C., Overeem, A., Leijnse, H., & Uijlenhoet, R. (2018). A measurement campaign to assess sources of error in microwave link rainfall estimation. Atmospheric Measurement Techniques, 11(8), 4645-4669.
Widrow, B., & Kollár, I. (2008). Quantization noise: Roundoff error in digital computation, signal processing, control, and communications. Cambridge University Press.
Zinevich, A., Messer, H., & Alpert, P. (2010). Prediction of rainfall intensity measurement errors using commercial microwave communication links. Atmospheric Measurement Techniques, 3(5), 1385-1402.
::end
