::title
基于跨链路回归与量化重放的短距毫米波链路测雨误差预算

!! 中文核验修订稿　2026年11月5日。实测分析固定于提交3174afd。链路筛选、回归及主要留月指标已由归档中间表复算。量化分项的共同有效支持已修复，相关历史统计和图3右图仍需原始时序重跑核验。

::abstract
短距离毫米波链路的雨致衰减较小，端点附加损耗和共享晴空参考误差会随路径缩短而放大为雨强误差。本文结合跨链路截距回归、功率量化重放和导频参考噪声分析，构建湿期偏移、参考漂移、报告量化及热噪声的衰减域误差预算。OpenMRG中{{n_use}}条有向子链路对应{{n_phys}}条物理路径，实际频率{{f_min}}–{{f_max}} GHz，长度{{L_min}}–{{L_max}} km。全季湿期回归截距为{{W}} dB，95%物理路径簇自助区间为{{W_ci}} dB。长度不超过0.5 km的50条子链路中，固定截距与湿期平均衰减之比的中位数为{{share_05}}%。按雷达雨强分档的截距由{{W_low}}增至{{W_high}} dB。干期参考误差在1 h后的标准差中位数为{{drift_1h}} dB。在792导频、15 dB逐导频信噪比及10 ms相干参考的理想模型下，参考热噪声标准差约为{{pilot_sd_15_10ms}} dB。采用6–7月标定和8月同网络留出，在雷达辅助干湿分类条件下，雨强依赖偏移校正使最短链路组小时RMSE中位数由{{ev_rmse_none_05}}降至{{ev_rmse_rain_05}} mm/h，累计雨量比中位数仍为{{ev_ratio_rain_05}}。结果支持将湿期偏移与参考稳定性作为短链路误差控制的主要对象，并以独立降雨观测检验校正后的剩余偏差。

::keywords
毫米波链路、降雨反演、误差预算、湿期偏移、共享参考、量化重放

# 引言

地面微波链路的雨致衰减提供路径平均的降雨信息。利用商用微波链路（CML）的接收电平测雨，已在城市和国家尺度得到验证（Messer et al., 2006，Leijnse et al., 2007，Overeem et al., 2011, 2013），并形成了干湿判别、基线估计、湿天线修正和雨强转换的处理流程（Overeem et al., 2016，Chwala & Kunstmann, 2019，Graf et al., 2020）。

毫米波小站回传和固定无线接入提供了用于降雨感知的短距离传播路径，高频链路对降雨具有较高灵敏度（Han et al., 2020）。通感一体化使基站可以直接利用解调参考信号（DM-RS）等导频的信道估计（Liu et al., 2022，3GPP, 2024），绝对L1 SS-RSRP/CSI-RSRP与L3 SS-RSRP报告采用1 dB分辨率，差分L1报告采用2 dB映射（3GPP, 2026）。这些接口为比较报告步长和细粒度幅度观测提供了测量背景。

路径变短改变了误差的构成。雨致衰减与路径长度成正比，而湿天线衰减、晴空基线误差和功率量化误差与路径长度无关（Kharadly & Ross, 2001，Leijnse et al., 2008，Zinevich et al., 2010）。换算为雨强后，这些误差随路径缩短而放大。已有湿天线研究（Schleiss et al., 2013，Fencl et al., 2019，Pastorek et al., 2022）和误差来源试验（van Leth et al., 2018）揭示了设备与环境对反演的影响。本文聚焦百米至公里级链路，在同一网络的处理支持下比较各误差分量，并检验偏移校正的跨月效果。

本文考察短距毫米波链路中各类长度无关误差的相对贡献，以及网络长度差异能够支持的偏移标定效果。此前研究已系统分析长度、频率、功率分辨率及湿天线误差（Leijnse et al., 2008），利用长度回归截距研究雾条件下的天线损耗（David et al., 2013），并实测检验过极值记录的量化偏差（Ostrometzky & Eshel, 2018）。在这些工作基础上，本文的贡献集中于三个方面。

（1）在相同衰减坐标中组织湿期统计截距、参考漂移、事件量化误差与理想导频参考噪声，并区分各量的偏差、标准差和均方根口径。

（2）利用{{n_use}}条有向子链路的长度差异及雨强分档，刻画湿期截距的条件变化，给出短链路偏移占比及其统计支持。

（3）在6–7月标定偏移函数、8月留出的设置下评价校正收益，同时报告剩余水量偏差与检测率，为通信分辨率改进和传播偏移校正的优先次序提供量化依据。

# 误差模型

## 2.1 短路径对衰减误差的放大

降雨事件中的衰减相对于事件前的晴空参考计算。令第$t$分钟的同步损耗为$\hat D_t=\mathrm{TSL}_t-\mathrm{RSL}_t$，事件前干期均值为$\hat B$，则附加衰减为

$$\hat{A}_t = \hat D_t - \hat B = kLR_t^{\alpha} + W_t + d_t + \eta_{\mathrm{ref}} - \eta_t$$ {#eq:decomp}

^ 式中，$kLR_t^{\alpha}$为长度$L$（km）的路径上雨强$R_t$（mm/h）引起的衰减，$k$和$\alpha$由频率和极化按ITU-R P.838-3确定（ITU-R, 2005），$W_t$为既定基线规则下的湿期统计偏移，天线水膜、天线罩湿润及湿期设备变化均可贡献该分量，$d_t$为晴空同步损耗相对参考时段的漂移。$\eta$定义为RSL测量误差减去TSL测量误差，$\eta_{\mathrm{ref}}$和$\eta_t$分别对应参考均值和当前观测，因此对衰减的贡献为参考项减当前项。

同一事件中，$\eta_{\mathrm{ref}}$对所有分钟相同。$W_t$和$d_t$可具有较强时间相关性，其平均效果取决于变化尺度，需保留相关项进行误差传播。一个衰减误差$\delta$（dB）相当于雨强

$$R_{\mathrm{eq}}(\delta) = \left(\frac{\delta}{kL}\right)^{1/\alpha}$$ {#eq:req}

^ 该量表示把正衰减尺度完全解释为降雨时得到的无雨等效雨强。在非零雨强$R$处，确定性正偏移的增量为$\Delta R=(R^{\alpha}+\delta/(kL))^{1/\alpha}-R$。28 GHz垂直极化时$k$ = {{k28}}，$\alpha$ = {{a28}}，0.5 km链路上0.1 dB的误差即相当于{{req_05_01}} mm/h的雨强。

## 2.2 量化报告的共享参考偏差

功率报告把电平取整到步长为$\Delta$的网格。

$$Q(x) = u + \Delta\,\operatorname{round}\!\left(\frac{x-u}{\Delta}\right)$$ {#eq:quant}

^ 式中$u$为每条子链路固定的网格偏置。重放使用NumPy最近邻取整，半步取最近偶数。1 dB与2 dB仅用于受控分辨率比较，真实接口的解码区间、饱和和过滤需另行匹配。若窗口内电平服从均值$c$、标准差$\sigma$的高斯边际分布，取整误差的期望为（Widrow et al., 1996）

$$b(c;\sigma) = \sum_{m=1}^{\infty}\frac{(-1)^{m}\Delta}{\pi m}\,\mathrm{e}^{-2\pi^{2}m^{2}\sigma^{2}/\Delta^{2}}\sin\frac{2\pi m(c-u)}{\Delta}$$ {#eq:bias}

^ $b$随$c$在网格中的位置周期变化，起伏越小，衰减越慢。对随机的网格偏置，其均方根为

$$\sigma_{b}(\sigma) = \left[\sum_{m=1}^{\infty}\frac{\Delta^{2}}{2\pi^{2}m^{2}}\,\mathrm{e}^{-4\pi^{2}m^{2}\sigma^{2}/\Delta^{2}}\right]^{1/2}$$ {#eq:brms}

^ $\sigma$ = 0时$\sigma_b = \Delta/\sqrt{12}$，$\sigma > 0.42\Delta$时$\sigma_b < 0.01\Delta$。晴空参考窗口内电平平稳（本文数据中$\sigma$的中位数为{{sig_ref_med}} dB），平均后的取整误差期望可持续偏离零。仅考虑RSL量化时，参考量化项的期望由$b(c_{\mathrm{ref}};\sigma_{\mathrm{ref}})$给出，一次有限窗口实现还包含与有效样本量及时间相关性有关的波动，并在整个事件中共享。降雨期电平随雨强变化（$\sigma$的中位数为{{sig_ev_med}} dB），更宽的电平分布可减小量化误差的期望。数值实现使用60个谐波，零起伏时的闭式极限为$\Delta/\sqrt{12}$。

## 2.3 导频测量的共享参考误差

NR下行DM-RS配置类型1位于OFDM符号2和11，每隔一个子载波放置一个导频（3GPP, 2024），66个资源块时每个时隙有$J$ = {{J}}个导频。用$M$个晴空时隙平均得到参考$r_j = h_j + e_{0j}$，降雨时隙的最小二乘导频估计为$z_{tj} = a_t\exp(\mathrm{i}\phi_t)h_j + e_{tj}$，其中幅度$a_t$满足$A_t = -20\log_{10}a_t$，$\phi_t$为未知公共相位，$e_{0j}$和$e_{tj}$为方差$v/M$和$v$的复高斯噪声。在信道固定、参考相干对齐且各导频噪声独立的条件下，对共同衰减参数进行$n$次观测，局部无偏估计的Cramér–Rao方差下界为（Kay, 1993）

$$\operatorname{Var}(\bar{A}_n) \ge \frac{K^{2}}{2J}\left(\frac{1}{\gamma_{\mathrm{r}}} + \frac{1}{n\gamma_{\mathrm{t}}}\right),\quad K = \frac{20}{\ln 10}$$ {#eq:crb}

^ 式中$\gamma_{\mathrm{r}} = M\,\mathrm{SNR}$和$\gamma_{\mathrm{t}} = a^{2}\mathrm{SNR}$分别为单导频的参考信噪比和当前信噪比（信道平均功率归一化为1）。第一项为共享参考误差，不随$n$减小。子载波间隔120 kHz时每个时隙0.125 ms，参考时长为$T_{\mathrm{ref}}$时，共享误差的标准差为

$$\sigma_{\mathrm{ref}} = \frac{K}{\sqrt{2J\gamma_{\mathrm{r}}}},\quad \gamma_{\mathrm{r}} = \frac{T_{\mathrm{ref}}}{0.125\ \mathrm{ms}}\,\mathrm{SNR}$$ {#eq:pilotref}

^ 对剖面幅度估计器作一阶传播，其共享分量为$K^2/(2J\gamma_{\mathrm{r}}) + K^2/[J(\gamma_{\mathrm{r}} + \gamma_{\mathrm{t}})^2]$，在$\gamma_{\mathrm{r}} \gg 1$时与下界一致。该模型描述幅度保真的复导频输入，实际测量还取决于波束、功率与接收增益校正。

## 2.4 利用链路长度差异识别偏移

多长度回归利用路径雨衰与$kL$的关系估计组内截距。给定雨强档、频率和极化后，其物理解释要求未建模降雨暴露、设备及参考误差在回归横轴上的条件均值为零。按参考雨强（本文用雷达路径平均雨强）把湿期分钟分档，在第$b$档内对全部链路回归

$$\bar{a}_{i,b} = W_b + \beta_b\,\bar{e}_{i,b}$$ {#eq:regress}

^ 式中$\bar{a}_{i,b}$为链路$i$在该档的平均测得衰减，$\bar{e}_{i,b}$为按链路$i$的$k$、$\alpha$和$L$由雷达路径平均雨强计算的平均雨致衰减$kL\bar R^{\alpha}$。截距$W_b$表示该统计模型中的长度无关分量。仅当横轴整体受到共同的常数比例缩放、有效样本和分档均保持固定时，OLS截距才具有尺度不变性。实际雷达雨强还用于干湿判别和分档，改变其标定可能同时改变样本支持。把各档的$W_b$按档内平均雨强插值为单调函数$W(R)$，雨强由下式求解。

$$W(\hat{R}) + kL\hat{R}^{\alpha} = \hat{A}$$ {#eq:invert}

^ 左端随$\hat{R}$单调增加，解唯一，$\hat{A} < W(0)$时取$\hat{R}$ = 0。

# 数据与方法

## 3.1 OpenMRG数据

OpenMRG数据集（Andersson et al., 2022）包含瑞典哥德堡728条子链路2015年6–8月的发射和接收电平（TSL、RSL），采样间隔10 s，分辨率约1/3 dB，并配有SMHI天气雷达和城市雨量计资料。本文选取24–40 GHz的{{n_sel}}条子链路（实际频率{{f_min}}–{{f_max}} GHz），其中{{n_use}}条（{{n_phys}}条物理链路）的湿期分钟数不少于600，用于分析。链路长度为{{L_min}}–{{L_max}} km，中位数{{L_med}} km，其中长度不超过1 km的{{n_short}}条，与5G-A毫米波链路的长度相当。主要设置见表@tab:data。

::table tab:data | 数据与处理设置 | 3.0,11.4
项目 | 设置
链路 | OpenMRG，2015年6–8月，24–40 GHz子链路{{n_sel}}条，分析{{n_use}}条（{{n_phys}}条物理链路）
电平 | TSL、RSL，10 s采样，约1/3 dB分辨率，按1 min平均，衰减取TSL − RSL
参考降雨 | SMHI雷达，2 km、5 min，$Z = 200R^{1.5}$，沿路径20点平均
湿期 | 雷达路径雨强在前后30 min内超过0.1 mm/h，连续湿期构成事件（不短于10 min）
晴空参考 | 事件开始前6 h内干期分钟的平均衰减，不少于60 min，在事件内保持不变
雨强分档 | 0、0.1、0.5、1、2、5、10、200 mm/h
报告量化 | 10 s RSL按1 dB和2 dB取整，各子链路的网格偏置随机
标定与检验 | 6–7月标定，8月检验
区间估计 | 以物理链路为单元的集群自助法，2000次
::end

## 3.2 处理流程

信号电平按1 min平均，衰减取TSL − RSL。雷达反射率按$Z = 200R^{1.5}$换算为雨强，沿路径等间隔取20个点平均。湿期由雷达判别，使各试验的干湿分类完全相同，从而只比较测量误差和偏移。每个事件以开始前6 h内干期分钟的平均衰减作为晴空参考，事件内各分钟减去同一参考得到降雨期衰减。共得到{{n_events}}条子链路湿段记录，持续时间的中位数为{{ev_med_min}} min。共享同一天气过程的多条链路分别计数，独立天气事件数需另行聚类确定。实际实现先对雷达路径采样点的雨强取均值，再作幂律换算，路径雨强非均匀性带来的近似误差需通过空间积分对照评价。

## 3.3 试验设计

（1）长度无关偏移。按式(@eq:regress)对全部湿期分钟和7个雨强档分别回归，以物理链路为单元做2000次集群自助，双向子链路共同重采样，得到95%区间。该区间反映当前季节记录条件下的路径抽样变化，跨天气过程的不确定性需增加事件层重采样。另以$k_iL_i$代替$\bar{e}_{i,b}$作长度自变量回归，沿用雷达辅助的湿段、参考及雨强档，以Theil–Sen回归检验离群链路的影响，并检验湿期边界（15、30、60 min）和参考窗口（3、6、12 h）的影响。

（2）参考漂移。在干期中重现一次降雨事件。以6 h参考预测其后一个完整干期小时的衰减（真值为0），参考窗口结束至检验小时结束的时间取1、2、4、6 h，共{{n_dry_hours}}个检验小时。

（3）报告量化。把10 s RSL按$\Delta$ = 1 dB（RSRP的上报步长）和2 dB（差分L1-RSRP的步长）取整（3GPP, 2026），各子链路的网格偏置$u$随机抽取，量化后的报告用相同的干期分钟重建参考。逐事件比较参考和降雨期的量化误差与式(@eq:bias)的预测值。

（4）导频仿真。在NR DM-RS时域链路中（子载波间隔120 kHz，1024点FFT，66个资源块，固定多径），取信噪比5、15、25 dB和参考长度1、8、80个时隙，每种条件生成400个参考，每个参考复用于128个恒定衰减时隙（0.5 dB），用同一参考下估计值的协方差估计共享分量。

（5）偏移标定的效果。用6–7月数据标定常数偏移和$W(R)$，在8月比较三种处理的小时雨强与雷达路径平均雨强。不扣除偏移、扣除常数偏移、按式(@eq:invert)扣除$W(R)$。原生、1 dB和2 dB方案分别用6–7月资料重新标定常数与$W(R)$。雷达干窗的输出直接设为零，干湿掩码使用前后30 min的居中窗口。因此，该检验评价雷达辅助分类条件下的离线幅度恢复。

## 3.4 参数依据与敏感性

24–40 GHz为固定筛选区间，600个可用湿分钟为子链路最低数据支持，每个雨强档至少保留30个有效分钟。0.1 mm/h阈值、前后30 min湿期扩展及6 h参考窗构成主处理配置。15、30、60 min扩展和3、6、12 h参考窗的对照产生{{W_sens}} dB的截距范围，反映处理参数对结果的影响。该敏感性范围与固定配置的抽样区间分别报告。6–7月和8月按时间顺序划分，用于检验同网络跨月迁移。

1 dB与2 dB步长参照绝对与差分RSRP报告的分辨率，网格偏置使用固定随机种子20261105，各分辨率分别重标偏移函数。后续通过多种网格偏置评估这一随机设置的影响。导频数$J=66\times2\times6=792$来自所选资源块与DM-RS映射。5、15、25 dB和0.125、1、10 ms参考时长用于比较噪声与参考积累尺度，每个条件使用400次独立参考。这些数值定义仿真对照的范围，实际接收机配置需由硬件测量对应。

# 结果

## 4.1 长度无关偏移主导短链路误差

图@fig:example为一条{{ex_L}} km、{{ex_f}} GHz子链路在一次{{ex_hours}} h降雨中的衰减。该子链路按固定规则选取。其8月小时雨强均方根误差在扣除$W(R)$后的降幅最接近长度不超过0.5 km的链路的中位数（{{ex_red_med}}%），事件为其8月雷达累计雨量最大的一次。测得衰减始终明显高于雷达对应的雨致衰减，差值随雨强增大，加上$W(R)$后两者基本吻合。该事件雷达路径累计雨量为{{ex_radar_mm}} mm，链路在不扣除偏移、扣除常数偏移和扣除$W(R)$时分别为{{ex_link_none_mm}}、{{ex_link_const_mm}}和{{ex_link_rain_mm}} mm。

跨链路回归给出湿期平均的长度无关偏移$W$ = {{W}} dB（95%区间{{W_ci}} dB），斜率$\beta$ = {{beta}}（{{beta_ci}}）（图@fig:offset(a)）。不依赖雷达雨强数值的回归给出{{W_rf}} dB（{{W_rf_ci}} dB），Theil–Sen回归给出{{W_ts}} dB。26–30、32–34和37–40 GHz三个频段分别为{{W_26_30}}、{{W_32_34}}和{{W_37_40}} dB，37–40 GHz最大。改变湿期边界和参考窗口时，$W$在{{W_sens}} dB之间变化。对每条子链路计算$W/\bar a_i$后，长度组中位数在长度不超过0.5 km的链路为{{share_05}}%，0.5–1 km占{{share_1}}%，4 km以上仍占{{share_4}}%。

![一条短链路的降雨事件（子链路{{ex_sublink}}，{{ex_date}}）。(a) 衰减，(b) 雨强。灰底为湿期。](fig1_example.png){#fig:example width=13.2}

偏移随雨强显著增大（图@fig:offset(b)）。雷达雨强低于{{bin_low}} mm/h时为{{W_low}} dB，超过{{bin_high}} mm/h时达{{W_high}} dB（95%区间{{W_high_ci}} dB）。1 mm/h时$W$ = {{W_at_1}} dB，是长度不超过0.5 km的链路雨致衰减的{{W_over_rain_1_05}}倍。只有雨强超过一定值，雨致衰减才大于偏移，这一临界雨强的中位数在长度不超过0.5 km的链路为{{Rc_05}} mm/h，0.5–1 km为{{Rc_1}} mm/h，1–2 km为{{Rc_2}} mm/h，4 km以上仅{{Rc_20}} mm/h。在除最高档以外的各档，按$k_iL_i$回归和Theil–Sen回归的截距与雷达回归相差不超过{{W_alt_maxdiff_low}} dB。

![长度无关偏移。(a) 各子链路湿期平均衰减与雷达对应雨致衰减，直线为跨链路回归，截距即$W$，(b) 各雨强档的截距$W(R)$（阴影为95%区间），细线为28 GHz垂直极化下0.3、1、3 km路径的雨致衰减$kLR^{\alpha}$。](fig2_offset.png){#fig:offset width=14.5}

## 4.2 参考漂移与报告量化

干期检验中，参考误差的小时平均标准差（各子链路的中位数）在参考结束1 h后为{{drift_1h}} dB，6 h后增至{{drift_6h}} dB（图@fig:drift(a)），1 h时各子链路平均偏差的中位数约为−{{drift_bias}} dB。逐分钟均方根误差为{{drift_min}} dB，与小时平均值接近，提示存在较强的缓慢变化分量，其时间相关性需在误差传播中保留。

按原归档脚本的RSL有效样本支持，1 dB报告下参考电平量化增量的事件均值均方根为{{q1_ref_rms}} dB，降雨期为{{q1_ev_rms}} dB，小雨时两者落在同一量化间隔内，相互部分抵消，两端增量之差的事件级均方根为{{q1_tot_rms}} dB。用式(@eq:bias)逐事件预测参考误差，回归斜率为{{q1_slope}}，$R^2$ = {{q1_r2}}（图@fig:drift(b)）。2 dB报告的对应值为{{q2_ref_rms}}、{{q2_ev_rms}}和{{q2_tot_rms}} dB，斜率{{q2_slope}}，$R^2$ = {{q2_r2}}。1 dB时预测的均方根（{{q1_ref_pred}} dB）略小于实测值，该差异可能同时受到电平分布、有限有效样本、时间相关性和原生约1/3 dB分辨率影响。

!! 核验发现，原量化分项按RSL有效性统计，实际雨衰同时要求TSL与RSL有效。已在分析代码中统一共同支持，并用缺测反例完成回归测试。本节0.19、0.15、0.14 dB及图3右图保留为修复前的归档统计，完成原始时序重跑后更新。

![参考漂移与报告量化。(a) 干期检验中参考误差随参考时效的变化，阴影为各子链路的四分位范围，(b) 1 dB报告下各事件参考量化误差的实测值与式(@eq:bias)预测值，红线为1比1线。右图沿用修复前的RSL有效支持，待共同支持重算。](fig3_drift_quant.png){#fig:drift width=14.5}

## 4.3 导频测量的共享参考误差

图@fig:pilot比较了仿真得到的共享参考误差与闭式结果。9种条件下，仿真值与闭式值相差{{pilot_rel_range}}，闭式值均落在仿真的95%区间内，区间针对共享方差分量。15 dB时，1 ms晴空导频的共享误差标准差为{{pilot_sd_15_1ms}} dB，10 ms时为{{pilot_sd_15_10ms}} dB，5 dB、10 ms时为{{pilot_sd_5_10ms}} dB。毫秒级参考下，导频测量的共享误差比参考漂移和1 dB量化误差小一到两个数量级，与湿期偏移相比处于更低量级。当前观测的热噪声需另行计入。零附加衰减下，单个当前时隙的标准差近似为0.03881 dB，与10 ms参考合成为0.03905 dB。参考与当前各积累80个独立时隙时，联合近似为0.00614 dB。

![导频测量的共享参考误差随晴空参考时长的变化。线为闭式结果，点为NR DM-RS链路仿真（误差线为95%区间），水平线为OpenMRG中的长度无关偏移（1 mm/h）、1 h参考漂移和1 dB报告的事件级量化误差。量化水平线沿用修复前归档值，待共同支持重算。](fig4_pilot.png){#fig:pilot width=10.5}

## 4.4 误差预算与偏移标定的效果

表@tab:budget把各正衰减尺度换算为无雨等效雨强，湿期偏移项使用$W$(1 mm/h)。表中同时列示偏移、标准差和均方根的量级，各项需按其统计身份解释。长度不超过0.5 km的链路上，雨致衰减只有{{b_L05_signal_db}} dB，长度无关偏移相当于{{b_L05_W1_mmh}} mm/h，扣除平均偏移后，链路间残余偏移的离散仍相当于{{b_L05_resid_mmh}} mm/h，参考漂移和1 dB量化分别相当于{{b_L05_drift_mmh}}和{{b_L05_quant_mmh}} mm/h，导频噪声只相当于{{b_L05_pilot_mmh}} mm/h。路径越长，各项越小。

::table tab:budget | 正衰减尺度的无雨等效雨强，偏移取W(1 mm/h)，组中位数，雨强单位mm/h | 1.4,1.45,1.8,2.3,1.8,1.8,1.8,1.8
长度 (km) | 链路数 | 雨致衰减 (dB) | 长度无关偏移 | 残余偏移 | 参考漂移 | 1 dB量化* | 参考热噪声
≤0.5 | {{b_L05_n}} | {{b_L05_signal_db}} | {{b_L05_W1_mmh}} | {{b_L05_resid_mmh}} | {{b_L05_drift_mmh}} | {{b_L05_quant_mmh}} | {{b_L05_pilot_mmh}}
0.5–1 | {{b_L1_n}} | {{b_L1_signal_db}} | {{b_L1_W1_mmh}} | {{b_L1_resid_mmh}} | {{b_L1_drift_mmh}} | {{b_L1_quant_mmh}} | {{b_L1_pilot_mmh}}
1–2 | {{b_L2_n}} | {{b_L2_signal_db}} | {{b_L2_W1_mmh}} | {{b_L2_resid_mmh}} | {{b_L2_drift_mmh}} | {{b_L2_quant_mmh}} | {{b_L2_pilot_mmh}}
2–4 | {{b_L4_n}} | {{b_L4_signal_db}} | {{b_L4_W1_mmh}} | {{b_L4_resid_mmh}} | {{b_L4_drift_mmh}} | {{b_L4_quant_mmh}} | {{b_L4_pilot_mmh}}
>4 | {{b_L20_n}} | {{b_L20_signal_db}} | {{b_L20_W1_mmh}} | {{b_L20_resid_mmh}} | {{b_L20_drift_mmh}} | {{b_L20_quant_mmh}} | {{b_L20_pilot_mmh}}
::end

^ 注。长度无关偏移取$W$(1 mm/h) = {{W_at_1}} dB，残余偏移取1–2 mm/h档内链路间的离散{{resid_1}} dB，参考漂移为1 h干期检验的标准差，1 dB量化为事件级误差的均方根，参考热噪声为15 dB、10 ms参考下的共享项{{pilot_sd_15_10ms}} dB。带星号的量化列来自修复前归档统计，须随共同支持重算更新。

8月同网络的时间留出检验中（图@fig:eval），扣除$W(R)$使长度不超过0.5 km的链路的逐子链路小时RMSE的组中位数由{{ev_rmse_none_05}} mm/h降至{{ev_rmse_rain_05}} mm/h，两个组中位数之比对应降幅{{ev_red_05}}%。0.5–1 km和1–2 km分别降低{{ev_red_1}}%和{{ev_red_2}}%。扣除常数偏移只能降至{{ev_rmse_constant_05}} mm/h。长度不超过0.5 km的链路的累计雨量比（链路/雷达）由{{ev_ratio_none_05}}降至{{ev_ratio_rain_05}}，4 km以上由{{ev_ratio_none_4}}降至{{ev_ratio_rain_4}}。最短组小时相关系数中位数由{{ev_r_none_05}}变为{{ev_r_rain_05}}，检出率中位数由1.00变为0.80。RMSE改善与检出率变化需共同报告。各分辨率分别重标$W(R)$后，1 dB方案相对原生方案的组RMSE中位数变化{{q1_inc_range}}，2 dB报告变化{{q2_inc_range}}。短链路的保留小时均值折算的累计雨量仍然偏高（长度不超过0.5 km的为{{ev_ratio_rain_05}}倍），说明残余偏移和参考漂移（表@tab:budget）在短链路上仍不可忽略。

![雷达辅助干湿分类条件下的8月小时雨强检验。(a) 均方根误差，(b) 链路与雷达的累计雨量比。点为各长度组子链路的中位数，各分辨率分别在6–7月标定。累计量由至少具有45个有效分钟的小时均值相加得到。](fig5_evaluation.png){#fig:eval width=14.5}

# 讨论

## 5.1 对5G-A降雨感知的含义

当前OpenMRG样本和处理支持下，湿期偏移与晴空参考稳定性对短链路误差具有较大影响。不同误差项包含均值偏移、跨链路离散、跨事件均方根和理想热噪声标准差，排序用于资源配置的量级判断。雨强依赖校正改善了同网络留月RMSE，同时最短组仍有约3.1倍水量偏差，说明进一步提升需要针对校正后的剩余误差开展验证。

导频级幅度测量可以提供更细粒度的参考观测。在固定信道、相干参考、15 dB逐导频信噪比和10 ms积累的条件下，参考热噪声约为0.0043 dB。当前时隙噪声、功率与增益校正、波束变化以及实际有效样本量共同决定完整测量精度。将该理论结果迁移到5G-A接收机，需要同设备报告与复导频的配对观测。

长度自变量回归得到{{W_rf}} dB，与雷达雨衰自变量的{{W}} dB接近。两种估计共用雷达辅助的湿段和晴空参考，该一致性支持继续检验长度信息的标定价值。自主运行还需要因果的干湿判别和参考更新，并在相同路径、相同天气支持上检验截距的可迁移性。

## 5.2 局限与后续工作

（1）湿期截距的物理归属与可识别性。后续在相近频率、极化与硬件条件下设置多长度观测，同时记录天线湿润、温度和增益。通过受控天线湿润与自然降雨对照，分离端点水膜、空间降雨不均匀和基线漂移。雷达乘性标定应连同分档、干湿掩码与有效支持一起重算，并采用路径和天气事件两个层级报告区间。

（2）独立降雨验证与共同支持。先完成量化分项的TSL/RSL共同有效支持重跑，再冻结所有参数。增加因果干湿判别，并按天气事件和物理路径留出，对照直接反演、常数偏移及已有湿天线修正方法。以未参与标定的雨量计或雨滴谱仪核验事件水量，雷达用于空间代表性分析。同时报告RMSE、偏差、检出率、干期虚报、有效覆盖率及配对改善区间。当前资料仅覆盖一个城市的一个夏季，跨季节与跨城市检验将直接检验偏移迁移。

（3）真实NR接收机的配对测量。采用可导出复导频和RSRP报告的固定链路，记录发射功率、增益和波束。使用已知衰减台阶与重复晴空参考，分开测定参考误差、当前观测误差和慢漂移，再增加天线湿润与自然降雨。对分辨率方案分别给出固定校正器和重新标定两组对照，以量化观测接口与传播校正各自带来的收益。

# 结论

（1）在本次OpenMRG样本中，湿期长度无关截距为{{W}} dB，物理路径簇自助95%区间为{{W_ci}} dB。50条长度不超过0.5 km的有向子链路中，截距与湿期平均衰减之比的中位数为{{share_05}}%。雨强分档截距由{{W_low}}增至{{W_high}} dB，显示偏移的条件变化需要在校正中保留。

（2）干期参考误差在1 h后具有{{drift_1h}} dB的标准差中位数。共享参考使事件内误差具有相关结构。给定理想导频配置的参考热噪声约为{{pilot_sd_15_10ms}} dB，完整噪声预算还包含当前观测项。量化分项已明确共同支持的计算要求，归档数值需随修复后的原始时序重跑更新。

（3）采用6–7月标定与8月同网络留出，在雷达辅助干湿分类条件下，最短链路组小时RMSE中位数从{{ev_rmse_none_05}}降至{{ev_rmse_rain_05}} mm/h，相对降幅约{{ev_red_05}}%。水量比中位数仍为{{ev_ratio_rain_05}}。这些结果支持优先估计湿期偏移与参考漂移，并以独立水量观测和真实NR配对实验确认剩余误差及测量分辨率收益。

#* 数据与代码可用性

OpenMRG公开数据的DOI为10.5281/zenodo.7107689，采用CC BY-SA 4.0许可。原分析固定提交为3174afd04993fc829858cf5fca0228dd6a351ed1。本轮修订由已归档链路汇总、参考漂移和8月评价表复算主要统计量，并新增量化共同支持回归测试。雨强分档曲线及其区间沿用原分析归档，待进一步复算。原始OpenMRG预处理、修复后的事件量化统计及NR时域仿真需按仓库脚本运行。代码和核验记录位于KJTWangyw/Rain-5G-A的paper/openmrg-error-budget-zh分支。
#* 参考文献

::refs
3GPP. (2024). NR. Physical channels and modulation (3GPP TS 38.211 version 18.3.0 Release 18). ETSI TS 138 211 V18.3.0.
3GPP. (2026). NR. Requirements for support of radio resource management (3GPP TS 38.133 version 15.32.0 Release 15). ETSI TS 138 133 V15.32.0.
Andersson, J. C. M., Olsson, J., van de Beek, R., & Hansryd, J. (2022). OpenMRG: Open data from Microwave links, Radar, and Gauges for rainfall quantification in Gothenburg, Sweden. Earth System Science Data, 14(12), 5411-5426.
Chwala, C., & Kunstmann, H. (2019). Commercial microwave link networks for rainfall observation: Assessment of the current status and future challenges. WIREs Water, 6(2), e1337.
David, N., Alpert, P., & Messer, H. (2013). The potential of commercial microwave networks to monitor dense fog—feasibility study. Journal of Geophysical Research. Atmospheres, 118, 11750–11761. DOI 10.1002/2013JD020346.
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
Ostrometzky, J., & Eshel, A. (2018). Empirical study of the quantization induced bias in commercial microwave links’ min/max attenuation measurements for rain monitoring. Environments, 5(7), 80. DOI 10.3390/environments5070080.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2011). Measuring urban rainfall using microwave links from commercial cellular communication networks. Water Resources Research, 47(12), W12505.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2013). Country-wide rainfall maps from cellular communication networks. Proceedings of the National Academy of Sciences, 110(8), 2741-2745.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2016). Retrieval algorithm for rainfall mapping from microwave links in a cellular communication network. Atmospheric Measurement Techniques, 9(5), 2425-2444.
Pastorek, J., Fencl, M., Rieckermann, J., & Bareš, V. (2022). Precipitation estimates from commercial microwave links: Practical approaches to wet-antenna correction. IEEE Transactions on Geoscience and Remote Sensing, 60, 1-9.
Schleiss, M., Rieckermann, J., & Berne, A. (2013). Quantification and modeling of wet-antenna attenuation for commercial microwave links. IEEE Geoscience and Remote Sensing Letters, 10(5), 1195-1199.
van Leth, T. C., Overeem, A., Leijnse, H., & Uijlenhoet, R. (2018). A measurement campaign to assess sources of error in microwave link rainfall estimation. Atmospheric Measurement Techniques, 11(8), 4645-4669.
Widrow, B., Kollár, I., & Liu, M. C. (1996). Statistical theory of quantization. IEEE Transactions on Instrumentation and Measurement, 45(2), 353–361. DOI 10.1109/19.492748.
Zinevich, A., Messer, H., & Alpert, P. (2010). Prediction of rainfall intensity measurement errors using commercial microwave communication links. Atmospheric Measurement Techniques, 3(5), 1385-1402.
::end

