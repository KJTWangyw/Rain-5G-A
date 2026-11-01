::title
5G NR导频降雨衰减估计中的共享晴空参考误差：<br>信息下界与闭式不确定性传播

::abstract
5G新空口（NR）的解调参考信号（DM-RS）为毫米波链路提供了导频级的降雨衰减观测，但导频需要与晴空参考比较，而该参考由有限个时隙估计，并在其后的每次观测中重复使用。参考误差为所有衰减估计所共有，不随时间平均而减小。本文定量刻画这一共享误差，并给出其闭式传播。对于剖面幅度估计器，平均衰减的Cramér–Rao下界可分解为不可平均消除的共享项与随平均减小的项，共享项所占比例只取决于衰减以及平均观测次数与参考时隙数之比。将参考噪声与当前噪声经估计器的驻点方程传播，得到完整的闭式方差，其中保留了常规固定权重线性化所忽略的噪声能量项；后者低于信息下界，相差因子$1+1/\gamma_{\mathrm{r}}$，其中$\gamma_{\mathrm{r}}$为单导频参考信噪比。基于28 GHz NR DM-RS OFDM链路的蒙特卡洛仿真验证了上述理论。在0 dB、2个参考时隙条件下，固定权重线性化给出的名义95%区间仅覆盖{{cov_loc}}%的情形，而闭式区间覆盖{{cov_comp}}%，并在全部16种信噪比与参考长度组合下保持在{{grid_comp_range}}，与参数化自助法结果一致，而计算成本仅为其很小一部分。研究结果明确了复用参考何时主导平均和累计降雨估计的不确定性，以及参考所需的长度。

::keywords
5G NR；降雨衰减；解调参考信号；晴空参考；Cramér–Rao下界；不确定性传播。

# 引言

地面微波与毫米波链路上的雨致衰减提供路径平均的降雨信息，可与雨量计和天气雷达互为补充。已有研究利用商用微波链路的接收功率记录，在城市和国家尺度上绘制降雨分布，并建立了干湿判别、基线估计、湿天线修正和雨强转换等处理流程（Messer et al., 2006；Leijnse et al., 2007；Overeem et al., 2013, 2016；Chwala & Kunstmann, 2019）。密集城区的毫米波回传链路对降雨更加敏感（Han et al., 2020）。

5G新空口（NR）提供了比接收功率更精细的观测。其解调参考信号（DM-RS）是放置在每个调度时隙固定子载波上的已知符号（3GPP, 2024），接收机因此能在每个导频上得到复信道观测。在通感一体化中，这类导频是探测传播介质的天然手段（Liu et al., 2022）。但导频并不直接测量衰减：接收机必须把它们与晴空参考信道比较，而该参考由有限个无雨时隙估计，随后在每次观测中重复使用。因此，参考误差为所有衰减估计所共有。时间平均可以减小独立的接收噪声，却不会改变共享分量。

已有链路降雨的不确定性研究讨论了功率量化、湿天线衰减以及降雨图的测量与插值误差（Zinevich et al., 2010；Schleiss et al., 2013；Rios Gaona et al., 2015；Pastorek et al., 2022），但尚未涉及有限导频参考的误差。导频衰减估计的不确定性通常采用一阶线性化传播（Kay, 1993），对本文的估计器而言，即在估计值处固定权重和相位。该方法在短参考和低信噪比（SNR）条件下的准确性尚未得到检验，而蒙特卡洛传播（Efron, 1979；JCGM, 2008）虽能回答这一问题，计算代价却很高。

本文只研究一个问题：复用晴空参考带来的不确定性有多大，以及如何用闭式公式正确地传播。主要贡献有三点：（1）推导了共享参考下平均衰减的Cramér–Rao下界，给出不可平均消除的下限，以及判断共享误差何时占主导的闭式判据；（2）得到剖面幅度估计器的完整闭式方差，并指出固定权重线性化低于下界的两个机制；（3）通过NR DM-RS OFDM链路的蒙特卡洛仿真，在16种信噪比与参考长度组合下验证理论，内容包括区间覆盖率以及与参数化自助法的计算成本比较。

本文其余部分安排如下：第2节建立信号模型，推导信息下界与方差传播；第3节给出仿真研究；第4节和第5节分别为讨论与结论。

# 信号模型与不确定性传播

## 2.1 NR DM-RS导频观测模型

对于长度为$L$（km）的链路，雨致衰减（dB）服从ITU-R P.838-3建议书的幂律关系（ITU-R, 2005）：

$$A_t = k L R_t^{\alpha}$$ {#eq:powerlaw}

^ 式中，$R_t$为第$t$次观测时的路径平均雨强（mm/h），$k$和$\alpha$取决于频率与极化方式。接收机观测承载PDSCH DM-RS的NR下行时隙，DM-RS采用配置类型1，位于OFDM符号2和11、分配带宽内每隔一个子载波处（3GPP, 2024），每个时隙共有$J$个导频（12个资源块时$J = 144$）。去除循环前缀、作FFT并除以已知DM-RS符号后，最小二乘导频估计建模为

$$r_j = h_j + e_{0j},\qquad z_{tj} = a_t\,e^{\mathrm{i}\phi_t} h_j + e_{tj},\qquad a_t = 10^{-A_t/20}$$ {#eq:model}

^ 式中，$h_j$为第$j$个导频上未知的晴空信道；$r_j$为对$M$个晴空时隙相干平均得到的参考；$z_{tj}$为第$t$次降雨观测的导频，其公共相位为$\phi_t$；$e_{0j}$和$e_{tj}$为相互独立的圆对称复高斯误差，方差分别为$v_0 = v/M$和$v$。噪声方差$v$由接收机标定得到。记$H = \sum_j |h_j|^2$，单导频参考信噪比和当前信噪比分别为$\gamma_{\mathrm{r}} = H/(Jv_0)$和$\gamma_{\mathrm{t}} = a^2 H/(Jv)$；导频平均增益为1时，$\gamma_{\mathrm{r}} = M\cdot\mathrm{SNR}$，$\gamma_{\mathrm{t}} = a^2\,\mathrm{SNR}$。图@fig:schematic给出了导频结构和参考的复用方式。

![导频结构与晴空参考的复用](fig1_schematic.png){#fig:schematic width=13.6}

## 2.2 剖面幅度估计器

将信道视为冗余参数，对式(@eq:model)的联合似然关于$h_j$和$\phi_t$取最大，得到剖面幅度估计器

$$\hat{a}_t=\arg\min_{a>0}\left\{\sum_{j=1}^{J}\frac{|z_{tj}|^2+a^2|r_j|^2}{v+a^2v_0}-2a\left|\sum_{j=1}^{J}\frac{r_j^{*}z_{tj}}{v+a^2v_0}\right|\right\},\qquad \hat{A}_t=-20\log_{10}\hat{a}_t$$ {#eq:profile}

记$P_0 = \sum_j|r_j|^2$、$P_t = \sum_j|z_{tj}|^2$、$C_t = |\sum_j r_j^{*}z_{tj}|$，驻点满足

$$C_t v_0\,a^2+(P_0 v-P_t v_0)\,a-C_t v=0$$ {#eq:root}

^ 其唯一正根即为估计值，无需数值搜索。本文关注$n$个估计值的均值$\bar{A}_n = n^{-1}\sum_{t=1}^{n}\hat{A}_t$，这些估计值均使用同一参考计算。

## 2.3 共享参考下的信息下界

考虑公共幅度为$a$的$n$次观测。式(@eq:model)关于$a$、相位$\phi_t$以及$h_j$实部和虚部的Fisher信息中，$a$与$\phi_t$之间不存在交叉项。对信道取Schur补，得

$$I_n(a)=\sum_{j=1}^{J}\frac{2n|h_j|^2}{v+n a^2 v_0}$$ {#eq:fisher}

^ 再由$K = 20/\ln 10$及$\mathrm{Var}(\hat{A}) \ge K^2/[a^2 I_n(a)]$，平均衰减的任意无偏估计满足

$$\operatorname{Var}(\bar{A}_n)\ \ge\ B_n=\frac{K^2}{2J}\left(\frac{1}{\gamma_{\mathrm r}}+\frac{1}{n\gamma_{\mathrm t}}\right)$$ {#eq:crb}

第一项是共享下限，任何平均都无法消除；第二项按$1/n$减小。共享项在下界中所占比例为

$$\rho_n=\frac{n\gamma_{\mathrm t}}{\gamma_{\mathrm r}+n\gamma_{\mathrm t}}=\frac{n a^2}{M+n a^2}$$ {#eq:share}

^ 后一等式在导频平均增益为1时成立，且与信噪比无关。因此，当$n > M/a^2$时共享参考占主导：延长参考与增加平均次数作用于不确定性的不同分量。

## 2.4 方差传播

**固定权重线性化。**常规的局部传播在估计值处固定权重和相位，对式(@eq:profile)作线性化。白噪声下，它给出逐次观测的方差和所有观测共享的协方差：

$$\nu_t=\frac{K^2 v}{2\hat{a}_t^2 P_0},\qquad c_0=\frac{K^2 v_0}{2P_0},\qquad \operatorname{Var}_{\mathrm{loc}}(\bar{A}_n)=c_0+\frac{1}{n^2}\sum_{t=1}^{n}\nu_t$$ {#eq:local}

由于作为插件量的参考能量期望为$\mathrm{E}P_0 = H + Jv_0 = H(1+1/\gamma_{\mathrm{r}})$，上述两项均小于式(@eq:crb)中的对应项：

$$\mathrm{E}\left[\operatorname{Var}_{\mathrm{loc}}(\bar{A}_n)\right]\approx\frac{B_n}{1+1/\gamma_{\mathrm r}}$$ {#eq:localbias}

因此，固定权重线性化预测的方差低于任何无偏估计器所能达到的下限；当$\gamma_{\mathrm{r}}$较小时，据此构造的区间过窄。

**完整传播。**估计值是式(@eq:root)的根，即统计量$P_0$、$P_t$和$C_t$的光滑函数，可在无噪声点附近对其作一阶展开。为此引入以下随机量（$Y_t$和$W_t$在相位对齐后定义）：

$$X=\mathrm{Re}\sum_j h_j^{*}e_{0j},\qquad Y_t=\mathrm{Re}\sum_j h_j^{*}e_{tj},\qquad W_t=\mathrm{Re}\sum_j e_{0j}^{*}e_{tj}$$

$$Q_0=\sum_j|e_{0j}|^2-Jv_0,\qquad Q_t=\sum_j|e_{tj}|^2-Jv$$

^ 它们互不相关，方差分别为$Hv_0/2$、$Hv/2$、$Jv_0v/2$、$Jv_0^2$和$Jv^2$。展开后有

$$\delta P_0=2X+Q_0,\qquad \delta C_t\approx aX+Y_t+W_t,\qquad \delta P_t=2aY_t+Q_t$$ {#eq:perturb}

^ 对式(@eq:root)应用隐函数定理，并记$D = H(v+a^2v_0)$，得

$$\delta a_t=-\frac{a}{H}X-\frac{va}{D}Q_0+\frac{1}{H}Y_t+\frac{v-a^2v_0}{D}W_t+\frac{v_0 a}{D}Q_t$$ {#eq:implicit}

前两项只含参考误差，为所有观测所共享；其余三项在各次观测之间相互独立。由$\delta A_t = -(K/a)\,\delta a_t$，平均衰减的方差为

$$\operatorname{Var}(\bar{A}_n)\approx S+\frac{U}{n},\qquad S=\frac{K^2}{2J\gamma_{\mathrm r}}+\frac{K^2}{J(\gamma_{\mathrm r}+\gamma_{\mathrm t})^2}$$ {#eq:complete}

$$U=\frac{K^2}{2J\gamma_{\mathrm t}}+\frac{K^2}{2J\gamma_{\mathrm t}\gamma_{\mathrm r}}\left(\frac{\gamma_{\mathrm r}-\gamma_{\mathrm t}}{\gamma_{\mathrm r}+\gamma_{\mathrm t}}\right)^2+\frac{K^2}{J(\gamma_{\mathrm r}+\gamma_{\mathrm t})^2}$$ {#eq:independent}

$X$和$Y_t$对应的项恰好重现式(@eq:crb)下界的两项：在噪声的一阶近似下，剖面估计器是有效的。其余各项来自噪声能量$Q_0$、$Q_t$和交叉项$W_t$，量级为$1/(J\gamma^2)$；当$\gamma_{\mathrm{r}}$和$\gamma_{\mathrm{t}}$接近1时与一阶项相当，高信噪比下可以忽略。共享项超出下界的相对量为$2\gamma_{\mathrm{r}}/(\gamma_{\mathrm{r}}+\gamma_{\mathrm{t}})^2$，当$\gamma_{\mathrm{r}}$远大于$\gamma_{\mathrm{t}}$时约为$2/\gamma_{\mathrm{r}}$。

## 2.5 区间构造

式(@eq:complete)和式(@eq:independent)只需可观测量：将信道能量替换为$\hat{H} = \max(P_0 - Jv_0,\,0.05P_0)$，将$a_t$替换为$\hat{a}_t$；当衰减在窗口内变化时，独立项对各次观测求和，共享噪声能量项的系数取平均。名义水平为$1-\beta$的区间为

$$\bar{A}_n\pm z_{1-\beta/2}\sqrt{\hat{S}+\hat{U}/n}$$ {#eq:interval}

由同一组估计值构造三种对照区间：式(@eq:local)的固定权重线性化；将共享协方差错误地视为各次观测间独立的同一线性化；以及参数化自助法（Efron, 1979），其$N_{\mathrm{b}}$次重抽样中的每一次只抽取一份参考误差，并由全部观测复用：

$$r^{*}_j=\tilde{h}_j+e^{*}_{0j},\qquad z^{*}_{tj}=\hat{a}_t\tilde{h}_j+e^{*}_{tj},\qquad I_{\mathrm{basic}}=\left[2\bar{A}_n-q^{*}_{1-\beta/2},\ 2\bar{A}_n-q^{*}_{\beta/2}\right]$$ {#eq:bootstrap}

^ 式中，$\tilde{h}_j = r_j(\hat{H}/P_0)^{1/2}$，$q^{*}_u$为自助均值的$u$分位数。

# 仿真研究

## 3.1 仿真配置

表@tab:config列出了仿真设置。每个时隙在时域生成：将TS 38.211的DM-RS序列映射到相应子载波，经IFFT并加循环前缀，通过时延均位于循环前缀内的三径信道，再按降雨幅度缩放、乘以随机公共相位并叠加高斯白噪声。接收机去除循环前缀、作FFT并除以已知导频。{{A_db}} dB的固定衰减对应{{length_km}} km链路上{{rain_rate}} mm/h的雨强。四种信噪比与四种参考长度组成16种条件。每种条件生成{{refs}}个独立参考，每个参考复用于32次降雨观测；对前$n$个估计值的均值（$n$取1、4、16和32）评估区间，名义水平取50%、80%、90%和95%。主条件为0 dB、$M = 2$（$\gamma_{\mathrm{r}} = {{g_r}}$，$\gamma_{\mathrm{t}} = {{g_t}}$），代表低信噪比下的短参考。

::table tab:config | 仿真设置 | 4.6,10.0
参数 | 设置
链路 | 28 GHz，水平极化，0.5 km；$k$ = {{k_h}}，$\alpha$ = {{alpha_h}}
NR参数集 | 子载波间隔120 kHz，1024点FFT，循环前缀72个采样点
导频 | 12个RB，类型1 DM-RS位于符号2和11，每时隙$J$ = 144
信道 | 三径，时延0、73、187 ns，相对功率0、−12.0、−18.4 dB
单导频晴空信噪比 | 0、5、15、25 dB
参考时隙数$M$ | 2、8、16、64
降雨观测 | 每个参考32次；$A$ = {{A_db}} dB；随机公共相位
平均次数与名义水平 | $n$ = 1、4、16、32；名义水平50%、80%、90%、95%
蒙特卡洛 | 每种条件{{refs}}个参考；$N_{\mathrm{b}}$ = {{draws}}
::end

## 3.2 平均衰减的方差

表@tab:components比较了主条件下的方差分量：式(@eq:crb)的下界、式(@eq:localbias)的固定权重线性化、式(@eq:complete)–(@eq:independent)的完整传播以及仿真结果（含95%自助置信区间）。信息下界的共享项为{{crb_S}} dB²，$n$ = 32时共享参考占其{{share_main}}%（式(@eq:share)）。完整传播预测的共享项为{{comp_S}} dB²，其中{{comp_Sq}} dB²来自参考噪声能量；仿真得到的共享协方差为{{sim_S}} dB²（95%区间{{sim_S_lo}}–{{sim_S_hi}}）。固定权重线性化的预测为{{loc_S}} dB²，比完整传播低{{loc_under_pct}}%，甚至低于下界本身。对于32个估计值的均值，仿真方差为{{emp_V}} dB²（{{emp_V_lo}}–{{emp_V_hi}}），完整传播给出{{comp_V}} dB²，固定权重线性化给出{{loc_V}} dB²。

::table tab:components | 主条件下的方差分量（dB²；0 dB，$M$ = 2，$J$ = 144） | 3.6,3.6,3.4,4.0
量 | 共享项$S$ | 独立项$U$ | $\mathrm{Var}(\bar{A}_{32})$
信息下界 | {{crb_S}} | {{crb_U}} | {{crb_V}}
固定权重线性化 | {{loc_S}} | {{loc_U}} | {{loc_V}}
完整传播 | {{comp_S}} | {{comp_U}} | {{comp_V}}
仿真 | {{sim_S}} ({{sim_S_lo}}–{{sim_S_hi}}) | {{sim_U}} | {{emp_V}} ({{emp_V_lo}}–{{emp_V_hi}})
::end

图@fig:variance给出方差随$n$的变化。$M = 2$时，仿真方差趋于式(@eq:complete)预测的共享项，固定权重的预测始终低于下界，而独立假设的预测按$1/n$持续下降。$M = 64$时共享项很小，除独立假设外各预测一致。在16种条件中，预测的共享项落在仿真共享协方差95%区间内的有{{S_cells_inside}}种；在{{n_combos}}个条件与$n$的组合中，$\bar{A}_n$的预测方差落在仿真方差95%区间内的有{{V_combos_inside}}个。

![0 dB时平均衰减的方差随平均观测次数的变化：(a) $M$ = 2；(b) $M$ = 64。点为仿真值及其95%区间。](fig2_variance.png){#fig:variance width=14.6}

## 3.3 共享分量与适用范围

图@fig:shared(a)将所有条件置于同一坐标下。剖面估计器的共享方差高于下界，倍数为$1 + 2\gamma_{\mathrm{r}}/(\gamma_{\mathrm{r}}+\gamma_{\mathrm{t}})^2$；固定权重线性化则低于下界，倍数为$\gamma_{\mathrm{r}}/(\gamma_{\mathrm{r}}+1)$。主条件下两者分别为{{comp_over_crb}}和{{loc_over_crb}}。两种偏差都随$\gamma_{\mathrm{r}}$增大而消失，因此噪声能量项只在短参考或低信噪比时起作用。图@fig:shared(b)给出式(@eq:share)的共享占比：主条件下$n$ = 32时达到{{share_main}}%，$M$ = 64时降至{{share_m64}}%。

![(a) 共享方差与信息下界之比，点为仿真值及其95%区间；(b) {{A_db}} dB衰减下共享项在下界中的占比，即式(@eq:share)。](fig3_shared.png){#fig:shared width=14.6}

## 3.4 覆盖率、宽度与计算成本

图@fig:coverage给出$n$ = 32时名义95%区间的实际覆盖率。式(@eq:interval)的闭式区间在全部16种条件下的覆盖率为{{grid_comp_range}}，在所有$n$下为{{gridall_comp_range}}。参考为2个时隙时，固定权重区间的覆盖率仅为{{loc_m2_range}}；忽略共享协方差时，覆盖率低至{{grid_ind_min}}%。参数化自助法的覆盖率为{{grid_boot_range}}。

![16种条件下$n$ = 32时名义95%区间的实际覆盖率（%），每种条件{{refs}}个参考。](fig4_coverage.png){#fig:coverage width=14.6}

表@tab:main和图@fig:calibration详细给出了主条件的结果。在{{refs}}个情形中，闭式区间覆盖{{cnt_comp}}个（{{cov_comp}}%），固定权重线性化覆盖{{cnt_loc}}个（{{cov_loc}}%）；配对比较中，仅被闭式区间覆盖的有{{mcn_comp_only}}个，仅被线性化区间覆盖的有{{mcn_loc_only}}个（精确McNemar检验，{{mcn_p}}）。闭式区间与自助法区间无显著差异（{{mcnb_p}}），平均宽度几乎相同（分别为{{w_comp}} dB和{{w_boot}} dB）。在50%至95%的名义水平范围内，闭式区间覆盖率与名义水平的偏差不超过{{levels_dev_max}}个百分点（图@fig:calibration(a)）。闭式区间比固定权重区间宽{{w_increase_pct}}%，多出的部分正是线性化所遗漏的方差。当衰减在32次观测中由{{ramp_start}} dB线性增至{{ramp_end}} dB时，闭式区间在所有$n$下的覆盖率保持在{{ramp_comp_range}}，而固定权重线性化在$n$ = 32时仅为{{ramp_loc_32}}%。闭式区间每个窗口的计算耗时为{{ms_complete}} ms，低于固定权重线性化（{{ms_local}} ms）；自助法耗时{{ms_boot}} ms，约为闭式方法的{{speedup}}倍。

::table tab:main | 主条件下的名义95%区间（0 dB，$M$ = 2，$n$ = 32，{{refs}}个参考） | 4.0,2.0,4.4,2.1,2.1
方法 | 覆盖次数 | 覆盖率/%（95%区间） | 宽度/dB | 耗时/ms
独立假设（错误） | {{cnt_ind}} | {{cov_ind}} ({{cov_ind_lo}}–{{cov_ind_hi}}) | {{w_ind}} | {{ms_local}}
固定权重线性化 | {{cnt_loc}} | {{cov_loc}} ({{cov_loc_lo}}–{{cov_loc_hi}}) | {{w_loc}} | {{ms_local}}
完整闭式传播 | {{cnt_comp}} | {{cov_comp}} ({{cov_comp_lo}}–{{cov_comp_hi}}) | {{w_comp}} | {{ms_complete}}
参数化自助法 | {{cnt_boot}} | {{cov_boot}} ({{cov_boot_lo}}–{{cov_boot_hi}}) | {{w_boot}} | {{ms_boot}}
::end

![主条件（0 dB，$M$ = 2，$n$ = 32）：(a) 实际覆盖率与名义水平的关系及Wilson 95%区间；(b) 95%区间的平均宽度。](fig5_calibration.png){#fig:calibration width=14.6}

# 讨论

## 4.1 共享参考何时起主导作用

式(@eq:share)将参考的作用归结为比值$n a^2/M$。平均观测次数超过$M/a^2$时，共享项即占主导：对于{{A_db}} dB衰减，2个参考时隙时约为{{n_cross_main}}次，64个时隙时约为{{n_cross_m64}}次。这在慢时间尺度上稀疏采样时尤为重要：若每分钟取一个时隙，$M$即为晴空参考的分钟数，而小时或降雨事件的累计估计需要平均数十至数百次观测。对于累计雨量，各时段的共享误差相干叠加，因此参考长度应随所需的平均窗口增长，而不应固定不变。按$\delta R/R \approx \delta A/(\alpha A)$换算到雨强，主条件下的共享下限为{{floor_main_sd}} dB，对应{{rain_rate}} mm/h雨强下{{floor_main_rel}}%的相对误差，且无法通过平均消除；参考时隙增至64个时，该误差降至{{floor_m64_rel}}%。

## 4.2 固定权重线性化失效的原因

固定权重线性化覆盖不足有两个可分离的原因。其一，作为插件量的含噪参考能量使其预测低于信息下界，相差因子$1+1/\gamma_{\mathrm{r}}$（式(@eq:localbias)）。其二，固定权重舍弃了式(@eq:implicit)中的噪声能量项，而剖面估计器的实际方差比下界高约$2/\gamma_{\mathrm{r}}$。$\gamma_{\mathrm{r}} = {{g_r}}$时，两者叠加使共享项被低估{{loc_under_pct}}%。参数化自助法通过重新计算估计器同时捕获了这两种效应，而闭式方法以低于线性化的计算成本达到了相同的覆盖率。对于没有闭式根的估计器，例如频率相关噪声的情形，自助法仍可作为检验手段。

## 4.3 局限与展望

推导假设导频间为白噪声，且在参考复用期间信道保持不变。信道或接收增益的缓慢漂移会引入随参考使用时长增长的共享误差，在已知其谱特性后，可作为附加的协方差项纳入。本文考虑的是总附加衰减，从中分离雨衰与湿天线衰减需要额外的模型（Schleiss et al., 2013；Pastorek et al., 2022）。验证基于符合标准DM-RS结构的仿真NR波形，下一步将把该方法应用于实际运营链路记录的DM-RS。

# 结论

本文定量刻画了复用晴空参考给5G NR导频降雨衰减估计带来的不确定性。Cramér–Rao下界将不可消除的共享下限与随平均减小的项分开，共享占比$n a^2/(M+n a^2)$可直接判断平均估计受参考误差还是接收噪声限制。常规固定权重线性化低于该下界，相差因子$1+1/\gamma_{\mathrm{r}}$，并遗漏了剖面估计器的噪声能量项；本文推导的完整闭式传播保留了这些项。在NR DM-RS链路仿真中，闭式区间在全部条件下对名义95%水平的覆盖率为{{grid_comp_range}}，包括固定权重线性化仅覆盖{{cov_loc}}%的0 dB、2时隙参考条件，并以很小的计算代价得到与参数化自助法一致的结果。因此，参考长度应与平均窗口一同选择，闭式方差则可以几乎零成本地随每个导频衰减估计一并给出。

#* 数据与代码可用性

复现全部图表所需的代码、配置文件和结果已公开于https://github.com/KJTWangyw/Rain-5G-A。

#* 参考文献

::refs
3GPP. (2024). NR; Physical channels and modulation (3GPP TS 38.211 version 18.3.0 Release 18). ETSI TS 138 211 V18.3.0.
Chwala, C., & Kunstmann, H. (2019). Commercial microwave link networks for rainfall observation: Assessment of the current status and future challenges. WIREs Water, 6(2), e1337.
Efron, B. (1979). Bootstrap methods: Another look at the jackknife. The Annals of Statistics, 7(1), 1-26.
Han, C., Huo, J., Gao, Q., Su, G., & Wang, H. (2020). Rainfall monitoring based on next-generation millimeter-wave backhaul technologies in a dense urban environment. Remote Sensing, 12(6), 1045.
ITU-R. (2005). Specific attenuation model for rain for use in prediction methods (Recommendation ITU-R P.838-3). International Telecommunication Union.
JCGM. (2008). Evaluation of measurement data—Supplement 1 to the "Guide to the expression of uncertainty in measurement"—Propagation of distributions using a Monte Carlo method (JCGM 101:2008). Joint Committee for Guides in Metrology.
Kay, S. M. (1993). Fundamentals of statistical signal processing: Estimation theory. Prentice Hall.
Leijnse, H., Uijlenhoet, R., & Stricker, J. N. M. (2007). Rainfall measurement using radio links from cellular communication networks. Water Resources Research, 43(3), W03201.
Liu, F., Cui, Y., Masouros, C., Xu, J., Han, T. X., Eldar, Y. C., & Buzzi, S. (2022). Integrated sensing and communications: Toward dual-functional wireless networks for 6G and beyond. IEEE Journal on Selected Areas in Communications, 40(6), 1728-1767.
Messer, H., Zinevich, A., & Alpert, P. (2006). Environmental monitoring by wireless communication networks. Science, 312(5774), 713.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2013). Country-wide rainfall maps from cellular communication networks. Proceedings of the National Academy of Sciences, 110(8), 2741-2745.
Overeem, A., Leijnse, H., & Uijlenhoet, R. (2016). Retrieval algorithm for rainfall mapping from microwave links in a cellular communication network. Atmospheric Measurement Techniques, 9(5), 2425-2444.
Pastorek, J., Fencl, M., Rieckermann, J., & Bareš, V. (2022). Precipitation estimates from commercial microwave links: Practical approaches to wet-antenna correction. IEEE Transactions on Geoscience and Remote Sensing, 60, 4104409.
Rios Gaona, M. F., Overeem, A., Leijnse, H., & Uijlenhoet, R. (2015). Measurement and interpolation uncertainties in rainfall maps from cellular communication networks. Hydrology and Earth System Sciences, 19(8), 3571-3584.
Schleiss, M., Rieckermann, J., & Berne, A. (2013). Quantification and modeling of wet-antenna attenuation for commercial microwave links. IEEE Geoscience and Remote Sensing Letters, 10(5), 1195-1199.
Zinevich, A., Messer, H., & Alpert, P. (2010). Prediction of rainfall intensity measurement errors using commercial microwave communication links. Atmospheric Measurement Techniques, 3(5), 1385-1402.
::end
