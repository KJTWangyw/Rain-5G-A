::title
Shared clear-sky reference errors in 5G NR pilot-based rain attenuation estimation: information bound and closed-form uncertainty propagation

::abstract
Demodulation reference signals (DM-RS) of 5G New Radio (NR) give millimetre-wave links a pilot-level measurement of rain attenuation, but the clear-sky reference against which the pilots are compared is estimated from a finite number of slots and reused for every subsequent observation. Its error is common to all attenuation estimates and does not decrease with temporal averaging. This study quantifies this shared error and propagates it in closed form. For the profile amplitude estimator, the Cramér–Rao bound of the averaged attenuation separates a shared floor from an averaging term, and the shared fraction depends only on the attenuation and on the number of averaged observations relative to the number of reference slots. Propagating the reference and current noise through the stationary-point equation of the estimator gives a complete closed-form variance that retains the noise-energy terms omitted by the conventional fixed-weight linearisation; the latter falls below the information bound by the factor $1+1/\gamma_{\mathrm{r}}$, where $\gamma_{\mathrm{r}}$ is the per-pilot reference signal-to-noise ratio. Monte Carlo simulations through an NR DM-RS OFDM chain at 28 GHz confirm the theory. With a two-slot reference at 0 dB, nominal 95% intervals from the fixed-weight linearisation cover {{cov_loc}}% of the cases, whereas the closed-form intervals cover {{cov_comp}}% and stay within {{grid_comp_range}} in all 16 signal-to-noise and reference conditions, matching a parametric bootstrap at a small fraction of its cost. The results define when the reused reference dominates the uncertainty of averaged and accumulated rain estimates and how long the reference must be.

::keywords
5G NR; rain attenuation; demodulation reference signal; clear-sky reference; Cramér–Rao bound; uncertainty propagation.

# Introduction

Rain attenuation along terrestrial microwave and millimetre-wave links provides path-averaged rainfall information that complements rain gauges and weather radar. Commercial microwave links have been used to map rainfall from received-power records at city and country scales, with processing chains for dry–wet classification, baseline estimation, wet-antenna correction and conversion to rain rate (Messer et al., 2006; Leijnse et al., 2007; Overeem et al., 2013, 2016; Chwala & Kunstmann, 2019). Millimetre-wave backhaul links in dense urban areas further increase the sensitivity to rain (Han et al., 2020).

The 5G New Radio (NR) air interface offers a finer measurement than received power. Its demodulation reference signals (DM-RS) are known symbols placed on fixed subcarriers of every scheduled slot (3GPP, 2024), so the receiver obtains a complex channel observation on each pilot. In integrated sensing and communication, such pilots are the natural probe of the propagation medium (Liu et al., 2022). The pilots do not, however, measure attenuation directly: the receiver must compare them with a clear-sky reference channel, which is estimated from a finite number of rain-free slots and then reused for every subsequent observation. The error of this reference is therefore common to all attenuation estimates. Averaging over time reduces the independent receiver noise but leaves the shared component unchanged.

Uncertainty studies of link-based rainfall have addressed power quantisation, wet-antenna attenuation and the measurement and interpolation errors of rainfall maps (Zinevich et al., 2010; Schleiss et al., 2013; Rios Gaona et al., 2015; Pastorek et al., 2022), but not the error of a finite pilot reference. The uncertainty of a pilot-based attenuation estimate is commonly propagated by first-order linearisation (Kay, 1993), which for the estimator considered here fixes the weights and phase at the estimate. Its accuracy for short references and low signal-to-noise ratio (SNR) has not been established, and Monte Carlo propagation (Efron, 1979; JCGM, 2008) resolves the question only at a high computational cost.

This study addresses a single question: how large is the uncertainty created by a reused clear-sky reference, and how can it be propagated correctly in closed form? The contributions are threefold. (1) The Cramér–Rao bound of the averaged attenuation with a shared reference is derived; it yields an irreducible floor and a closed-form criterion for when the shared error dominates. (2) A complete closed-form variance of the profile amplitude estimator is obtained, and the two mechanisms by which the fixed-weight linearisation falls below the bound are identified. (3) The theory is validated by Monte Carlo simulation through an NR DM-RS OFDM chain over 16 SNR and reference-length conditions, including interval coverage and computational cost against a parametric bootstrap.

The following sections are arranged as follows. Section 2 introduces the signal model and derives the bound and the variance propagation. Section 3 presents the simulation study. Discussions and conclusions are given in Sections 4 and 5, respectively.

# Signal model and uncertainty propagation

## 2.1 NR DM-RS pilot observation model

For a link of length $L$ (km), the rain-induced attenuation (dB) follows the power law of Recommendation ITU-R P.838-3 (ITU-R, 2005):

$$A_t = k L R_t^{\alpha}$$ {#eq:powerlaw}

where $R_t$ is the path-averaged rain rate (mm/h) at observation $t$ and $k$ and $\alpha$ depend on frequency and polarisation. The receiver observes NR downlink slots carrying PDSCH DM-RS of configuration type 1 on OFDM symbols 2 and 11, on every other subcarrier of the allocation (3GPP, 2024), which gives $J$ pilots per slot ($J = 144$ for 12 resource blocks). After cyclic-prefix removal, FFT and division by the known DM-RS symbols, the least-squares pilot estimates are modelled as

$$r_j = h_j + e_{0j},\qquad z_{tj} = a_t\,e^{\mathrm{i}\phi_t} h_j + e_{tj},\qquad a_t = 10^{-A_t/20}$$ {#eq:model}

where $h_j$ is the unknown clear-sky channel on pilot $j$, $r_j$ is the reference obtained by coherently averaging $M$ clear-sky slots, $z_{tj}$ is the pilot of rain observation $t$ with common phase $\phi_t$, and $e_{0j}$ and $e_{tj}$ are independent circular complex Gaussian errors with variances $v_0 = v/M$ and $v$. The noise variance $v$ is known from the receiver calibration. With $H = \sum_j |h_j|^2$, the per-pilot reference and current SNRs are $\gamma_{\mathrm{r}} = H/(Jv_0)$ and $\gamma_{\mathrm{t}} = a^2 H/(Jv)$; for unit mean pilot gain, $\gamma_{\mathrm{r}} = M\cdot\mathrm{SNR}$ and $\gamma_{\mathrm{t}} = a^2\,\mathrm{SNR}$. Figure @fig:schematic shows the pilot structure and the reuse of the reference.

![Pilot structure and reference reuse: (a) DM-RS type 1 on symbols 2 and 11 (two of 12 RB shown); (b) one reference from $M$ clear-sky slots serves $n$ rain observations.](fig1_schematic.png){#fig:schematic width=14.6}

## 2.2 Profile amplitude estimator

Treating the channel as a nuisance parameter and maximising the joint likelihood of Equation (@eq:model) over $h_j$ and $\phi_t$ gives the profile amplitude estimator

$$\hat{a}_t=\arg\min_{a>0}\left\{\sum_{j=1}^{J}\frac{|z_{tj}|^2+a^2|r_j|^2}{v+a^2v_0}-2a\left|\sum_{j=1}^{J}\frac{r_j^{*}z_{tj}}{v+a^2v_0}\right|\right\},\qquad \hat{A}_t=-20\log_{10}\hat{a}_t$$ {#eq:profile}

With $P_0 = \sum_j|r_j|^2$, $P_t = \sum_j|z_{tj}|^2$ and $C_t = |\sum_j r_j^{*}z_{tj}|$, the stationary point satisfies

$$C_t v_0\,a^2+(P_0 v-P_t v_0)\,a-C_t v=0$$ {#eq:root}

whose unique positive root is the estimate, so no numerical search is required. The quantity of interest is the mean of $n$ estimates, $\bar{A}_n = n^{-1}\sum_{t=1}^{n}\hat{A}_t$, all computed with the same reference.

## 2.3 Information bound with a shared reference

Consider $n$ observations with a common amplitude $a$. The Fisher information of Equation (@eq:model) with respect to $a$, the phases $\phi_t$ and the real and imaginary parts of $h_j$ has no cross term between $a$ and $\phi_t$. Taking the Schur complement with respect to the channel gives

$$I_n(a)=\sum_{j=1}^{J}\frac{2n|h_j|^2}{v+n a^2 v_0}$$ {#eq:fisher}

and, with $K = 20/\ln 10$ and $\mathrm{Var}(\hat{A}) \ge K^2/[a^2 I_n(a)]$, the bound for any unbiased estimate of the averaged attenuation is

$$\operatorname{Var}(\bar{A}_n)\ \ge\ B_n=\frac{K^2}{2J}\left(\frac{1}{\gamma_{\mathrm r}}+\frac{1}{n\gamma_{\mathrm t}}\right)$$ {#eq:crb}

The first term is a shared floor that no amount of averaging removes; the second decreases as $1/n$. The shared fraction of the bound is

$$\rho_n=\frac{n\gamma_{\mathrm t}}{\gamma_{\mathrm r}+n\gamma_{\mathrm t}}=\frac{n a^2}{M+n a^2}$$ {#eq:share}

where the second form holds for unit mean pilot gain and does not depend on the SNR. The shared reference therefore dominates once $n > M/a^2$: lengthening the reference and averaging more observations act on different components of the uncertainty.

## 2.4 Variance propagation

**Fixed-weight linearisation.** The conventional local propagation linearises Equation (@eq:profile) with the weights and phase fixed at the estimate. For white noise it gives a per-observation variance and a covariance shared by all observations,

$$\nu_t=\frac{K^2 v}{2\hat{a}_t^2 P_0},\qquad c_0=\frac{K^2 v_0}{2P_0},\qquad \operatorname{Var}_{\mathrm{loc}}(\bar{A}_n)=c_0+\frac{1}{n^2}\sum_{t=1}^{n}\nu_t$$ {#eq:local}

Because the plug-in reference energy has expectation $\mathrm{E}P_0 = H + Jv_0 = H(1+1/\gamma_{\mathrm{r}})$, both terms are smaller than the corresponding terms of Equation (@eq:crb):

$$\mathrm{E}\left[\operatorname{Var}_{\mathrm{loc}}(\bar{A}_n)\right]\approx\frac{B_n}{1+1/\gamma_{\mathrm r}}$$ {#eq:localbias}

The fixed-weight linearisation therefore predicts a variance below what any unbiased estimator can attain, and intervals built on it are too narrow whenever $\gamma_{\mathrm{r}}$ is small.

**Complete propagation.** The estimate is the root of Equation (@eq:root), a smooth function of the statistics $P_0$, $P_t$ and $C_t$. Write $X = \mathrm{Re}\sum_j h_j^{*}e_{0j}$, $Y_t = \mathrm{Re}\sum_j h_j^{*}e_{tj}$ (after phase alignment), $W_t = \mathrm{Re}\sum_j e_{0j}^{*}e_{tj}$, $Q_0 = \sum_j|e_{0j}|^2 - Jv_0$ and $Q_t = \sum_j|e_{tj}|^2 - Jv$; these terms are mutually uncorrelated with variances $Hv_0/2$, $Hv/2$, $Jv_0v/2$, $Jv_0^2$ and $Jv^2$. Around the noise-free point,

$$\delta P_0=2X+Q_0,\qquad \delta C_t\approx aX+Y_t+W_t,\qquad \delta P_t=2aY_t+Q_t$$ {#eq:perturb}

and the implicit-function theorem applied to Equation (@eq:root) gives, with $D = H(v+a^2v_0)$,

$$\delta a_t=-\frac{a}{H}X-\frac{va}{D}Q_0+\frac{1}{H}Y_t+\frac{v-a^2v_0}{D}W_t+\frac{v_0 a}{D}Q_t$$ {#eq:implicit}

The first two terms contain only the reference error and are shared by all observations; the remaining three are independent across observations. Since $\delta A_t = -(K/a)\,\delta a_t$, the variance of the averaged attenuation is

$$\operatorname{Var}(\bar{A}_n)\approx S+\frac{U}{n},\qquad S=\frac{K^2}{2J\gamma_{\mathrm r}}+\frac{K^2}{J(\gamma_{\mathrm r}+\gamma_{\mathrm t})^2}$$ {#eq:complete}

$$U=\frac{K^2}{2J\gamma_{\mathrm t}}+\frac{K^2}{2J\gamma_{\mathrm t}\gamma_{\mathrm r}}\left(\frac{\gamma_{\mathrm r}-\gamma_{\mathrm t}}{\gamma_{\mathrm r}+\gamma_{\mathrm t}}\right)^2+\frac{K^2}{J(\gamma_{\mathrm r}+\gamma_{\mathrm t})^2}$$ {#eq:independent}

The terms in $X$ and $Y_t$ reproduce the two terms of the bound in Equation (@eq:crb): to first order in the noise, the profile estimator is efficient. The remaining terms come from the noise energies $Q_0$, $Q_t$ and the cross term $W_t$; they are of order $1/(J\gamma^2)$, comparable to the first-order terms when $\gamma_{\mathrm{r}}$ and $\gamma_{\mathrm{t}}$ are close to one, and negligible at high SNR. The shared excess over the bound is $2\gamma_{\mathrm{r}}/(\gamma_{\mathrm{r}}+\gamma_{\mathrm{t}})^2$, or about $2/\gamma_{\mathrm{r}}$ when $\gamma_{\mathrm{r}}$ is much larger than $\gamma_{\mathrm{t}}$.

## 2.5 Interval construction

Equations (@eq:complete) and (@eq:independent) require only observable quantities. The channel energy is replaced by $\hat{H} = \max(P_0 - Jv_0,\,0.05P_0)$ and $a_t$ by $\hat{a}_t$; when the attenuation varies over the window, the independent terms are summed over the observations and the shared noise-energy coefficient is averaged. The interval at nominal level $1-\beta$ is

$$\bar{A}_n\pm z_{1-\beta/2}\sqrt{\hat{S}+\hat{U}/n}$$ {#eq:interval}

Three reference intervals are formed from the same estimates: the fixed-weight linearisation of Equation (@eq:local); the same linearisation with the shared covariance wrongly treated as independent across observations; and a parametric bootstrap (Efron, 1979) in which each of $N_{\mathrm{b}}$ replicates draws one reference error, reused by all observations,

$$r^{*}_j=\tilde{h}_j+e^{*}_{0j},\qquad z^{*}_{tj}=\hat{a}_t\tilde{h}_j+e^{*}_{tj},\qquad I_{\mathrm{basic}}=\left[2\bar{A}_n-q^{*}_{1-\beta/2},\ 2\bar{A}_n-q^{*}_{\beta/2}\right]$$ {#eq:bootstrap}

where $\tilde{h}_j = r_j(\hat{H}/P_0)^{1/2}$ and $q^{*}_u$ is the $u$-quantile of the bootstrap means.

# Simulation study

## 3.1 Configuration

Table @tab:config lists the simulation settings. Each slot is generated in the time domain: the DM-RS sequence of TS 38.211 is mapped to its subcarriers, transformed by the IFFT, extended by the cyclic prefix, passed through a three-path channel whose delays lie inside the cyclic prefix, scaled by the rain amplitude, rotated by a random common phase and corrupted by white Gaussian noise. The receiver removes the cyclic prefix, applies the FFT and divides by the known pilots. A fixed attenuation of {{A_db}} dB corresponds to {{rain_rate}} mm/h on the {{length_km}} km link. Four SNRs and four reference lengths form 16 conditions; in each, {{refs}} independent references are generated and each is reused for 32 rain observations, and intervals are evaluated for the mean of the first $n$ = 1, 4, 16 and 32 estimates at nominal levels of 50%, 80%, 90% and 95%. The main condition, 0 dB with $M = 2$ ($\gamma_{\mathrm{r}} = {{g_r}}$, $\gamma_{\mathrm{t}} = {{g_t}}$), represents a short reference at low SNR.

::table tab:config | Simulation settings. | 4.6,10.0
Parameter | Setting
Link | 28 GHz, horizontal polarisation, 0.5 km; $k$ = {{k_h}}, $\alpha$ = {{alpha_h}}
NR numerology | 120 kHz subcarrier spacing, 1024-point FFT, 72-sample cyclic prefix
Pilots | 12 RB, DM-RS type 1 on symbols 2 and 11, $J$ = 144 per slot
Channel | Three paths at 0, 73 and 187 ns; 0, −12.0 and −18.4 dB
Clear-sky SNR per pilot | 0, 5, 15 and 25 dB
Reference slots $M$ | 2, 8, 16 and 64
Rain observations | 32 per reference; $A$ = {{A_db}} dB; random common phase
Averaging and levels | $n$ = 1, 4, 16, 32; nominal 50%, 80%, 90%, 95%
Monte Carlo | {{refs}} references per condition; $N_{\mathrm{b}}$ = {{draws}}
::end

## 3.2 Variance of the averaged attenuation

Table @tab:components compares the variance components at the main condition: the bound of Equation (@eq:crb), the fixed-weight linearisation of Equation (@eq:localbias), the complete propagation of Equations (@eq:complete)–(@eq:independent) and the simulation, with 95% bootstrap intervals. The information bound has a shared term of {{crb_S}} dB², and at $n$ = 32 the shared reference accounts for {{share_main}}% of it (Equation @eq:share). The complete propagation predicts a shared term of {{comp_S}} dB², of which {{comp_Sq}} dB² comes from the reference noise energy, and the simulated shared covariance is {{sim_S}} dB² (95% interval {{sim_S_lo}}–{{sim_S_hi}}). The fixed-weight linearisation predicts {{loc_S}} dB², {{loc_under_pct}}% below the complete value and below the bound itself. For the mean of 32 estimates, the simulated variance is {{emp_V}} dB² ({{emp_V_lo}}–{{emp_V_hi}}), the complete propagation gives {{comp_V}} dB², and the fixed-weight linearisation gives {{loc_V}} dB².

::table tab:components | Variance components at the main condition (dB²; 0 dB, $M$ = 2, $J$ = 144). | 3.6,3.6,3.4,4.0
Quantity | Shared $S$ | Independent $U$ | $\mathrm{Var}(\bar{A}_{32})$
Information bound | {{crb_S}} | {{crb_U}} | {{crb_V}}
Fixed-weight | {{loc_S}} | {{loc_U}} | {{loc_V}}
Complete | {{comp_S}} | {{comp_U}} | {{comp_V}}
Simulation | {{sim_S}} ({{sim_S_lo}}–{{sim_S_hi}}) | {{sim_U}} | {{emp_V}} ({{emp_V_lo}}–{{emp_V_hi}})
::end

Figure @fig:variance shows the variance as a function of $n$. At $M = 2$ the simulated variance levels off at the shared term predicted by Equation (@eq:complete), while the fixed-weight prediction stays below the bound and the independence assumption keeps decreasing as $1/n$. At $M = 64$ the shared term is small and all predictions except the independence assumption agree. Over the 16 conditions, the predicted shared term lies inside the 95% interval of the simulated shared covariance in {{S_cells_inside}} conditions, and the predicted variance of $\bar{A}_n$ lies inside the 95% interval of the simulated variance in {{V_combos_inside}} of {{n_combos}} combinations of condition and $n$.

![Variance of the averaged attenuation at 0 dB versus the number of averaged observations: (a) $M$ = 2; (b) $M$ = 64. Points: simulation with 95% intervals.](fig2_variance.png){#fig:variance width=14.6}

## 3.3 Shared component and operating regime

Figure @fig:shared(a) places all conditions on a common axis. The shared variance of the profile estimator exceeds the bound by the factor $1 + 2\gamma_{\mathrm{r}}/(\gamma_{\mathrm{r}}+\gamma_{\mathrm{t}})^2$, and the fixed-weight linearisation falls below it by the factor $\gamma_{\mathrm{r}}/(\gamma_{\mathrm{r}}+1)$; at the main condition the two ratios are {{comp_over_crb}} and {{loc_over_crb}}. Both deviations vanish as $\gamma_{\mathrm{r}}$ grows, so the noise-energy terms matter only for short references or low SNR. Figure @fig:shared(b) maps the shared fraction of Equation (@eq:share): it reaches {{share_main}}% at the main condition with $n$ = 32 and drops to {{share_m64}}% for $M$ = 64.

![(a) Shared variance relative to the information bound; points: simulation with 95% intervals. (b) Shared fraction of the bound, Equation (@eq:share), at {{A_db}} dB.](fig3_shared.png){#fig:shared width=14.6}

## 3.4 Coverage, width and cost

Figure @fig:coverage shows the empirical coverage of nominal 95% intervals for $n$ = 32. The closed-form intervals of Equation (@eq:interval) cover {{grid_comp_range}} in all 16 conditions and {{gridall_comp_range}} over all values of $n$. The fixed-weight intervals reach only {{loc_m2_range}} with a two-slot reference, and ignoring the shared covariance leaves the coverage as low as {{grid_ind_min}}%. The parametric bootstrap covers {{grid_boot_range}}.

![Empirical coverage (%) of nominal 95% intervals for $n$ = 32 in the 16 conditions ({{refs}} references each).](fig4_coverage.png){#fig:coverage width=14.6}

Table @tab:main and Figure @fig:calibration detail the main condition. The closed-form interval covers {{cnt_comp}} of {{refs}} cases ({{cov_comp}}%), against {{cnt_loc}} ({{cov_loc}}%) for the fixed-weight linearisation; in paired comparison, {{mcn_comp_only}} cases are covered only by the closed form and {{mcn_loc_only}} only by the linearisation (exact McNemar test, {{mcn_p}}). The closed form and the bootstrap do not differ significantly ({{mcnb_p}}), and their mean widths are nearly equal ({{w_comp}} and {{w_boot}} dB). The closed-form coverage deviates from the nominal level by at most {{levels_dev_max}} percentage points between 50% and 95% (Figure @fig:calibration(a)). The width of the closed-form interval exceeds that of the fixed-weight interval by {{w_increase_pct}}%, which is the variance the linearisation omits. When the attenuation increases linearly from {{ramp_start}} to {{ramp_end}} dB over the 32 observations, the closed-form coverage stays within {{ramp_comp_range}} for all $n$, against {{ramp_loc_32}}% for the fixed-weight linearisation at $n$ = 32. Computing the closed-form interval takes {{ms_complete}} ms per window, less than the fixed-weight linearisation ({{ms_local}} ms), whereas the bootstrap takes {{ms_boot}} ms, about {{speedup}} times longer.

::table tab:main | Nominal 95% intervals at the main condition (0 dB, $M$ = 2, $n$ = 32, {{refs}} references). | 4.2,2.2,3.8,2.2,2.2
Method | Covered | Coverage, % (95% CI) | Width (dB) | Time (ms)
Independent (wrong) | {{cnt_ind}} | {{cov_ind}} ({{cov_ind_lo}}–{{cov_ind_hi}}) | {{w_ind}} | {{ms_local}}
Fixed-weight | {{cnt_loc}} | {{cov_loc}} ({{cov_loc_lo}}–{{cov_loc_hi}}) | {{w_loc}} | {{ms_local}}
Complete (closed form) | {{cnt_comp}} | {{cov_comp}} ({{cov_comp_lo}}–{{cov_comp_hi}}) | {{w_comp}} | {{ms_complete}}
Parametric bootstrap | {{cnt_boot}} | {{cov_boot}} ({{cov_boot_lo}}–{{cov_boot_hi}}) | {{w_boot}} | {{ms_boot}}
::end

![Main condition (0 dB, $M$ = 2, $n$ = 32): (a) empirical versus nominal coverage with Wilson 95% intervals; (b) mean width of the 95% intervals.](fig5_calibration.png){#fig:calibration width=14.6}

# Discussions

## 4.1 When the shared reference matters

Equation (@eq:share) reduces the role of the reference to the ratio $n a^2/M$. The shared term dominates once more than $M/a^2$ observations are averaged, about {{n_cross_main}} observations for a two-slot reference at {{A_db}} dB and {{n_cross_m64}} for 64 slots. This matters when slots are sampled sparsely in slow time: with one slot per minute, $M$ counts minutes of clear sky, and an hourly or event-accumulated estimate averages tens to hundreds of observations. For accumulated rainfall the shared error of every interval adds coherently, so the reference length should grow with the intended averaging window rather than be fixed. Mapped to rain rate through $\delta R/R \approx \delta A/(\alpha A)$, the shared floor at the main condition is {{floor_main_sd}} dB, a relative rain-rate error of {{floor_main_rel}}% at {{rain_rate}} mm/h that no averaging removes; with 64 reference slots it falls to {{floor_m64_rel}}%.

## 4.2 Why the fixed-weight linearisation fails

The fixed-weight linearisation under-covers for two separable reasons. First, the noisy reference energy used as plug-in places its prediction below the information bound by the factor $1+1/\gamma_{\mathrm{r}}$ (Equation @eq:localbias). Second, freezing the weights discards the noise-energy terms of Equation (@eq:implicit), so the actual variance of the profile estimator exceeds the bound by about $2/\gamma_{\mathrm{r}}$. At $\gamma_{\mathrm{r}} = {{g_r}}$ the two effects add up to an underestimate of the shared term by {{loc_under_pct}}%. The parametric bootstrap captures both effects because it re-evaluates the estimator, but the closed form reaches the same coverage at a cost below that of the linearisation. The bootstrap remains useful as a check for estimators without a closed-form root, for instance with frequency-dependent noise.

## 4.3 Limitations and future study

The derivation assumes white noise across pilots and a channel that stays fixed while the reference is reused. A slowly drifting channel or receiver gain adds a shared error that grows with the age of the reference; it can be included as an additional covariance term once its spectrum is known. The attenuation considered here is the total excess attenuation; separating rain from wet-antenna attenuation requires additional models (Schleiss et al., 2013; Pastorek et al., 2022). The validation uses simulated NR waveforms with a standard-compliant DM-RS structure; applying the method to DM-RS recorded on operational links is the next step.

# Conclusions

This study quantified the uncertainty that a reused clear-sky reference introduces into 5G NR pilot-based rain attenuation estimation. The Cramér–Rao bound separates an irreducible shared floor from an averaging term, and the shared fraction $n a^2/(M+n a^2)$ indicates directly whether the reference or the receiver noise limits an averaged estimate. The conventional fixed-weight linearisation falls below this bound by the factor $1+1/\gamma_{\mathrm{r}}$ and omits the noise-energy terms of the profile estimator; the complete closed-form propagation derived here retains them. In simulations through an NR DM-RS chain, the closed-form intervals covered {{grid_comp_range}} at the nominal 95% level in all conditions, including a two-slot reference at 0 dB where the linearisation covered {{cov_loc}}%, and matched a parametric bootstrap at a small fraction of its cost. The reference length should therefore be chosen together with the averaging window, and the closed-form variance can accompany every pilot-based attenuation estimate at negligible cost.

#* Data and code availability

The code, configuration files and results that reproduce all figures and tables are openly available at https://github.com/KJTWangyw/Rain-5G-A.

#* References

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
