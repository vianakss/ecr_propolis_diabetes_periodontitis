********************************************************************************
*** LONGITUDINAL ANCOVA MIXED MODEL - CYTOKINES ***
********************************************************************************
capture log close
log using "log_table_mixed_model_cytokines.log", text replace
xtset id time

use "dataset_propolis_analysis.dta", clear

capture drop lab_bac miss90_bac
bysort id: egen lab_bac = max(time==0 & !missing(log_un))
bysort id: egen miss90_bac = max(time==2 & missing(log_un))
replace miss90_bac = 0 if !lab_bac
count if miss90_bac & time==2 & !missing(il6)

*** 0. DATA PREPARATION AND MISSING DATA PATTERN ***
* 0.1. Baseline values as covariates
capture drop il6_base
bysort id: egen il6_base = max(cond(time==0, il6, .))
capture drop log_tnf
gen log_tnf = log10(tnf)
capture drop log_tnf_base
bysort id: egen log_tnf_base = max(cond(time==0, log_tnf, .))

* 0.2. Missing data pattern
misstable summarize il6 log_tnf
table group time, statistic(count il6) statistic(count log_tnf)

* 0.3. Baseline availability (participants retained by the analysis of covariance)
count if time==0 & !missing(il6)
count if time==0 & !missing(log_tnf)


*** 1. CYTOKINES ***
* 1.1. INTERLEUKIN-6
* 1.1.1. Observed descriptive values
table group time, statistic(count il6) statistic(mean il6) statistic(sd il6)
table group time, statistic(median il6) statistic(p25 il6) statistic(p75 il6)

* 1.1.2. Longitudinal analysis of covariance
mixed il6 i.group##i.time c.il6_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.1.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.1.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.1.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.1.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_il6, replace)

* 1.1.7. Model diagnostics
capture drop fit_il6 res_il6
predict fit_il6, fitted
predict res_il6, rstandard
scatter res_il6 fit_il6, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - IL-6", size(medium)) name(resfit_il6, replace)
qnorm res_il6 if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids IL-6 - 60 days", size(medium)) name(qnorm_il6_t1, replace)
qnorm res_il6 if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids IL-6 - 90 days", size(medium)) name(qnorm_il6_t2, replace)
graph combine qnorm_il6_t1 qnorm_il6_t2, name(il6_resid_panel, replace)
graph export "il6_resid_panel.png", name(il6_resid_panel) width(3000) replace



* 1.2. TUMOUR NECROSIS FACTOR-ALPHA (log10 scale)
* 1.2.1. Observed descriptive values
table group time, statistic(count tnf) statistic(mean tnf) statistic(sd tnf)
table group time, statistic(median tnf) statistic(p25 tnf) statistic(p75 tnf)

* 1.2.2. Longitudinal analysis of covariance
mixed log_tnf i.group##i.time c.log_tnf_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.2.3. Estimated treatment effect at 60 days: Test / Control
lincom 1.group, small
display "Ratio of geometric means at 60 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.2.4. Estimated treatment effect at 90 days: Test / Control
lincom 1.group + 1.group#2.time, small
display "Ratio of geometric means at 90 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.2.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.2.6. Model-estimated means and profile
margins group#time
matrix M = r(table)
display "Geometric mean, Control 60 days: " %5.2f 10^M[1,1] " [" %5.2f 10^M[5,1] "; " %5.2f 10^M[6,1] "]"
display "Geometric mean, Control 90 days: " %5.2f 10^M[1,2] " [" %5.2f 10^M[5,2] "; " %5.2f 10^M[6,2] "]"
display "Geometric mean, Test 60 days: " %5.2f 10^M[1,3] " [" %5.2f 10^M[5,3] "; " %5.2f 10^M[6,3] "]"
display "Geometric mean, Test 90 days: " %5.2f 10^M[1,4] " [" %5.2f 10^M[5,4] "; " %5.2f 10^M[6,4] "]"
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_log_tnf, replace)

* 1.2.7. Model diagnostics
capture drop fit_log_tnf res_log_tnf
predict fit_log_tnf, fitted
predict res_log_tnf, rstandard
scatter res_log_tnf fit_log_tnf, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - TNF-alpha (log10)", size(medium)) name(resfit_log_tnf, replace)
qnorm res_log_tnf if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids TNF-alpha (log10) - 60 days", size(medium)) name(qnorm_log_tnf_t1, replace)
qnorm res_log_tnf if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids TNF-alpha (log10) - 90 days", size(medium)) name(qnorm_log_tnf_t2, replace)
graph combine qnorm_log_tnf_t1 qnorm_log_tnf_t2, name(log_tnf_resid_panel, replace)
graph export "log_tnf_resid_panel.png", name(log_tnf_resid_panel) width(3000) replace


graph combine resfit_il6 resfit_log_tnf, rows(1) name(scatter_cytokines_panel, replace)
graph export "scatter_cytokines_panel.png", name(scatter_cytokines_panel) width(3000) replace

log close


* 1.3. Sensitivity analysis
capture log close
log using "log_leave_one_out_salivary_parameters.log", replace
* 1.3.1. Leave-one-out analysis
* 1.3.1.1. Participant identifiers with cytokine data and labels (0 is reserved for the full sample)
assert id != 0
levelsof id if !missing(il6), local(ids)
local idlbl : value label id
label save `idlbl' using "id_labels.do", replace

* 1.3.1.2. Storage of estimates
capture postutil clear
tempname loo
postfile `loo' str14 outcome excl_id est60 lb60 ub60 p60 est90 lb90 ub90 p90 F p_global conv nobs ngroups using "loo_cytokines.dta", replace

* 1.3.1.3. Full-sample model (excl_id = 0) and one model per excluded participant
foreach y in il6 log_tnf {
foreach i in 0 `ids' {
capture quietly mixed `y' i.group##i.time c.`y'_base if time>0 & id!=`i' || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)
if _rc==0 {
local conv = e(converged)
local nobs = e(N)
local ngr = el(e(N_g),1,1)
quietly lincom 1.group, small
local e60 = r(estimate)
local l60 = r(lb)
local u60 = r(ub)
local p60 = r(p)
quietly lincom 1.group + 1.group#2.time, small
local e90 = r(estimate)
local l90 = r(lb)
local u90 = r(ub)
local p90 = r(p)
quietly test 1.group 1.group#2.time, small
local F = r(F)
local pg = r(p)
post `loo' ("`y'") (`i') (`e60') (`l60') (`u60') (`p60') (`e90') (`l90') (`u90') (`p90') (`F') (`pg') (`conv') (`nobs') (`ngr')
}
else post `loo' ("`y'") (`i') (.) (.) (.) (.) (.) (.) (.) (.) (.) (.) (0) (.) (.)
}
}
postclose `loo'

* 1.3.1.4. Convergence check
preserve
use "loo_cytokines.dta", clear
do "id_labels.do"
label values excl_id `idlbl'
gen full = (excl_id==0)
count if conv != 1

* 1.3.1.5. Full-sample estimates (must match the main analysis)
list outcome est60 est90 p_global if full, noobs

* 1.3.1.6. Shift of the 90-day treatment effect relative to the full sample (analysis scale)
bysort outcome: egen full90 = max(cond(full, est90, .))
gen shift90 = abs(est90 - full90)

* 1.3.1.7. Back-transformation of tumour necrosis factor-alpha to ratios of geometric means
foreach v in est60 lb60 ub60 est90 lb90 ub90 {
replace `v' = 10^`v' if outcome=="log_tnf"
}

* 1.3.1.8. Range of estimates and of the global P value across exclusions
table outcome if !full, statistic(min est60) statistic(max est60) statistic(min est90) statistic(max est90) statistic(min p_global) statistic(max p_global)

* 1.3.1.9. Number of exclusions with a global P value below 0.05
gen sig = (p_global < 0.05) if !missing(p_global)
table outcome if !full, statistic(count sig) statistic(sum sig)

* 1.3.1.10. The three most influential participants per outcome
gsort outcome -shift90
by outcome: list excl_id est90 lb90 ub90 p90 p_global ngroups if _n<=3 & !full, noobs

* 1.3.1.11. Exclusions that change the conclusion of the global test relative to the full sample
bysort outcome: egen full_sig = max(cond(full, sig, .))
list outcome excl_id est90 lb90 ub90 p90 p_global if !full & sig != full_sig, noobs

* 1.3.1.12. Graphical summary: 90-day treatment effect by excluded participant
foreach y in il6 log_tnf {
quietly summarize est90 if full & outcome=="`y'", meanonly
local f = r(mean)
local null = cond("`y'"=="log_tnf", 1, 0)
local ttl = cond("`y'"=="il6", "Interleukin-6", "Tumour necrosis factor-alpha")
local ytl = cond("`y'"=="il6", "Treatment effect at 90 days (pg/mL)", "Ratio of geometric means at 90 days")
twoway (rcap lb90 ub90 excl_id if !full & outcome=="`y'", lcolor("#5B6E8C")) (scatter est90 excl_id if !full & outcome=="`y'", mcolor("#384358") msize(medsmall)), yline(`f', lcolor("#541A2E") lwidth(medthick)) yline(`null', lcolor(gs10) lpattern(dash)) xlabel(`ids', valuelabel angle(90) labsize(vsmall) nogrid) xtitle("Excluded participant") ytitle("`ytl'") title("`ttl'", size(medium)) legend(off) graphregion(color(white)) plotregion(color(white)) name(loo_`y', replace)
}
graph combine loo_il6 loo_log_tnf, rows(1) graphregion(color(white)) name(loo_cytokines_panel, replace)
graph export "loo_cytokines_panel.png", name(loo_cytokines_panel) width(3000) replace
restore


* 1.3.1.2. Sensitivity analysis: exclusion of 90-day samples without valid bacterial quantification
mixed il6 i.group##i.time c.il6_base if time>0 & !(miss90_bac & time==2) || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)
lincom 1.group, small
lincom 1.group + 1.group#2.time, small
test 1.group 1.group#2.time, small


log close

