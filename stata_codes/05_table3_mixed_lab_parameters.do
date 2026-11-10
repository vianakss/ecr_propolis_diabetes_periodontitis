
********************************************************************************
*** LONGITUDINAL ANCOVA MIXED MODEL - LABORATORY PARAMETERS ***
********************************************************************************

xtset id time


capture log close
log using "log_table_mixed_model_laboratory_parameters.log", replace

use "dataset_propolis_analysis.dta", clear

*** 0. DATA PREPARATION AND MISSING DATA PATTERN ***
* 0.1. Baseline values as covariates
capture drop hba1c_base
bysort id: egen hba1c_base = max(cond(time==0, hba1c, .))
capture drop glucose_base
bysort id: egen glucose_base = max(cond(time==0, glucose, .))
capture drop log_crp
gen log_crp = log10(crp)
capture drop log_crp_base
bysort id: egen log_crp_base = max(cond(time==0, log_crp, .))


* 0.2. Missing data pattern
misstable summarize hba1c glucose log_crp
table group time, statistic(count hba1c) statistic(count glucose) statistic(count log_crp)

* 0.3. Baseline availability (participants retained by the analysis of covariance)
count if time==0 & !missing(hba1c)
count if time==0 & !missing(glucose)
count if time==0 & !missing(log_crp)

* 0.4. Balance of the laboratory complete-case selection between groups
capture drop has_lab
bysort id: egen has_lab = max(!missing(il6))
preserve
keep if time==0
tabulate group has_lab, row exact
restore


*** 1. LABORATORY PARAMETERS ***
* 1.1. HBA1C
* 1.1.1. Observed descriptive values
table group time, statistic(count hba1c) statistic(mean hba1c) statistic(sd hba1c)
table group time, statistic(median hba1c) statistic(p25 hba1c) statistic(p75 hba1c)

* 1.1.2. Longitudinal analysis of covariance
mixed hba1c i.group##i.time c.hba1c_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.1.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.1.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.1.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.1.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_hba1c, replace)

* 1.1.7. Model diagnostics
capture drop fit_hba1c res_hba1c
predict fit_hba1c, fitted
predict res_hba1c, rstandard
scatter res_hba1c fit_hba1c, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - HbA1c", size(medium)) name(resfit_hba1c, replace)
qnorm res_hba1c if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids HbA1c - 60 days", size(medium)) name(qnorm_hba1c_t1, replace)
qnorm res_hba1c if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids HbA1c - 90 days", size(medium)) name(qnorm_hba1c_t2, replace)
graph combine qnorm_hba1c_t1 qnorm_hba1c_t2, name(hba1c_resid_panel, replace)
graph export "hba1c_resid_panel.png", name(hba1c_resid_panel) width(3000) replace



* 1.2. FASTING PLASMA GLUCOSE
* 1.2.1. Observed descriptive values
table group time, statistic(count glucose) statistic(mean glucose) statistic(sd glucose)
table group time, statistic(median glucose) statistic(p25 glucose) statistic(p75 glucose)

* 1.2.2. Longitudinal analysis of covariance
mixed glucose i.group##i.time c.glucose_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.2.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.2.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.2.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.2.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_glucose, replace)

* 1.2.7. Model diagnostics
capture drop fit_glucose res_glucose
predict fit_glucose, fitted
predict res_glucose, rstandard
scatter res_glucose fit_glucose, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Glucose", size(medium)) name(resfit_glucose, replace)
qnorm res_glucose if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids glucose - 60 days", size(medium)) name(qnorm_glucose_t1, replace)
qnorm res_glucose if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids glucose - 90 days", size(medium)) name(qnorm_glucose_t2, replace)
graph combine qnorm_glucose_t1 qnorm_glucose_t2, name(glucose_resid_panel, replace)
graph export "glucose_resid_panel.png", name(glucose_resid_panel) width(3000) replace


* 1.3. C-REACTIVE PROTEIN (log10 scale)
* 1.3.1. Observed descriptive values
table group time, statistic(count crp) statistic(mean crp) statistic(sd crp) statistic(median crp) statistic(p25 crp) statistic(p75 crp)
table group time, statistic(count log_crp) statistic(mean log_crp) statistic(sd log_crp)

* 1.3.2. Longitudinal analysis of covariance
mixed log_crp i.group##i.time c.log_crp_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.3.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small
display "Ratio of geometric means at 60 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.3.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small
display "Ratio of geometric means at 90 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.3.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.3.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_log_crp, replace)

* 1.3.7. Model diagnostics
capture drop fit_log_crp res_log_crp
predict fit_log_crp, fitted
predict res_log_crp, rstandard
scatter res_log_crp fit_log_crp, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - CRP (log10)", size(medium)) name(resfit_log_crp, replace)
qnorm res_log_crp if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids CRP (log10) - 60 days", size(medium)) name(qnorm_log_crp_t1, replace)
qnorm res_log_crp if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids CRP (log10) - 90 days", size(medium)) name(qnorm_log_crp_t2, replace)
graph combine qnorm_log_crp_t1 qnorm_log_crp_t2, name(log_crp_resid_panel, replace)
graph export "log_crp_resid_panel.png", name(log_crp_resid_panel) width(3000) replace

graph combine resfit_hba1c resfit_glucose resfit_log_crp, rows(1) name(scatter_laboratory_panel, replace)
graph export "scatter_laboratory_panel.png", name(scatter_laboratory_panel) width(3000) replace

log close


* 1.4. Sensitivity analysis
capture log close
log using "log_leave_one_out_laboratory_parameters.log", replace
* 1.4.1. Leave-one-out analysis
* 1.4.1.1. Participant identifiers (0 is reserved for the full sample)
assert id != 0
levelsof id, local(ids)

* 1.4.1.2. Storage of estimates
capture postutil clear
tempname loo
postfile `loo' str12 outcome excl_id est60 lb60 ub60 p60 est90 lb90 ub90 p90 F p_global conv nobs ngroups using "loo_laboratory.dta", replace

* 1.4.1.3. Full-sample model (excl_id = 0) and one model per excluded participant
foreach y in hba1c glucose log_crp {
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

* 1.4.1.4. Convergence check
preserve
use "loo_laboratory.dta", clear
gen full = (excl_id==0)
count if conv != 1

* 1.4.1.5. Full-sample estimates (must match the main analysis)
list outcome est60 est90 p_global if full, noobs

* 1.4.1.6. Shift of the 90-day treatment effect relative to the full sample (analysis scale)
bysort outcome: egen full90 = max(cond(full, est90, .))
gen shift90 = abs(est90 - full90)

* 1.4.1.7. Back-transformation of C-reactive protein to ratios of geometric means
foreach v in est60 lb60 ub60 est90 lb90 ub90 {
replace `v' = 10^`v' if outcome=="log_crp"
}

* 1.4.1.8. Range of estimates and of the global P value across exclusions
table outcome if !full, statistic(min est60) statistic(max est60) statistic(min est90) statistic(max est90) statistic(min p_global) statistic(max p_global)

* 1.4.1.9. Number of exclusions with a global P value below 0.05
gen sig = (p_global < 0.05) if !missing(p_global)
table outcome if !full, statistic(count sig) statistic(sum sig)

* 1.4.1.10. The three most influential participants per outcome
gsort outcome -shift90
by outcome: list excl_id est90 lb90 ub90 p90 p_global ngroups if _n<=3 & !full, noobs

* 1.4.1.11. Graphical summary: 90-day treatment effect by excluded participant
foreach y in hba1c glucose log_crp {
quietly summarize est90 if full & outcome=="`y'", meanonly
local f = r(mean)
local null = cond("`y'"=="log_crp", 1, 0)
local ttl = cond("`y'"=="hba1c", "HbA1c", cond("`y'"=="glucose", "Fasting plasma glucose", "C-reactive protein"))
local ytl = cond("`y'"=="hba1c", "Treatment effect at 90 days (%)", cond("`y'"=="glucose", "Treatment effect at 90 days (mg/dL)", "Ratio of geometric means at 90 days"))
twoway (rcap lb90 ub90 excl_id if !full & outcome=="`y'", lcolor("#5B6E8C")) (scatter est90 excl_id if !full & outcome=="`y'", mcolor("#384358") msize(medsmall)), yline(`f', lcolor("#541A2E") lwidth(medthick)) yline(`null', lcolor(gs10) lpattern(dash)) xlabel(1 5 10 15 20 24) xtitle("Excluded participant") ytitle("`ytl'") title("`ttl'", size(medium)) legend(off) graphregion(color(white)) plotregion(color(white)) name(loo_`y', replace)
}
graph combine loo_hba1c loo_glucose loo_log_crp, rows(1) graphregion(color(white)) name(loo_laboratory_panel, replace)
graph export "loo_laboratory_panel.png", name(loo_laboratory_panel) width(3000) replace
restore

log close
