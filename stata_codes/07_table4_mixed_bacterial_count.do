
********************************************************************************
*** LONGITUDINAL ANCOVA MIXED MODEL - BACTERIAL COUNTS ***
********************************************************************************
capture log close
log using "log_table_mixed_model_bacterial_counts.log", text replace
xtset id time

*** 0. DATA PREPARATION AND MISSING DATA PATTERN ***
* 0.1. Baseline values as covariates (log10 scale)
capture drop log_un_base
bysort id: egen log_un_base = max(cond(time==0, log_un, .))
capture drop log_pg_base
bysort id: egen log_pg_base = max(cond(time==0, log_pg, .))
capture drop log_fn_base
bysort id: egen log_fn_base = max(cond(time==0, log_fn, .))
capture drop log_so_base
bysort id: egen log_so_base = max(cond(time==0, log_so, .))

* 0.2. Missing data pattern
misstable summarize log_un log_pg log_fn log_so
table group time, statistic(count log_un) statistic(count log_pg) statistic(count log_fn) statistic(count log_so)

* 0.3. Baseline availability (participants retained by the analysis of covariance)
count if time==0 & !missing(log_un)
count if time==0 & !missing(log_pg)
count if time==0 & !missing(log_fn)
count if time==0 & !missing(log_so)


*** 1. BACTERIAL COUNTS ***
* 1.1. TOTAL BACTERIAL LOAD (log10 scale)
* 1.1.1. Observed descriptive values
table group time, statistic(count log_un) statistic(mean log_un) statistic(sd log_un)
table group time, statistic(median log_un) statistic(p25 log_un) statistic(p75 log_un)

* 1.1.2. Longitudinal analysis of covariance
mixed log_un i.group##i.time c.log_un_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.1.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small
display "Ratio of geometric means at 60 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.1.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small
display "Ratio of geometric means at 90 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.1.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.1.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_log_un, replace)

* 1.1.7. Model diagnostics
capture drop fit_log_un res_log_un
predict fit_log_un, fitted
predict res_log_un, rstandard
scatter res_log_un fit_log_un, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Total bacterial load (log10)", size(medium)) name(resfit_log_un, replace)
qnorm res_log_un if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids total load (log10) - 60 days", size(medium)) name(qnorm_log_un_t1, replace)
qnorm res_log_un if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids total load (log10) - 90 days", size(medium)) name(qnorm_log_un_t2, replace)
graph combine qnorm_log_un_t1 qnorm_log_un_t2, name(log_un_resid_panel, replace)
graph export "log_un_resid_panel.png", name(log_un_resid_panel) width(3000) replace



* 1.2. PORPHYROMONAS GINGIVALIS (log10 scale)
* 1.2.1. Observed descriptive values
table group time, statistic(count log_pg) statistic(mean log_pg) statistic(sd log_pg)
table group time, statistic(median log_pg) statistic(p25 log_pg) statistic(p75 log_pg)

* 1.2.2. Longitudinal analysis of covariance
mixed log_pg i.group##i.time c.log_pg_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.2.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small
display "Ratio of geometric means at 60 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.2.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small
display "Ratio of geometric means at 90 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.2.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.2.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_log_pg, replace)

* 1.2.7. Model diagnostics
capture drop fit_log_pg res_log_pg
predict fit_log_pg, fitted
predict res_log_pg, rstandard
scatter res_log_pg fit_log_pg, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - P. gingivalis (log10)", size(medium)) name(resfit_log_pg, replace)
qnorm res_log_pg if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids P. gingivalis (log10) - 60 days", size(medium)) name(qnorm_log_pg_t1, replace)
qnorm res_log_pg if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids P. gingivalis (log10) - 90 days", size(medium)) name(qnorm_log_pg_t2, replace)
graph combine qnorm_log_pg_t1 qnorm_log_pg_t2, name(log_pg_resid_panel, replace)
graph export "log_pg_resid_panel.png", name(log_pg_resid_panel) width(3000) replace



* 1.3. FUSOBACTERIUM NUCLEATUM (log10 scale)
* 1.3.1. Observed descriptive values
table group time, statistic(count log_fn) statistic(mean log_fn) statistic(sd log_fn)
table group time, statistic(median log_fn) statistic(p25 log_fn) statistic(p75 log_fn)

* 1.3.2. Longitudinal analysis of covariance
mixed log_fn i.group##i.time c.log_fn_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

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
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_log_fn, replace)

* 1.3.7. Model diagnostics
capture drop fit_log_fn res_log_fn
predict fit_log_fn, fitted
predict res_log_fn, rstandard
scatter res_log_fn fit_log_fn, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - F. nucleatum (log10)", size(medium)) name(resfit_log_fn, replace)
qnorm res_log_fn if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids F. nucleatum (log10) - 60 days", size(medium)) name(qnorm_log_fn_t1, replace)
qnorm res_log_fn if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids F. nucleatum (log10) - 90 days", size(medium)) name(qnorm_log_fn_t2, replace)
graph combine qnorm_log_fn_t1 qnorm_log_fn_t2, name(log_fn_resid_panel, replace)
graph export "log_fn_resid_panel.png", name(log_fn_resid_panel) width(3000) replace



* 1.4. STREPTOCOCCUS ORALIS (log10 scale)
* 1.4.1. Observed descriptive values
table group time, statistic(count log_so) statistic(mean log_so) statistic(sd log_so)
table group time, statistic(median log_so) statistic(p25 log_so) statistic(p75 log_so)

* 1.4.2. Longitudinal analysis of covariance
mixed log_so i.group##i.time c.log_so_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.4.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small
display "Ratio of geometric means at 60 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.4.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small
display "Ratio of geometric means at 90 days: " %5.2f 10^r(estimate) " [" %5.2f 10^r(lb) "; " %5.2f 10^r(ub) "]"

* 1.4.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small
display %9.6f r(p)

* 1.4.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_log_so, replace)

* 1.4.7. Model diagnostics
capture drop fit_log_so res_log_so
predict fit_log_so, fitted
predict res_log_so, rstandard
scatter res_log_so fit_log_so, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - S. oralis (log10)", size(medium)) name(resfit_log_so, replace)
qnorm res_log_so if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids S. oralis (log10) - 60 days", size(medium)) name(qnorm_log_so_t1, replace)
qnorm res_log_so if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids S. oralis (log10) - 90 days", size(medium)) name(qnorm_log_so_t2, replace)
graph combine qnorm_log_so_t1 qnorm_log_so_t2, name(log_so_resid_panel, replace)
graph export "log_so_resid_panel.png", name(log_so_resid_panel) width(3000) replace


graph combine resfit_log_un resfit_log_pg resfit_log_fn resfit_log_so, cols(2) name(scatter_bacterial_panel, replace)
graph export "scatter_bacterial_panel.png", name(scatter_bacterial_panel) width(3000) replace

log close






* 1.5. Sensitivity analysis
capture log close
log using "log_leave_one_out_bacterial_counts.log", replace
* 1.5.1. Leave-one-out analysis
* 1.5.1.1. Participant identifiers with bacterial data and labels (0 is reserved for the full sample)
assert id != 0
levelsof id if !missing(log_un), local(ids)
local idlbl : value label id
label save `idlbl' using "id_labels.do", replace

* 1.5.1.2. Storage of estimates
capture postutil clear
tempname loo
postfile `loo' str14 outcome excl_id est60 lb60 ub60 p60 est90 lb90 ub90 p90 F p_global conv nobs ngroups using "loo_bacterial.dta", replace

* 1.5.1.3. Full-sample model (excl_id = 0) and one model per excluded participant
foreach y in log_un log_pg log_fn log_so {
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

* 1.5.1.4. Convergence check
preserve
use "loo_bacterial.dta", clear
do "id_labels.do"
label values excl_id `idlbl'
gen full = (excl_id==0)
count if conv != 1

* 1.5.1.5. Full-sample estimates (must match the main analysis)
list outcome est60 est90 p_global if full, noobs

* 1.5.1.6. Shift of the 90-day treatment effect relative to the full sample
bysort outcome: egen full90 = max(cond(full, est90, .))
gen shift90 = abs(est90 - full90)

* 1.5.1.7. Range of estimates and of the global P value across exclusions
table outcome if !full, statistic(min est60) statistic(max est60) statistic(min est90) statistic(max est90) statistic(min p_global) statistic(max p_global)

* 1.5.1.8. Number of exclusions with a global P value below 0.05
gen sig = (p_global < 0.05) if !missing(p_global)
table outcome if !full, statistic(count sig) statistic(sum sig)

* 1.5.1.9. The three most influential participants per outcome
gsort outcome -shift90
by outcome: list excl_id est90 lb90 ub90 p90 p_global ngroups if _n<=3 & !full, noobs

* 1.5.1.10. Exclusions that change the conclusion of the global test relative to the full sample
bysort outcome: egen full_sig = max(cond(full, sig, .))
list outcome excl_id est90 lb90 ub90 p90 p_global if !full & sig != full_sig, noobs

* 1.5.1.11. Graphical summary: 90-day treatment effect by excluded participant
foreach y in log_un log_pg log_fn log_so {
quietly summarize est90 if full & outcome=="`y'", meanonly
local f = r(mean)
if "`y'"=="log_un" local ttl "Total bacterial load"
if "`y'"=="log_pg" local ttl "{it:Porphyromonas gingivalis}"
if "`y'"=="log_fn" local ttl "{it:Fusobacterium nucleatum}"
if "`y'"=="log_so" local ttl "{it:Streptococcus oralis}"
twoway (rcap lb90 ub90 excl_id if !full & outcome=="`y'", lcolor("#5B6E8C")) (scatter est90 excl_id if !full & outcome=="`y'", mcolor("#384358") msize(medsmall)), yline(`f', lcolor("#541A2E") lwidth(medthick)) yline(0, lcolor(gs10) lpattern(dash)) xlabel(`ids', valuelabel angle(90) labsize(vsmall) nogrid) xtitle("Excluded participant") ytitle("Treatment effect at 90 days (log10)") title("`ttl'", size(medium)) legend(off) graphregion(color(white)) plotregion(color(white)) name(loo_`y', replace)
}
graph combine loo_log_un loo_log_pg loo_log_fn loo_log_so, cols(2) graphregion(color(white)) name(loo_bacterial_panel, replace)
graph export "loo_bacterial_panel.png", name(loo_bacterial_panel) width(3000) replace
restore

foreach y in hba1c glucose log_crp il6 log_tnf log_un log_pg log_fn log_so {
count if abs(res_`y')>3 & !missing(res_`y')
list id group time `y' res_`y' if abs(res_`y')>3 & !missing(res_`y'), noobs
}

log close

