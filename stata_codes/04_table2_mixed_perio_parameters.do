********************************************************************************
*** LONGITUDINAL ANCOVA MIXED MODEL - PERIODONTAL PARAMETERS ***
********************************************************************************

xtset id time


capture log close
log using "log_table_mixed_model_periodontal_parameters.log", replace

use "dataset_propolis_analysis.dta", clear

*** 0. DATA PREPARATION AND MISSING DATA PATTERN ***
* 0.1. Baseline values as covariates
capture drop plaque_base
bysort id: egen plaque_base = max(cond(time==0, plaque, .))
capture drop pctbop_base
bysort id: egen pctbop_base = max(cond(time==0, pctbop, .))
capture drop meanpd_base
bysort id: egen meanpd_base = max(cond(time==0, meanpd, .))
capture drop meancal_base
bysort id: egen meancal_base = max(cond(time==0, meancal, .))
capture drop pctsite4mm_base
bysort id: egen pctsite4mm_base = max(cond(time==0, pctsite4mm, .))
capture drop pctsite5_6mm_base
bysort id: egen pctsite5_6mm_base = max(cond(time==0, pctsite5_6mm, .))
capture drop pctsite7plus_base
bysort id: egen pctsite7plus_base = max(cond(time==0, pctsite7plus, .))

* 0.2. Missing data pattern
misstable summarize plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus
table group time, statistic(count plaque) statistic(count pctbop) statistic(count meanpd) statistic(count meancal) statistic(count pctsite4mm) statistic(count pctsite5_6mm) statistic(count pctsite7plus)

* 0.3. Baseline availability (participants retained by the analysis of covariance)
count if time==0 & !missing(plaque)
count if time==0 & !missing(pctbop)
count if time==0 & !missing(meanpd)
count if time==0 & !missing(meancal)
count if time==0 & !missing(pctsite4mm)
count if time==0 & !missing(pctsite5_6mm)
count if time==0 & !missing(pctsite7plus)


*** 1. PERIODONTAL PARAMETERS ***
* 1.1. PLAQUE INDEX
* 1.1.1. Observed descriptive values
table group time, statistic(count plaque) statistic(mean plaque) statistic(sd plaque)
table group time, statistic(median plaque) statistic(p25 plaque) statistic(p75 plaque)

* 1.1.2. Longitudinal analysis of covariance
mixed plaque i.group##i.time c.plaque_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.1.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.1.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.1.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.1.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_plaque, replace)

* 1.1.7. Model diagnostics
capture drop fit_plaque res_plaque
predict fit_plaque, fitted
predict res_plaque, rstandard
scatter res_plaque fit_plaque, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Plaque index", size(medium)) name(resfit_plaque, replace)
qnorm res_plaque if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids plaque index - 60 days", size(medium)) name(qnorm_plaque_t1, replace)
qnorm res_plaque if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids plaque index - 90 days", size(medium)) name(qnorm_plaque_t2, replace)
graph combine qnorm_plaque_t1 qnorm_plaque_t2, name(plaque_resid_panel, replace)
graph export "plaque_resid_panel.png", name(plaque_resid_panel) width(3000) replace



* 1.2. BLEEDING ON PROBING (%)
* 1.2.1. Observed descriptive values
table group time, statistic(count pctbop) statistic(mean pctbop) statistic(sd pctbop)
table group time, statistic(median pctbop) statistic(p25 pctbop) statistic(p75 pctbop)

* 1.2.2. Longitudinal analysis of covariance
mixed pctbop i.group##i.time c.pctbop_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.2.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.2.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.2.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.2.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_pctbop, replace)

* 1.2.7. Model diagnostics
capture drop fit_pctbop res_pctbop
predict fit_pctbop, fitted
predict res_pctbop, rstandard
scatter res_pctbop fit_pctbop, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - BOP", size(medium)) name(resfit_pctbop, replace)
qnorm res_pctbop if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids BOP - 60 days", size(medium)) name(qnorm_pctbop_t1, replace)
qnorm res_pctbop if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids BOP - 90 days", size(medium)) name(qnorm_pctbop_t2, replace)
graph combine qnorm_pctbop_t1 qnorm_pctbop_t2, name(pctbop_resid_panel, replace)
graph export "pctbop_resid_panel.png", name(pctbop_resid_panel) width(3000) replace



* 1.3. MEAN PROBING DEPTH
* 1.3.1. Observed descriptive values
table group time, statistic(count meanpd) statistic(mean meanpd) statistic(sd meanpd)
table group time, statistic(median meanpd) statistic(p25 meanpd) statistic(p75 meanpd)

* 1.3.2. Longitudinal analysis of covariance
mixed meanpd i.group##i.time c.meanpd_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.3.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.3.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.3.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.3.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_meanpd, replace)

* 1.3.7. Model diagnostics
capture drop fit_meanpd res_meanpd
predict fit_meanpd, fitted
predict res_meanpd, rstandard
scatter res_meanpd fit_meanpd, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Mean PD", size(medium)) name(resfit_meanpd, replace)
qnorm res_meanpd if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids mean PD - 60 days", size(medium)) name(qnorm_meanpd_t1, replace)
qnorm res_meanpd if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids mean PD - 90 days", size(medium)) name(qnorm_meanpd_t2, replace)
graph combine qnorm_meanpd_t1 qnorm_meanpd_t2, name(meanpd_resid_panel, replace)
graph export "meanpd_resid_panel.png", name(meanpd_resid_panel) width(3000) replace



* 1.4. MEAN CLINICAL ATTACHMENT LEVEL
* 1.4.1. Observed descriptive values
table group time, statistic(count meancal) statistic(mean meancal) statistic(sd meancal)
table group time, statistic(median meancal) statistic(p25 meancal) statistic(p75 meancal)

* 1.4.2. Longitudinal analysis of covariance
mixed meancal i.group##i.time c.meancal_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.4.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.4.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.4.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.4.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_meancal, replace)

* 1.4.7. Model diagnostics
capture drop fit_meancal res_meancal
predict fit_meancal, fitted
predict res_meancal, rstandard
scatter res_meancal fit_meancal, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Mean CAL", size(medium)) name(resfit_meancal, replace)
qnorm res_meancal if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids mean CAL - 60 days", size(medium)) name(qnorm_meancal_t1, replace)
qnorm res_meancal if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids mean CAL - 90 days", size(medium)) name(qnorm_meancal_t2, replace)
graph combine qnorm_meancal_t1 qnorm_meancal_t2, name(meancal_resid_panel, replace)
graph export "meancal_resid_panel.png", name(meancal_resid_panel) width(3000) replace



* 1.5. SITES WITH PD 4 MM (%)
* 1.5.1. Observed descriptive values
table group time, statistic(count pctsite4mm) statistic(mean pctsite4mm) statistic(sd pctsite4mm)
table group time, statistic(median pctsite4mm) statistic(p25 pctsite4mm) statistic(p75 pctsite4mm)

* 1.5.2. Longitudinal analysis of covariance
mixed pctsite4mm i.group##i.time c.pctsite4mm_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.5.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.5.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.5.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.5.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_pctsite4mm, replace)

* 1.5.7. Model diagnostics
capture drop fit_pctsite4mm res_pctsite4mm
predict fit_pctsite4mm, fitted
predict res_pctsite4mm, rstandard
scatter res_pctsite4mm fit_pctsite4mm, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Sites PD 4 mm", size(medium)) name(resfit_pctsite4mm, replace)
qnorm res_pctsite4mm if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids sites PD 4 mm - 60 days", size(medium)) name(qnorm_pctsite4mm_t1, replace)
qnorm res_pctsite4mm if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids sites PD 4 mm - 90 days", size(medium)) name(qnorm_pctsite4mm_t2, replace)
graph combine qnorm_pctsite4mm_t1 qnorm_pctsite4mm_t2, name(pctsite4mm_resid_panel, replace)
graph export "pctsite4mm_resid_panel.png", name(pctsite4mm_resid_panel) width(3000) replace



* 1.6. SITES WITH PD 5-6 MM (%)
* 1.6.1. Observed descriptive values
table group time, statistic(count pctsite5_6mm) statistic(mean pctsite5_6mm) statistic(sd pctsite5_6mm)
table group time, statistic(median pctsite5_6mm) statistic(p25 pctsite5_6mm) statistic(p75 pctsite5_6mm)

* 1.6.2. Longitudinal analysis of covariance
mixed pctsite5_6mm i.group##i.time c.pctsite5_6mm_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.6.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.6.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.6.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.6.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_pctsite5_6mm, replace)

* 1.6.7. Model diagnostics
capture drop fit_pctsite5_6mm res_pctsite5_6mm
predict fit_pctsite5_6mm, fitted
predict res_pctsite5_6mm, rstandard
scatter res_pctsite5_6mm fit_pctsite5_6mm, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Sites PD 5-6 mm", size(medium)) name(resfit_pctsite5_6mm, replace)
qnorm res_pctsite5_6mm if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids sites PD 5-6 mm - 60 days", size(medium)) name(qnorm_pctsite5_6mm_t1, replace)
qnorm res_pctsite5_6mm if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids sites PD 5-6 mm - 90 days", size(medium)) name(qnorm_pctsite5_6mm_t2, replace)
graph combine qnorm_pctsite5_6mm_t1 qnorm_pctsite5_6mm_t2, name(pctsite5_6mm_resid_panel, replace)
graph export "pctsite5_6mm_resid_panel.png", name(pctsite5_6mm_resid_panel) width(3000) replace



* 1.7. SITES WITH PD >=7 MM (%)
* 1.7.1. Observed descriptive values
table group time, statistic(count pctsite7plus) statistic(mean pctsite7plus) statistic(sd pctsite7plus)
table group time, statistic(median pctsite7plus) statistic(p25 pctsite7plus) statistic(p75 pctsite7plus)

* 1.7.2. Longitudinal analysis of covariance
mixed pctsite7plus i.group##i.time c.pctsite7plus_base if time>0 || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)

* 1.7.3. Estimated treatment effect at 60 days: Test - Control
lincom 1.group, small

* 1.7.4. Estimated treatment effect at 90 days: Test - Control
lincom 1.group + 1.group#2.time, small

* 1.7.5. Global test of treatment effect across both follow-up time points
test 1.group 1.group#2.time, small

* 1.7.6. Model-estimated means and profile
margins group#time
marginsplot, plot1opts(lcolor("#384358") mcolor("#384358")) plot2opts(lcolor("#541A2E") mcolor("#541A2E")) ci1opts(lcolor("#384358")) ci2opts(lcolor("#541A2E")) graphregion(color(white)) plotregion(color(white)) name(margins_pctsite7plus, replace)

* 1.7.7. Model diagnostics
capture drop fit_pctsite7plus res_pctsite7plus
predict fit_pctsite7plus, fitted
predict res_pctsite7plus, rstandard
scatter res_pctsite7plus fit_pctsite7plus, mcolor("#5B6E8C") msize(medsmall) yline(0, lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Scatter plot - Sites PD >=7 mm", size(medium)) name(resfit_pctsite7plus, replace)
qnorm res_pctsite7plus if time==1, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids sites PD >=7 mm - 60 days", size(medium)) name(qnorm_pctsite7plus_t1, replace)
qnorm res_pctsite7plus if time==2, mcolor("#5B6E8C") rlopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) plotregion(color(white)) title("Resids sites PD >=7 mm - 90 days", size(medium)) name(qnorm_pctsite7plus_t2, replace)
graph combine qnorm_pctsite7plus_t1 qnorm_pctsite7plus_t2, name(pctsite7plus_resid_panel, replace)
graph export "pctsite7plus_resid_panel.png", name(pctsite7plus_resid_panel) width(3000) replace

graph combine resfit_plaque resfit_pctbop resfit_meanpd resfit_meancal resfit_pctsite4mm resfit_pctsite5_6mm resfit_pctsite7plus, cols(4) name(scatter_periodontal_panel, replace)
graph export "scatter_periodontal_panel.png", name(scatter_periodontal_panel) width(3000) replace

log close
 

 
capture log close
log using "log_leave_one_out_periodontal_parameters.log", replace
* 1.8. Sensitivity analysis - influential residuals
* 1.8.1. Identification of observations with absolute standardized residual > 3
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
list id group time `y' res_`y' if abs(res_`y') > 3 & !missing(res_`y'), noobs
}

* 1.8.2. Flags and full trajectories of the participants with influential observations
capture drop infl_plaque pt_infl_plaque
gen infl_plaque = abs(res_plaque) > 3 & !missing(res_plaque)
bysort id: egen pt_infl_plaque = max(infl_plaque)
list id group time plaque plaque_base if pt_infl_plaque, noobs sepby(id)
capture drop infl_meancal pt_infl_meancal
gen infl_meancal = abs(res_meancal) > 3 & !missing(res_meancal)
bysort id: egen pt_infl_meancal = max(infl_meancal)
list id group time meancal meancal_base if pt_infl_meancal, noobs sepby(id)

* 1.8.3. Plaque index without the influential observation
mixed plaque i.group##i.time c.plaque_base if time>0 & !infl_plaque || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)
lincom 1.group, small
lincom 1.group + 1.group#2.time, small
test 1.group 1.group#2.time, small

* 1.8.4. Mean clinical attachment level without the influential observation
mixed meancal i.group##i.time c.meancal_base if time>0 & !infl_meancal || id:, nocons residuals(unstructured, t(time)) reml dfmethod(kroger)
lincom 1.group, small
lincom 1.group + 1.group#2.time, small
test 1.group 1.group#2.time, small



* 1.9. Leave-one-out analysis
* 1.9.1. Participant identifiers and labels (0 is reserved for the full sample)
assert id != 0
levelsof id, local(ids)
local idlbl : value label id
label save `idlbl' using "id_labels.do", replace

* 1.9.2. Storage of estimates
capture postutil clear
tempname loo
postfile `loo' str14 outcome excl_id est60 lb60 ub60 p60 est90 lb90 ub90 p90 F p_global conv nobs ngroups using "loo_periodontal.dta", replace

* 1.9.3. Full-sample model (excl_id = 0) and one model per excluded participant
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
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

* 1.9.4. Convergence check
preserve
use "loo_periodontal.dta", clear
do "id_labels.do"
label values excl_id `idlbl'
gen full = (excl_id==0)
count if conv != 1

* 1.9.5. Full-sample estimates (must match the main analysis)
list outcome est60 est90 p_global if full, noobs

* 1.9.6. Shift of the 90-day treatment effect relative to the full sample
bysort outcome: egen full90 = max(cond(full, est90, .))
gen shift90 = abs(est90 - full90)

* 1.9.7. Range of estimates and of the global P value across exclusions
table outcome if !full, statistic(min est60) statistic(max est60) statistic(min est90) statistic(max est90) statistic(min p_global) statistic(max p_global)

* 1.9.8. Number of exclusions with a global P value below 0.05
gen sig = (p_global < 0.05) if !missing(p_global)
table outcome if !full, statistic(count sig) statistic(sum sig)

* 1.9.9. The three most influential participants per outcome
gsort outcome -shift90
by outcome: list excl_id est90 lb90 ub90 p90 p_global ngroups if _n<=3 & !full, noobs

* 1.9.10. Graphical summary: 90-day treatment effect by excluded participant
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
quietly summarize est90 if full & outcome=="`y'", meanonly
local f = r(mean)
if "`y'"=="plaque" local ttl "Plaque index"
if "`y'"=="pctbop" local ttl "Bleeding on probing"
if "`y'"=="meanpd" local ttl "Mean probing depth"
if "`y'"=="meancal" local ttl "Mean clinical attachment level"
if "`y'"=="pctsite4mm" local ttl "Sites with PD 4 mm"
if "`y'"=="pctsite5_6mm" local ttl "Sites with PD 5-6 mm"
if "`y'"=="pctsite7plus" local ttl "Sites with PD >=7 mm"
local ytl = cond(inlist("`y'","meanpd","meancal"), "Treatment effect at 90 days (mm)", "Treatment effect at 90 days (%)")
twoway (rcap lb90 ub90 excl_id if !full & outcome=="`y'", lcolor("#5B6E8C")) (scatter est90 excl_id if !full & outcome=="`y'", mcolor("#384358") msize(medsmall)), yline(`f', lcolor("#541A2E") lwidth(medthick)) yline(0, lcolor(gs10) lpattern(dash)) xlabel(1(1)24, valuelabel angle(90) labsize(vsmall)) xtitle("Excluded participant") ytitle("`ytl'") title("`ttl'", size(medium)) legend(off) graphregion(color(white)) plotregion(color(white)) name(loo_`y', replace)
}
graph combine loo_plaque loo_pctbop loo_meanpd loo_meancal loo_pctsite4mm loo_pctsite5_6mm loo_pctsite7plus, cols(4) graphregion(color(white)) name(loo_periodontal_panel, replace)
graph export "loo_periodontal_panel.png", name(loo_periodontal_panel) width(3000) replace
restore

* 1.9.11. Exclusions with a global P value below 0.05
preserve
use "loo_periodontal.dta", clear
do "id_labels.do"
label values excl_id id
list outcome excl_id est90 lb90 ub90 p90 p_global if p_global<0.05 & excl_id!=0, noobs
restore
log close




log using "log_endpoint_feres.log", text replace

*** 1. PROPORTION ACHIEVING THE ENDPOINT BY GROUP AND TIME
foreach t in 0 1 2 {
display as text _newline "Time `t'"
tabulate group endpoint if time==`t', row exact
}

*** 2. RISK DIFFERENCE, RISK RATIO AND FISHER'S EXACT TEST AT 60 AND 90 DAYS (Test versus Control)
foreach t in 1 2 {
display as text _newline "Time `t'"
cs endpoint group if time==`t', exact
}

*** 3. SENSITIVITY: WORST- AND BEST-CASE SCENARIOS FOR MISSING 90-DAY VALUES
preserve
keep if time==2
count if missing(endpoint)
generate byte worst = endpoint
replace worst = 0 if missing(endpoint) & group==1
replace worst = 1 if missing(endpoint) & group==0
generate byte best = endpoint
replace best = 1 if missing(endpoint) & group==1
replace best = 0 if missing(endpoint) & group==0
display as text _newline "Worst case for the Test group"
cs worst group, exact
display as text _newline "Best case for the Test group"
cs best group, exact
restore

*** 4. RISK DIFFERENCE WITH NEWCOMBE HYBRID SCORE CONFIDENCE INTERVALS (method 10)
foreach t in 1 2 {
quietly count if time==`t' & group==1 & !missing(endpoint)
local n1 = r(N)
quietly count if time==`t' & group==1 & endpoint==1
local x1 = r(N)
quietly count if time==`t' & group==0 & !missing(endpoint)
local n0 = r(N)
quietly count if time==`t' & group==0 & endpoint==1
local x0 = r(N)
quietly cii proportions `n1' `x1', wilson
local p1 = r(proportion)
local l1 = r(lb)
local u1 = r(ub)
quietly cii proportions `n0' `x0', wilson
local p0 = r(proportion)
local l0 = r(lb)
local u0 = r(ub)
local d = `p1' - `p0'
local lo = `d' - sqrt((`p1'-`l1')^2 + (`u0'-`p0')^2)
local hi = `d' + sqrt((`u1'-`p1')^2 + (`p0'-`l0')^2)
display "Time `t': risk difference (Test-Control) = " %6.2f 100*`d' " [" %6.2f 100*`lo' "; " %6.2f 100*`hi' "] percentage points"
}

log close
