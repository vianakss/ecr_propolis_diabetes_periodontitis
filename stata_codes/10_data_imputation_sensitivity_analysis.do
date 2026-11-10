********************************************************************************
*** MULTIPLE IMPUTATION - BLOOD CHEMISTRY ***
********************************************************************************
version 19
set more off
capture log close
log using "log_mi_blood_chemistry.log", text replace
use "dataset_-_ecr_propolis_long.dta", clear
xtset, clear

*** 0. DATA PREPARATION ***
* 0.1. Analysis scale (C-reactive protein on the log10 scale)
assert crp>0 if !missing(crp)
capture drop log_crp
gen log_crp = log10(crp)
keep id group time hba1c glucose log_crp

* 0.2. Wide format: 0 = baseline, 1 = 60 days, 2 = 90 days
reshape wide hba1c glucose log_crp, i(id) j(time)

* 0.3. Missing data pattern and confirmation that baseline values are complete
misstable patterns hba1c1 hba1c2 glucose1 glucose2 log_crp1 log_crp2, frequency
foreach y in hba1c glucose log_crp {
assert !missing(`y'0)
}

* 0.4. Indicators of originally missing follow-up values and their distribution by group
foreach y in hba1c glucose log_crp {
gen miss_`y'1 = missing(`y'1)
gen miss_`y'2 = missing(`y'2)
tabulate group miss_`y'1
tabulate group miss_`y'2
}


*** 1. MULTIPLE IMPUTATION UNDER MISSING AT RANDOM ***
* 1.1. Declaration of the imputation data
mi set flong
mi register imputed hba1c1 hba1c2 glucose1 glucose2 log_crp1 log_crp2
mi register regular hba1c0 glucose0 log_crp0

* 1.2. Convergence of the data augmentation: worst linear function over 1,000 iterations
foreach y in hba1c glucose log_crp {
mi impute mvn `y'1 `y'2 = i.group `y'0, mcmconly burnin(1000) rseed(20261100) savewlf(wlf_`y', replace)
preserve
use wlf_`y', clear
tsset iter
tsline wlf, ytitle("Worst linear function") xtitle("Burn-in period") title("`y'", size(medium)) graphregion(color(white)) plotregion(color(white)) name(wlf_`y', replace)
restore
}
graph combine wlf_hba1c wlf_glucose wlf_log_crp, rows(1) graphregion(color(white)) name(wlf_blood_panel, replace)
graph export "wlf_blood_panel.png", name(wlf_blood_panel) width(3000) replace

* 1.3. Multivariate normal imputation of 60- and 90-day values on group and baseline value (100 imputations)
mi impute mvn hba1c1 hba1c2 = i.group hba1c0, add(100) rseed(20261101)
mi impute mvn glucose1 glucose2 = i.group glucose0, replace rseed(20261102)
mi impute mvn log_crp1 log_crp2 = i.group log_crp0, replace rseed(20261103)

* 1.4. Diagnostics: observed versus imputed follow-up values by group
foreach y in hba1c glucose log_crp {
foreach t in 1 2 {
tabstat `y'`t' if _mi_m==0, by(group) statistics(n mean sd min max) nototal
capture noisily tabstat `y'`t' if _mi_m>0 & miss_`y'`t', by(group) statistics(n mean sd min max) nototal
}
}

* 1.5. Number of negative imputed values (raw-scale outcomes)
foreach y in hba1c glucose {
foreach t in 1 2 {
quietly count if _mi_m>0 & miss_`y'`t'
local n = r(N)
quietly count if _mi_m>0 & miss_`y'`t' & `y'`t'<0
display "`y'`t': " r(N) " negative values among " `n' " imputed values"
}
}

* 1.6. Storage of the imputed dataset
save "mi_blood_chemistry.dta", replace


*** 2. ANALYSIS OF THE IMPUTED DATA (MISSING AT RANDOM) ***
* 2.1. HbA1c and fasting plasma glucose: treatment effects at 60 and 90 days and global test
foreach y in hba1c glucose {
mi estimate: mvreg `y'1 `y'2 = i.group `y'0
display "Complete-data degrees of freedom: " e(df_c_mi)
mi test [`y'1]1.group [`y'2]1.group
display "Global test `y': F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

* 2.2. C-reactive protein: ratios of geometric means at 60 and 90 days and global test
mi estimate: mvreg log_crp1 log_crp2 = i.group log_crp0
display "Complete-data degrees of freedom: " e(df_c_mi)
matrix T = r(table)
foreach t in 1 2 {
local c = colnumb(T, "log_crp`t':1.group")
display "Ratio of geometric means, time `t': " %5.2f 10^T[1,`c'] " [" %5.2f 10^T[5,`c'] "; " %5.2f 10^T[6,`c'] "]"
}
mi test [log_crp1]1.group [log_crp2]1.group
display "Global test log_crp: F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)

* 2.3. Sensitivity of the global test to the equal fraction of missing information assumption
foreach y in hba1c glucose log_crp {
quietly mi estimate: mvreg `y'1 `y'2 = i.group `y'0
quietly mi test [`y'1]1.group [`y'2]1.group
display "`y' equal FMI:        F = " %5.2f r(F) ", P = " %5.3f r(p)
quietly mi test [`y'1]1.group [`y'2]1.group, ufmitest
display "`y' unrestricted FMI: F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

* 2.4. Adjusted means by group at 60 and 90 days, evaluated at the mean baseline value (Rubin's rules)
* C-reactive protein: back-transform the estimates and confidence limits with 10^x to obtain geometric means
foreach y in hba1c glucose log_crp {
quietly summarize `y'0 if _mi_m==0
local mb = r(mean)
display as text _newline "Adjusted means: `y' (baseline mean = " %9.4f `mb' ")"
mi estimate (c60: [`y'1]_cons + [`y'1]`y'0*`mb') (t60: [`y'1]_cons + [`y'1]1.group + [`y'1]`y'0*`mb') (c90: [`y'2]_cons + [`y'2]`y'0*`mb') (t90: [`y'2]_cons + [`y'2]1.group + [`y'2]`y'0*`mb'), nocoef: mvreg `y'1 `y'2 = i.group `y'0
}

log close




********************************************************************************
*** MULTIPLE IMPUTATION AND TIPPING-POINT ANALYSIS - PERIODONTAL PARAMETERS ***
********************************************************************************
version 19
set more off
capture log close
log using "log_mi_periodontal.log", text replace
use "dataset_-_ecr_propolis_long.dta", clear
xtset, clear

*** 0. DATA PREPARATION ***
* 0.1. Periodontal outcomes in wide format: 0 = baseline, 1 = 60 days, 2 = 90 days
keep id group time plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus
reshape wide plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus, i(id) j(time)

* 0.2. Missing data pattern and confirmation that baseline values are complete
misstable patterns plaque1 plaque2 pctbop1 pctbop2 meanpd1 meanpd2 meancal1 meancal2 pctsite4mm1 pctsite4mm2 pctsite5_6mm1 pctsite5_6mm2 pctsite7plus1 pctsite7plus2, frequency
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
assert !missing(`y'0)
}

* 0.3. Indicators of originally missing follow-up values and their distribution by group
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
gen miss_`y'1 = missing(`y'1)
gen miss_`y'2 = missing(`y'2)
tabulate group miss_`y'2
}


*** 1. MULTIPLE IMPUTATION UNDER MISSING AT RANDOM ***
* 1.1. Declaration of the imputation data
mi set flong
mi register imputed plaque1 plaque2 pctbop1 pctbop2 meanpd1 meanpd2 meancal1 meancal2 pctsite4mm1 pctsite4mm2 pctsite5_6mm1 pctsite5_6mm2 pctsite7plus1 pctsite7plus2
mi register regular plaque0 pctbop0 meanpd0 meancal0 pctsite4mm0 pctsite5_6mm0 pctsite7plus0

* 1.2. Convergence of the data augmentation: worst linear function over 1,000 iterations
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
mi impute mvn `y'1 `y'2 = i.group `y'0, mcmconly burnin(1000) rseed(20261200) savewlf(wlf_`y', replace)
preserve
use wlf_`y', clear
tsset iter
tsline wlf, ytitle("Worst linear function") xtitle("Burn-in period") title("`y'", size(medium)) graphregion(color(white)) plotregion(color(white)) name(wlf_`y', replace)
restore
}
graph combine wlf_plaque wlf_pctbop wlf_meanpd wlf_meancal wlf_pctsite4mm wlf_pctsite5_6mm wlf_pctsite7plus, cols(4) graphregion(color(white)) name(wlf_periodontal_panel, replace)
graph export "wlf_periodontal_panel.png", name(wlf_periodontal_panel) width(3000) replace

* 1.3. Multivariate normal imputation of 60- and 90-day values on group and baseline value (100 imputations)
local opt "add(100)"
local s = 20261201
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
mi impute mvn `y'1 `y'2 = i.group `y'0, `opt' rseed(`s')
local opt "replace"
local s = `s' + 1
}

* 1.4. Diagnostics: observed versus imputed 90-day values by group
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
tabstat `y'2 if _mi_m==0, by(group) statistics(n mean sd min max) nototal
capture noisily tabstat `y'2 if _mi_m>0 & miss_`y'2, by(group) statistics(n mean sd min max) nototal
}

* 1.5. Number of negative imputed values
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
foreach t in 1 2 {
quietly count if _mi_m>0 & miss_`y'`t'
local n = r(N)
quietly count if _mi_m>0 & miss_`y'`t' & `y'`t'<0
display "`y'`t': " r(N) " negative values among " `n' " imputed values"
}
}

* 1.6. Storage of the imputed dataset
save "mi_periodontal.dta", replace


*** 2. ANALYSIS OF THE IMPUTED DATA (MISSING AT RANDOM) ***
* 2.1. Multivariate analysis of covariance: treatment effects at 60 and 90 days and global test
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
mi estimate: mvreg `y'1 `y'2 = i.group `y'0
display "Complete-data degrees of freedom: " e(df_c_mi)
mi test [`y'1]1.group [`y'2]1.group
display "Global test `y': F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

* 2.2. Sensitivity of the global test to the equal fraction of missing information assumption
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
quietly mi estimate: mvreg `y'1 `y'2 = i.group `y'0
quietly mi test [`y'1]1.group [`y'2]1.group
display "`y' equal FMI:        F = " %5.2f r(F) ", P = " %5.3f r(p)
quietly mi test [`y'1]1.group [`y'2]1.group, ufmitest
display "`y' unrestricted FMI: F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

* 2.3. Adjusted means by group at 60 and 90 days, evaluated at the mean baseline value (Rubin's rules)
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
quietly summarize `y'0 if _mi_m==0
local mb = r(mean)
display as text _newline "Adjusted means: `y' (baseline mean = " %9.4f `mb' ")"
mi estimate (c60: [`y'1]_cons + [`y'1]`y'0*`mb') (t60: [`y'1]_cons + [`y'1]1.group + [`y'1]`y'0*`mb') (c90: [`y'2]_cons + [`y'2]`y'0*`mb') (t90: [`y'2]_cons + [`y'2]1.group + [`y'2]`y'0*`mb'), nocoef: mvreg `y'1 `y'2 = i.group `y'0
}


*** 3. TIPPING-POINT ANALYSIS OF THE PRIMARY OUTCOME (MEAN PROBING DEPTH) ***
* 3.1. Shift of imputed 90-day values by delta (mm) in each group, from -1 to 1 in steps of 0.25
capture postutil clear
tempname tp
postfile `tp' dc dt est90 lb90 ub90 p90 F p_global using "tipping_periodontal.dta", replace
forvalues ic = -4/4 {
forvalues it = -4/4 {
local dc = `ic'/4
local dt = `it'/4
preserve
quietly replace meanpd2 = meanpd2 + `dc' if _mi_m>0 & miss_meanpd2 & group==0
quietly replace meanpd2 = meanpd2 + `dt' if _mi_m>0 & miss_meanpd2 & group==1
quietly mi estimate: mvreg meanpd1 meanpd2 = i.group meanpd0
matrix T = r(table)
local c = colnumb(T, "meanpd2:1.group")
local e90 = T[1,`c']
local p90 = T[4,`c']
local l90 = T[5,`c']
local u90 = T[6,`c']
quietly mi test [meanpd1]1.group [meanpd2]1.group
post `tp' (`dc') (`dt') (`e90') (`l90') (`u90') (`p90') (r(F)) (r(p))
restore
}
}
postclose `tp'

* 3.2. Results: scenario without adjustment, grid of global P values and heat map
use "tipping_periodontal.dta", clear
list est90 lb90 ub90 p90 F p_global if dc==0 & dt==0, noobs
table (dc) (dt), statistic(mean p_global) nformat(%5.3f)
summarize p_global
twoway (contour p_global dt dc, ccuts(0.01 0.05 0.10 0.20) ccolors("#541A2E" "#A0607A" "#E6E9EE" "#8A9BB5" "#384358") interp(none)) (scatteri 0 0, msymbol(X) mcolor(black) msize(large)), xlabel(-1(0.25)1) ylabel(-1(0.25)1, angle(0)) xtitle("Delta in the Control group (mm)") ytitle("Delta in the Test group (mm)") ztitle("Overall P value") title("Mean probing depth", size(medium)) legend(off) graphregion(color(white)) plotregion(color(white)) name(tipping_meanpd, replace)
graph export "tipping_meanpd.png", name(tipping_meanpd) width(3000) replace
graph save tipping_meanpd "tipping_meanpd.gph", replace


*** 4. CONFIRMATION WITH PREDICTIVE MEAN MATCHING FOR PERCENTAGE OUTCOMES (BOUNDED AT ZERO) ***
* The imputed dataset is reloaded so that the stored multivariate normal imputations are not overwritten
use "mi_periodontal.dta", clear
local s = 20261301
foreach y in plaque pctbop pctsite4mm pctsite5_6mm pctsite7plus {
mi impute pmm `y'2 = i.group `y'0 `y'1, knn(5) replace rseed(`s')
local s = `s' + 1
mi estimate: mvreg `y'1 `y'2 = i.group `y'0
mi test [`y'1]1.group [`y'2]1.group
display "Global test `y' (predictive mean matching): F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

log close




********************************************************************************
*** MULTIPLE IMPUTATION AND TIPPING-POINT ANALYSIS - BACTERIAL COUNTS ***
********************************************************************************
version 19
set more off
capture log close
log using "log_mi_bacterial_counts.log", text replace
use "dataset_-_ecr_propolis_long.dta", clear
xtset, clear

*** 0. DATA PREPARATION ***
* 0.1. Required variables (log10 scale) and participants with bacterial quantification (the remaining four were not assayed by design)
foreach v in log_un log_pg log_fn log_so {
confirm numeric variable `v'
}
capture drop lab_bac
bysort id: egen lab_bac = max(time==0 & !missing(log_un))
keep if lab_bac
keep id group time log_un log_pg log_fn log_so

* 0.2. Wide format: 0 = baseline, 1 = 60 days, 2 = 90 days
reshape wide log_un log_pg log_fn log_so, i(id) j(time)

* 0.3. Missing data pattern and confirmation that only 90-day values are missing
misstable patterns log_un0 log_un1 log_un2 log_pg2 log_fn2 log_so2, frequency
foreach y in log_un log_pg log_fn log_so {
assert !missing(`y'0, `y'1)
}

* 0.4. Indicators of originally missing 90-day values (used in the tipping-point analysis)
foreach y in log_un log_pg log_fn log_so {
gen miss_`y'2 = missing(`y'2)
}
assert miss_log_un2==miss_log_pg2 & miss_log_un2==miss_log_fn2 & miss_log_un2==miss_log_so2
tabulate group miss_log_un2


*** 1. MULTIPLE IMPUTATION UNDER MISSING AT RANDOM ***
* 1.1. Declaration of the imputation data
mi set flong
mi register imputed log_un2 log_pg2 log_fn2 log_so2
mi register regular log_un0 log_un1 log_pg0 log_pg1 log_fn0 log_fn1 log_so0 log_so1

* 1.2. Bayesian linear regression imputation of 90-day values on group, baseline and 60-day values (100 imputations)
mi impute regress log_un2 = i.group log_un0 log_un1, add(100) rseed(20261029)
mi impute regress log_pg2 = i.group log_pg0 log_pg1, replace rseed(20261030)
mi impute regress log_fn2 = i.group log_fn0 log_fn1, replace rseed(20261031)
mi impute regress log_so2 = i.group log_so0 log_so1, replace rseed(20261032)

* 1.3. Diagnostics: observed versus imputed 90-day values by group
foreach y in log_un log_pg log_fn log_so {
tabstat `y'2 if _mi_m==0, by(group) statistics(n mean sd min max) nototal
tabstat `y'2 if _mi_m>0 & miss_`y'2, by(group) statistics(n mean sd min max) nototal
}

* 1.4. Storage of the imputed dataset
save "mi_bacterial_counts.dta", replace


*** 2. ANALYSIS OF THE IMPUTED DATA (MISSING AT RANDOM) ***
* 2.1. Multivariate analysis of covariance: treatment effects at 60 and 90 days and global test
foreach y in log_un log_pg log_fn log_so {
mi estimate: mvreg `y'1 `y'2 = i.group `y'0
display "Complete-data degrees of freedom: " e(df_c_mi)
mi test [`y'1]1.group [`y'2]1.group
display "Global test `y': F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

* 2.2. Sensitivity of the global test to the equal fraction of missing information assumption
foreach y in log_un log_pg log_fn log_so {
quietly mi estimate: mvreg `y'1 `y'2 = i.group `y'0
quietly mi test [`y'1]1.group [`y'2]1.group
display "`y' equal FMI:        F = " %5.2f r(F) ", P = " %5.3f r(p)
quietly mi test [`y'1]1.group [`y'2]1.group, ufmitest
display "`y' unrestricted FMI: F(" r(df) "; " %6.2f r(df_r) ") = " %5.2f r(F) ", P = " %5.3f r(p)
}

* 2.3. Adjusted means by group at 60 and 90 days on the log10 scale, evaluated at the mean baseline value (Rubin's rules)
foreach y in log_un log_pg log_fn log_so {
quietly summarize `y'0 if _mi_m==0
local mb = r(mean)
display as text _newline "Adjusted means: `y' (baseline mean = " %9.4f `mb' ")"
mi estimate (c60: [`y'1]_cons + [`y'1]`y'0*`mb') (t60: [`y'1]_cons + [`y'1]1.group + [`y'1]`y'0*`mb') (c90: [`y'2]_cons + [`y'2]`y'0*`mb') (t90: [`y'2]_cons + [`y'2]1.group + [`y'2]`y'0*`mb'), nocoef: mvreg `y'1 `y'2 = i.group `y'0
}


*** 3. TIPPING-POINT ANALYSIS (DELTA ADJUSTMENT) ***
* 3.1. Shift of imputed 90-day values by delta (log10 units) in each group, from -2 to 2 in steps of 0.5
capture postutil clear
tempname tp
postfile `tp' str8 outcome dc dt est90 lb90 ub90 p90 F p_global using "tipping_bacterial.dta", replace
foreach y in log_fn log_so {
forvalues ic = -4/4 {
forvalues it = -4/4 {
local dc = `ic'/2
local dt = `it'/2
preserve
quietly replace `y'2 = `y'2 + `dc' if _mi_m>0 & miss_`y'2 & group==0
quietly replace `y'2 = `y'2 + `dt' if _mi_m>0 & miss_`y'2 & group==1
quietly mi estimate: mvreg `y'1 `y'2 = i.group `y'0
matrix T = r(table)
local c = colnumb(T, "`y'2:1.group")
local e90 = T[1,`c']
local p90 = T[4,`c']
local l90 = T[5,`c']
local u90 = T[6,`c']
quietly mi test [`y'1]1.group [`y'2]1.group
post `tp' ("`y'") (`dc') (`dt') (`e90') (`l90') (`u90') (`p90') (r(F)) (r(p))
restore
}
}
}
postclose `tp'

* 3.2. Results: scenario without adjustment and grids of P values (global test and 90-day contrast)
preserve
use "tipping_bacterial.dta", clear
list outcome est90 lb90 ub90 p90 F p_global if dc==0 & dt==0, noobs
table (dc) (dt) if outcome=="log_fn", statistic(mean p_global) nformat(%5.3f)
table (dc) (dt) if outcome=="log_so", statistic(mean p_global) nformat(%5.3f)
table (dc) (dt) if outcome=="log_fn", statistic(mean p90) nformat(%5.3f)
table (dc) (dt) if outcome=="log_so", statistic(mean p90) nformat(%5.3f)
restore

* 3.3. Refined tipping point for F. nucleatum: delta in the Control group from 0 to 1 in steps of 0.05, no adjustment in the Test group
* Reported for the 90-day contrast and for the global test with and without the equal fraction of missing information assumption
use "mi_bacterial_counts.dta", clear
forvalues ic = 0/20 {
local dc = `ic'/20
preserve
quietly replace log_fn2 = log_fn2 + `dc' if _mi_m>0 & miss_log_fn2 & group==0
quietly mi estimate: mvreg log_fn1 log_fn2 = i.group log_fn0
matrix T = r(table)
local c = colnumb(T, "log_fn2:1.group")
local p90 = T[4,`c']
quietly mi test [log_fn1]1.group [log_fn2]1.group
local pe = r(p)
quietly mi test [log_fn1]1.group [log_fn2]1.group, ufmitest
display "delta Control = " %4.2f `dc' "   90-day P = " %5.3f `p90' "   global P (equal FMI) = " %5.3f `pe' "   global P (unrestricted FMI) = " %5.3f r(p)
restore
}

* 3.4. Heat maps of the P value for the 90-day treatment effect across the grid of delta values (not affected by the equal fraction of missing information assumption)
use "tipping_bacterial.dta", clear
foreach y in log_fn log_so {
local ttl = cond("`y'"=="log_fn", "{it:Fusobacterium nucleatum}", "{it:Streptococcus oralis}")
twoway (contour p90 dt dc if outcome=="`y'", ccuts(0.01 0.05 0.10 0.20) ccolors("#541A2E" "#A0607A" "#E6E9EE" "#8A9BB5" "#384358") interp(none)) (scatteri 0 0, msymbol(X) mcolor(black) msize(large)), xlabel(-2(0.5)2) ylabel(-2(0.5)2, angle(0)) xtitle("Delta in the Control group (log10)") ytitle("Delta in the Test group (log10)") ztitle("P value, 90-day treatment effect") title("`ttl'", size(medium)) legend(off) graphregion(color(white)) plotregion(color(white)) name(tp_`y', replace)
}
graph combine tp_log_fn tp_log_so, rows(1) graphregion(color(white)) name(tipping_bacterial_panel, replace)
graph export "tipping_bacterial_panel_p90.png", name(tipping_bacterial_panel) width(3000) replace

log close
