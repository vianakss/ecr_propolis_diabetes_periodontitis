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

* 1.2. Multivariate normal imputation of 60- and 90-day values on group and baseline value (100 imputations)
mi impute mvn hba1c1 hba1c2 = i.group hba1c0, add(100) rseed(20261101)
mi impute mvn glucose1 glucose2 = i.group glucose0, replace rseed(20261102)
mi impute mvn log_crp1 log_crp2 = i.group log_crp0, replace rseed(20261103)

* 1.3. Diagnostics: observed versus imputed follow-up values by group
foreach y in hba1c glucose log_crp {
foreach t in 1 2 {
tabstat `y'`t' if _mi_m==0, by(group) statistics(n mean sd min max) nototal
capture noisily tabstat `y'`t' if _mi_m>0 & miss_`y'`t', by(group) statistics(n mean sd min max) nototal
}
}

* 1.4. Number of negative imputed values (raw-scale outcomes)
foreach y in hba1c glucose {
foreach t in 1 2 {
quietly count if _mi_m>0 & miss_`y'`t'
local n = r(N)
quietly count if _mi_m>0 & miss_`y'`t' & `y'`t'<0
display "`y'`t': " r(N) " negative values among " `n' " imputed values"
}
}

* 1.5. Storage of the imputed dataset
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
*** MULTIPLE IMPUTATION - PERIODONTAL PARAMETERS ***
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

* 1.2. Multivariate normal imputation of 60- and 90-day values on group and baseline value (100 imputations)
local opt "add(100)"
local s = 20261201
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
mi impute mvn `y'1 `y'2 = i.group `y'0, `opt' rseed(`s')
local opt "replace"
local s = `s' + 1
}

* 1.3. Diagnostics: observed versus imputed 90-day values by group
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
tabstat `y'2 if _mi_m==0, by(group) statistics(n mean sd min max) nototal
capture noisily tabstat `y'2 if _mi_m>0 & miss_`y'2, by(group) statistics(n mean sd min max) nototal
}

* 1.4. Number of negative imputed values
foreach y in plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus {
foreach t in 1 2 {
quietly count if _mi_m>0 & miss_`y'`t'
local n = r(N)
quietly count if _mi_m>0 & miss_`y'`t' & `y'`t'<0
display "`y'`t': " r(N) " negative values among " `n' " imputed values"
}
}

* 1.5. Storage of the imputed dataset
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

* 0.4. Indicators of originally missing 90-day values (used in the diagnostics and in the tipping-point analysis)
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


*** 3. TIPPING-POINT ANALYSIS FOR F. NUCLEATUM (DELTA ADJUSTMENT) ***
* 3.1. Shift of the imputed 90-day values in the Control group by delta, from 0 to 1 log10 in steps of 0.05, with no adjustment in the Test group
* Reported for the 90-day treatment effect and for the global test with and without the equal fraction of missing information assumption
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

log close
