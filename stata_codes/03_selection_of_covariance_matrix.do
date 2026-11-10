version 19.0

capture log close
log using "covariance_structure_selection.log", text replace

use "dataset_propolis_analysis.dta", clear

* 0. Design and variable checks
isid id time
assert inlist(time,0,1,2)
local outcomes "hba1c log_crp glucose plaque pctbop meanpd meancal pctsite4mm pctsite5_6mm pctsite7plus il6 log_tnf log_un log_pg log_fn log_so"
foreach y of local outcomes {
confirm numeric variable `y'
}

* 1. Time in days: used only in t() of the exponential structure; i.time remains the fixed-effect factor
capture drop days
generate days = 0 if time==0
replace days = 60 if time==1
replace days = 90 if time==2
label variable days "Days since baseline"
assert !missing(days)

* 2. Detailed inspection for the primary outcome (mean probing depth)
mixed meanpd i.group##i.time || id:, noconstant residuals(unstructured, t(time)) reml
estimates store un
mixed meanpd i.group##i.time || id:, noconstant residuals(exchangeable) reml
estimates store cs
mixed meanpd i.group##i.time || id:, noconstant residuals(exponential, t(days)) reml
estimates store expo
mixed meanpd i.group##i.time || id:, residuals(exponential, t(days)) reml
estimates store ri_expo
foreach m in un cs expo ri_expo {
estimates restore `m'
display as text "Model `m': converged = " e(converged)
estat wcorrelation, format(%5.3f)
}
* Likelihood-ratio tests (all nested in unstructured; same fixed effects, so REML comparison is valid)
* ri_expo vs un: the null may lie on the boundary (random-intercept variance = 0), so the test is conservative
lrtest un cs
lrtest un expo
lrtest un ri_expo
estimates stats un cs expo ri_expo

* 3. All outcomes: information criteria, convergence and random-intercept standard deviation
local s1 "noconstant residuals(unstructured, t(time))"
local s2 "noconstant residuals(exchangeable)"
local s3 "noconstant residuals(exponential, t(days))"
local s4 "residuals(exponential, t(days))"
local n1 "Unstructured"
local n2 "Exchangeable"
local n3 "Exponential"
local n4 "RI+Exponential"
capture postutil clear
tempname cov
postfile `cov' str16 outcome str16 structure int rc byte converged double N ll df aic bic sd_ri using "covariance_comparison.dta", replace
foreach y of local outcomes {
forvalues k = 1/4 {
capture quietly mixed `y' i.group##i.time || id:, `s`k'' reml
local rc = _rc
if `rc'==0 {
quietly estat ic
matrix S = r(S)
local sdri = .
if `k'==4 {
capture local sdri = exp(_b[lns1_1_1:_cons])
}
post `cov' ("`y'") ("`n`k''") (`rc') (e(converged)) (S[1,1]) (S[1,3]) (S[1,4]) (S[1,5]) (S[1,6]) (`sdri')
}
else {
post `cov' ("`y'") ("`n`k''") (`rc') (.) (.) (.) (.) (.) (.) (.)
}
}
}
postclose `cov'

* 4. Summary table
preserve
use "covariance_comparison.dta", clear
bysort outcome: egen double sd_N = sd(N)
assert sd_N==0 | missing(sd_N)
drop sd_N
bysort outcome (aic): generate byte best_aic = (_n==1 & !missing(aic))
bysort outcome (bic): generate byte best_bic = (_n==1 & !missing(bic))
format ll aic bic %10.2f
format sd_ri %9.2e
sort outcome aic
list outcome structure rc converged N ll df aic bic sd_ri best_aic best_bic, sepby(outcome) noobs abbreviate(16)
tabulate structure best_aic if best_aic==1
export delimited using "covariance_comparison.csv", replace
restore
log close
