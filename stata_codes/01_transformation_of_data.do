********************************************************************************
*** DISTRIBUTION DIAGNOSTICS AND LOG10 TRANSFORMATION OF OUTCOMES ***
********************************************************************************
version 19
capture log close
log using "log_distribution_diagnostics.log", text replace
use "dataset_-_ecr_propolis_long.dta", clear


********************************************************************************
*** 1. BACTERIAL COUNTS ***
********************************************************************************
* 1.1. Check for zero or negative counts (log10 is undefined for values <= 0)
foreach b in un pg fn so {
count if `b'_count <= 0
}

* 1.2. Logarithmic transformation and check of losses
foreach b in un pg fn so {
capture drop log_`b'
gen log_`b' = log10(`b'_count) if `b'_count > 0 & !missing(`b'_count)
count if missing(log_`b') & !missing(`b'_count)
}
label variable log_un "Total bacterial count, log10"
label variable log_pg "P. gingivalis, log10"
label variable log_fn "F. nucleatum, log10"
label variable log_so "S. oralis, log10"

* 1.3. Descriptive diagnosis of raw vs. log10 scale
summarize un_count pg_count fn_count so_count
summarize log_un log_pg log_fn log_so
display _newline "Skewness and kurtosis"
display "Variable" _col(20) "Skewness" _col(36) "Kurtosis"
foreach v of varlist un_count pg_count fn_count so_count log_un log_pg log_fn log_so {
quietly summarize `v', detail
display "`v'" _col(20) %10.2f r(skewness) _col(36) %10.2f r(kurtosis)
}

* 1.4. Histogram distributions
foreach b in un pg fn so {
histogram `b'_count, by(time, cols(3) title("`b' - raw scale")) fcolor("#5B6E8C") lcolor("#5B6E8C") normal normopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) name(histogram_`b'_raw, replace)
histogram log_`b', by(time, cols(3) title("`b' - log10")) fcolor("#5B6E8C") lcolor("#5B6E8C") normal normopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) name(histogram_`b'_log, replace)
graph combine histogram_`b'_raw histogram_`b'_log, cols(1) graphregion(color(white)) name(histogram_`b', replace)
graph export "histogram_`b'.png", name(histogram_`b') width(3000) replace
}

* 1.5. Box distributions
foreach b in un pg fn so {
graph box `b'_count, over(time) box(1, color("#5B6E8C") lcolor("#5B6E8C")) medtype(cline) medline(lcolor("#5B6E8C") lwidth(medthick)) cwhiskers lines(lcolor("#5B6E8C")) marker(1, mcolor("#541A2E")) title("`b' - raw scale") graphregion(color(white)) name(box_`b'_raw, replace)
graph box log_`b', over(time) box(1, color("#5B6E8C") lcolor("#5B6E8C")) medtype(cline) medline(lcolor("#5B6E8C") lwidth(medthick)) cwhiskers lines(lcolor("#5B6E8C")) marker(1, mcolor("#541A2E")) title("`b' - log10") graphregion(color(white)) name(box_`b'_log, replace)
graph combine box_`b'_raw box_`b'_log, cols(2) graphregion(color(white)) name(box_`b', replace)
graph export "box_`b'.png", name(box_`b') width(3000) replace
}

* 1.6. Scatterplot matrices of the repeated measurements across assessment time points (raw and log10 scales)
foreach b in un pg fn so {
preserve
keep id group time `b'_count log_`b'
reshape wide `b'_count log_`b', i(id group) j(time)
graph matrix `b'_count0 `b'_count1 `b'_count2, half msymbol(O) mfcolor("#5B6E8C") mlcolor("#5B6E8C") title("`b' - raw scale") graphregion(color(white)) name(matrix_`b'_raw, replace)
graph matrix log_`b'0 log_`b'1 log_`b'2, half msymbol(O) mfcolor("#5B6E8C") mlcolor("#5B6E8C") title("`b' - log10") graphregion(color(white)) name(matrix_`b'_log, replace)
graph combine matrix_`b'_raw matrix_`b'_log, cols(2) graphregion(color(white)) name(matrix_`b', replace)
graph export "matrix_`b'.png", name(matrix_`b') width(3000) replace
restore
}


********************************************************************************
*** 2. SALIVARY CYTOKINES ***
********************************************************************************
* 2.1. Logarithmic transformation and check of losses
capture drop log_il6
capture drop log_tnf
gen log_il6 = log10(il6)
gen log_tnf = log10(tnf)
count if missing(log_il6) & !missing(il6)
count if missing(log_tnf) & !missing(tnf)
label variable log_il6 "Interleukin-6, log10"
label variable log_tnf "Tumour necrosis factor-alpha, log10"

* 2.2. Descriptive diagnosis of raw vs. log10 scale
summarize il6 tnf
summarize log_il6 log_tnf
display _newline "Skewness and kurtosis"
display "Variable" _col(20) "Skewness" _col(36) "Kurtosis"
foreach v of varlist il6 tnf log_il6 log_tnf {
quietly summarize `v', detail
display "`v'" _col(20) %10.2f r(skewness) _col(36) %10.2f r(kurtosis)
}

* 2.3. Histogram distributions
foreach b in il6 tnf {
histogram `b', by(time, cols(3) title("`b' - raw scale")) fcolor("#5B6E8C") lcolor("#5B6E8C") normal normopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) name(histogram_`b'_raw, replace)
histogram log_`b', by(time, cols(3) title("`b' - log10")) fcolor("#5B6E8C") lcolor("#5B6E8C") normal normopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) name(histogram_`b'_log, replace)
graph combine histogram_`b'_raw histogram_`b'_log, cols(1) graphregion(color(white)) name(histogram_`b', replace)
graph export "histogram_`b'.png", name(histogram_`b') width(3000) replace
}

* 2.4. Box distributions
foreach b in il6 tnf {
graph box `b', over(time) box(1, color("#5B6E8C") lcolor("#5B6E8C")) medtype(cline) medline(lcolor("#5B6E8C") lwidth(medthick)) cwhiskers lines(lcolor("#5B6E8C")) marker(1, mcolor("#541A2E")) title("`b' - raw scale") graphregion(color(white)) name(box_`b'_raw, replace)
graph box log_`b', over(time) box(1, color("#5B6E8C") lcolor("#5B6E8C")) medtype(cline) medline(lcolor("#5B6E8C") lwidth(medthick)) cwhiskers lines(lcolor("#5B6E8C")) marker(1, mcolor("#541A2E")) title("`b' - log10") graphregion(color(white)) name(box_`b'_log, replace)
graph combine box_`b'_raw box_`b'_log, cols(2) graphregion(color(white)) name(box_`b', replace)
graph export "box_`b'.png", name(box_`b') width(3000) replace
}

* 2.6. Scatterplot matrices of the repeated measurements across assessment time points (raw and log10 scales)
foreach b in il6 tnf {
preserve
keep id group time `b' log_`b'
reshape wide `b' log_`b', i(id group) j(time)
graph matrix `b'0 `b'1 `b'2, half msymbol(O) mfcolor("#5B6E8C") mlcolor("#5B6E8C") title("`b' - raw scale") graphregion(color(white)) name(matrix_`b'_raw, replace)
graph matrix log_`b'0 log_`b'1 log_`b'2, half msymbol(O) mfcolor("#5B6E8C") mlcolor("#5B6E8C") title("`b' - log10") graphregion(color(white)) name(matrix_`b'_log, replace)
graph combine matrix_`b'_raw matrix_`b'_log, cols(2) graphregion(color(white)) name(matrix_`b', replace)
graph export "matrix_`b'.png", name(matrix_`b') width(3000) replace
restore
}


********************************************************************************
*** 3. BLOOD CHEMISTRY ***
********************************************************************************
* 3.1. Logarithmic transformation and check of losses
capture drop log_hba1c
capture drop log_crp
capture drop log_glucose
gen log_hba1c = log10(hba1c)
gen log_crp = log10(crp)
gen log_glucose = log10(glucose)
count if missing(log_hba1c) & !missing(hba1c)
count if missing(log_crp) & !missing(crp)
count if missing(log_glucose) & !missing(glucose)
label variable log_hba1c "HbA1c, log10"
label variable log_crp "C-reactive protein, log10"
label variable log_glucose "Fasting plasma glucose, log10"

* 3.2. Descriptive diagnosis of raw vs. log10 scale
summarize hba1c crp glucose
summarize log_hba1c log_crp log_glucose
display _newline "Skewness and kurtosis"
display "Variable" _col(20) "Skewness" _col(36) "Kurtosis"
foreach v of varlist hba1c crp glucose log_hba1c log_crp log_glucose {
quietly summarize `v', detail
display "`v'" _col(20) %10.2f r(skewness) _col(36) %10.2f r(kurtosis)
}

* 3.3. Histogram distributions
foreach b in hba1c crp glucose {
histogram `b', by(time, cols(3) title("`b' - raw scale")) fcolor("#5B6E8C") lcolor("#5B6E8C") normal normopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) name(histogram_`b'_raw, replace)
histogram log_`b', by(time, cols(3) title("`b' - log10")) fcolor("#5B6E8C") lcolor("#5B6E8C") normal normopts(lcolor("#541A2E") lwidth(medthick)) graphregion(color(white)) name(histogram_`b'_log, replace)
graph combine histogram_`b'_raw histogram_`b'_log, cols(1) graphregion(color(white)) name(histogram_`b', replace)
graph export "histogram_`b'.png", name(histogram_`b') width(3000) replace
}

* 3.4. Box distributions
foreach b in hba1c crp glucose {
graph box `b', over(time) box(1, color("#5B6E8C") lcolor("#5B6E8C")) medtype(cline) medline(lcolor("#5B6E8C") lwidth(medthick)) cwhiskers lines(lcolor("#5B6E8C")) marker(1, mcolor("#541A2E")) title("`b' - raw scale") graphregion(color(white)) name(box_`b'_raw, replace)
graph box log_`b', over(time) box(1, color("#5B6E8C") lcolor("#5B6E8C")) medtype(cline) medline(lcolor("#5B6E8C") lwidth(medthick)) cwhiskers lines(lcolor("#5B6E8C")) marker(1, mcolor("#541A2E")) title("`b' - log10") graphregion(color(white)) name(box_`b'_log, replace)
graph combine box_`b'_raw box_`b'_log, cols(2) graphregion(color(white)) name(box_`b', replace)
graph export "box_`b'.png", name(box_`b') width(3000) replace
}

* 3.5. Scatterplot matrices of the repeated measurements across assessment time points (raw and log10 scales)
foreach b in hba1c crp glucose {
preserve
keep id group time `b' log_`b'
reshape wide `b' log_`b', i(id group) j(time)
graph matrix `b'0 `b'1 `b'2, half msymbol(O) mfcolor("#5B6E8C") mlcolor("#5B6E8C") title("`b' - raw scale") graphregion(color(white)) name(matrix_`b'_raw, replace)
graph matrix log_`b'0 log_`b'1 log_`b'2, half msymbol(O) mfcolor("#5B6E8C") mlcolor("#5B6E8C") title("`b' - log10") graphregion(color(white)) name(matrix_`b'_log, replace)
graph combine matrix_`b'_raw matrix_`b'_log, cols(2) graphregion(color(white)) name(matrix_`b', replace)
graph export "matrix_`b'.png", name(matrix_`b') width(3000) replace
restore
}

log close
