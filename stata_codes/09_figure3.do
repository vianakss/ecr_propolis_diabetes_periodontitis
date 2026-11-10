
********************************************************************************
*** FIGURE 2 - FOREST PLOT OF THE 90-DAY TREATMENT EFFECTS ***
********************************************************************************
capture log close
log using "log_figure3_forest_plot.log", text replace

* 1. Treatment effects at 90 days from the main analysis (panel 6 on the log10 scale)
clear
input byte panel byte row str60 label double(est lb ub) byte sig
1 1 "Mean PD (primary outcome)" 0.1052611 -0.1236268 0.334149 0
1 2 "Mean CAL" 0.0444731 -0.2822213 0.3711675 0
2 1 "Plaque index" -11.35407 -23.29056 0.5824269 0
2 2 "Bleeding on probing" -4.271435 -11.9392 3.396327 0
2 3 "Sites with PD 4 mm" -2.6797 -5.350694 -0.0087055 0
2 4 "Sites with PD 5–6 mm" 1.602913 -0.6193285 3.825154 0
2 5 "Sites with PD ≥7 mm" 0.1015435 -0.5717529 0.77484 0
3 1 "Glycated haemoglobin" 0.4715974 -0.0199981 0.9631929 0
4 1 "Fasting plasma glucose" 10.1623 -19.54698 39.87159 0
5 1 "Salivary interleukin-6*†" 6.981553 2.446481 11.51662 1
6 1 "C-reactive protein" -0.2727218 -0.5994444 0.0540008 0
6 2 "Salivary tumour necrosis factor-α" 0.1504902 -0.1531642 0.4541446 0
6 3 "Total bacterial load" -0.1700173 -0.9070056 0.5669709 0
6 4 "{it:Porphyromonas gingivalis}" 0.1875253 -0.6289928 1.004043 0
6 5 "{it:Fusobacterium nucleatum}*†" 1.424093 0.4372221 2.410965 1
6 6 "{it:Streptococcus oralis}*" -1.359436 -2.13245 -0.5864218 1
end

* 2. Back-transformation of panel 6 to ratios of geometric means
foreach v in est lb ub {
replace `v' = 10^`v' if panel==6
}

* 3. Panel-specific titles, axis titles and heights
local t1 "A. Periodontal outcomes (mm)"
local t2 "B. Periodontal outcomes (%)"
local t3 "C. HbA1c (%)"
local t4 "D. Fasting plasma glucose (mg/dL)"
local t5 "E. Salivary interleukin-6 (pg/mL)"
local t6 "F. Outcomes analysed on the log10 scale"
local x1 "Difference (Test − Control), mm"
local x2 "Difference (Test − Control), percentage points"
local x3 "Difference (Test − Control), %"
local x4 "Difference (Test − Control), mg/dL"
local x5 "Difference (Test − Control), pg/mL"
local x6 "Ratio of geometric means (Test/Control), log scale"
local h1 = 26
local h2 = 52
local h3 = 17
local h4 = 17
local h5 = 17
local h6 = 58

* 4. One graph per panel (the null value is always included in the horizontal axis)
forvalues p = 1/6 {
local yl ""
quietly levelsof row if panel==`p', local(rows)
foreach r of local rows {
quietly levelsof label if panel==`p' & row==`r', local(lab) clean
local yl `"`yl' `r' "`lab'""'
}
quietly summarize row if panel==`p', meanonly
local nr = r(max)
local null = cond(`p'==6, 1, 0)
local xopt ""
if `p'==6 local xopt "xscale(log range(0.005 500)) xlabel(0.01 "0.01" 0.1 "0.1" 1 "1" 10 "10" 100 "100", labsize(small))"
if `p'!=6 local xopt "xscale(range(`null')) xlabel(#5, labsize(small))"
twoway (rcap lb ub row if panel==`p', horizontal lcolor("#5B6E8C") lwidth(medthick)) (scatter row est if panel==`p' & sig==0, mcolor("#384358") msymbol(D) msize(medium)) (scatter row est if panel==`p' & sig==1, mcolor("#541A2E") msymbol(D) msize(medium)), xline(`null', lcolor(gs10) lpattern(dash)) ylabel(`yl', angle(0) noticks labsize(small) nogrid) yscale(reverse range(0.5 `=`nr'+0.5')) ytitle("") `xopt' xtitle("`x`p''", size(small)) title("`t`p''", size(medsmall) position(11) justification(left)) legend(off) graphregion(color(white)) plotregion(color(white)) fysize(`h`p'') name(fp`p', replace)
}

* 5. Combined figure: periodontal panels on the left, laboratory panels on the right
graph combine fp1 fp2, cols(1) imargin(small) graphregion(color(white)) name(fp_left, replace)
graph combine fp3 fp4 fp5 fp6, cols(1) imargin(small) graphregion(color(white)) name(fp_right, replace)
graph combine fp_left fp_right, cols(2) xsize(12) ysize(8) imargin(small) graphregion(color(white)) name(figure2_forest, replace)
graph export "figure2_forest_plot.png", name(figure2_forest) width(3600) replace
graph export "figure2_forest_plot.tif", name(figure2_forest) width(3600) replace
graph save figure2_forest "figure2_forest_plot.gph", replace

log close
