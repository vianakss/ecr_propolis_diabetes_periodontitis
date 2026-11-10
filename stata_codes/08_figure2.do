
********************************************************************************
*** FIGURE 2 - INDIVIDUAL TRAJECTORIES OF BACTERIAL COUNTS ***
********************************************************************************
capture log close
log using "log_figure2_bacterial_trajectories.log", text replace
use "dataset_-_ecr_propolis_long.dta", clear
xtset, clear

* 1. Participants with bacterial quantification
capture drop lab_bac
bysort id: egen lab_bac = max(time==0 & !missing(log_un))
keep if lab_bac
keep id group time log_un log_pg log_fn log_so

* 2. Participants without a valid 90-day sample and group means of observed values
bysort id: egen miss90 = max(time==2 & missing(log_un))
foreach y in log_un log_pg log_fn log_so {
bysort group time: egen m_`y' = mean(`y')
}
egen tag_gt = tag(group time)
sort id time

* 3. Label positions for participants without a valid 90-day sample (clock positions)
foreach y in log_un log_pg log_fn log_so {
gen pos_`y' = 3
}
replace pos_log_un = 2 if id=="B7":id
replace pos_log_un = 4 if id=="B3":id
replace pos_log_fn = 2 if id=="B3":id
replace pos_log_fn = 4 if id=="B7":id
replace pos_log_so = 2 if id=="B7":id
replace pos_log_so = 4 if id=="B4":id

* 4. One panel per bacterium and group
foreach y in log_un log_pg log_fn log_so {
foreach g in 0 1 {
local col = cond(`g'==0, "#384358", "#541A2E")
local grp = cond(`g'==0, "Control", "Test")
if "`y'"=="log_un" local sp "Total bacterial load"
if "`y'"=="log_pg" local sp "{it:Porphyromonas gingivalis}"
if "`y'"=="log_fn" local sp "{it:Fusobacterium nucleatum}"
if "`y'"=="log_so" local sp "{it:Streptococcus oralis}"
if "`y'"=="log_un" local yr "11(1)15"
if "`y'"=="log_pg" local yr "3(1)8"
if "`y'"=="log_fn" local yr "6(1)11"
if "`y'"=="log_so" local yr "3(1)8"
twoway (line `y' time if group==`g' & !miss90, connect(L) lcolor(gs11) lwidth(thin)) (scatter `y' time if group==`g' & !miss90, mcolor(gs11) msize(small)) (line `y' time if group==`g' & miss90, connect(L) lcolor("#A0607A") lpattern(dash) lwidth(medthin)) (scatter `y' time if group==`g' & miss90, mcolor("#A0607A") msymbol(Oh) msize(small)) (scatter `y' time if group==`g' & miss90 & time==1, msymbol(none) mlabel(id) mlabcolor("#A0607A") mlabsize(vsmall) mlabvposition(pos_`y')) (connected m_`y' time if group==`g' & tag_gt, sort lcolor("`col'") lwidth(thick) mcolor("`col'") msymbol(D) msize(medium)), xlabel(0 "Baseline" 1 "60 days" 2 "90 days", labsize(small)) xscale(range(-0.2 2.3)) ylabel(`yr', angle(0) labsize(small)) xtitle("") ytitle("log10 count", size(small)) title("`sp' – `grp'", size(medsmall)) legend(off) graphregion(color(white)) plotregion(color(white)) name(tr_`y'_`g', replace)
}
}

* 5. Combined figure: bacteria in rows, groups in columns
graph combine tr_log_un_0 tr_log_un_1 tr_log_pg_0 tr_log_pg_1 tr_log_fn_0 tr_log_fn_1 tr_log_so_0 tr_log_so_1, cols(2) xsize(10) ysize(14) imargin(small) graphregion(color(white)) name(figure3_trajectories, replace)
graph export "figure2_bacterial_trajectories.png", name(figure2_trajectories) width(3000) replace
graph export "figure2_bacterial_trajectories.tif", name(figure2_trajectories) width(3000) replace
graph save figure2_trajectories "figure2_bacterial_trajectories.gph", replace






log close

