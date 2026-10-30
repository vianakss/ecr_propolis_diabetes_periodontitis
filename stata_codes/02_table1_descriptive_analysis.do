********************************************************************************
 *** Descriptive analysis - Table 1 ***
********************************************************************************
version 19
capture log close
log using "log_table1_sample_characteristics.log", text replace

preserve

keep if time == 0
isid id
assert _N == 24
assert inlist(group, 0, 1)

* 0. COMPLETENESS OF BASELINE CHARACTERISTICS
misstable summarize age sex education income diabetes_duration family_hx_diabetes bmi bmi_cat nteeth subging_instr_time

* 1. AGE
table group, statistic(count age) statistic(mean age) statistic(sd age) statistic(p25 age) statistic(median age) statistic(p75 age)
ttest age, by(group) unequal
ranksum age, by(group) exact

* 2. SEX
tabulate sex group, column exact

* 3. EDUCATION
tabulate education group, column exact

* 4. FAMILY INCOME
tabulate income group, column exact

* 5. DIABETES DURATION
table group, statistic(count diabetes_duration) statistic(mean diabetes_duration) statistic(sd diabetes_duration) statistic(p25 diabetes_duration) statistic(median diabetes_duration) statistic(p75 diabetes_duration)
ttest diabetes_duration, by(group) unequal
ranksum diabetes_duration, by(group) exact

* 6. FAMILY HISTORY OF DIABETES
tabulate family_hx_diabetes group, column exact

* 7. BODY MASS INDEX
table group, statistic(count bmi) statistic(mean bmi) statistic(sd bmi) statistic(p25 bmi) statistic(median bmi) statistic(p75 bmi)
ttest bmi, by(group) unequal
tabulate bmi_cat group, column exact
ranksum bmi, by(group) exact

* 8. NUMBER OF TEETH
table group, statistic(count nteeth) statistic(mean nteeth) statistic(sd nteeth) statistic(p25 nteeth) statistic(median nteeth) statistic(p75 nteeth)
ttest nteeth, by(group) unequal
ranksum nteeth, by(group) exact

* 9. SUBGINGIVAL INSTRUMENTATION TIME
table group, statistic(count subging_instr_time) statistic(mean subging_instr_time) statistic(sd subging_instr_time) statistic(p25 subging_instr_time) statistic(median subging_instr_time) statistic(p75 subging_instr_time)
ttest subging_instr_time, by(group) unequal
ranksum subging_instr_time, by(group) exact

restore
log close
