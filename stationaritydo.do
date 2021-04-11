* Stationarity tests 
xtset mfiid year

* loop stationarity test 
foreach var of varlist liabilities_and_equity - finsoc_no_breadththeta_hat{
set more off
display "**********************************"
display `var'
display "**********************************"
xtunitroot fisher `var', dfuller trend lags(0)
}

* loop stationarity test 
foreach var of varlist gross_loan_portfolio_to_total_as - education{
set more off
display "**********************************"
display `var'
display "**********************************"
xtunitroot fisher `var', dfuller trend lags(0)
}




