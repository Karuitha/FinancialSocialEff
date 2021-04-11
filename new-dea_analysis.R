# DEA analysis afresh ----
## Load the data ----
setwd("C://Users//John Karuitha//OneDrive - University of Witwatersrand//Documents//My Thesis//Karuitha and Ojah Data//THESIS//Dissertation//Objective 1 Binary//Chapters in Progress//FinancialSocialEff")
load("C://Users//John Karuitha//OneDrive - University of Witwatersrand//Documents//My Thesis//Karuitha and Ojah Data//THESIS//Dissertation//Objective 1 Binary//Chapters in Progress//FinancialSocialEff//new_dea_data.RData")
## Load required packages ----
library(pacman)
library(tidyverse)
library(ggthemes)
library(gridExtra)
library(GGally)
library(VGAM)
library(plm)
library(rDEA)
library(censReg)
library(stargazer)
library(caret)
library(car)
library(bbplot)
library(skimr)
library(broom)
library(corrplot)
library(stargazer)
library(lmtest)
library(bbplot)
library(robustHD)
library(gplots)
library(pcse)
library(car)
library(gt)
library(GGally)
library(tidyquant)
library(plotly)
library(gghalves)
library(gganimate)
library(lubridate)
library(gifski)
library(av)
library(psych)

###########################################################################
# load the data ----
amelia <- read.csv("amelia.csv") ## full dataset

# Convert columns to factors ----
amelia$currentlegalstatus <- factor(amelia$currentlegalstatus, 
                                    levels = c("NGO", "Bank", "NBFI", 
                                    "Credit Union/ Cooperative", "Rural Bank"))

amelia$age <- factor(amelia$age, levels = c("New", "Young", "Mature"))

## Rename MENA to North Africa region and convert to factor
amelia$region <- ifelse(amelia$region == "Middle East and North Africa", "North_Africa", 
                        amelia$region)

amelia$region <- factor(amelia$region, levels = 
                          c("North_Africa", "Africa"))

## lag some variables ----
amelia <- amelia %>% 
  mutate(capital_.asset_ratio_lag = dplyr::lag(capital_.asset_ratio, 2), 
         donations_assets_ratio_lag = dplyr::lag(donations_assets_ratio, 2),
         debt_to_equity_ratio_lag = dplyr::lag(debt_to_equity_ratio, 2))

## Subset data for DEA analysis variables only plus identifiers ----
dea <- amelia %>% select(mfiid, year, liabilities_and_equity,
                          operating_expense_._assets, 
                          percent_of_female_borrowers, 
                          average_loan_balance_per_borrower,
                          gross_loan_portfolio_to_total_assets,
                          operational_self_sufficiency)

## Omit missing data for DEA analysis----
dea <- dea %>% na.omit()

# Make required functions ----
## Remove negatives and zeros function ---
addTwo <- function(x){x+3}

## The DEA function ----
###This function takes matrices of inputs and outputs and 
### Returns DEA bootstrap scores at 95% Confidence Interval
envelope <- function(input, output){
                      library(rDEA)
                      dea.robust (input, 
                      output, W = NULL, model = "output", 
                      RTS = "variable", B = 1000, 
                      alpha = 0.05, bw = "bw.ucv", 
                      bw_mult = 1)}

## Plotting function----
##Use paste0 as suggested in stack overflow 
median_n <- function(x){median(x, na.rm = TRUE)}

plotter <- function(data, x , y, z, xlabel, ylabel, title){
  library(tidyverse)
  library(ggthemes)
  library(gghalves)
  ggplot(data = data, mapping = aes(x = reorder({{x}}, {{y}}, median_n), 
  y = {{y}}, fill = {{z}})) + 
  geom_half_violin(draw_quantiles = c(0.25, 0.5, 0.75), 
  scale = "count") + 
  scale_y_log10() + labs(y = ylabel, x = xlabel, 
  title = title, 
  caption =
  paste0("Source: Authors' construction from the MIX data\n",
        "*The red point is the mean. Horizontal lines show the ",
        "first, second, and third quartiles\n",
        "*The size of the violins is proportional ",
        "to the data points")) + 
  tidyquant::theme_tq() + 
  theme(legend.position = "none") + 
  stat_summary(fun = mean, geom = "point", 
                size = 1, color = "red")}

## Hausmann Test function ----
hausmann_test <- function(data, depvar){
  
  library(plm)
  library(broom)
  library(pcse)
  library(car)
  library(zoo)
  library(lmtest)
  library(broom)
  
  fixed <- plm(depvar ~ currentlegalstatus + age + 
                 operating_expense_._assets + 
                 debt_to_equity_ratio + 
                 donations_assets_ratio + 
                 capital_.asset_ratio +
                 asset_structure + kkm + education + year, 
               data = data, effect = "individual", model = "within", 
      index = c("mfiid", "year")) 
  
  random <- plm(depvar ~ currentlegalstatus + age + 
          operating_expense_._assets + debt_to_equity_ratio + 
          donations_assets_ratio + capital_.asset_ratio +
          asset_structure + kkm + education + year, 
          data = data, effect = "individual", model = "random", 
      index = c("mfiid", "year"))
 
  broom::tidy(phtest(fixed, random))
}

## Regression function - fixed, random and pooling ---- 
modelling <- function(data, depvar,
                      effect = c("individual", "twoways", "time", "nested"),
                      model = c("within", "random", "ht", "between", "pooling", "fd"),
                      inst.method = c("bvk", "baltagi", "am", "bms"),
                      random.method = c("swar", "amemiya", "walhus", "nerlove"),
                      digits = 4, index = c("mfiid", "year")){
  
  ## Load libraries
  library(plm)
  library(broom)
  library(pcse)
  library(car)
  library(zoo)
  library(lmtest)
  library(broom)
  library(stargazer)
  
  ## match arguments
  effect <- match.arg(effect)
  model <- match.arg(model)
  inst.method <- match.arg(inst.method)
  random.method <- match.arg(random.method)
  
  ## Run the Model
  unadjusted <- plm(depvar ~ currentlegalstatus + 
                      age + region + 
                      operating_expense_._assets +
                      donations_assets_ratio + 
                      debt_to_equity_ratio + 
                      capital_.asset_ratio +
                      asset_structure + kkm + 
                      education + year, 
                    data = data, effect = effect,
                    model = model, inst.method = inst.method, index = index, 
                    digits = digits)
  
  ## Output results
  list(broom::glance(unadjusted),
       
       ## Correct standard errors for heteroscedasticity 
       ## And cross-sectional dependence
       
       coeftest(unadjusted, vcov. = function(x) {
         vcovBK(x, method = "arellano", type="HC1", cluster = "group")
       }))
}

## Winsorizing data function ----
winsorizeR <- function(x){
           library(robustHD)
            winsorize(x, minval = NULL, 
            maxval = NULL, probs = c(0.10, 0.90),
            na.rm = TRUE, type = 7)}

## QQlines function ----
qqliner <- function(model, main){
qqnorm(residuals(model), ylab = 'Residuals', col = "red", 
       main = main)

qqline(residuals(model))}

#################################################################################
# Clean data for DEA analysis ----
## Removing zeros and negatives
head(dea)
dea %>% skim()

## Add 3 to each variable to eliminate zeros and negatives then log to scale 

dea$liabilities_and_equity <- addTwo(dea$liabilities_and_equity) %>% log()
dea$operating_expense_._assets <- addTwo(dea$operating_expense_._assets) %>% log()
dea$percent_of_female_borrowers <- addTwo(dea$percent_of_female_borrowers) %>% log()
dea$average_loan_balance_per_borrower <- addTwo(dea$average_loan_balance_per_borrower) 
dea$average_loan_balance_per_borrower <- 1/dea$average_loan_balance_per_borrower
dea$gross_loan_portfolio_to_total_assets <- addTwo(dea$gross_loan_portfolio_to_total_assets) %>% log()
dea$operational_self_sufficiency <- addTwo(dea$operational_self_sufficiency) %>% log()
####################################################################################
####################################################################################
# View transformed data
dea %>% skim()
##################################################################################
## Running the DEA models ----
## Run sample DEA
### Social efficiency no breadth 
#se_no_breadth <- envelope(dea[,3:4], dea[,5:6])

###Social efficiency with breadth
#se_with_breadth <- envelope(dea[,3:4], dea[,5:7])

###Financial efficiency 
#fineff <- envelope(dea[,3:4], dea[,8])

### Social and financial efficiency no breadth
#finsoc_no_breadth <- envelope(dea[,3:4], dea[,c(5,6,8)])

### Social and financial efficiency with breadth
#finsoc_with_breadth <- envelope(dea[,3:4], dea[,c(5,6,7,8)])

## join efficiency scores to the main dataset dea----
dea <- cbind(dea$mfiid, dea$year, fineff$theta_hat, finsoc_no_breadth$theta_hat, 
             finsoc_with_breadth$theta_hat, se_no_breadth$theta_hat, 
             se_with_breadth$theta_hat) %>% data.frame()

names(dea) <- c("mfiid", "year", "finefftheta_hat", "finsoc_no_breadththeta_hat", 
                "finsoc_with_breadththeta_hat", "se_no_breadththeta_hat", 
                "se_with_breadththeta_hat")

## Left join DEA with original amelia data----
amelia <- left_join(amelia, dea, by = c("mfiid", "year"))
names(amelia) <- str_replace_all(names(amelia), "\\$", "")

## Generate data for 3 years or more and for 5 years or more of data ----
amelia3 <- amelia %>% dplyr::group_by(mfiid) %>% dplyr::filter(n() >= 3) %>% dplyr::ungroup()
glimpse(amelia3)
amelia5 <- amelia %>% dplyr::group_by(mfiid) %>% dplyr::filter(n() >= 5) %>% dplyr::ungroup()

## Plot raw data ----
amelia %>% plotter(x = currentlegalstatus, y = operational_self_sufficiency,
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "Operational Self-Sufficiency", 
                   title = "Operational Self-Sufficiency by MFI Legal Status")


amelia %>% plotter(x = currentlegalstatus, y = percent_of_female_borrowers,
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "% of Female Borrowers", 
                   title = "Percent of Female Borrowers by MFI Legal Status")

amelia %>% plotter(x = currentlegalstatus, y = average_loan_balance_per_borrower,
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "Average Loan Balance/ Borrower", 
                   title = "Average Loan Balance per Borrower by MFI Legal Status")

amelia %>% plotter(x = currentlegalstatus, y = gross_loan_portfolio_to_total_assets,
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "Gross Loan Portfolio to Total Assets", 
                   title = "Gross Loan Portfolio to Total Assets by MFI Legal Status")

## Plot efficiency scores ----
amelia %>% plotter(x = currentlegalstatus , y = finefftheta_hat, 
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "Financial Efficiency", 
                   title = "Financial Efficiency of MFIs by Legal Status")

amelia %>% plotter(x = currentlegalstatus , y = finsoc_no_breadththeta_hat, 
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "Financial and Social Efficiency", 
                   title = "Financial and Social Efficiency (No Breadth) of MFIs by Legal Status")

amelia %>% plotter(x = currentlegalstatus , y = finsoc_with_breadththeta_hat, 
                   z = currentlegalstatus, xlabel = "", 
                   ylabel = "Financial and Social Efficiency", 
                   title = "Financial and Social Efficiency (With Breadth) of MFIs by Legal Status")

amelia %>% plotter(x = currentlegalstatus , y = se_no_breadththeta_hat, 
                   z = currentlegalstatus, xlabel = "", ylabel = "Social Efficiency", 
                   title = "Social Efficiency (No Breadth) of MFIs by Legal Status")

amelia %>% plotter(x = currentlegalstatus ,y = se_with_breadththeta_hat, 
                   z = currentlegalstatus, xlabel = "", ylabel = "Social Efficiency", 
                   title = "Social Efficiency (With Breadth) of MFIs by Legal Status")

## Do the Hausmann tests ----
hausmann_test_results <- sapply(amelia[, c("finefftheta_hat", "finsoc_no_breadththeta_hat", 
                                           "finsoc_with_breadththeta_hat", "se_no_breadththeta_hat", 
                                           "se_with_breadththeta_hat")], 
                                hausmann_test, data = amelia)

out <- hausmann_test_results

out

## Running the models ----
## Try sapply to loop over all variables at once
## Fixed effects models ----
## regression using the full dataset 
my_models_fixed <- sapply(amelia[, c("finefftheta_hat", 
                                     "finsoc_no_breadththeta_hat", 
                                     "finsoc_with_breadththeta_hat", 
                                     "se_no_breadththeta_hat", 
                                     "se_with_breadththeta_hat")], 
                                      modelling, data = amelia)

write_csv(cbind(my_models_fixed[[2]], my_models_fixed[[4]], my_models_fixed[[6]], 
     my_models_fixed[[8]], my_models_fixed[[10]]) %>% 
       data.frame(), "fixed_full.csv")

write_csv(rbind(my_models_fixed[[1]], my_models_fixed[[3]], my_models_fixed[[5]], 
      my_models_fixed[[7]], my_models_fixed[[9]]) %>% t() %>% 
        data.frame(), "fixed_full_one.csv")

##regression using MFIs with 3 or more years of data 
my_models_fixed3 <- sapply(amelia3[, c("finefftheta_hat", 
                                       "finsoc_no_breadththeta_hat", 
                                       "finsoc_with_breadththeta_hat", 
                                       "se_no_breadththeta_hat", 
                                       "se_with_breadththeta_hat")], 
                                        modelling, data = amelia3)

write_csv(cbind(my_models_fixed3[[2]], my_models_fixed3[[4]], my_models_fixed3[[6]], 
      my_models_fixed3[[8]], my_models_fixed3[[10]]) %>% data.frame(), "fixed31.csv")

write_csv(rbind(my_models_fixed3[[1]], my_models_fixed3[[3]], my_models_fixed3[[5]], 
      my_models_fixed3[[7]], my_models_fixed3[[9]]) %>% t() %>% 
        data.frame(), "fixed32.csv")

##regression using MFIs with 5 or more years of data
my_models_fixed5 <- sapply(amelia5[, c("finefftheta_hat", 
                                       "finsoc_no_breadththeta_hat", 
                                      "finsoc_with_breadththeta_hat", 
                                      "se_no_breadththeta_hat", 
                                      "se_with_breadththeta_hat")], 
                           modelling, data = amelia5)

write_csv(cbind(my_models_fixed5[[2]], my_models_fixed5[[4]], my_models_fixed5[[6]], 
      my_models_fixed5[[8]], my_models_fixed5[[10]]) %>% data.frame(), "fixed51.csv")

write_csv(rbind(my_models_fixed5[[1]], my_models_fixed5[[3]], my_models_fixed5[[5]], 
      my_models_fixed5[[7]], my_models_fixed5[[9]]) %>% t() %>% 
        as.data.frame(), "fixed52.csv")


## Random effects model ----
my_models_random <- sapply(amelia[, c("finefftheta_hat", 
                                      "finsoc_no_breadththeta_hat", 
                                      "finsoc_with_breadththeta_hat", 
                                      "se_no_breadththeta_hat", 
                                      "se_with_breadththeta_hat")], 
                                       modelling, data = amelia, 
                                       model = "random")

write_csv(cbind(my_models_random[[2]], my_models_random[[4]], my_models_random[[6]], 
                my_models_random[[8]], my_models_random[[10]]) %>% 
            data.frame(), "random_full1.csv")

write_csv(rbind(my_models_random[[1]], my_models_random[[3]], my_models_random[[5]], 
                my_models_random[[7]], my_models_random[[9]]) %>% t() %>% 
            as.data.frame(), "random_full2.csv")


my_models_random3 <- sapply(amelia3[, c("finefftheta_hat", "finsoc_no_breadththeta_hat", 
                                      "finsoc_with_breadththeta_hat", "se_no_breadththeta_hat", 
                                      "se_with_breadththeta_hat")], modelling, data = amelia3, 
                                      model = "random")

write_csv(cbind(my_models_random3[[2]], my_models_random3[[4]], my_models_random3[[6]], 
                my_models_random3[[8]], my_models_random3[[10]]) %>% 
            data.frame(), "random31.csv")

write_csv(rbind(my_models_random3[[1]], my_models_random3[[3]], my_models_random3[[5]], 
                my_models_random3[[7]], my_models_random3[[9]]) %>% 
            t() %>% as.data.frame(), "random32.csv")


my_models_random5 <- sapply(amelia5[, c("finefftheta_hat", 
                                        "finsoc_no_breadththeta_hat", 
                                      "finsoc_with_breadththeta_hat", 
                                      "se_no_breadththeta_hat", 
                                      "se_with_breadththeta_hat")], 
                            modelling, data = amelia5, 
                           model = "random")

write_csv(cbind(my_models_random5[[2]], my_models_random5[[4]], 
                my_models_random5[[6]], 
                my_models_random5[[8]], my_models_random5[[10]]) %>% 
                data.frame(), "random51.csv")

write_csv(rbind(my_models_random5[[1]], my_models_random5[[3]], 
                my_models_random5[[5]], 
                my_models_random5[[7]], my_models_random5[[9]]) %>% t() %>% 
                as.data.frame(), "random52.csv")

## Output the models ----
random_output1 <- cbind(my_models_random[[2]], my_models_random[[4]], 
                        my_models_random[[6]], 
      my_models_random[[8]], my_models_random[[10]])

random_output1

write.csv(random_output1, "output1.csv")

random_output2 <- rbind(my_models_random[[1]], my_models_random[[3]], 
                        my_models_random[[5]], 
      my_models_random[[7]], my_models_random[[9]]) %>% t()

random_output2

write.csv(random_output2, "output2.csv")


ggpairs(amelia[, c("finefftheta_hat", "finsoc_no_breadththeta_hat", 
                   "finsoc_with_breadththeta_hat", "se_no_breadththeta_hat", 
                   "se_with_breadththeta_hat")], ggplot2::aes(color = amelia$currentlegalstatus),
        lower = list(continuous = wrap("points", alpha = 0.2, size=0.1)),
        upper = list(continuous= wrap("cor",params=c(size = 0.2))),
        diag = list(continuous = "barDiag"),
        title = "Correlations Between MFI Efficiency Scores (Grouped by MFI Legal Status)",
        axisLabels = "show", 
        columnLabels = c("Financial", "Financial/Social without Breadth", 
              "Financial/Social with Breadth", "Social without Breadth", 
              "Social with Breadth"), legend = c(1,1)) + theme(legend.position = "bottom") + 
  theme(legend.title=element_blank()) + 
  labs(caption = "Breadth refers to the extent of outreach measured using gross loans 
       to assets advanced by Microfinance Institutions")

psych::pairs.panels(amelia[, c("finefftheta_hat", "finsoc_no_breadththeta_hat", 
                               "finsoc_with_breadththeta_hat", "se_no_breadththeta_hat", 
                               "se_with_breadththeta_hat")], stars = TRUE, digits = 3, pch = ".", 
                    smooth = FALSE, 
                    main = "Correlation Matrix for MFI Efficiency Scores")
  

#list(continuous = "points", alpha = 0.3, size=0.1)
## Generate meanplots ----
## Financial Efficiency ----
amelia %>% group_by(year, currentlegalstatus) %>% 
  summarize(mean_fin_eff = mean(finefftheta_hat, na.rm = TRUE), 
  median_fin_eff = median(finefftheta_hat, na.rm = TRUE)) %>% 
  ggplot(aes(x = year, y = mean_fin_eff)) + 
  geom_label(label="Mean", 
            x=2018,
            y=0.44,
            label.padding = unit(0.35, "lines"), # Rectangle size around label
            label.size = 0.05,
            color = "red",
            fill="white") + 
  geom_line(col = "red") + geom_line(aes(x = year, 
  y = median_fin_eff), col = "blue") + 
  tidyquant::theme_tq() + facet_wrap(~currentlegalstatus) + 
  labs(x = "Year", y = "Financial Efficiency") + 
  geom_label(label="Median", 
             x=2018,
             y=0.43,
             label.padding = unit(0.35, "lines"), # Rectangle size around label
             label.size = 0.05,
             color = "blue",
             fill="white") + 
  geom_line(col = "red") + ggtitle("Trends in Financial Efficiency")

## Social Efficiency with Breadth ----
amelia %>% group_by(year, currentlegalstatus) %>% 
  summarize(mean_soc_eff = mean(se_with_breadththeta_hat, na.rm = TRUE), 
  median_soc_eff = median(se_with_breadththeta_hat, na.rm = TRUE)) %>% 
  ggplot(aes(x = year, y = mean_soc_eff)) + geom_line(col = "red") + 
  geom_line(aes(x = year, y = median_soc_eff), col = "blue") + 
  tidyquant::theme_tq() + facet_wrap(~currentlegalstatus) + 
  labs(x = "Year", y = "Social Efficiency", title = "Trends in Social Efficiency") + 
  geom_label(label="Median", 
             x=2018,
             y=0.89,
             label.padding = unit(0.35, "lines"), # Rectangle size around label
             label.size = 0.05,
             color = "blue",
             fill="white") + 
  geom_line(col = "red") + 
  geom_label(label="Mean", 
             x=2018,
             y=0.87,
             label.padding = unit(0.35, "lines"), # Rectangle size around label
             label.size = 0.05,
             color = "red",
             fill="white") 

## Financial and social efficiency ---- 
amelia %>% group_by(year, currentlegalstatus) %>% 
  summarize(mean_finsoc_eff = mean(finsoc_with_breadththeta_hat, na.rm = TRUE), 
  median_finsoc_eff = median(finsoc_with_breadththeta_hat, na.rm = TRUE)) %>% 
  ggplot(aes(x = year, y = mean_finsoc_eff)) + geom_line(col = "red") + 
  geom_line(aes(x = year, y = median_finsoc_eff), col = "blue") + 
  tidyquant::theme_tq() + facet_wrap(~currentlegalstatus) + 
  labs(x = "Year", y = "Financial and Social Efficiency")+ 
  geom_label(label="Median", 
             x=2018,
             y=0.89,
             label.padding = unit(0.35, "lines"), # Rectangle size around label
             label.size = 0.05,
             color = "blue",
             fill="white") + 
  geom_line(col = "red") + 
  geom_label(label="Mean", 
             x=2018,
             y=0.87,
             label.padding = unit(0.35, "lines"), # Rectangle size around label
             label.size = 0.05,
             color = "red",
             fill="white") + ggtitle("Trends in Financial and Social Efficiency")

## Do summary statistics of DEA inputs and Outputs ----
summary <- amelia %>% select(liabilities_and_equity,
                         operating_expense_._assets, 
                         percent_of_female_borrowers, 
                         average_loan_balance_per_borrower,
                         gross_loan_portfolio_to_total_assets,
                         operational_self_sufficiency) %>% skim()

summary %>% write_csv(., "summary_inputs.csv")

correlation  <- amelia %>% select(liabilities_and_equity,
                             operating_expense_._assets, 
                             percent_of_female_borrowers, 
                             average_loan_balance_per_borrower,
                             gross_loan_portfolio_to_total_assets,
                             operational_self_sufficiency) %>% na.omit() %>% cor()
correlation

input_outputs_summary <- cbind(summary, correlation) %>% 
  write.csv(., "input_output_summary.csv")

## Do summary of efficiency scores 
cbind(amelia %>% select(finefftheta_hat, se_no_breadththeta_hat, 
                  se_with_breadththeta_hat, finsoc_with_breadththeta_hat, 
                  finsoc_no_breadththeta_hat) %>% skim(), 

amelia %>% select(finefftheta_hat, se_no_breadththeta_hat, 
                  se_with_breadththeta_hat, finsoc_with_breadththeta_hat, 
                  finsoc_no_breadththeta_hat) %>% na.omit() %>% cor()) %>% 
  write_csv(., "summary_eff_scores.csv")

## Summarize dependent and independent variables ----
cbind(
amelia %>% select(finefftheta_hat, se_no_breadththeta_hat, 
  se_with_breadththeta_hat, finsoc_no_breadththeta_hat, 
  finsoc_with_breadththeta_hat,
  operating_expense_._assets, debt_to_equity_ratio, 
donations_assets_ratio, capital_.asset_ratio, 
asset_structure, kkm, education) %>% skim(),

amelia %>% select(finefftheta_hat, se_no_breadththeta_hat, 
    se_with_breadththeta_hat, finsoc_no_breadththeta_hat, 
    finsoc_with_breadththeta_hat,
    operating_expense_._assets, debt_to_equity_ratio, 
    donations_assets_ratio, capital_.asset_ratio, 
    asset_structure, kkm, education) %>% 
  na.omit() %>% cor()
) %>% write_csv(., "regression_variables_summary.csv")

## Plot the median social and financial efficiency scores ----
amelia %>% group_by(currentlegalstatus, age) %>% 
  summarize(social_eff = median(se_with_breadththeta_hat, na.rm = TRUE), 
            fin_eff = median(finefftheta_hat, na.rm = TRUE)) %>% 
  tibble() %>% ggplot(aes(x = social_eff, y = fin_eff, col = currentlegalstatus)) + 
  geom_label(aes(label = currentlegalstatus)) + 
  labs(x = "Median Social Efficiency", y = "Median Financial Efficiency", 
       title = "Financial Versus Social Efficiency of MFIs in Africa Faceted by Age") + 
  tidyquant::theme_tq() + theme(legend.position = "none") + facet_wrap(~age)

amelia %>% group_by(currentlegalstatus, outreach) %>% 
  summarize(social_eff = median(se_with_breadththeta_hat, na.rm = TRUE), 
            fin_eff = median(finefftheta_hat, na.rm = TRUE)) %>% 
  tibble() %>% ggplot(aes(x = social_eff, y = fin_eff, col = currentlegalstatus)) + 
  geom_label(aes(label = currentlegalstatus)) + 
  labs(x = "Median Social Efficiency", y = "Median Financial Efficiency", 
       title = "Financial Versus Social Efficiency of MFIs in Africa Faceted by Outreach") + 
  tidyquant::theme_tq() + theme(legend.position = "none") + facet_wrap(~ outreach)

amelia %>% group_by(currentlegalstatus, legal_tradition) %>% 
  summarize(social_eff = median(se_with_breadththeta_hat, na.rm = TRUE), 
            fin_eff = median(finefftheta_hat, na.rm = TRUE)) %>% 
  tibble() %>% ggplot(aes(x = social_eff, y = fin_eff, col = currentlegalstatus)) + 
  geom_label(aes(label = currentlegalstatus)) + 
  labs(x = "Median Social Efficiency", y = "Median Financial Efficiency", 
  title = "Financial Versus Social Efficiency of MFIs in Africa Faceted by Legal Tradition") + 
  tidyquant::theme_tq() + theme(legend.position = "none") + 
  facet_wrap(~ legal_tradition)

## data for rerunning stationarity tests in stata ----
amelia %>% select(mfiid, year,
liabilities_and_equity,
operating_expense_._assets, 
percent_of_female_borrowers, 
average_loan_balance_per_borrower,
gross_loan_portfolio_to_total_assets,
operational_self_sufficiency, lassets,
operating_expense_._assets, 
debt_to_equity_ratio, 
donations_assets_ratio, capital_.asset_ratio,
asset_structure, 
kkm, education,
finefftheta_hat, se_no_breadththeta_hat, 
se_with_breadththeta_hat, finsoc_with_breadththeta_hat, 
finsoc_no_breadththeta_hat) %>% write_csv(., "stationarity.csv")

## Winsorize data and run regressions -----
data_to_winsorize <- amelia %>% select(mfiid, year,
                  liabilities_and_equity,
                  operating_expense_._assets, 
                  percent_of_female_borrowers, 
                  average_loan_balance_per_borrower,
                  gross_loan_portfolio_to_total_assets,
                  operational_self_sufficiency, lassets,
                  operating_expense_._assets, 
                  debt_to_equity_ratio, 
                  donations_assets_ratio, 
                  capital_.asset_ratio, 
                  asset_structure, 
                  kkm, education, finefftheta_hat, se_no_breadththeta_hat, 
                  se_with_breadththeta_hat, finsoc_with_breadththeta_hat, 
                  finsoc_no_breadththeta_hat) %>% na.omit()

  
Winsorized_data <- sapply(data_to_winsorize, winsorizeR)

bio <- amelia %>% select(mfiid, year, currentlegalstatus, age, region)

Winsorized_data <- left_join(bio, Winsorized_data, 
                    by = c("mfiid", "year"), copy = TRUE)

head(Winsorized_data)

## Run hausmann test of data ----
hausmann_test_wins <- sapply(Winsorized_data[, c("finefftheta_hat", 
                                                 "finsoc_no_breadththeta_hat", 
                                                 "finsoc_with_breadththeta_hat", 
                                                 "se_no_breadththeta_hat", 
                                                 "se_with_breadththeta_hat")], 
                             hausmann_test, 
                             data = Winsorized_data)

#write_csv(hausmann_test_wins, "Hausmann_wins.csv")

## Run regression on winsorized datasets ----
my_models_random_wins <- sapply(Winsorized_data[, c("finefftheta_hat", 
                                            "finsoc_no_breadththeta_hat", 
                                            "finsoc_with_breadththeta_hat", 
                                            "se_no_breadththeta_hat", 
                                            "se_with_breadththeta_hat")], modelling, 
                                data = Winsorized_data, model = "random")


my_models_random_wins[[2]]

my_models_fixed_wins <- sapply(Winsorized_data[, c("finefftheta_hat", 
                                                   "finsoc_no_breadththeta_hat", 
                                                   "finsoc_with_breadththeta_hat", 
                                                   "se_no_breadththeta_hat", 
                                                   "se_with_breadththeta_hat")], modelling, 
                                data = Winsorized_data, model = "within")

my_models_fixed_wins[[10]]


## Output fixed effects models as excel/ csv ----
fixed_wins1 <- cbind(my_models_fixed_wins[[2]], my_models_fixed_wins[[4]], 
                        my_models_fixed_wins[[6]], my_models_fixed_wins[[8]], 
                        my_models_fixed_wins[[10]]) %>% data.frame()

fixed_wins1 %>% write.csv(., "fixed_wins1.csv")

fixed_wins2 <- rbind(my_models_fixed_wins[[1]], my_models_fixed_wins[[3]], 
                     my_models_fixed_wins[[5]], my_models_fixed_wins[[7]], 
                     my_models_fixed_wins[[9]]) %>% t() %>% data.frame()

fixed_wins2 %>% write.csv(., "fixed_wins2.csv")

## Output the random effects models as excel/csv ----
random_wins1 <- cbind(my_models_random_wins[[2]], my_models_random_wins[[4]], 
                     my_models_random_wins[[6]], my_models_random_wins[[8]], 
                     my_models_random_wins[[10]]) %>% data.frame()

random_wins1 %>% write.csv(., "random_wins1.csv")

random_wins2 <- rbind(my_models_random_wins[[1]], my_models_random_wins[[3]], 
                      my_models_random_wins[[5]], my_models_random_wins[[7]], 
                      my_models_random_wins[[9]]) %>% data.frame()

random_wins2 %>% write.csv(., "random_wins2.csv")

## Do a transtition plot ----
amelia %>% ggplot(aes(x = year, 
  y = se_with_breadththeta_hat, 
  col = currentlegalstatus)) + 
  geom_point(aes(size = gross_loan_portfolio)) + 
  theme(legend.position = "none") + scale_x_log10() + 
  scale_y_log10() + transition_time(year) + shadow_wake(0.2)

#####################################################################################
## Draw QQ plots ----
fixedqq <- plm(finefftheta_hat ~ currentlegalstatus + age + region +
      operating_expense_._assets + debt_to_equity_ratio + 
      donations_assets_ratio + capital_.asset_ratio + 
      asset_structure + kkm + 
      education + year, data = amelia, effect = "individual", 
    model = "within", index = c("mfiid", "year"))

randomqq <- plm(finefftheta_hat ~ currentlegalstatus + age + region +
      operating_expense_._assets + debt_to_equity_ratio + 
      donations_assets_ratio + capital_.asset_ratio +
      asset_structure + kkm + 
      education + year, data = amelia, effect = "individual", 
    model = "random", index = c("mfiid", "year"))

fixedqq1 <- plm(se_with_breadththeta_hat ~ currentlegalstatus + age + region +
      operating_expense_._assets + debt_to_equity_ratio + 
      donations_assets_ratio + capital_.asset_ratio + 
      asset_structure + kkm + 
      education + year, data = amelia, effect = "individual", 
    model = "within", index = c("mfiid", "year"))

randomqq1 <- plm(se_with_breadththeta_hat ~ currentlegalstatus + age + region +
      operating_expense_._assets + debt_to_equity_ratio + 
      donations_assets_ratio + capital_.asset_ratio + 
      asset_structure + kkm + education + 
      year, data = amelia, effect = "individual", 
    model = "random", index = c("mfiid", "year"))

summary(fixedqq)
summary(fixedqq1)
summary(randomqq)
summary(randomqq1)

par(mfrow = c(2,2))
qqliner(fixedqq, main = "Financial Efficiency- Fixed Effects Model")
qqliner(fixedqq1, main = "Social Efficiency - Fixed Effects Model")
qqliner(randomqq, main = "Financial Efficiency- Random Effects Model")
qqliner(randomqq1, main = "Social Efficiency - Random Effects Model")

##Save the session ----
save.image("C:/Users/John Karuitha/OneDrive - University of Witwatersrand/Documents/My Thesis/Karuitha and Ojah Data/THESIS/Dissertation/Objective 1 Binary/Chapters in Progress/FinancialSocialEff/new_dea_data.RData")

##############################################################################
linkTest<-function(model){
  e <- pryr::where(deparse(substitute(model)))
  cc<-class(model)[[1]]
  # Define variables for predicted values using names not likely to
  # conflict with other names in the environment of model
  .predicted <- predict(model)
  `.predicted^2`<-.predicted^2
  # Check to see that the predicted and predicted^2 variable actually
  # vary.
  if(round(var(.predicted), digits=2) == 0){
    stop("No parameters that vary. Cannot perform test.")
  } 
  # The variables need to be found in the environment of model for 
  # the update() call to work.
  assign('.predicted', .predicted, envir = e)
  assign('.predicted^2', `.predicted^2`, envir = e)
  model<-update(model, . ~ .predicted + `.predicted^2`)
  # Clean up
  rm(.predicted, envir = e)
  rm(`.predicted^2`, envir = e)
  model
}

linkTest(my_models_random)
