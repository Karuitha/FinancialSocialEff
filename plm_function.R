## Plm function ----
plm_function <- function(data, depvar,
        effect = c("individual", "time", "twoways", "nested"),
        model = c("within", "random", "ht", "between", "pooling", "fd"),
        inst.method = c("bvk", "baltagi", "am", "bms"),
        random.method = c("swar", "amemiya", "walhus", "nerlove"),
        digits = 4, index){
  
  ## Load libraries
  library(plm)
  library(broom)
  library(pcse)
  library(car)
  library(zoo)
  library(lmtest)
  library(broom)
  
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



vars <- data.frame(amelia$finefftheta_hat, 
  amelia$se_no_breadththeta_hat, 
  amelia$se_with_breadththeta_hat, 
  amelia$finsoc_no_breadththeta_hat, 
  amelia$finsoc_with_breadththeta_hat)


plm_function(data = amelia, depvar = amelia$se_with_breadththeta_hat, 
             index = c("mfiid", "year"), model = "within")


sapply(amelia[,65:69], plm_function, data = amelia, 
       index = c("mfiid", "year"), model = "random")




### Old Plm function 
modelling <- function(data, indepvar, effect = "twoways", model = "within"){
  library(plm)
  library(broom)
  library(pcse)
  library(car)
  library(zoo)
  library(lmtest)
  library(broom)
  
  unadjusted <- plm(indepvar ~ currentlegalstatus + age + region +
                      operating_expense_._assets + debt_to_equity_ratio + 
                      debt_to_equity_ratio ^ 2 + donations_assets_ratio + 
                      asset_structure + kkm + 
                      education, data = data, effect = "twoways", 
                    model = model, index = c("mfiid", "year"))
  
  list(broom::glance(unadjusted),
       
       # Correct standard errors for heteroscedasticity 
       # and cross-sectional dependence
       coeftest(unadjusted, vcov. = function(x) vcovBK(x, method = "arellano", 
                                                       type="HC1", cluster = "group")))
}

