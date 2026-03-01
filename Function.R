#' Ensure a Directory Exists
#'
#' Checks whether the directory at the given path exists; if not, creates it
#' recursively (including any missing parent directories).
#'
#' @param file_path Character scalar. Path to the target directory.
#'
#' @returns Invisibly returns `NULL`. Called for its side effect of ensuring the
#'   directory exists.
#' @export
#'
#' @examples
#' # Create a nested directory under a temporary folder
#' d <- file.path(tempdir(), "a", "b", "c")
#' ensure_dir(d)
#' dir.exists(d)
ensure_dir <- function(file_path){
  if(!dir.exists(file_path)) dir.create(file_path,recursive = T)
}


#' Format AUC with 95% CI from a pROC ROC Object
#'
#' Computes the area under the ROC curve (AUC) and its confidence interval using
#' `pROC::ci.auc()`, then returns a single formatted string like
#' `"0.823 (0.755-0.891)"`.
#'
#' @param roc_object A `pROC::roc` object (as returned by `pROC::roc()`).
#'
#' @returns A character scalar containing the formatted AUC and confidence
#'   interval: `"AUC (lower-upper)"`, rounded to three decimals.
#' @export
#'
#' @examples
#' if (requireNamespace("pROC", quietly = TRUE)) {
#'   roc_obj <- pROC::roc(response = c(0, 0, 1, 1), predictor = c(0.1, 0.4, 0.35, 0.8))
#'   print_auc_and_ci(roc_obj)
#' }
print_auc_and_ci <- function(roc_object){
  AUC_roc <- signif(pROC::ci.auc(roc_object)[2],2)
  ci_roc <- paste0(signif(pROC::ci.auc(roc_object)[1],2),"-",
                   signif(pROC::ci.auc(roc_object)[3],2))
  AUC_formatted <- str_c(AUC_roc," (", ci_roc, ")")
  return(AUC_formatted)
}

#' Calculate Discrimination and Calibration Metrics (AUC, Intercept, Slope)
#'
#' Computes model performance for binary outcomes using predicted probabilities:
#' discrimination (AUC with confidence interval) and calibration (intercept and
#' slope from a logistic regression of observed outcomes on predicted values).
#'
#' @param observed Vector of observed binary outcomes (e.g., 0/1 or a two-level
#'   factor).
#' @param predicted Numeric vector of predicted probabilities or risk scores
#'   (higher values indicate higher risk).
#'
#' @returns A named vector containing:
#' \describe{
#'   \item{auc}{AUC formatted to 3 decimals (character).}
#'   \item{auc_CI}{AUC confidence interval formatted as "lower-upper" (character).}
#'   \item{calibration intercept}{Calibration intercept estimate (numeric).}
#'   \item{calibration.intercept.SE}{Standard error of the intercept (numeric).}
#'   \item{calibration slope}{Calibration slope estimate (numeric).}
#'   \item{calibration.slope.SE}{Standard error of the slope (numeric).}
#' }
#' @export
#'
#' @examples
#' if (requireNamespace("pROC", quietly = TRUE)) {
#'   set.seed(1)
#'   observed <- rbinom(50, 1, 0.4)
#'   predicted <- runif(50)
#'   calculate_performance2(observed, predicted)
#' }
calculate_performance2 <- function(observed, predicted){
  ## AUC
  roc_object <- pROC::roc(observed~predicted)
  auc <- signif(pROC::ci.auc(roc_object)[2],2)
  ci_auc <- paste0(signif(pROC::ci.auc(roc_object)[1],2), "-",
                   signif(pROC::ci.auc(roc_object)[3],2))
  
  glm.fit <- summary(glm(observed~predicted, family = binomial))
  ## calibration.intercept
  calibration.intercept <- glm.fit$coefficients[1,"Estimate"]
  calibration.intercept.SE <- glm.fit$coefficients[1,"Std. Error"]
  
  ## calibration.slope
  calibration.slope <- glm.fit$coefficients[2,"Estimate"]
  calibration.slope.SE <- glm.fit$coefficients[2,"Std. Error"]
  
  vec <- c(auc, ci_auc, calibration.intercept, calibration.intercept.SE,
           calibration.slope, calibration.slope.SE)
  names(vec) <- c("auc", "auc_CI","calibration intercept", "calibration.intercept.SE",
                  "calibration slope","calibration.slope.SE")
  return(vec)
}


#' Calculate Diagnostic Performance Indices (with CIs) for Multiple Indicators
#'
#' Calculates common diagnostic accuracy indices for a set of binary indicators
#' against a binary reference outcome, including sensitivity, specificity, PPV,
#' NPV, positive/negative likelihood ratios, and accuracy. Confidence intervals
#' are computed using `bootLR::diagCI()` (and Wilson CI for accuracy), and results
#' are returned in a tidy data frame with formatted strings.
#'
#' @param data A data frame containing the indicators and the reference outcome.
#' @param indicators Character vector of column names in `data` representing
#'   diagnostic indicators (assumed coded as 0/1). Default:
#'   `c("P85","VCA_IgA","EA_IgA","NA1_IgA","EBVDNA")`.
#' @param outcome Character scalar. Column name of the binary reference outcome
#'   in `data` (e.g., 0/1). Default is `"group_ref"`.
#' @param calcLRCI_use Character scalar passed to `bootLR::diagCI(calcLRCI = ...)`
#'   to control the likelihood ratio CI method. Common options include
#'   `"analytic"` and `"BayesianLR.test"`.
#' @param binomMethod_use Character scalar passed to `bootLR::diagCI(binomMethod = ...)`
#'   for binomial proportion CIs. Options include `"wilson"`, `"exact"`, `"ac"`,
#'   `"asymptotic"`, `"prop.test"`, `"bayes"`, `"logit"`, `"cloglog"`, `"probit"`.
#'
#' @returns A data frame with one row per indicator and the following columns:
#' \describe{
#'   \item{Indicator}{Indicator name.}
#'   \item{TP, FP, FN, TN}{Confusion-matrix counts.}
#'   \item{Sensitivity, Specificity, PPV, NPV}{Percent estimates formatted as
#'     `"xx.x\n(lower to upper)"`.}
#'   \item{posLR, negLR}{Likelihood ratios formatted as `"x.xx\n(lower to upper)"`.}
#'   \item{accuracy}{Accuracy formatted as `"xx.x\n(lower to upper)"` (Wilson CI).}
#' }
#'
#' @export
#'
#' @examples
#' if (requireNamespace("bootLR", quietly = TRUE) &&
#'     requireNamespace("cutpointr", quietly = TRUE) &&
#'     requireNamespace("Hmisc", quietly = TRUE) &&
#'     requireNamespace("stringr", quietly = TRUE)) {
#'   dat <- data.frame(
#'     group_ref = c(1, 1, 0, 0, 1, 0),
#'     P85       = c(1, 0, 0, 1, 1, 0),
#'     EBVDNA    = c(1, 1, 0, 0, 0, 0)
#'   )
#'   calcuate_diagnosis_index2(dat, indicators = c("P85", "EBVDNA"), outcome = "group_ref")
#' }
calcuate_diagnosis_index2 <- function(data,
                                      indicators=c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA","EBVDNA"),
                                      outcome="group_ref",
                                      calcLRCI_use="analytic", #"BayesianLR.test" or "analytic"
                                      binomMethod_use = "wilson") # c("all","exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit")
{
  if(!requireNamespace("cutpointr"))install.packages("cutpointr")
  if(!requireNamespace("tidyverse"))install.packages("tidyverse")
  if(!requireNamespace("bootLR"))install.packages("bootLR")
  suppressWarnings(library(bootLR))
  suppressWarnings(library(cutpointr))
  suppressWarnings(library(tidyverse))
  
  result_df <- data.frame(
    Indicator = character(0),
    TP=numeric(0),
    FP=numeric(0),
    FN=numeric(0),
    TN=numeric(0),
    Sensitivity = character(0),
    Specificity = character(0),
    PPV = character(0),
    NPV = character(0),
    posLR = character(0),
    negLR = character(0),
    accuracy = character(0)
  )
  
  for (indicator in indicators) {
    
    pred <- data[[indicator]]
    
    valid_rows <- !is.na(pred)
    pred <- pred[valid_rows]
    truth <- data[[outcome]][valid_rows]
    
    data_withoutNA <- data.frame(predictiveIndex=pred,
                                 outcomeIndex=truth)
    
    TP = sum(pred[truth==1]==1)
    FP = sum(pred[truth==0]==1)
    FN = sum(pred[truth==1]==0)
    TN = sum(pred[truth==0]==0)
    
    diagCI_result <- bootLR::diagCI(truePos = TP,
                                    totalDzPos = (TP+FN),
                                    trueNeg = TN,
                                    totalDzNeg = (TN+FP),
                                    calcLRCI = calcLRCI_use,
                                    binomMethod = binomMethod_use
    ) %>% as.data.frame()
    diagCI_result$Indicator <- indicator
    
    diagCI_result_tidy <- diagCI_result %>%
      dplyr::mutate(TP=TP,FP=FP,FN=FN,TN=TN) %>%
      dplyr::select(-c(truePos:totalDzNeg)) %>%
      dplyr::select(Indicator,TP:TN,everything()) %>%
      dplyr::rename("Sensitivity"="sens",
                    "Specificity"="spec") %>%
      dplyr::mutate(across(c(Sensitivity:NPV_UB),~sprintf("%.1f", .x*100))) %>%
      dplyr::mutate(across(c(posLR:negLR_UB),~sprintf("%.2f", .x))) %>%
      dplyr::mutate(
        Sensitivity = str_c(Sensitivity,
                            "\n(",sens_LB, " to ", sens_UB, ")"),
        Specificity = str_c(Specificity,
                            "\n(",spec_LB, " to ", spec_UB, ")"),
        PPV = str_c(PPV,
                    "\n(",PPV_LB, " to ", PPV_UB, ")"),
        NPV = str_c(NPV,
                    "\n(",NPV_LB, " to ", NPV_UB, ")"),
        posLR = str_c(posLR,
                      "\n(",posLR_LB, " to ", posLR_UB, ")"),
        negLR = str_c(negLR,
                      "\n(",negLR_LB, " to ", negLR_UB, ")"),
      ) %>%
      dplyr::select(Indicator,
                    TP:TN,
                    Sensitivity,
                    Specificity,
                    PPV,
                    NPV,
                    posLR,
                    negLR)
    
    #### accuracy
    accuracy_raw <- cutpointr::accuracy(tp=diagCI_result_tidy$TP,fp = diagCI_result_tidy$FP,
                                        tn = diagCI_result_tidy$TN,fn = diagCI_result_tidy$FN)
    accuracy <- sprintf("%.1f", accuracy_raw*100)
    wilson_accuracy <- Hmisc::binconf(x = (diagCI_result_tidy$TP+diagCI_result_tidy$TN),
                                      n = (diagCI_result_tidy$TP+diagCI_result_tidy$TN+
                                             diagCI_result_tidy$FP+diagCI_result_tidy$FN),
                                      method = "wilson")
    accuracy_upper <- sprintf("%.1f",wilson_accuracy[3]*100)
    accuracy_lower <- sprintf("%.1f",wilson_accuracy[2]*100)
    
    diagCI_result_tidy$accuracy <- str_c(accuracy,
                                         "\n(",accuracy_lower, " to ", accuracy_upper, ")")
    
    result_df <- rbind(result_df, diagCI_result_tidy)
  }
  
  colnames(result_df) <- c("Indicator",
                           "TP","FP","FN","TN",
                           "Sensitivity",
                           "Specificity",
                           "PPV",
                           "NPV",
                           "posLR",
                           "negLR",
                           "accuracy")
  
  return(result_df)
}

#' Paired Comparison of Diagnostic Tests (Sensitivity/Specificity/PPV/NPV) with Raw P-values
#'
#' Compares one binary diagnostic indicator against multiple other indicators
#' using paired (within-subject) test comparison methods. Sensitivity and
#' specificity are compared using McNemar's test when the number of discordant
#' pairs is sufficiently large; otherwise a McNemar mid-P test is used. PPV and
#' NPV are compared using paired predictive value methods from `DTComPair`.
#' The function returns both raw p-values and formatted p-value strings.
#'
#' @param data A data frame containing the reference outcome and binary indicator
#'   variables (coded as 0/1; factors are converted to numeric).
#' @param indicator1 Character scalar. Column name of the primary indicator to be
#'   compared (default `"P85"`).
#' @param compare_indicator Character vector. Column names of indicators to be
#'   compared against `indicator1`.
#' @param outcome Character scalar. Column name of the binary reference outcome
#'   (default `"group_ref"`).
#'
#' @returns A data frame with one row per comparison (`indicator1` vs each
#'   `compare_indicator`), including sensitivity/specificity/PPV/NPV for both
#'   tests and their paired-comparison p-values (raw and formatted).
#' @export
#'
#' @examples
#' if (requireNamespace("DTComPair", quietly = TRUE) &&
#'     requireNamespace("dplyr", quietly = TRUE) &&
#'     requireNamespace("stringr", quietly = TRUE)) {
#'   dat <- data.frame(
#'     group_ref = c(1,1,0,0,1,0,1,0),
#'     P85       = c(1,0,0,1,1,0,1,0),
#'     VCA_IgA   = c(1,1,0,0,1,0,0,0),
#'     EA_IgA    = c(0,0,0,1,1,0,1,0)
#'   )
#'   DTC_pair(dat, indicator1 = "P85",
#'            compare_indicator = c("VCA_IgA", "EA_IgA"),
#'            outcome = "group_ref")
#' }
DTC_pair <- function(data,
                     indicator1="P85",
                     compare_indicator=c("VCA_IgA", "EA_IgA", "NA1_IgA"),
                     outcome="group_ref") {
  #### Ensure required packages are available ####
  if(!requireNamespace("tidyverse")) install.packages("tidyverse")
  if(!requireNamespace("DTComPair")) install.packages("DTComPair")
  if(!requireNamespace("contingencytables")) remotes::install_github("ocbe-uio/contingencytables", ref="develop")
  library(DTComPair)
  library(tidyverse)
  library(contingencytables)
  
  #### Initialize a results data frame ####
  result_compare <- data.frame(
    compare = character(0),
    sen1=character(0),
    sen2=character(0),
    sen_compM=character(0),
    sen_method=character(0),
    sen.p.raw=numeric(),
    sen.p.value=character(0),
    spe1=character(0),
    spe2=character(0),
    spe_compM=character(0),
    spe_method=character(0),
    spe.p.raw=numeric(),
    spe.p.value=character(0),
    ppv1=character(0),
    ppv2=character(0),
    ppv.p.value=character(0),
    npv1=character(0),
    npv2=character(0),
    npv.p.value=character(0)
  )
  
  for (i in 1:length(compare_indicator)) {
    indicator2 <- compare_indicator[i]
    data_sub <- data %>% dplyr::select(!!sym(outcome),!!sym(indicator1),!!sym(indicator2)) %>%
      dplyr::mutate(across(where(is.factor), ~ as.numeric(as.character(.x))))
    
    colnames(data_sub) <- c("group","indicator1","indicator2")
    
    result_pair <- DTComPair::tab.paired(d=group, y1=indicator1, y2=indicator2, data = data_sub)
    accuray_compare <- DTComPair::sesp.mcnemar(result_pair)
    PV_compare <- DTComPair::pv.wgs(result_pair)
    
    #### Sensitivity ----
    matrix_sens <- result_pair$diseased[1:2, 1:2]
    
    ## 2x2 matrix elements for sensitivity (row-major) stored as "a,b,c,d"
    sen_compM <- str_c(matrix_sens[1,1], matrix_sens[1,2],
                       matrix_sens[2,1], matrix_sens[2,2], sep = ",")
    
    sens_b <- matrix_sens[1,2]
    sens_c <- matrix_sens[2,1]
    
    if ((sens_b + sens_c) >= 25) {
      sen.p.raw <- mcnemar.test(
        matrix(c(matrix_sens[1,1], matrix_sens[1,2],
                 matrix_sens[2,1], matrix_sens[2,2]),
               nrow = 2, byrow = TRUE),
        correct = FALSE
      )$p.value
      sen.p.value <- signif(accuray_compare$sensitivity[5], 2)
      sen_method_used <- "McNemar"
    } else {
      midp_value_sens <- contingencytables::McNemar_midP_test_paired_2x2(matrix_sens)
      sen.p.raw <- midp_value_sens$midP
      sen.p.value <- signif(midp_value_sens$midP, 2)
      sen_method_used <- "McNemar mid-P"
    }
    
    #### Specificity ----
    matrix_spec <- result_pair$non.diseased[1:2, 1:2]
    
    ## 2x2 matrix elements for specificity (row-major) stored as "a,b,c,d"
    spe_compM <- str_c(matrix_spec[1,1], matrix_spec[1,2],
                       matrix_spec[2,1], matrix_spec[2,2], sep = ",")
    
    spec_b <- matrix_spec[1,2]
    spec_c <- matrix_spec[2,1]
    
    if ((spec_b + spec_c) >= 25) {
      spe.p.raw <- mcnemar.test(
        matrix(c(matrix_spec[1,1], matrix_spec[1,2],
                 matrix_spec[2,1], matrix_spec[2,2]),
               nrow = 2, byrow = TRUE),
        correct = FALSE
      )$p.value
      spe.p.value <- signif(accuray_compare$specificity[5], 2)
      spe_method_used <- "McNemar"
    } else {
      midp_value_spe <- contingencytables::McNemar_midP_test_paired_2x2(matrix_spec)
      spe.p.raw <- midp_value_spe$midP
      spe.p.value <- signif(midp_value_spe$midP, 2)
      spe_method_used <- "McNemar mid-P"
    }
    
    ppv.p.value <- signif(PV_compare$ppv["p.value"], 2)
    npv.p.value <- signif(PV_compare$npv["p.value"], 2)
    
    result_compare1 <- data.frame(
      compare=stringr::str_c(indicator1, " vs ", indicator2),
      sen1=sprintf("%0.1f", accuray_compare$sensitivity["test1"]*100),
      sen2=sprintf("%0.1f", accuray_compare$sensitivity["test2"]*100),
      sen_compM=sen_compM,
      sen_method=sen_method_used,
      sen.p.raw=sen.p.raw,
      sen.p.value=ifelse(sen.p.value < 0.0001, "p<0.0001", str_c("p=", sen.p.value)),
      spe1=sprintf("%0.1f", accuray_compare$specificity["test1"]*100),
      spe2=sprintf("%0.1f", accuray_compare$specificity["test2"]*100),
      spe_compM=spe_compM,
      spe_method=spe_method_used,
      spe.p.raw=spe.p.raw,
      spe.p.value=ifelse(spe.p.value < 0.0001, "p<0.0001", str_c("p=", spe.p.value)),
      ppv1=sprintf("%0.1f", PV_compare$ppv["test1"]*100),
      ppv2=sprintf("%0.1f", PV_compare$ppv["test2"]*100),
      ppv.p.raw=PV_compare$ppv["p.value"],
      ppv.p.value=ifelse(ppv.p.value < 0.0001, "p<0.0001", str_c("p=", ppv.p.value)),
      npv1=sprintf("%0.1f", PV_compare$npv["test1"]*100),
      npv2=sprintf("%0.1f", PV_compare$npv["test2"]*100),
      npv.p.raw=PV_compare$npv["p.value"],
      npv.p.value=ifelse(npv.p.value < 0.0001, "p<0.0001", str_c("p=", npv.p.value))
    )
    result_compare <- rbind(result_compare, result_compare1)
  }
  
  return(result_compare)
}



#' Fit a GLM/GAM Model and Append Predictions to the Data
#'
#' Fits either a generalized linear model (GLM) or a generalized additive model
#' (GAM) to the provided data, then adds two new columns containing the linear
#' predictor (link scale) and the predicted values on the response scale
#' (e.g., predicted probabilities for binomial models).
#'
#' @param data A data frame used for model fitting and prediction.
#' @param formula A model formula used to fit the model.
#' @param model Character scalar specifying the model type: `"glm"` or `"gam"`.
#' @param pred_link_col Character scalar. Name of the column to store the linear
#'   predictor (link-scale predictions).
#' @param pred_prob_col Character scalar. Name of the column to store response-
#'   scale predictions (e.g., probabilities for binomial family).
#' @param family A `stats::family()` object passed to `glm()`/`mgcv::gam()`.
#'   Defaults to `binomial()`.
#' @param ... Additional arguments passed to `glm()` or `mgcv::gam()`.
#'
#' @returns A list with components:
#' \describe{
#'   \item{data}{The input data with two added prediction columns.}
#'   \item{model}{The fitted model object (`glm` or `gam`).}
#'   \item{pred}{A numeric vector of response-scale predictions (same as `pred_prob_col`).}
#' }
#'
#' @export
#'
#' @examples
#' set.seed(1)
#' dat <- data.frame(
#'   y = rbinom(100, 1, 0.4),
#'   x = rnorm(100)
#' )
#'
#' res <- run_model_add_pred(
#'   data = dat,
#'   formula = y ~ x,
#'   model = "glm",
#'   pred_link_col = "lp",
#'   pred_prob_col = "p"
#' )
#'
#' head(res$data[, c("y", "x", "lp", "p")])
run_model_add_pred <- function(data, formula,
                               model = c("glm","gam"),
                               pred_link_col,
                               pred_prob_col,
                               family = binomial(),
                               ...) {
  model <- match.arg(model)
  
  fitted_model <- if(model == "glm"){
    glm(formula, data = data, family = family, ...)
  } else {
    mgcv::gam(formula, data = data, family = family, optimizer = "efs", ...)
  }
  
  if(model == "glm"){
    pred_link <- stats::predict.glm(object = fitted_model, newdata = data, type = "link")
    pred_prob <- stats::predict.glm(object = fitted_model, newdata = data, type = "response")
  } else {
    pred_link <- mgcv::predict.gam(object = fitted_model, newdata = data)
    pred_prob <- mgcv::predict.gam(object = fitted_model, newdata = data, type = "response")
  }
  
  data[[pred_link_col]] <- as.numeric(pred_link)
  data[[pred_prob_col]] <- as.numeric(pred_prob)
  
  list(data = data, model = fitted_model, pred = pred_prob)
}


#' Derive ROC "Best" Threshold and Add a Binary Cutoff Column
#'
#' Builds a ROC curve using `pROC::roc()`, selects the "best" cutoff via
#' `pROC::coords(..., "best")` (by default maximizing Youden's index), and adds a
#' new binary column to the data based on that threshold.
#'
#' @param dat A data frame containing the predictor column.
#' @param y Vector of observed class labels for ROC analysis (typically 0/1 or a
#'   two-level factor).
#' @param pred_name Character scalar. Column name in `dat` containing the numeric
#'   prediction scores/probabilities.
#' @param cut_name Character scalar. Name of the new column to create (0/1),
#'   indicating whether `pred_name` is above the selected threshold.
#' @param direction Passed to `pROC::roc()`. Controls the comparison direction;
#'   default is `"<"`.
#' @param levels Passed to `pROC::roc()`. The order of outcome levels (controls
#'   which class is treated as "control" vs "case"); default is `c(0, 1)`.
#' @param plot_it Logical. If `TRUE`, plots the ROC curve and marks the "best"
#'   threshold.
#' @param ... Additional arguments passed to `pROC::roc()`.
#'
#' @returns A list with components:
#' \describe{
#'   \item{data}{`dat` with an added binary cutoff column `cut_name`.}
#'   \item{roc}{The `pROC::roc` object.}
#'   \item{best}{The result returned by `pROC::coords(roc_obj, "best")`.}
#'   \item{thr}{Numeric value of the selected threshold.}
#' }
#'
#' @export
#'
#' @examples
#' if (requireNamespace("pROC", quietly = TRUE)) {
#'   dat <- data.frame(pred = c(0.1, 0.4, 0.35, 0.8))
#'   y <- c(0, 0, 1, 1)
#'   res <- roc_best_addcut(dat, y, pred_name = "pred", cut_name = "pred_cut")
#'   res$thr
#'   res$data
#' }
roc_best_addcut <- function(dat, y, pred_name, cut_name,
                            direction = "<", levels = c(0,1),
                            plot_it = FALSE, ...) {
  
  roc_obj <- pROC::roc(response = y,
                       predictor = dat[[pred_name]],
                       direction = direction,
                       levels = levels,
                       ...)
  
  best <- pROC::coords(roc_obj, "best")
  thr  <- as.numeric(best["threshold"]$threshold)
  
  dat[[cut_name]] <- ifelse(dat[[pred_name]] >= thr, 1, 0)
  
  if(plot_it){
    plot(roc_obj, legacy.axes = TRUE,
         thresholds="best",
         print.thres="best")
  }
  
  list(data = dat, roc = roc_obj, best = best, thr = thr)
}

format_p_special <- function(p, inclusive = FALSE, ...) {
  ifelse(
    is.na(p), NA_character_,
    ifelse(
      p < 0.001,
      "P<.001",
      ifelse(
        p < 0.01,  # 新增：0.001 ≤ p < 0.01 用三位小数
        stringr::str_c("P=", sprintf("%0.3f", p)),
        {
          use3 <- if (inclusive) (p >= 0.045 & p <= 0.05) else (p > 0.045 & p < 0.05)
          fmt <- ifelse(use3, "%0.3f", "%0.2f")
          stringr::str_c("P=", sprintf(fmt, p))
        }
      )
    )
  )
}

