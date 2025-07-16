ensure_dir <- function(file_path){
  if(!dir.exists(file_path)) dir.create(file_path,recursive = T)
}

print_auc_and_ci <- function(roc_object){
  AUC_roc <- sprintf("%0.3f",
                     pROC::ci.auc(roc_object)[2])
  ci_roc <- paste0(sprintf("%0.3f",
                           pROC::ci.auc(roc_object)[1]),"-",
                   sprintf("%0.3f",
                           pROC::ci.auc(roc_object)[3]))
  AUC_formatted <- str_c(AUC_roc," (", ci_roc, ")")
  return(AUC_formatted)
}

calculate_performance2 <- function(observed, predicted){ 
  auc <- pROC::auc(observed~predicted) 
  glm.fit <- summary(glm(observed~predicted, family = binomial)) 
  calibration.intercept <- glm.fit$coef[1,1] 
  calibration.slope <- glm.fit$coef[2,1] 
  vec <- c(auc, calibration.intercept, calibration.slope) 
  names(vec) <- c("auc", "calibration intercept", "calibration slope") 
  return(vec) 
}

calcuate_diagnosis_index2 <- function(data,
                                      indicators=c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA","EBVDNA"),
                                      outcome="group_ref",
                                      calcLRCI_use="analytic", #"BayesianLR.test" or "analytic"
                                      binomMethod_use = "wilson") # c("all","exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit")
{
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
    negLR = character(0)
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
    
    ## Sen, spe, PPV, NPV, posLR, negLR
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
      # dplyr::mutate(across(c(posLR:negLR_UB), ~ sprintf("%0.3f", .x))) %>% 
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
    
    result_df <- rbind(result_df, diagCI_result_tidy)
  }
  
  colnames(result_df) <- c("Indicator",
                           "TP","FP","FN","TN",
                           "Sensitivity",
                           "Specificity",
                           "PPV",
                           "NPV",
                           "posLR",
                           "negLR")
  
  return(result_df)
}


DTC_pair <- function(data,
                     indicator1="P85",
                     compare_indicator=c("VCA_IgA", "EA_IgA", "NA1_IgA"),
                     outcome="group_ref") {
  if(!requireNamespace("tidyverse"))install.packages("tidyverse")
  if(!requireNamespace("DTComPair"))install.packages("DTComPair")  
  if(!requireNamespace("contingencytables"))remotes::install_github("ocbe-uio/contingencytables", ref="develop")
  library(DTComPair)
  library(tidyverse)
  library(contingencytables)
  
  result_compare <- data.frame(
    compare = character(0),
    sen1=character(0),
    sen2=character(0),
    sen_compM=character(0),
    sen_method=character(0),
    sen.p.value=character(0),
    spe1=character(0),
    spe2=character(0),
    spe_compM=character(0),
    spe_method=character(0),
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
      dplyr::mutate(across(where(is.factor),  ~ as.numeric(as.character(.x)))) 
    
    colnames(data_sub) <- c("group","indicator1","indicator2")
    
    result_pair <- DTComPair::tab.paired(d=group, y1=indicator1, y2=indicator2, data = data_sub)
    accuray_compare <- DTComPair::sesp.mcnemar(result_pair)
    PV_compare <- DTComPair::pv.wgs(result_pair)
    
    #### sensitivity
    matrix_sens <- result_pair$diseased[1:2,1:2]
    
    ## matrix for sensitivity
    sen_compM <- str_c(matrix_sens[1,1],matrix_sens[1,2],
                       matrix_sens[2,1],matrix_sens[2,2],sep = ",")
    
    sens_b <- matrix_sens[1,2]
    sens_c <- matrix_sens[2,1]
    
    if((sens_b+sens_c)>=25){
      sen.p.value <- signif(accuray_compare$sensitivity[5],2)
      sen_method_used <- "McNemar"
    } else if ((sens_b+sens_c)<25) {
      midp_value_sens <- contingencytables::McNemar_midP_test_paired_2x2(matrix_sens)
      sen.p.value <- signif(midp_value_sens$midP,2)
      sen_method_used <- "McNemar mid-P"
    }
    
    #### specificity
    matrix_spec <- result_pair$non.diseased[1:2,1:2]
    
    ## matrix for specificity
    spe_compM <- str_c(matrix_spec[1,1],matrix_spec[1,2],
                       matrix_spec[2,1],matrix_spec[2,2],sep = ",")
    
    spec_b <- matrix_spec[1,2]
    spec_c <- matrix_spec[2,1]
    
    if((spec_b+spec_c)>=25){
      spe.p.value <- signif(accuray_compare$specificity[5],2)
      spe_method_used <- "McNemar"
    } else if ((spec_b+spec_c)<25) {
      midp_value_spe <- contingencytables::McNemar_midP_test_paired_2x2(matrix_spec)
      spe.p.value <- signif(midp_value_spe$midP,2)
      spe_method_used <- "McNemar mid-P"
    } 
    
    ppv.p.value <- signif(PV_compare$ppv["p.value"],2)
    npv.p.value <- signif(PV_compare$npv["p.value"],2)
    
    result_compare1 <- data.frame(
      compare=stringr::str_c(indicator1, " vs ", indicator2),
      sen1=sprintf("%0.1f", accuray_compare$sensitivity["test1"]*100),
      sen2=sprintf("%0.1f", accuray_compare$sensitivity["test2"]*100),
      sen_compM=sen_compM,
      sen_method=sen_method_used,
      sen.p.value=ifelse(sen.p.value<0.0001,"p<0.0001",str_c("p=",sen.p.value)),
      spe1=sprintf("%0.1f", accuray_compare$specificity["test1"]*100),
      spe2=sprintf("%0.1f", accuray_compare$specificity["test2"]*100),
      spe_compM=spe_compM,
      spe_method=spe_method_used,
      spe.p.value=ifelse(spe.p.value<0.0001,"p<0.0001",str_c("p=",spe.p.value)),
      ppv1=sprintf("%0.1f", PV_compare$ppv["test1"]*100),
      ppv2=sprintf("%0.1f", PV_compare$ppv["test2"]*100),
      ppv.p.value=ifelse(ppv.p.value<0.0001,"p<0.0001",str_c("p=",ppv.p.value)),
      npv1=sprintf("%0.1f", PV_compare$npv["test1"]*100),
      npv2=sprintf("%0.1f", PV_compare$npv["test2"]*100),
      npv.p.value=ifelse(npv.p.value<0.0001,"p<0.0001",str_c("p=",npv.p.value))
    )
    result_compare <- rbind(result_compare, result_compare1)
  }
  
  return(result_compare)
}
