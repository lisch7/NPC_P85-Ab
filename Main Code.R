#### R code for BNLF2b cohort study
rm(list=ls())

# Environmental preparation ----
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/")) 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

if (!require("pak", quietly = TRUE)) install.packages("pak")

CRAN_packages_install <- c("remotes","devtools",
                           "readxl","writexl",
                           "tableone","compareGroups", # baseline table
                           "tidyverse","tidyplots","pROC","patchwork",# plot
                           "DTComPair","contingencytables","rstatix","statpsych",
                           "this.path","ggVennDiagram","autoReg","mgcv"
                           )
for (package in CRAN_packages_install) {
  pak::pak(pkg = package)
}

#### load function
R_script_path <- this.path::this.dir()
source(sprintf("%s/Function.R",R_script_path))

#### data prepare
suppressMessages(library(tidyverse))
suppressMessages(library(tidyplots))
suppressMessages(library(pROC))

current_path <- getwd()
current_path

save.path <- sprintf("%s/results/",current_path)
ensure_dir(save.path)

# data <- readxl::read_xlsx("combined_data.xlsx")

data$age <- as.numeric(data$age)
data$gender <- factor(data$gender,levels = c("Male","Female"))
data$NPC_stage <- factor(data$NPC_stage,levels=c("I","II","III","IVA","IVB","Unknown","/"))
data$group_ref <- factor(data$group,levels = c("NPC","nonNPC"),labels = c("1","0"))
data$group <- factor(data$group,levels = c("NPC","nonNPC"),labels = c("NPC","nonNPC"))
data$P85 <- factor(data$P85,levels = c("(+)","(-)"),labels = c("1","0"))
data$VCA_IgA <- factor(data$VCA_IgA,levels = c("(+)","(-)"),labels = c("1","0"))
data$EA_IgA <- factor(data$EA_IgA,levels = c("(+)","(-)"),labels = c("1","0"))
data$NA1_IgA <- factor(data$NA1_IgA,levels = c("(+)","(-)"),labels = c("1","0"))
data$EBV_DNA_absolute <- ifelse(data$EBV_DNA_absolute=="NA",NA,data$EBV_DNA_absolute)
data$EBV_DNA_absolute <- as.numeric(data$EBV_DNA_absolute)
data$EBV_DNA0 <- ifelse(data$EBV_DNA0=="NA",NA,data$EBV_DNA0)
data$EBV_DNA0 <- factor(data$EBV_DNA0,levels = c("(+)","(-)"),labels = c("1","0"))
data$center <- factor(data$center,levels = c("SYSUCC","ZSCPH","TJH-HUST","WZRCH","FJCH")) 

#### Parameter Definition 
correct_method <- "holm"

# Data Preparation ----
## ColorSet ----
if(F){
  my_color <- list()
  my_color$single <- c("#374E55ff","#DF8F44FF","#968b79","#94b2be")
  names(my_color$single) <- c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")
  my_color$model <- c("#374E55ff","#6A6599FF","#DF8F44FF","#00A1D5FF","#79AF97FF")
  names(my_color$model) <- c("P85-Ab","VCA_EBNA1","P85_VCA","P85_EBNA1","P85_VCA_EBNA1")
  save(my_color,file = "./my_color.rdata")
}

## Model D&V ----
if(T){
  if (!require("autoReg")) install.packages("autoReg")
  if(!requireNamespace("mgcv",quietly = TRUE))install.packages("mgcv")
  suppressWarnings(library(mgcv))
  set.seed(100)
  
  data_logistic <- data
  data_logistic$group_ref <- ifelse(data_logistic$group=="NPC",1,0)
  
  models <- list(
    list(
      name = "dual_method",
      type = "glm",
      formula = group_ref ~ VCA_IgA_SCO + NA1_IgA_SCO
    ),
    list(
      name = "gam_P85_VCA",
      type = "gam",
      formula = group_ref == 1 ~ s(P85_SCO) + s(VCA_IgA_SCO)
    ),
    list(
      name = "gam_P85_EBNA1",
      type = "gam",
      formula = group_ref == 1 ~ s(P85_SCO) + s(NA1_IgA_SCO)
    ),
    list(
      name = "gam_P85_dual",
      type = "gam",
      formula = group_ref == 1 ~ s(P85_SCO) + s(NA1_IgA_SCO) + s(VCA_IgA_SCO)
    )
  )
  
  fits <- list()
  rocs <- list()
  perf_list <- list()
  thr_list  <- list()
  
  for(m in models){
    
    tmp <- run_model_add_pred(
      data = data_logistic,
      formula = m$formula,
      model = m$type,
      pred_link_col = sprintf("%s",m$name),
      pred_prob_col = sprintf("%s_prob",m$name)
    )
    
    data_logistic <- tmp$data
    fits[[m$name]] <- tmp$model
    
    cat("model done:", m$name, "\n")
    
    tmp2 <- roc_best_addcut(
      dat = data_logistic,
      y = data_logistic$group_ref,
      pred_name = sprintf("%s_prob",m$name),
      cut_name = sprintf("%s_cut",m$name),
      plot_it = FALSE                  # 想每个都画就 TRUE
    )
    
    data_logistic <- tmp2$data
    rocs[[m$name]] <- tmp2
    
    thr_list[[m$name]] <- tmp2$thr
    
    cat("ROC done:", m$name, "\n")
    
    perf_list[[m$name]] <- calculate_performance2(
      observed  = data_logistic$group_ref == 1,
      predicted = data_logistic[[m$name]]   # link列
    )
    
  }
  
  ## ** save model
  thr_list %>% names()
  model_data <- list(
    "P85-Ab+VCA-IgA" = list(
      model =  fits[["gam_P85_VCA"]],
      cutoff = thr_list[["gam_P85_VCA"]] 
    ),
    "P85-Ab+EBNA1-IgA" = list(
      model = fits[["gam_P85_EBNA1"]],
      cutoff = thr_list[["gam_P85_EBNA1"]]  
    ),
    "Triplet-Antibody Strategy" = list(
      model = fits[["gam_P85_dual"]],
      cutoff = thr_list[["gam_P85_dual"]]  
    )
  )
  
  saveRDS(model_data,sprintf("%s/Model.rds",save.path))
}

# Main Table ----
## Table2 ----
saveDir <- "Table2"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 
data_pre <- data %>% as.data.frame()
result_df <- calcuate_diagnosis_index2(data=data_pre,
                                       indicators = indicators,
                                       binomMethod_use = "wilson",
                                       outcome = "group_ref")

write.csv(result_df,
          sprintf("%s/%s_Raw_Data.csv",save_path,saveDir),row.names = F)

result_compare <- DTC_pair(data=data_pre,
                           indicator1="P85",
                           compare_indicator=indicators[-1],
                           outcome="group_ref")
result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
result_compare$sens.p.adjust <- ifelse(result_compare$sen_method=="McNemar",
                                       result_compare$sens.p.adjust,str_c(result_compare$sens.p.adjust,"b"))
result_compare$spec.p.adjust <- ifelse(result_compare$spe_method=="McNemar",
                                       result_compare$spec.p.adjust,str_c(result_compare$spec.p.adjust,"b"))
write.csv(result_compare,
          sprintf("%s/%s_Raw_Data_compare.csv",save_path,saveDir),row.names = F)


#### Combined for publication & verification 
merged.data <- read.csv(sprintf("%s/%s_Raw_Data.csv",save_path,saveDir)) %>% 
  dplyr::select(-c(PPV:negLR))
merged.statistical.data <- read.csv(sprintf("%s/%s_Raw_Data_compare.csv",save_path,saveDir)) %>% 
  dplyr::select(compare,sens.p.adjust,spec.p.adjust)

value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  tidyr::separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged.final <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator") %>% 
  dplyr::select(-refGroup) %>% 
  dplyr::mutate(across(everything(), ~ replace_na(., "ref")))

write.csv(merged.final,
          sprintf("%s/%s_Combined_for_publication.csv",save_path,saveDir),
          row.names = F)

## Table3 ----
saveDir <- "Table3"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

symptom_data <- data_logistic

# diagnostic performance
indicators <- c("P85","gam_P85_dual_cut")
candidate_symptom <- c('Asymptomatic','NPC-specific Symptoms Cluster',
                       'NPC-nonspecific Symptoms Cluster')

result_list <- list()
result_statistical_list <- list()
for (i in 1:length(candidate_symptom)) {
  used_symptom <- candidate_symptom[i]
  print(used_symptom)
  data_sub <- symptom_data %>% dplyr::filter(!!sym(used_symptom)==1) %>% as.data.frame()
  
  result_df <- calcuate_diagnosis_index2(data=data_sub,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df$symptom <- used_symptom
  result_list[[i]] <- result_df
  
  ## compare_result
  result_compare <- DTC_pair(data=data_sub,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$symptom <-  used_symptom
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_statistical_list[[i]] <- result_compare
}

merged.data <- map_dfr(result_list,bind_rows)
writexl::write_xlsx(merged.data,
                    sprintf("%s/%s_Five_centers_Model_symptom_combined.xlsx",save_path,saveDir))

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
writexl::write_xlsx(merged.statistical.data.combined,
                    sprintf("%s/%s_Five_centers_Model_symptom_compare_combined.xlsx",save_path,saveDir))

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(sens_p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b")),
    spec_p_value=case_when(
      spe_method=="McNemar"~spec.p.adjust,
      spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
    )) 

write.csv(merged.statistical.data,
          sprintf("%s/%s_Five_centers_Model_symptom_compare.csv",save_path,saveDir),
          row.names = F)

# Main Figure ----
## Figure 2 ----
saveDir <- "Figure2"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_combine <- data

centers <- unique(data_combine$center)
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 

result_list <- list()
result_statistical_list <- list()

for (i in 1:length(centers)) {
  center_name <- centers[i]
  print(center_name)
  
  data_pre <- data_combine %>% 
    filter(.$center == center_name) %>%
    as.data.frame()
  
  if(nrow(data_pre) == 0) next 
  
  result_df <- calcuate_diagnosis_index2(data=data_pre,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df$Center <- center_name
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=data_pre,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Center <- center_name
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
  result_compare$ppv.p.adjust.raw <- p.adjust(result_compare$ppv.p.raw, method = correct_method)
  result_compare$npv.p.adjust.raw <- p.adjust(result_compare$npv.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_compare$ppv.p.adjust <- format_p_special(result_compare$ppv.p.adjust.raw)
  result_compare$npv.p.adjust <- format_p_special(result_compare$npv.p.adjust.raw)
  result_compare
  result_statistical_list[[i]] <- result_compare
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  dplyr::select(Indicator,Center,everything())

write.csv(merged.data,
          sprintf("%s/%s_Combined_Center_Raw_Data.csv",save_path,saveDir),
          row.names = F)

merged.data.statistics <- map_dfr(result_statistical_list,bind_rows) %>% 
  dplyr::select(compare,Center,everything()) %>% 
  dplyr::mutate(sens_p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b")),
    spec_p_value=case_when(
      spe_method=="McNemar"~spec.p.adjust,
      spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
    )) 
merged.data.statistics <- merged.data.statistics[order(merged.data.statistics$Center),]
merged.data.statistics
write.csv(merged.data.statistics,
          sprintf("%s/%s_Combined_Center_Raw_Data_statistics.csv",
                  save_path,saveDir),
          row.names = F)

data_forest <- merged.data %>% as.data.frame()
indicators_used <- c("P85","VCA_IgA","EA_IgA","NA1_IgA")

data_forest$Center <- factor(data_forest$Center,
                             levels = c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH"))

matrics <- c("Sensitivity","Specificity")
plot_list <- list()

for (i in 1:length(matrics)){
  used_matrics <- matrics[i]
  print(used_matrics)
  
  p <- data_forest %>%
    dplyr::select(Indicator,Center,!!sym(used_matrics)) %>%
    dplyr::rename(est=!!sym(used_matrics)) %>%
    mutate(
      est=str_extract(.$est,".*(?=\n)"),
      low_ci=str_extract(.$est,"(?<=\\().*(?= to)"),
      high_ci=str_extract(.$est,"(?<=to ).*(?=\\))")) %>%
    dplyr::filter(Indicator%in%c("P85","VCA_IgA","EA_IgA","NA1_IgA")) %>%
    mutate(
      est = as.numeric(est),
      low_ci = as.numeric(low_ci),
      high_ci = as.numeric(high_ci)
    ) %>%
    mutate(Indicator=factor(.$Indicator,levels = rev(c("P85","VCA_IgA","EA_IgA","NA1_IgA")),
                            labels = rev(c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")))) %>%
    ggplot(aes(x=est,y=Indicator))+
    geom_point(aes(color=Indicator))+
    geom_errorbar(aes(xmin =low_ci,xmax = high_ci,color=Indicator,group=Indicator),
                  width=0.5,position = position_dodge(0.6)) +
    geom_text(aes(x = est,y=as.numeric(Indicator)+0.5,
                  label = sprintf("%.1f",est),
                  color=Indicator))+
    geom_vline(xintercept=90,linetype="dashed",color="grey50")+
    facet_wrap(~Center,ncol=1, strip.position = "left")+
    ggpubr::theme_classic2()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),
      strip.placement = "outside",
      panel.spacing.y = unit(0.5, "lines"),
      axis.text = element_text(color="black")
    ) +
    scale_color_manual(values = my_color$single)+
    xlab(str_c(used_matrics," (%)"))+
    ylab("")
  
  plot_list[[i]] <- p
  ggsave(sprintf("%s/%s_PanelA_%s.pdf",save_path,saveDir,used_matrics),plot = p,
         width = 6,height = 8)
}

# eFigure ----
## eFigure 3 ----
saveDir <- "eFigure3"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_pre <- data  %>% as.data.frame()

### eFigure 3A ----
library(ggVennDiagram)

#### NPC cases
data_plot <- data %>% dplyr::filter(group=="NPC") %>% 
  mutate(`P85-Ab`=as.character(.$P85),
         `VCA-IgA`=as.character(.$VCA_IgA),
         `EA-IgA`=as.character(.$EA_IgA),
         `EBNA1-IgA`=as.character(.$NA1_IgA))

sets <- list(
  `P85-Ab`= which(data_plot$`P85-Ab` == 1),
  `VCA-IgA`= which(data_plot$`VCA-IgA` == 1),
  `EA-IgA`=which(data_plot$`EA-IgA` == 1),
  `EBNA1-IgA`=which(data_plot$`EBNA1-IgA` == 1))

ggVennDiagram(sets,set_color =  my_color$single[1:4],label="count")+
  scale_fill_gradient(low="grey90",high = "#B24745FF")+
  ggtitle("NPC group (n=1,680)")

ggsave(sprintf("%s/%s_panelA_Veen_NPC.pdf",save_path,saveDir),width = 5,height = 5)

#### non-NPC
data_plot <- data %>% dplyr::filter(group=="nonNPC") %>% 
  mutate(`P85-Ab`=as.character(.$P85),
         `VCA-IgA`=as.character(.$VCA_IgA),
         `EA-IgA`=as.character(.$EA_IgA),
         `EBNA1-IgA`=as.character(.$NA1_IgA))

sets <- list(
  `P85-Ab`= which(data_plot$`P85-Ab` == 0),
  `VCA-IgA`= which(data_plot$`VCA-IgA` == 0),
  `EA-IgA`=which(data_plot$`EA-IgA` == 0),
  `EBNA1-IgA`=which(data_plot$`EBNA1-IgA` == 0))

ggVennDiagram(sets,set_color = my_color$single[1:4],label="count")+
  scale_fill_gradient(low="grey90",high = "#B24745FF")+
  ggtitle("nonNPC group (n=2,097)")

ggsave(sprintf("%s/%s_panelA_Veen_nonNPC.pdf",save_path,saveDir),width = 5,height = 5)


### eFigure 2B ----

indicators <- c("P85_SCO", "VCA_IgA_SCO","EA_IgA_SCO","NA1_IgA_SCO") 

result_roc_df <- data.frame(
  Indicator = character(0),
  AUC = character(0),
  P_value.raw = numeric(),
  P_value = character(0),
  stringsAsFactors = FALSE)

roc_p85 <- pROC::roc(data_pre$group_ref, 
                     as.numeric(data_pre$P85_SCO), 
                     direction="<", levels = c(0, 1))

for(i in 1:length(indicators)){
  indicator <- indicators[i]
  data_pre[[indicators[i]]] <- as.numeric(data_pre[[indicators[i]]])
  roc_object <- pROC::roc(data_pre$group_ref, data_pre[[indicators[i]]],
                          direction="<",
                          levels=c(0, 1))
  AUC_value_with_DeLongCI <- ci.auc(roc_object)
  AUC_value_with_3_digits <- signif(AUC_value_with_DeLongCI[2],2)
  AUC_CI_lower <- signif(AUC_value_with_DeLongCI[1],2)
  AUC_CI_upper <- signif(AUC_value_with_DeLongCI[3],2)
  AUC_formatted <- str_c(AUC_value_with_3_digits, 
                         "\n(", AUC_CI_lower, "-", AUC_CI_upper, ")")
  
  if (indicator != "P85_SCO") {
    roc_test_result <- roc.test(roc_p85, roc_object)
    P_value.raw <- roc_test_result$p.value
    P_value <- format_p_special(roc_test_result$p.value)
  } else {
    P_value.raw <- NA
    P_value <- "ref"
  }

  result_roc_df <- rbind(result_roc_df, c(indicator, AUC_formatted, P_value.raw,P_value))
}

colnames(result_roc_df) <- c("Indicator", "AUC", "P_value.raw","P_value")
result_roc_df

result_roc_df$p.adjust.raw <- p.adjust(result_roc_df$P_value.raw, method = correct_method,n=3)
result_roc_df$p.adjust <- format_p_special(result_roc_df$p.adjust.raw)

write.csv(result_roc_df,
          sprintf("%s/%s_panelB_AUROC_Raw_Data.csv",save_path,saveDir),
          row.names = F)

#### plot
roc_p85 <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$P85_SCO),direction="<",
                     levels = c(0, 1))
value_p85 <- print_auc_and_ci(roc_p85)

roc_vca <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$VCA_IgA_SCO),direction="<",
                     levels=c(0, 1))
value_vca <- print_auc_and_ci(roc_vca)

roc_ea <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$EA_IgA_SCO),direction="<",
                    levels=c(0, 1))
value_ea <- print_auc_and_ci(roc_ea)

roc_na1 <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$NA1_IgA_SCO),direction="<",
                     levels=c(0, 1))
value_na1 <- print_auc_and_ci(roc_na1)

pdf(sprintf("%s/%s_ROC.pdf",save_path,saveDir),width = 6,height = 6)
{
  plot(roc_p85, 
       col=my_color$single["P85-Ab"],
       legacy.axes=TRUE)
  
  plot.roc(roc_vca,
           add=TRUE, 
           col=my_color$single["VCA-IgA"], 
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  
  plot.roc(roc_ea,
           add=TRUE,
           col=my_color$single["EA-IgA"],
           print.auc.x=0.35,
           print.auc.y=0.50) 
  
  plot.roc(roc_na1,
           add=TRUE,
           col=my_color$single["EBNA1-IgA"],
           print.auc.x=0.35,
           print.auc.y=0.45)
  
  text <- c(paste0("P85-Ab\n",value_p85),
            paste0("VCA-IgA\n",value_vca),
            paste0("EA-IgA\n",value_ea),
            paste0("EBNA1-IgA\n",value_na1))
  
  legend(0.5, 0.5, 
         bty = "n", 
         legend=text,
         text.col = c(my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")]),
         col= c(my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")]),
         lwd=3)
}
dev.off()

## eFigure 4 ----
saveDir <- "eFigure4"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_combine <- data

centers <- unique(data_combine$center)
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 

result_list <- list()
result_statistical_list <- list()
AUC_list <- list()
delong_list <- list()

for (i in 1:length(centers)) {
  center_name <- centers[i]
  print(center_name)

  data_pre <- data_combine %>% 
    filter(.$center == center_name) %>%
    as.data.frame()
  
  if(nrow(data_pre) == 0) next 
  
  result_df <- calcuate_diagnosis_index2(data=data_pre,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df$Center <- center_name
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=data_pre,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Center <- center_name
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
  result_compare$ppv.p.adjust.raw <- p.adjust(result_compare$ppv.p.raw, method = correct_method)
  result_compare$npv.p.adjust.raw <- p.adjust(result_compare$npv.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_compare$ppv.p.adjust <- format_p_special(result_compare$ppv.p.adjust.raw)
  result_compare$npv.p.adjust <- format_p_special(result_compare$npv.p.adjust.raw)
  result_compare
  result_statistical_list[[i]] <- result_compare
  
  if(T){
    roc_p85 <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$P85_SCO),direction="<",
                         levels = c(0, 1))
    value_p85 <- print_auc_and_ci(roc_p85)
    
    roc_vca <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$VCA_IgA_SCO),direction="<",
                         levels=c(0, 1))
    value_vca <- print_auc_and_ci(roc_vca)
    
    roc_ea <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$EA_IgA_SCO),direction="<",
                        levels=c(0, 1))
    value_ea <- print_auc_and_ci(roc_ea)
    
    roc_na1 <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$NA1_IgA_SCO),direction="<",
                         levels=c(0, 1))
    value_na1 <- print_auc_and_ci(roc_na1)
    
    pdf(sprintf("%s/%s_ROC_%s.pdf",save_path,saveDir,center_name),width = 6,height = 6)
    {
      plot(roc_p85, 
           col=my_color$single["P85-Ab"],
           legacy.axes=TRUE)
      
      plot.roc(roc_vca,
               add=TRUE,
               col=my_color$single["VCA-IgA"],
               print.auc.x=0.35,
               print.auc.y=0.55)

      plot.roc(roc_ea,
               add=TRUE,
               col=my_color$single["EA-IgA"],
               print.auc.x=0.35,
               print.auc.y=0.50) 
      
      plot.roc(roc_na1,
               add=TRUE,
               col=my_color$single["EBNA1-IgA"],
               print.auc.x=0.35,
               print.auc.y=0.45)
      
      text<-c(paste0("P85-Ab\n",value_p85),
              paste0("VCA-IgA\n",value_vca),
              paste0("EA-IgA\n",value_ea),
              paste0("EBNA1-IgA\n",value_na1))
      
      legend(0.5, 0.5,
             bty = "n",
             legend=text,
             text.col = my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")],
             col= my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")],
             lwd=3)
    }
    dev.off()
    
    compare1 <- roc.test(roc_p85,roc_vca)
    compare1_pvalue.raw <- compare1$p.value
    compare1_pvalue <- format_p_special(compare1$p.value)
    
    compare2 <- roc.test(roc_p85,roc_ea)
    compare2_pvalue.raw <- compare2$p.value
    compare2_pvalue <- format_p_special(compare2$p.value)
    
    compare3 <- roc.test(roc_p85,roc_na1)
    compare3_pvalue.raw <- compare3$p.value
    compare3_pvalue <- format_p_special(compare3$p.value)
    
    delong_df <- data.frame(
      center=center_name,
      Indicator = c("P85 vs VCA","P85 vs EA","P85 vs NA1"),
      P_value.raw = c(compare1_pvalue.raw,compare2_pvalue.raw,compare3_pvalue.raw),
      P_value = c(compare1_pvalue,compare2_pvalue,compare3_pvalue))
    delong_df$p.adjust.raw <- p.adjust(delong_df$P_value.raw, method = correct_method,n=3)
    delong_df$p.adjust <- format_p_special(delong_df$p.adjust.raw)
    delong_list[[i]] <- delong_df
    
    AUC_df <- data.frame(
      center=center_name,
      Indicator = c("P85","VCA","EA","NA1"), 
      AUC = c(value_p85,value_vca,value_ea,value_na1) 
    )
    
    AUC_list[[i]] <- AUC_df
  }
}

AUC.merge <- map_dfr(AUC_list,bind_rows)
ensure_dir(save_path)
write.csv(AUC.merge,
          sprintf("%s/%s_AUROC_Raw_Data.csv",save_path,saveDir),
          row.names = F)

delong.merge <- map_dfr(delong_list,bind_rows)
delong.merge <- delong.merge[order(delong.merge$center),]
write.csv(delong.merge,
          sprintf("%s/%s_AUROC_delong_Raw_Data.csv",save_path,saveDir),
          row.names = F)

## eFigure 5 ----
saveDir <- "eFigure5"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_combine <- data

centers <- unique(data_combine$center)
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 

result_list <- list()
result_statistical_list <- list()

for (i in 1:length(centers)) {
  center_name <- centers[i]
  print(center_name)
  
  data_pre <- data_combine %>% 
    filter(.$center == center_name) %>%
    as.data.frame()
  
  if(nrow(data_pre) == 0) next 
  
  result_df <- calcuate_diagnosis_index2(data=data_pre,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df$Center <- center_name
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=data_pre,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Center <- center_name
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
  result_compare$ppv.p.adjust.raw <- p.adjust(result_compare$ppv.p.raw, method = correct_method)
  result_compare$npv.p.adjust.raw <- p.adjust(result_compare$npv.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_compare$ppv.p.adjust <- format_p_special(result_compare$ppv.p.adjust.raw)
  result_compare$npv.p.adjust <- format_p_special(result_compare$npv.p.adjust.raw)
  result_compare
  result_statistical_list[[i]] <- result_compare
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  dplyr::select(Indicator,Center,everything())

write.csv(merged.data,
          sprintf("%s/%s_Combined_Center_Raw_Data.csv",save_path,saveDir),
          row.names = F)

merged.data.statistics <- map_dfr(result_statistical_list,bind_rows) %>% 
  dplyr::select(compare,Center,everything()) %>% 
  dplyr::mutate(sens_p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b")),
    spec_p_value=case_when(
      spe_method=="McNemar"~spec.p.adjust,
      spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
    )) 
merged.data.statistics <- merged.data.statistics[order(merged.data.statistics$Center),]
merged.data.statistics
write.csv(merged.data.statistics,
          sprintf("%s/%s_Combined_Center_Raw_Data_statistics.csv",
                  save_path,saveDir),
          row.names = F)

data_forest <- merged.data %>% as.data.frame()
indicators_used <- c("P85","VCA_IgA","EA_IgA","NA1_IgA")

data_forest$Center <- factor(data_forest$Center,
                             levels = c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH"))

matrics <- c("PPV","NPV")
plot_list <- list()

for (i in 1:length(matrics)){
  used_matrics <- matrics[i]
  print(used_matrics)
  
  p <- data_forest %>%
    dplyr::select(Indicator,Center,!!sym(used_matrics)) %>%
    dplyr::rename(est=!!sym(used_matrics)) %>%
    mutate(
      est=str_extract(.$est,".*(?=\n)"),
      low_ci=str_extract(.$est,"(?<=\\().*(?= to)"),
      high_ci=str_extract(.$est,"(?<=to ).*(?=\\))")) %>%
    dplyr::filter(Indicator%in%c("P85","VCA_IgA","EA_IgA","NA1_IgA")) %>%
    mutate(
      est = as.numeric(est),
      low_ci = as.numeric(low_ci),
      high_ci = as.numeric(high_ci)
    ) %>%
    # mutate(text=sprintf("%.1f (%.1f to %.1f)",est,low_ci,high_ci)) %>% 
    mutate(Indicator=factor(.$Indicator,levels = rev(c("P85","VCA_IgA","EA_IgA","NA1_IgA")),
                            labels = rev(c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")))) %>%
    ggplot(aes(x=est,y=Indicator))+
    geom_point(aes(color=Indicator),shape=15)+
    geom_errorbar(aes(xmin =low_ci,xmax = high_ci,color=Indicator,group=Indicator),
                  width=0.5,position = position_dodge(0.6)) +
    geom_text(aes(x = est,y=as.numeric(Indicator)+0.5,
                  label = sprintf("%.1f",est),
                  color=Indicator))+
    geom_vline(xintercept=90,linetype="dashed",color="grey50")+
    facet_wrap(~Center,ncol=1, strip.position = "left")+
    ggpubr::theme_classic2()+
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.background = element_blank(),  # 去掉背景
      strip.placement = "outside",
      panel.spacing.y = unit(0.5, "lines"),
      axis.text = element_text(color="black")
    ) +
    scale_color_manual(values = my_color$single)+
    xlab(str_c(used_matrics," (%)"))+
    ylab("")
  
  plot_list[[i]] <- p
  ggsave(sprintf("%s/%s_%s_across_centers.pdf",save_path,saveDir,used_matrics),plot = p,
         width = 6,height = 5)
}

## eFigure 6 & eFigure 7 ----
saveDir <- "eFigure6"
saveDir_AUROC <- "eFigure7"
save_path <- file.path(save.path,saveDir)
save_path_AUROC <- file.path(save.path,saveDir_AUROC)
ensure_dir(save_path)
ensure_dir(save_path_AUROC)
load("./my_color.rdata")

### Sex ----
data_combine <- data
sex_class <- unique(data_combine$gender)
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA")

result_list <- list()
auc_list <- list()

for (i in 1:length(sex_class)) {
  sex_name <- sex_class[i]
  print(sex_name)
  
  data_used <- data_combine %>% 
    filter(.$gender == sex_name) %>%
    as.data.frame()
  
  if(nrow(data_used) == 0) next 
  
  result_df <- calcuate_diagnosis_index2(data=data_used,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df$group <- sex_name
  result_df$Total <- as.numeric(result_df$TP)+as.numeric(result_df$TN)+
    as.numeric(result_df$FP)+as.numeric(result_df$FN)
  
  result_df <- result_df %>%
    mutate(
      sen_est1 = as.numeric(str_extract(.$Sensitivity,".*(?=\n)"))/100,
      sen_low1 = as.numeric(str_extract(.$Sensitivity,"(?<=\\().*(?= to)")) / 100, 
      sen_hi1 = as.numeric(str_extract(.$Sensitivity,"(?<=to ).*(?=\\))")) / 100, 
      spe_est2 = as.numeric(str_extract(.$Specificity,".*(?=\n)"))/100,
      spe_low2 = as.numeric(str_extract(.$Specificity,"(?<=\\().*(?= to)")) / 100, 
      spe_hi2 = as.numeric(str_extract(.$Specificity,"(?<=to ).*(?=\\))")) / 100, 
    )
  
  result_list[[i]] <- result_df
  
  roc_p85 <- pROC::roc(data_used$group_ref, as.numeric(data_used$P85_SCO),direction="<",
                       levels = c(0, 1))
  value_p85 <- print_auc_and_ci(roc_p85)
  
  roc_vca <- pROC::roc(data_used$group_ref, as.numeric(data_used$VCA_IgA_SCO),direction="<",
                       levels=c(0, 1))
  value_vca <- print_auc_and_ci(roc_vca)
  
  roc_ea <- pROC::roc(data_used$group_ref, as.numeric(data_used$EA_IgA_SCO),direction="<",
                      levels=c(0, 1))
  value_ea <- print_auc_and_ci(roc_ea)
  
  roc_na1 <- pROC::roc(data_used$group_ref, as.numeric(data_used$NA1_IgA_SCO),direction="<",
                       levels=c(0, 1))
  value_na1 <- print_auc_and_ci(roc_na1)
  
  
  auc_df <- data.frame(
    Indicator = c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA"),
    AUC = c(signif(pROC::ci.auc(roc_p85)[2],2),
            signif(pROC::ci.auc(roc_vca)[2],2),
            signif(pROC::ci.auc(roc_ea)[2],2),
            signif(pROC::ci.auc(roc_na1)[2],2)),
    low = c(signif(pROC::ci.auc(roc_p85)[1],2),
            signif(pROC::ci.auc(roc_vca)[1],2),
            signif(pROC::ci.auc(roc_ea)[1],2),
            signif(pROC::ci.auc(roc_na1)[1],2)),
    hi = c(signif(pROC::ci.auc(roc_p85)[3],2),
           signif(pROC::ci.auc(roc_vca)[3],2),
           signif(pROC::ci.auc(roc_ea)[3],2),
           signif(pROC::ci.auc(roc_na1)[3],2)),
    format= c(value_p85,value_vca,value_ea,value_na1),
    group=sex_name)
  
  auc_list[[i]] <- auc_df
}

merged.data <- map_dfr(result_list,bind_rows) 
write.csv(merged.data,
          sprintf("%s/%s_sex_diagnostic_performance.csv",save_path,saveDir),row.names = F)

merged.auc_list <- map_dfr(auc_list,bind_rows) 
write.csv(merged.auc_list,
          sprintf("%s/%s_sex_AUROC.csv",save_path_AUROC,saveDir_AUROC),row.names = F)

### Age ----
data_combine <- data
data_age_calculate <- data_combine %>% 
  mutate(age_cut = case_when(
    age < 45 ~ "less than 45",
    age >= 45 ~ "more tham 45"
  ))

age_class <- unique(data_age_calculate$age_cut)
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 

result_list <- list()
auc_list <- list()

for (i in 1:length(age_class)) {
  age_name <- age_class[i]
  print(age_name)
  
  data_used <- data_age_calculate %>% 
    filter(.$age_cut == age_name) %>%
    as.data.frame()
  
  if(nrow(data_used) == 0) next 
  
  result_df <- calcuate_diagnosis_index2(data=data_used,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df$group <- age_name
  result_df$Total <- as.numeric(result_df$TP)+as.numeric(result_df$TN)+
    as.numeric(result_df$FP)+as.numeric(result_df$FN)
  
  result_df <- result_df %>%
    mutate(
      sen_est1 = as.numeric(str_extract(.$Sensitivity,".*(?=\n)"))/100,
      sen_low1 = as.numeric(str_extract(.$Sensitivity,"(?<=\\().*(?= to)")) / 100, 
      sen_hi1 = as.numeric(str_extract(.$Sensitivity,"(?<=to ).*(?=\\))")) / 100, 
      spe_est2 = as.numeric(str_extract(.$Specificity,".*(?=\n)"))/100,
      spe_low2 = as.numeric(str_extract(.$Specificity,"(?<=\\().*(?= to)")) / 100, 
      spe_hi2 = as.numeric(str_extract(.$Specificity,"(?<=to ).*(?=\\))")) / 100, 
    )
  
  result_list[[i]] <- result_df
  
  roc_p85 <- pROC::roc(data_used$group_ref, as.numeric(data_used$P85_SCO),direction="<",
                       levels = c(0, 1))
  value_p85 <- print_auc_and_ci(roc_p85)
  
  roc_vca <- pROC::roc(data_used$group_ref, as.numeric(data_used$VCA_IgA_SCO),direction="<",
                       levels=c(0, 1))
  value_vca <- print_auc_and_ci(roc_vca)
  
  roc_ea <- pROC::roc(data_used$group_ref, as.numeric(data_used$EA_IgA_SCO),direction="<",
                      levels=c(0, 1))
  value_ea <- print_auc_and_ci(roc_ea)
  
  roc_na1 <- pROC::roc(data_used$group_ref, as.numeric(data_used$NA1_IgA_SCO),direction="<",
                       levels=c(0, 1))
  value_na1 <- print_auc_and_ci(roc_na1)
  
  auc_df <- data.frame(
    Indicator = c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA"),
    AUC = c(signif(pROC::ci.auc(roc_p85)[2],2),
            signif(pROC::ci.auc(roc_vca)[2],2),
            signif(pROC::ci.auc(roc_ea)[2],2),
            signif(pROC::ci.auc(roc_na1)[2],2)),
    low = c(signif(pROC::ci.auc(roc_p85)[1],2),
            signif(pROC::ci.auc(roc_vca)[1],2),
            signif(pROC::ci.auc(roc_ea)[1],2),
            signif(pROC::ci.auc(roc_na1)[1],2)),
    hi = c(signif(pROC::ci.auc(roc_p85)[3],2),
           signif(pROC::ci.auc(roc_vca)[3],2),
           signif(pROC::ci.auc(roc_ea)[3],2),
           signif(pROC::ci.auc(roc_na1)[3],2)),
    format= c(value_p85,value_vca,value_ea,value_na1),
    group=age_name)
  
  auc_list[[i]] <- auc_df
}

merged.data <- map_dfr(result_list,bind_rows) 
write.csv(merged.data,
          sprintf("%s/%s_age_diagnostic_performance.csv",save_path,saveDir),row.names = F)

merged.auc_list <- map_dfr(auc_list,bind_rows) 
write.csv(merged.auc_list,
          sprintf("%s/%s_age_AUROC.csv",save_path_AUROC,saveDir_AUROC),row.names = F)


### eFigure 6 Plot ----
data_age <- read.csv(sprintf("%s/%s_age_diagnostic_performance.csv",save_path,saveDir))
data_sex <- read.csv(sprintf("%s/%s_sex_diagnostic_performance.csv",save_path,saveDir))
data_subgroup <- rbind(data_age,data_sex) 

suppressMessages(library(forestploter))

data_forest <- data_subgroup %>% 
  dplyr::select(Indicator,group,Total,TP:TN,sen_est1:spe_hi2)

custom_order <- c("Male", "Female", "less than 45", "more tham 45")
ordered_indices <- order(match(data_forest$group, custom_order))
sorted_data <- data_forest[ordered_indices, ]
data_forest <- sorted_data %>% as_tibble() %>% 
  add_row(group="Sex",.before = 1) %>% 
  add_row(group="Age,yr",.before = 10) 

writexl::write_xlsx(data_forest,
                    sprintf("%s/%s_combine_for_forest.xlsx",save_path,saveDir))

indictors_all <- c("P85","VCA_IgA","EA_IgA","NA1_IgA")

for (i in 1:length(indictors_all)){
  indicators_used <- indictors_all[i]
  print(indictors_all[i])
  data_forest_plot <- data_forest %>% 
    dplyr::filter(Indicator %in% indictors_all[i] | is.na(Indicator)) 
  data_forest_plot$group <- ifelse(is.na(data_forest_plot$Total), 
                                   data_forest_plot$group,
                                   paste0("   ", data_forest_plot$group))
  data_forest_plot <- data_forest_plot %>%
    mutate(across(c(Total, TP, FP, FN, TN), ~ifelse(is.na(.), "", .))) %>% 
    mutate(across(c(sen_est1, sen_low1, sen_hi1, spe_est2, spe_low2, spe_hi2), ~ .x * 100))
  
  data_forest_plot$Sensitivity <- paste(rep(" ", 10), collapse = " ")
  data_forest_plot$Specificity <- paste(rep(" ", 10), collapse = " ")
  
  data_forest_plot$`Sensitivity (95%CI)` <- ifelse(is.na(data_forest_plot$sen_est1),"",
                                                   sprintf("%.1f (%.1f to %.1f)",
                                                           data_forest_plot$sen_est1,
                                                           data_forest_plot$sen_low1,
                                                           data_forest_plot$sen_hi1))
  
  data_forest_plot$`Specificity (95%CI)` <- ifelse(is.na(data_forest_plot$spe_est2),"",
                                                   sprintf("%.1f (%.1f to %.1f)",
                                                           data_forest_plot$spe_est2,
                                                           data_forest_plot$spe_low2,
                                                           data_forest_plot$spe_hi2))
  data_forest_plot <- data_forest_plot %>% 
    dplyr::select(group:TN,Sensitivity,`Sensitivity (95%CI)`,Specificity,`Specificity (95%CI)`,sen_est1:spe_hi2)
  
  tm <- forestploter::forest_theme(base_size = 10,
                                   refline_gp = grid::gpar(lwd = 1, lty = "dashed", col = "red"),
                                   arrow_type = "closed",
                                   footnote_gp = grid::gpar(col = "skyblue"))
  
  p <- forestploter::forest(data_forest_plot[,c(1:6,7:10)], 
                            est=list(data_forest_plot$sen_est1, 
                                     data_forest_plot$spe_est2),
                            lower=list(data_forest_plot$sen_low1,
                                       data_forest_plot$spe_low2), 
                            upper=list(data_forest_plot$sen_hi1,
                                       data_forest_plot$spe_hi2),
                            sizes= list(as.numeric(data_forest_plot$Total)/3777,
                                        as.numeric(data_forest_plot$Total)/3777),
                            ci_column  = c(7,9),
                            ref_line  = c(0.9,0.95),
                            xlim = list(c(0.5, 1.0),
                                        c(0.6, 1.0)),
                            ticks_at=list(c(0.5, 0.8,1.0),
                                          c(0.6, 0.8,0.9,1.0)),
                            theme  =  tm) ;p
  
  pdf(sprintf("%s/%s_forest_plot_subgroup_%s.pdf",save_path,saveDir,indicators_used))
  print(p)
  dev.off()
}

### eFigure 7 plot ----
data_age <- read.csv(sprintf("%s/%s_age_AUROC.csv",save_path_AUROC,saveDir_AUROC))
data_sex <- read.csv(sprintf("%s/%s_sex_AUROC.csv",save_path_AUROC,saveDir_AUROC))
data_subgroup <- rbind(data_age,data_sex)

suppressMessages(library(forestploter))

data_forest <- data_subgroup %>% 
  select(Indicator,group,AUC:format)

custom_order <- c("Male", "Female", "less than 45", "more tham 45")
ordered_indices <- order(match(data_forest$group, custom_order))
sorted_data <- data_forest[ordered_indices, ]
data_forest <- sorted_data %>% as_tibble() %>% 
  add_row(group="Sex",.before = 1) %>% 
  add_row(group="Age,yr",.before = 10) %>% 
  dplyr::rename("AUC (95%CI)"="format") %>% 
  mutate(Total=case_when(
    group=="Male"~ dim(dplyr::filter(data_combine,gender=="Male"))[1],
    group=="Female"~ dim(dplyr::filter(data_combine,gender=="Female"))[1],
    group=="less than 45"~ dim(dplyr::filter(data_age_calculate,age_cut=="less than 45"))[1],
    group=="more tham 45"~ dim(dplyr::filter(data_age_calculate,age_cut=="more tham 45"))[1],
    TRUE~NA
  ))

writexl::write_xlsx(data_forest,
                    sprintf("%s/%s_Combine_for_forest_AUC.xlsx",save_path_AUROC,saveDir_AUROC))

data_forest <- data_forest %>% 
  # tidyr::drop_na() %>% 
  # distinct() %>% 
  mutate(Total=case_when(
    group=="Male"~ dim(dplyr::filter(data_combine,gender=="Male"))[1],
    group=="Female"~ dim(dplyr::filter(data_combine,gender=="Female"))[1],
    group=="less than 45"~ dim(dplyr::filter(data_age_calculate,age_cut=="less than 45"))[1],
    group=="more tham 45"~ dim(dplyr::filter(data_age_calculate,age_cut=="more tham 45"))[1],
    TRUE~NA
  )) %>% 
  mutate(connect_key=str_c(Indicator,group,sep = "_")) %>% 
  dplyr::select(connect_key,Total)

# data_forest <- data_forest %>% 
#   mutate(connect_key=str_c(Indicator,group,sep = "_")) %>% 
#   left_join(totol_num,by="connect_key") %>% 
#   dplyr::select(-connect_key)

data_forest$group <- ifelse(is.na(data_forest$Total), 
                            data_forest$group,
                            paste0("   ", data_forest$group))

data_forest <- data_forest %>%
  mutate(across(c(Total,`AUC (95%CI)`), ~ifelse(is.na(.), "", .)))

data_forest$AUC_ci  <- paste(rep(" ", 10), collapse = " ")

tm <- forestploter::forest_theme(base_size = 10,
                                 refline_gp = grid::gpar(lwd = 1, lty = "dashed", col = "red"),
                                 arrow_type = "closed",
                                 footnote_gp = grid::gpar(col = "skyblue"))

p <- forest(data_forest[,c(2,7,8,6)],
            est=data_forest$AUC,
            lower=data_forest$low, 
            upper=data_forest$hi, 
            ci_column  = 3, 
            ref_line  = 0.95,
            xlim = c(0.70, 1.0),
            ticks_at=c(0.70,0.80,0.9,1.0),
            theme  =  tm) ;p
pdf(sprintf("%s/%s_forest_plot_subgroup_AUC.pdf",save_path_AUROC,saveDir_AUROC))
print(p)
dev.off()

## eFigure 9 ----
saveDir <- "eFigure9"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_kappa <- data %>% 
  as.data.frame() %>% 
  tidyr::drop_na(P85_plasma) 

library(ggpubr)
ggscatter(data_kappa, x = "P85_SCO", 
          y = "P85_plasma_SCO",
          size = 1.5,
          color = "#0072B2",
          add = "reg.line", 
          add.params = list(color = "black", fill = "grey50", linewidth = 1),  
          conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = 300, label.y = 100, label.sep = "\n", size=3) +
  xlab("Serum P85-Ab COI") +
  ylab("Plasam P85-Ab COI")+
  theme_classic()+theme(aspect.ratio = 1)

ggsave(sprintf("%s/%s_Correlation_plasma_serum.pdf",save_path,saveDir),
       width = 6,height = 6)

suppressMessages(library(irr))
suppressMessages(library(statpsych))

data_kappa$P85_plasma <- factor(data_kappa$P85_plasma,levels = c("(+)","(-)"),labels = c(1,0))

#### with CI version 
indicators <- c("P85", "P85_plasma")

kappa_matrix <- data.frame(matrix(ncol = length(indicators), 
                                  nrow = length(indicators)))
colnames(kappa_matrix) <- c(str_c(indicators[1],"(+)"),str_c(indicators[1],"(-)"))
rownames(kappa_matrix) <- c(str_c(indicators[2],"(+)"),str_c(indicators[2],"(-)"))

rater1 <- data_kappa[[indicators[1]]]
rater2 <- data_kappa[[indicators[2]]]

f00 <- sum(rater1 == 0 & rater2 == 0)
f01 <- sum(rater1 == 0 & rater2 == 1)
f10 <- sum(rater1 == 1 & rater2 == 0)
f11 <- sum(rater1 == 1 & rater2 == 1)

ci_result <- statpsych::ci.kappa(0.05, f00, f01, f10, f11)
ci_result

kappa_with_ci <- sprintf("%s\n(%s - %s)", 
                         signif(ci_result[2, 1],2), 
                         signif(ci_result[2, 3],2), 
                         signif(ci_result[2, 4],2))
kappa_with_ci
write.csv(as.data.frame(kappa_with_ci),
          file = sprintf("%s/%s_Kappa_Raw_Data.csv",save_path,saveDir))

kappa_matrix[1,1] <- f11
kappa_matrix[1,2] <- f01
kappa_matrix[2,1] <- f10
kappa_matrix[2,2] <- f00
write.csv(as.data.frame(kappa_matrix),
          file = sprintf("%s/%s_Data_2x2_table.csv",save_path,saveDir))

## eFigure 10 ----
saveDir <- "eFigure10"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

perf_table <- do.call(rbind, lapply(names(perf_list), function(nm){
  cbind(
    model = nm,
    t(perf_list[[nm]]),
    threshold = thr_list[[nm]]
  )
}))

perf_table <- as.data.frame(perf_table)
write.csv(perf_table,sprintf("%s/%s_apparent_performance_5centres.csv",save_path,saveDir),
          row.names = F)

### eFigure 10A ----
roc_p85 <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$P85_SCO),
                     direction="<",
                     levels = c(0, 1))
value_p85 <- print_auc_and_ci(roc_p85);value_p85

## P85-Ab+vca-IgA
roc_vca <- rocs[["gam_P85_VCA"]]$roc
value_vca <- print_auc_and_ci(roc_vca);value_vca

## P85-Ab+EBNA1-IgA 
roc_na1 <- rocs[["gam_P85_EBNA1"]]$roc
value_na1 <- print_auc_and_ci(roc_na1);value_na1

# VCA-IgA+EBNA1-IgA 
roc_db <- rocs[["dual_method"]]$roc
value_db <- print_auc_and_ci(roc_db);value_db

## P85-Ab+VCA-IgA+EBNA1-IgA 
roc_P85_db <- rocs[["gam_P85_dual"]]$roc
value_P85_db <- print_auc_and_ci(roc_P85_db);value_P85_db

pdf(sprintf("%s/%sA_Model_combination.pdf",save_path,saveDir),
    width = 6,height = 6)
{
  plot(roc_p85, 
       col=my_color$model["P85-Ab"],
       legacy.axes=TRUE)
  
  plot.roc(roc_db,
           add=TRUE,
           col=my_color$model["VCA_EBNA1"],
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  plot.roc(roc_vca,
           add=TRUE,
           col=my_color$model["P85_VCA"],
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  plot.roc(roc_na1,
           add=TRUE,
           col=my_color$model["P85_EBNA1"],
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  
  plot.roc(roc_P85_db,
           add=TRUE,
           col=my_color$model["P85_VCA_EBNA1"],
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  text<-c(paste0("P85-Ab\n",value_p85),
          paste0("VCA-IgA+EBNA1-IgA\n",value_db),
          paste0("P85-Ab+VCA-IgA\n",value_vca),
          paste0("P85-Ab+EBNA1-IgA\n",value_na1),
          paste0("P85-Ab+VCA-IgA+EBNA1-IgA\n",value_P85_db)
  )
  
  legend(0.8, 0.5,
         bty = "n",
         legend=text, 
         col=my_color$model,
         lwd=3) 
}
dev.off()


roc_names <- c("roc_p85", "roc_vca", "roc_na1", "roc_db", "roc_P85_db")

results <- data.frame(
  comparisons = character(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

comparisons <- combn(roc_names, 2, simplify = FALSE)

for (i in seq_along(comparisons)) {
  
  roc1 <- get(comparisons[[i]][1])
  roc2 <- get(comparisons[[i]][2])
  
  test_result <- roc.test(roc1, roc2)
  
  results <- rbind(results, data.frame(
    comparisons = paste(comparisons[[i]][1], "vs", comparisons[[i]][2]),
    p_value = test_result$p.value
  ))
}

results$p.adjust.raw <- p.adjust(results$p_value, method = correct_method)
results$p.adjust <- format_p_special(results$p.adjust.raw)
custom_order <- c("roc_p85 vs roc_db","roc_p85 vs roc_vca","roc_p85 vs roc_na1","roc_p85 vs roc_P85_db",
                  "roc_vca vs roc_db","roc_na1 vs roc_db","roc_db vs roc_P85_db",
                  "roc_vca vs roc_na1","roc_vca vs roc_P85_db",
                  "roc_na1 vs roc_P85_db")
results$comparisons <- factor(results$comparisons, levels = custom_order)
results <- results[order(results$comparisons), ]
write.csv(results,sprintf("%s/%sA_Model_combination_delong_comprisons.csv",save_path,saveDir))

### eFigure 10B ----
pdf(sprintf("%s/%sB_CalibrationCurves_VCA_EBNA1.pdf",save_path,saveDir),
    width = 5.5,height = 6)
CalibrationCurves::val.prob.ci.2(p=data_logistic$dual_method_prob, 
                                 y=data_logistic$group_ref,
                                 allowPerfectPredictions = T,roundstats = 2)
dev.off()

pdf(sprintf("%s/%sB_CalibrationCurves_P85_VCA.pdf",save_path,saveDir),
    width = 5.5,height = 6)
CalibrationCurves::val.prob.ci.2(p=data_logistic$gam_P85_VCA_prob, 
                                 y=data_logistic$group_ref,
                                 allowPerfectPredictions = T,roundstats = 3)
dev.off()

pdf(sprintf("%s/%sB_CalibrationCurves_P85_EBNA1.pdf",save_path,saveDir),
    width = 5.5,height = 6)
CalibrationCurves::val.prob.ci.2(p=data_logistic$gam_P85_EBNA1_prob,
                                 y=data_logistic$group_ref,
                                 allowPerfectPredictions = T,roundstats = 3)
dev.off()

pdf(sprintf("%s/%sB_CalibrationCurves_P85_VCA_EBNA1.pdf",save_path,saveDir),
    width = 5.5,height = 6)
CalibrationCurves::val.prob.ci.2(p=data_logistic$gam_P85_dual_prob, 
                                 y=data_logistic$group_ref,
                                 allowPerfectPredictions = T,roundstats = 3)
dev.off()

## eFigure 11, eFigure 13-14 ----
saveDir <- "eFigure11_13_14"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 
data_combine <- data
centers <- unique(data_combine$center)

### Data calculation ----
data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(NPC_stage %in% c("I","II","III","IVA","IVB"))
dim(data_pre)

stages <- c("I","II","III","IVA","IVB")
unique_combinations <- unique(data_pre[, c("NPC_stage", "center")])

result_list <- list()
result_statistical_list <- list()

for (i in 1:nrow(unique_combinations)) {
  current_NPC_stage <- unique_combinations[i, "NPC_stage"]
  current_center <- unique_combinations[i, "center"]
  print(str_c(current_NPC_stage,"_",current_center))
  
  subset_data <- data_pre[data_pre$NPC_stage == current_NPC_stage & 
                            data_pre$center == current_center, ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df$Stage <- current_NPC_stage
  result_df$Center <- current_center
  
  result_df <- result_df %>% 
    select(Indicator,Center,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$Center <- current_center
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/%s_Individual_centers_NPC_5_stages_combined.csv",save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,I,II,III,IVA,IVB)

write.csv(merged.data,
          sprintf("%s/%s_Individual_centers_NPC_5_stages.csv",save_path,saveDir),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_Individual_centers_NPC_5_stages_compare_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols =-c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare,Center,I,II,III,IVA,IVB)

merged.statistical.data <- merged.statistical.data[order(merged.statistical.data$Center), ]

write.csv(merged.statistical.data,
          sprintf("%s/%s_Individual_centers_NPC_5_stages_compare.csv",save_path,saveDir),
          row.names = F)

#### 3 stages 
data_pre <- data %>% 
  as.data.frame() %>%
  filter(group_ref == 1) %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB")) %>% 
  mutate(NPC_stage = case_when(
    NPC_stage %in% c("I","II") ~ "Early stage",
    NPC_stage %in% c("III","IVA") ~ "Localregionally advanced stage",
    NPC_stage %in% c("IVB") ~ "Advanced stage"
  ))

unique_combinations <- unique(data_pre[, c("NPC_stage", "center")])

result_list <- list()
result_statistical_list <- list()

for (i in 1:nrow(unique_combinations)) {
  current_NPC_stage <- unique_combinations[i, "NPC_stage"]
  current_center <- unique_combinations[i, "center"]
  print(str_c(current_NPC_stage,"_",current_center))
  
  subset_data <- data_pre[data_pre$NPC_stage == current_NPC_stage & 
                            data_pre$center == current_center, ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df$Stage <- current_NPC_stage
  result_df$Center <- current_center
  
  result_df <- result_df %>% 
    select(Indicator,Center,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$Center <- current_center
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/%s_Individual_centers_NPC_3_stages_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.data,
          sprintf("%s/%s_Individual_centers_NPC_3_stages.csv",
                  save_path,saveDir),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_Individual_centers_NPC_3_stages_compare_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"a"),
  )) %>% 
  pivot_wider(id_cols =-c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare ,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

merged.statistical.data <- merged.statistical.data[order(merged.statistical.data$Center), ]

write.csv(merged.statistical.data,
          sprintf("%s/%s_Individual_centers_NPC_3_stages_compare.csv",
                  save_path,saveDir),
          row.names = F)

### plot ----
data_forest_original <- read.csv(sprintf("%s/%s_Individual_centers_NPC_3_stages.csv",
                                         save_path,saveDir)) %>% 
  as.data.frame()

data_forest_stat <- read.csv(sprintf("%s/%s_Individual_centers_NPC_3_stages_compare.csv",
                                     save_path,saveDir)) %>% 
  as.data.frame() %>% 
  mutate(Indicator=str_extract(.$compare,"(?<=vs ).*")) %>% 
  dplyr::select(Center:Indicator) %>% 
  dplyr::mutate(
    p_early=Early.stage,
    p_lr=Localregionally.advanced.stage,
    p_ad=Advanced.stage) 

data_forest <- left_join(data_forest_original,
                         dplyr::select(data_forest_stat,Center,Indicator:p_ad),
                         by=c("Indicator","Center")) %>% 
  mutate(
    p_early = ifelse(Indicator == "P85" & is.na(p_early), "ref", p_early),
    p_lr    = ifelse(Indicator == "P85" & is.na(p_lr),    "ref", p_lr),
    p_ad    = ifelse(Indicator == "P85" & is.na(p_ad),    "ref", p_ad)
  )

data_forest$Center <- factor(data_forest$Center,
                             levels = c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH"))

three_stages <- c("Early.stage","Localregionally.advanced.stage","Advanced.stage")
p_three <- c("p_early" ,"p_lr" ,"p_ad")

for(i in 1:length(three_stages)){
  print(three_stages[i])
  data_sen <- data_forest %>% 
    dplyr::select(Indicator,Center,!!sym(three_stages[i]),!!sym(p_three[i]))
  
  colnames(data_sen) <- c("Indicator","Center","sens","pvalue")
  
  data_sen <- data_sen %>% 
    mutate(
      est=str_extract(.$sens,".*(?=\n)"),
      low_ci=str_extract(.$sens,"(?<=\\().*(?= to)"),
      high_ci=str_extract(.$sens,"(?<=to ).*(?=\\))")) %>% 
    dplyr::filter(Indicator%in%c("P85","VCA_IgA","EA_IgA","NA1_IgA")) %>% 
    mutate(
      est = as.numeric(est), 
      low_ci = as.numeric(low_ci), 
      high_ci = as.numeric(high_ci) 
    )
  
  data_sen %>% 
    mutate(Indicator=factor(.$Indicator,levels = c("P85","VCA_IgA","EA_IgA","NA1_IgA"),
                            labels = c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA"))) %>% 
    mutate(label=str_c(ifelse(est==100,"100",sprintf("%0.1f",est)),
                       "\n",pvalue)) %>% 
    tidyplot(x = Indicator, y=est,color = Indicator,dodge_width = 0.6) %>%
    add_mean_dot(size = 1) %>% 
    add(ggplot2::geom_errorbar(aes(ymin =low_ci, ymax = high_ci,
                                   group=Indicator),width=0.5,
                               position = position_dodge(0.6))) %>% 
    remove_x_axis_title() %>% 
    add(ggplot2::theme(axis.text.x = element_text(angle=60,hjust = 1,vjust = 1))) %>% 
    adjust_colors(new_colors = my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")]) %>%  # "#4DBBD5FF","#00A087FF","#3C5488FF"
    adjust_y_axis(limits=c(0,100)) %>% 
    add(ggplot2::geom_label(aes(y=low_ci-15,label = label),fill="transparent",label.size=NA,size=3)) %>%
    adjust_y_axis_title("Sensitivity (%)") %>% 
    remove_x_axis_labels() %>% 
    split_plot(by = Center, ncol = 2, nrow = 3)
  
  ggsave(sprintf("%s/%s_Sensitivity_%s_across_centers.pdf",
                 save_path,saveDir,three_stages[i]),
         width = 4,height = 12)
}



# eTable ----
## eTable 2 ----
saveDir <- "eTable2"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA")
data_pre <- data %>% as.data.frame()
library(lme4)
library(broom.mixed)
fit_P85_GLMM <- glmer(group_ref ~ P85 + (1|center), data = data_pre, family = binomial)
fit_VCA_GLMM <- glmer(group_ref ~ VCA_IgA + (1|center), data = data_pre, family = binomial)
fit_EA_GLMM <- glmer(group_ref ~ EA_IgA + (1|center), data = data_pre, family = binomial)
fit_EBNA_GLMM <- glmer(group_ref ~ NA1_IgA + (1|center), data = data_pre, family = binomial)


ranef_P85 <- ranef(fit_P85_GLMM)$center
ranef_VCA <- ranef(fit_VCA_GLMM)$center
ranef_EA <- ranef(fit_EA_GLMM)$center
ranef_EBNA <- ranef(fit_EBNA_GLMM)$center

var_P85 <- VarCorr(fit_P85_GLMM)$center[1]
var_VCA <- VarCorr(fit_VCA_GLMM)$center[1]
var_EA <- VarCorr(fit_EA_GLMM)$center[1]
var_EBNA <- VarCorr(fit_EBNA_GLMM)$center[1]

results <- bind_rows(
  tidy(fit_P85_GLMM, effects = "fixed")  %>% mutate(Model = "P85", Random_Effect_Variance = var_P85),
  tidy(fit_VCA_GLMM, effects = "fixed")  %>% mutate(Model = "VCA_IgA", Random_Effect_Variance = var_VCA),
  tidy(fit_EA_GLMM, effects = "fixed")   %>% mutate(Model = "EA_IgA", Random_Effect_Variance = var_EA),
  tidy(fit_EBNA_GLMM, effects = "fixed") %>% mutate(Model = "EBNA1_IgA", Random_Effect_Variance = var_EBNA)
)

results_summary <- results %>%
  dplyr::select(Model, term, estimate, std.error, statistic, p.value, Random_Effect_Variance)


write.csv(results_summary,
          sprintf("%s/%s_GLMM_results_summary.csv",save_path,saveDir))


#### meta between-center Heterogeneity
library(mada)
library(broom)

data_forest_center <- read.csv(sprintf("%s/Figure2/Figure2_Combined_Center_Raw_Data.csv",save.path)) %>% 
  as.data.frame()

fit_p85   <- reitsma(data_forest_center[data_forest_center$Indicator=="P85", c("TP","FN","FP","TN")])
fit_vca   <- reitsma(data_forest_center[data_forest_center$Indicator=="VCA_IgA", c("TP","FN","FP","TN")])
fit_ea    <- reitsma(data_forest_center[data_forest_center$Indicator=="EA_IgA", c("TP","FN","FP","TN")])
fit_ebna1 <- reitsma(data_forest_center[data_forest_center$Indicator=="NA1_IgA", c("TP","FN","FP","TN")])

s <- summary(fit_p85);s
s1 <- summary(fit_vca);s1
s2 <- summary(fit_ea);s2
s3 <- summary(fit_ebna1);s3

extract_results <- function(fit, name){
  s <- summary(fit)
  tibble(
    Marker = name,
    Sensitivity  = paste0(signif(s$coefficients[3,1]*100,3), "% (", 
                          signif(s$coefficients[3,5]*100,3), " to ", 
                          signif(s$coefficients[3,6]*100,3), ")"),
    Specificity  = paste0(signif((100 - s$coefficients[4,1]*100),3), "% (", 
                          signif((100 - s$coefficients[4,6]*100),3), " to ", 
                          signif((100 - s$coefficients[4,5]*100),3), ")"),
    AUC          = signif(s$AUC$AUC,2),
    pAUC         = signif(s$AUC$pAUC,2),
    I2_Zhou      = paste0(signif(s$i2[1]*100,3), "%"),
    `I2(Holling sample size unadjusted approaches)` = paste0(round(min(s$i2[2],s$i2[3],s$i2[4])*100,1), "% to ", round(max(s$i2[2],s$i2[3],s$i2[4])*100,1), "%"),
    `I2(Holling sample size adjusted approaches)`= paste0(round(min(s$i2[5],s$i2[6],s$i2[7])*100,1), "% to ", round(max(s$i2[5],s$i2[6],s$i2[7])*100,1), "%")
  )
}

results_table <- bind_rows(
  extract_results(fit_p85, "P85"),
  extract_results(fit_vca, "VCA_IgA"),
  extract_results(fit_ea, "EA_IgA"),
  extract_results(fit_ebna1, "NA1_IgA")
)

results_table

write.csv(results_table,sprintf("%s/%s_Meta_combined_across_centers.csv",save_path,saveDir))


##  eTable 3 ----
saveDir <- "eTable3"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA")

data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(group_ref == 0)

stages <- unique(data_pre$category_green)
result_list <- list()
result_statistical_list <- list()

for (i in 1:length(stages)) {
  current_NPC_stage <- stages[i]
  print(current_NPC_stage)
  subset_data <- data_pre[data_pre$category_green ==  current_NPC_stage , ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df
  
  result_df$Class <- current_NPC_stage
  
  result_df <- result_df %>% 
    dplyr::select(Indicator,Class,Specificity,FP,TN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare <- result_compare %>% dplyr::select(compare,spe1:spe.p.value)
  
  result_compare$Class <- current_NPC_stage
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = "holm") #bonferroni
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Class,spe1:spe.p.value,spec.p.adjust.raw:spec.p.adjust)
  
  result_statistical_list[[i]] <- result_compare
}

merged.data.combined <- map_dfr(result_list,bind_rows) 
write.csv(merged.data.combined,
          sprintf("%s/%s_pooled_analysis_Control_combined.csv",save_path,saveDir),row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=FP+TN) %>% 
  pivot_wider(id_cols =-c(sum,FP:TN),
              names_from = Class,
              values_from = Specificity) %>% 
  dplyr::select(Indicator,Lymphoma,`EBV-related Benign Disease`,
                `Non-NPC Head and Neck Cancer`,`Other Malignancy`,`Other Benign Disease`)
write.csv(merged.data,
          sprintf("%s/%s_pooled_analysis_Control.csv",save_path,saveDir),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,
                                            bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_pooled_analysis_Control_compare_combined.csv",save_path,saveDir),row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    spe_method=="McNemar"~spec.p.adjust,
    spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
  )) %>%
  pivot_wider(id_cols = -c(spe1:spec.p.adjust),
              names_from = Class,
              values_from = p_value) %>% 
  dplyr::select(compare,Lymphoma,`EBV-related Benign Disease`,
                `Non-NPC Head and Neck Cancer`,`Other Malignancy`,`Other Benign Disease`)
write.csv(merged.statistical.data,
          sprintf("%s/%s_pooled_analysis_Control_compare.csv",save_path,saveDir),
          row.names = F)

#### combine EBV-related and unrelated diseases
data_pre <- data %>% 
  as.data.frame() %>% 
  filter(group_ref == 0)

stages <- unique(data_pre$category_orange)
result_list <- list()
result_statistical_list <- list()

for (i in 1:length(stages)) {
  current_NPC_stage <- stages[i]
  print(current_NPC_stage)
  
  subset_data <- data_pre[data_pre$category_orange ==  current_NPC_stage , ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df
  
  result_df$Class <- current_NPC_stage
  
  result_df <- result_df %>% 
    dplyr::select(Indicator,Class,Specificity,FP,TN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare <- result_compare %>% dplyr::select(compare,spe1:spe.p.value)
  
  result_compare$Class <- current_NPC_stage
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method) 
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Class,spe1:spe.p.value,spec.p.adjust.raw:spec.p.adjust)
  
  result_statistical_list[[i]] <- result_compare
}

merged.data.combined <- map_dfr(result_list,bind_rows) 
write.csv(merged.data.combined,
          sprintf("%s/%s_pooled_analysis_Control_combined_EBV.csv",save_path,saveDir),row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=FP+TN) %>% 
  pivot_wider(id_cols =-c(sum,FP:TN),
              names_from = Class,
              values_from = Specificity)
write.csv(merged.data,
          sprintf("%s/%s_pooled_analysis_Control_EBV.csv",save_path,saveDir),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,
                                            bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_pooled_analysis_Control_compare_combined_EBV.csv",save_path,saveDir),row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    spe_method=="McNemar"~spec.p.adjust,
    spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
  )) %>%
  pivot_wider(id_cols = -c(spe1:spec.p.adjust),
              names_from = Class,
              values_from = p_value) 

write.csv(merged.statistical.data,
          sprintf("%s/%s_pooled_analysis_Control_compare_EBV.csv",
                  save_path,saveDir),row.names = F)


##### Combined for publication & verification 
### minor
merged.statistical.data <- read.csv(sprintf("%s/%s_pooled_analysis_Control_compare.csv",
                                            save_path,saveDir))
merged.data <- read.csv(sprintf("%s/%s_pooled_analysis_Control.csv",
                                save_path,saveDir))

value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator")

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part1 <- merged2 %>%
  dplyr::select(Indicator, all_of(value_cols))

#### major
merged.statistical.data <- read.csv(sprintf("%s/%s_pooled_analysis_Control_compare_EBV.csv",
                                            save_path,saveDir))
merged.data <- read.csv(sprintf("%s/%s_pooled_analysis_Control_EBV.csv",
                                save_path,saveDir))

value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator")

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part2 <- merged2 %>%
  dplyr::select(Indicator, all_of(value_cols))

merged.final <- left_join(merged.final.part1,merged.final.part2) %>% 
  dplyr::select(Indicator,Lymphoma,EBV.related.Benign.Disease,EBV.related.Diseases,
                Non.NPC.Head.and.Neck.Cancer,Other.Malignancy,Other.Benign.Disease,EBV.unrelated.Diseases)
merged.final
write.csv(merged.final,sprintf("%s/%s_article_Combined.csv",save_path,saveDir))

## eTable 4 & eTable 5 ----
saveDir <- "eTable4_5"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_IECV <- data_logistic
clusters <- unique(data_IECV$center)
N.clust <- length(clusters) 
data.in <- data.leftout <- list()
results <- list()

for(i in 1:N.clust){
  data.in[[i]] <- data_IECV[data_IECV$center!=clusters[i],]
  data.leftout[[i]] <- data_IECV[data_IECV$center==clusters[i],]
}

### dual-antibody strategy ----
model <- "dual-antibody strategy"
model.fit <- fits[["dual_method"]]
leftout.prediction.logistic <- list()
leftout.performance.logistic <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- rms::lrm(group_ref~VCA_IgA_SCO + NA1_IgA_SCO,
                            data=data.in[[i]]) 
  
  leftout.prediction.logistic[[i]] <- predict(object = model.CV[[i]],
                                              data.leftout[[i]])
  roc <- pROC::roc(data.leftout[[i]]$group_ref, leftout.prediction.logistic[[i]])
  AUC_value <- pROC::auc(roc)[1]
  LCI_value <- pROC::ci(roc)[1]
  UCI_value <-pROC::ci(roc)[3]
  glm.model <- glm(data.leftout[[i]]$group_ref~leftout.prediction.logistic[[i]],
                   family = binomial)
  glm.fit <- summary(glm.model)
  calibration.intercept <- glm.fit$coef[1,1]
  calibration.intercept.SE <- glm.fit$coef[1,2]
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.SE <- glm.fit$coef[2,2]
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.logistic[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept.SE=calibration.intercept.SE,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
    calibration.slope.SE=calibration.slope.SE,
    calibration.slope_LCI=calibration.slope.LCI,
    calibration.slope_UCI=calibration.slope.UCI,
    cluster=unique(data.leftout[[i]]$center)
  )
}

desired_order <- c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH")
leftout.performance.logistic.agg <- map_dfr(leftout.performance.logistic,bind_rows) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)

leftout.prediction.logistic.agg <- do.call(c, leftout.prediction.logistic)
IECV.observed <- do.call(rbind, data.leftout)$group_ref
IECV.cluster <- do.call(rbind, data.leftout)$center
IECV.logistic <- calculate_performance2(IECV.observed, 
                                        leftout.prediction.logistic.agg)[-c(1:2)]
IECV.logistic

auc.clusters <- data.frame(auc=rep(NA,N.clust), SE=NA, LCI=NA,UCI=NA, cluster=NA)
for(i in 1:N.clust){
  d.cl <- data.leftout[[i]]
  roc1 <- pROC::roc(d.cl$group_ref, as.numeric(leftout.prediction.logistic[[i]]))
  auc.clusters$auc[i] <- pROC::auc(roc1)
  auc.clusters$LCI[i] <-pROC::ci(roc1)[1]
  auc.clusters$UCI[i] <-pROC::ci(roc1)[3]
  auc.clusters$SE[i] <-(pROC::ci(roc1)[3]-pROC::ci(roc1)[1])/3.92
  auc.clusters$cluster[i]<- as.character(clusters[i])
}

auc.clusters

library(meta) 
meta.AUC <- metagen(TE=auc, seTE=SE, studlab = cluster, data=auc.clusters) 
forestplot <- meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
                           rightcols = c("effect", "ci"),digits = 2,
                           leftlabs = c("Cluster", "AUC", "seTE"))

pdf(sprintf("%s/%s_%s_auc.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
             rightcols = c("effect", "ci"),digits = 2,
             leftlabs = c("Cluster", "AUC", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

meta.intercept <- metagen(TE=calibration.intercept, 
                          seTE=calibration.intercept.SE, 
                          studlab = cluster, 
                          data=leftout.performance.logistic.agg) 
forestplot <- meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration intercept", "seTE"))
pdf(sprintf("%s/%s_%s_intercept.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration intercept", "seTE"))
dev.off()

meta.slope <- metagen(TE=calibration.slope, 
                      seTE=calibration.slope.SE, 
                      studlab = cluster, 
                      data=leftout.performance.logistic.agg) 
forestplot <- meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration slope", "seTE"))
pdf(sprintf("%s/%s_%s_slope.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration slope", "seTE"))
dev.off()

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.2f (%0.2f to %0.2f)",
                     signif(auc,2),signif(LCI,2),signif(UCI,2))) %>% 
  left_join(leftout.performance.logistic.agg %>% 
              mutate(calibration.intercept=sprintf("%0.2f (%0.2f to %0.2f)",
                                                   signif(calibration.intercept,3),
                                                   signif(calibration.intercept_LCI,3),
                                                   signif(calibration.intercept_UCI,3)),
                     calibration.slope=sprintf("%0.2f (%0.2f to %0.2f)",
                                               signif(calibration.slope,3),
                                               signif(calibration.slope_LCI,3),
                                               signif(calibration.slope_UCI,3))) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

data.save
results[[model]] <- data.save

### ** P85-Ab+VCA-IgA ----
model <- "P85-Ab+VCA-IgA"
model.fit <- fits[["gam_P85_VCA"]]
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- mgcv::gam(model.fit$formula,optimizer="efs",
                             data = data.in[[i]], family = binomial)
  
  leftout.prediction.gam[[i]] <- mgcv::predict.gam(object = model.CV[[i]],
                                                   data.leftout[[i]])
  roc <- pROC::roc(data.leftout[[i]]$group_ref, leftout.prediction.gam[[i]])
  AUC_value <- pROC::auc(roc)[1]
  LCI_value <- pROC::ci(roc)[1]
  UCI_value <-pROC::ci(roc)[3]
  glm.model <- glm(data.leftout[[i]]$group_ref~leftout.prediction.gam[[i]],
                   family = binomial)
  glm.fit <- summary(glm.model)
  calibration.intercept <- glm.fit$coef[1,1]
  calibration.intercept.SE <- glm.fit$coef[1,2]
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.SE <- glm.fit$coef[2,2]
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept.SE=calibration.intercept.SE,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
    calibration.slope.SE=calibration.slope.SE,
    calibration.slope_LCI=calibration.slope.LCI,
    calibration.slope_UCI=calibration.slope.UCI,
    cluster=unique(data.leftout[[i]]$center)
  )
}

desired_order <- c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH")
leftout.performance.gam.agg <- map_dfr(leftout.performance.gam,bind_rows) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)

leftout.prediction.gam.agg <- do.call(c, leftout.prediction.gam)
IECV.observed <- do.call(rbind, data.leftout)$group_ref
IECV.cluster <- do.call(rbind, data.leftout)$center
IECV.gam <- calculate_performance2(IECV.observed, leftout.prediction.gam.agg)

auc.clusters <- data.frame(auc=rep(NA,N.clust), SE=NA, LCI=NA,UCI=NA, cluster=NA)
for(i in 1:N.clust){
  d.cl <- data.leftout[[i]]
  roc1 <- pROC::roc(d.cl$group_ref, as.numeric(leftout.prediction.gam[[i]]))
  auc.clusters$auc[i] <- pROC::auc(roc1)
  auc.clusters$LCI[i] <-pROC::ci(roc1)[1]
  auc.clusters$UCI[i] <-pROC::ci(roc1)[3]
  auc.clusters$SE[i] <-(pROC::ci(roc1)[3]-pROC::ci(roc1)[1])/3.92
  auc.clusters$cluster[i]<- as.character(clusters[i])
}

auc.clusters

library(meta) 
meta.AUC <- metagen(TE=auc, seTE=SE, studlab = cluster, data=auc.clusters) 
forestplot <- meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
                           rightcols = c("effect", "ci"),digits = 2,
                           leftlabs = c("Cluster", "AUC", "seTE"))

pdf(sprintf("%s/%s_%s_auc.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
             rightcols = c("effect", "ci"),digits = 2,
             leftlabs = c("Cluster", "AUC", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

meta.intercept <- metagen(TE=calibration.intercept, 
                          seTE=calibration.intercept.SE, 
                          studlab = cluster, 
                          data=leftout.performance.gam.agg) 
forestplot <- meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration intercept", "seTE"))
pdf(sprintf("%s/%s_%s_intercept.pdf",save_path,saveDir, model),width = 10,height=6)
meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration intercept", "seTE"))
dev.off()

meta.slope <- metagen(TE=calibration.slope, 
                      seTE=calibration.slope.SE, 
                      studlab = cluster, 
                      data=leftout.performance.gam.agg) 
forestplot <- meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration slope", "seTE"))
pdf(sprintf("%s/%s_%s_slope.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration slope", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.2f (%0.2f to %0.2f)",
                     signif(auc,2),signif(LCI,2),signif(UCI,2))) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.2f (%0.2f to %0.2f)",
                                                   signif(calibration.intercept,3),
                                                   signif(calibration.intercept_LCI,3),
                                                   signif(calibration.intercept_UCI,3)),
                     calibration.slope=sprintf("%0.2f (%0.2f to %0.2f)",
                                               signif(calibration.slope,3),
                                               signif(calibration.slope_LCI,3),
                                               signif(calibration.slope_UCI,3))) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)
data.save
results[[model]] <- data.save

### ** P85-Ab+EBNA1-IgA -----------------------------
model <- "P85-Ab+EBNA1-IgA"
model.fit <- fits[["gam_P85_EBNA1"]]
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,optimizer="efs",
                       data = data.in[[i]], family = binomial)
  
  leftout.prediction.gam[[i]] <- mgcv::predict.gam(object = model.CV[[i]],
                                                   data.leftout[[i]])
  roc <- pROC::roc(data.leftout[[i]]$group_ref, leftout.prediction.gam[[i]])
  AUC_value <- pROC::auc(roc)[1]
  LCI_value <- pROC::ci(roc)[1]
  UCI_value <-pROC::ci(roc)[3]
  glm.model <- glm(data.leftout[[i]]$group_ref~leftout.prediction.gam[[i]],
                   family = binomial)
  glm.fit <- summary(glm.model)
  calibration.intercept <- glm.fit$coef[1,1]
  calibration.intercept.SE <- glm.fit$coef[1,2]
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.SE <- glm.fit$coef[2,2]
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept.SE=calibration.intercept.SE,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
    calibration.slope.SE=calibration.slope.SE,
    calibration.slope_LCI=calibration.slope.LCI,
    calibration.slope_UCI=calibration.slope.UCI,
    cluster=unique(data.leftout[[i]]$center)
  )
  
}

desired_order <- c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH")
leftout.performance.gam.agg <- map_dfr(leftout.performance.gam,bind_rows) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)

leftout.prediction.gam.agg <- do.call(c, leftout.prediction.gam)
IECV.observed <- do.call(rbind, data.leftout)$group_ref
IECV.cluster <- do.call(rbind, data.leftout)$center
IECV.gam <- calculate_performance2(IECV.observed, leftout.prediction.gam.agg)
IECV.gam

auc.clusters <- data.frame(auc=rep(NA,N.clust), SE=NA, LCI=NA,UCI=NA, cluster=NA)
for(i in 1:N.clust){
  d.cl <- data.leftout[[i]]
  roc1 <- pROC::roc(d.cl$group_ref, as.numeric(leftout.prediction.gam[[i]]))
  auc.clusters$auc[i] <- pROC::auc(roc1)
  auc.clusters$LCI[i] <-pROC::ci(roc1)[1]
  auc.clusters$UCI[i] <-pROC::ci(roc1)[3]
  auc.clusters$SE[i] <-(pROC::ci(roc1)[3]-pROC::ci(roc1)[1])/3.92
  auc.clusters$cluster[i]<- as.character(clusters[i])
}

auc.clusters

library(meta) 
meta.AUC <- metagen(TE=auc, seTE=SE, studlab = cluster, data=auc.clusters) 
forestplot <- meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
                           rightcols = c("effect", "ci"),digits = 2,
                           leftlabs = c("Cluster", "AUC", "seTE"))
pdf(sprintf("%s/%s_%s_auc.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
             rightcols = c("effect", "ci"),digits = 2,
             leftlabs = c("Cluster", "AUC", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

meta.intercept <- metagen(TE=calibration.intercept, 
                          seTE=calibration.intercept.SE, 
                          studlab = cluster, 
                          data=leftout.performance.gam.agg) 
forestplot <- meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration intercept", "seTE"))
pdf(sprintf("%s/%s_%s_intercept.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration intercept", "seTE"))
dev.off()

meta.slope <- metagen(TE=calibration.slope, 
                      seTE=calibration.slope.SE, 
                      studlab = cluster, 
                      data=leftout.performance.gam.agg) 
forestplot <- meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration slope", "seTE"))
pdf(sprintf("%s/%s_%s_slope.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration slope", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.2f (%0.2f to %0.2f)",
                     signif(auc,2),signif(LCI,2),signif(UCI,2))) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.2f (%0.2f to %0.2f)",
                                                   signif(calibration.intercept,3),
                                                   signif(calibration.intercept_LCI,3),
                                                   signif(calibration.intercept_UCI,3)),
                     calibration.slope=sprintf("%0.2f (%0.2f to %0.2f)",
                                               signif(calibration.slope,3),
                                               signif(calibration.slope_LCI,3),
                                               signif(calibration.slope_UCI,3))) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)
data.save
results[[model]] <- data.save

### ** triplet-antibody strategy ----
model <- "triplet-antibody strategy"
model.fit <- fits[["gam_P85_dual"]]
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,optimizer="efs",
                       data = data.in[[i]], family = binomial)
  leftout.prediction.gam[[i]] <- mgcv::predict.gam(object = model.CV[[i]],
                                                   data.leftout[[i]])
  roc <- pROC::roc(data.leftout[[i]]$group_ref, leftout.prediction.gam[[i]])
  AUC_value <- pROC::auc(roc)[1]
  LCI_value <- pROC::ci(roc)[1]
  UCI_value <-pROC::ci(roc)[3]
  glm.model <- glm(data.leftout[[i]]$group_ref~leftout.prediction.gam[[i]],
                   family = binomial)
  glm.fit <- summary(glm.model)
  calibration.intercept <- glm.fit$coef[1,1]
  calibration.intercept.SE <- glm.fit$coef[1,2]
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.SE <- glm.fit$coef[2,2]
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept.SE=calibration.intercept.SE,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
    calibration.slope.SE=calibration.slope.SE,
    calibration.slope_LCI=calibration.slope.LCI,
    calibration.slope_UCI=calibration.slope.UCI,
    cluster=unique(data.leftout[[i]]$center)
  )
  
}

desired_order <- c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH")
leftout.performance.gam.agg <- map_dfr(leftout.performance.gam,bind_rows) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)

leftout.prediction.gam.agg <- do.call(c, leftout.prediction.gam)
IECV.observed <- do.call(rbind, data.leftout)$group_ref
IECV.cluster <- do.call(rbind, data.leftout)$center
IECV.gam <- calculate_performance2(IECV.observed, leftout.prediction.gam.agg)
IECV.gam

auc.clusters <- data.frame(auc=rep(NA,N.clust), SE=NA, LCI=NA,UCI=NA, cluster=NA)
for(i in 1:N.clust){
  d.cl <- data.leftout[[i]]
  roc1 <- pROC::roc(d.cl$group_ref, as.numeric(leftout.prediction.gam[[i]]))
  auc.clusters$auc[i] <- pROC::auc(roc1)
  auc.clusters$LCI[i] <-pROC::ci(roc1)[1]
  auc.clusters$UCI[i] <-pROC::ci(roc1)[3]
  auc.clusters$SE[i] <-(pROC::ci(roc1)[3]-pROC::ci(roc1)[1])/3.92
  auc.clusters$cluster[i]<- as.character(clusters[i])
}

auc.clusters

library(meta) 
meta.AUC <- metagen(TE=auc, seTE=SE, studlab = cluster, data=auc.clusters) 
forestplot <- meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
                           rightcols = c("effect", "ci"),digits = 2,
                           leftlabs = c("Cluster", "AUC", "seTE"))
pdf(sprintf("%s/%s_%s_auc.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.AUC, prediction = T, xlim=c(0.8,1), colgap.left="5mm", 
             rightcols = c("effect", "ci"),digits = 2,
             leftlabs = c("Cluster", "AUC", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

meta.intercept <- metagen(TE=calibration.intercept, 
                          seTE=calibration.intercept.SE, 
                          studlab = cluster, 
                          data=leftout.performance.gam.agg) 
forestplot <- meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration intercept", "seTE"))
pdf(sprintf("%s/%s_%s_intercept.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.intercept, prediction = T, xlim=c(-2,2), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration intercept", "seTE"))
dev.off()

meta.slope <- metagen(TE=calibration.slope, 
                      seTE=calibration.slope.SE, 
                      studlab = cluster, 
                      data=leftout.performance.gam.agg) 
forestplot <- meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
                           rightcols = c("effect", "ci"), digits = 2,
                           leftlabs = c("Cluster", "calibration slope", "seTE"))
pdf(sprintf("%s/%s_%s_slope.pdf",save_path,saveDir,model),width = 10,height=6)
meta::forest(meta.slope, prediction = T, xlim=c(-3,3), colgap.left="5mm", 
             rightcols = c("effect", "ci"), digits = 2,
             leftlabs = c("Cluster", "calibration slope", "seTE"))
dev.off()

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.2f (%0.2f to %0.2f)",
                     signif(auc,2),signif(LCI,2),signif(UCI,2))) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.2f (%0.2f to %0.2f)",
                                                   signif(calibration.intercept,3),
                                                   signif(calibration.intercept_LCI,3),
                                                   signif(calibration.intercept_UCI,3)),
                     calibration.slope=sprintf("%0.2f (%0.2f to %0.2f)",
                                               signif(calibration.slope,3),
                                               signif(calibration.slope_LCI,3),
                                               signif(calibration.slope_UCI,3))) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)
data.save
results[[model]] <- data.save

### ** Combined for publication & verification  ----
results_combine <- map_dfr(results,bind_rows)
write.csv(results_combine,sprintf("%s/%s_Article_combined.csv",save_path,saveDir))

## eTable 6 ----
saveDir <- "eTable6"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_pre <- data_logistic
indicators <- c("P85", "dual_method_cut", "gam_P85_VCA_cut","gam_P85_EBNA1_cut","gam_P85_dual_cut")
data_pre[indicators] <- lapply(data_pre[indicators], function(x) factor(x, levels = c("1", "0")))

result_df <- calcuate_diagnosis_index2(data=data_pre,
                                       indicators = indicators,
                                       binomMethod_use = "wilson",
                                       outcome = "group_ref")
result_df

write.csv(result_df,sprintf("%s/%s_Pooled_analysis_model.csv",save_path,saveDir),
          row.names = F)

result_compare <- DTC_pair(data=data_pre,
                           indicator1="P85",
                           compare_indicator=indicators[-1],
                           outcome="group_ref")
result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
result_compare$sens.p.adjust <- ifelse(result_compare$sen_method=="McNemar",
                                       result_compare$sens.p.adjust,str_c(result_compare$sens.p.adjust,"b"))
result_compare$spec.p.adjust <- ifelse(result_compare$spe_method=="McNemar",
                                       result_compare$spec.p.adjust,str_c(result_compare$spec.p.adjust,"b"))
result_compare
write.csv(result_compare,sprintf("%s/%s_Pooled_analysis_model_compare.csv",save_path,saveDir),
          row.names = F)


##### Combined for publication & verification 
merged.data <- read.csv(sprintf("%s/%s_Pooled_analysis_model.csv",save_path,saveDir)) %>% 
  dplyr::select(-c(PPV:negLR))
merged.statistical.data <- read.csv(sprintf("%s/%s_Pooled_analysis_model_compare.csv",save_path,saveDir)) %>% 
  dplyr::select(compare,sens.p.adjust,spec.p.adjust)

value_cols <- setdiff(colnames(merged.data), "Indicator")
value_cols

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged.final <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator") %>% 
  dplyr::select(-refGroup) %>% 
  mutate(across(everything(), ~ replace_na(., "ref")))

write.csv(merged.final,
          sprintf("%s/%s_Article_Combined.csv",save_path,saveDir),
          row.names = F)

## eTable 7 ----
saveDir <- "eTable7"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

symptom_data <- data_logistic

candidate_symptom <- c("Neck mass","Nasal Symptoms","Aural Symptoms",
                       "Cranial Symptoms","Ophthalmic Symptoms",
                       "Neurological Deficits")

result_list <- list()
result_statistical_list <- list()

for (i in 1:length(candidate_symptom)) {
  used_symptom <- candidate_symptom[i]
  print(used_symptom)
  data_sub <- symptom_data %>% dplyr::filter(!!sym(used_symptom)==1) %>% as.data.frame()
  
  result_df <- calcuate_diagnosis_index2(data=data_sub,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df$symptom <- used_symptom
  result_list[[i]] <- result_df
  
  ## compare_result
  result_compare <- DTC_pair(data=data_sub,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$symptom <-  used_symptom
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_statistical_list[[i]] <- result_compare
}

merged.data <- map_dfr(result_list,bind_rows)
writexl::write_xlsx(merged.data,
                    sprintf("%s/%s_Model_symptom_combined_individual_symptom.xlsx",
                            save_path,saveDir))

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
writexl::write_xlsx(merged.statistical.data.combined,
                    sprintf("%s/%s_Model_symptom_compare_combined_individual_symptom.xlsx",
                            save_path,saveDir))

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(sens_p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b")),
    spec_p_value=case_when(
      spe_method=="McNemar"~spec.p.adjust,
      spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
    )) 

writexl::write_xlsx(merged.statistical.data,
                    sprintf("%s/%s_Model_symptom_compare_individual_symptom.xlsx",
                            save_path,saveDir))

## eTable 8 ----
saveDir <- "eTable8"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA")

data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(NPC_stage %in% c("I","II","III","IVA","IVB","/"))
dim(data_pre)

stages <- c("I","II","III","IVA","IVB")
result_list <- list()
result_statistical_list <- list()

for (i in 1:length(stages)) {
  current_NPC_stage <- stages[i]
  
  subset_data <- data_pre[data_pre$NPC_stage ==  current_NPC_stage , ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df
  
  result_df$Stage <- current_NPC_stage
  
  result_df <- result_df %>% 
    dplyr::select(Indicator,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = "holm")
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=TP+FN) %>% 
  pivot_wider(id_cols =-c(sum,TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>%
  dplyr::select(Indicator,I,II,III,IVA,IVB)

write.csv(merged.data,
          sprintf("%s/%s_pooled_analysis_NPC_5_stages.csv",save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/%s_pooled_analysis_NPC_5_stages_compare.csv",save_path,saveDir),
          row.names = F)

#### for article
value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator")

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part1 <- merged2 %>%
  dplyr::select(Indicator, all_of(value_cols))


##### 3 stages 
data_pre <- data %>% 
  as.data.frame() %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB","/")) %>% 
  mutate(NPC_stage = case_when(
    NPC_stage %in% c("I","II") ~ "Early stage",
    NPC_stage %in% c("III","IVA") ~ "Localregionally advanced stage",
    NPC_stage %in% c("IVB") ~ "Advanced stage",
    NPC_stage %in% c("/") ~ "/"
  ))

stages <- c("Early stage","Localregionally advanced stage","Advanced stage")
result_list <- list()
result_statistical_list <- list()

for (i in 1:length(stages)) {
  current_NPC_stage <- stages[i]
  
  subset_data <- data_pre[data_pre$NPC_stage ==  current_NPC_stage , ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df
  
  result_df$Stage <- current_NPC_stage
  
  result_df <- result_df %>% 
    dplyr::select(Indicator,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = "holm")
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity)

write.csv(merged.data,
          sprintf("%s/%s_pooled_analysis_NPC_3_stages.csv",save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/%s_pooled_analysis_NPC_3_stages_compare.csv",save_path,saveDir),
          row.names = F)

#### for article
value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator")

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part2 <- merged2 %>%
  dplyr::select(Indicator, all_of(value_cols))

merged.final <- left_join(merged.final.part1,merged.final.part2) %>% 
  dplyr::select(Indicator,I:II,`Early stage`,III:IVA,`Localregionally advanced stage`,IVB,`Advanced stage`)

write.csv(merged.final, 
          sprintf("%s/%s_Article.csv",save_path,saveDir),
          row.names = F)

## eTable 9 ----
saveDir <- "eTable9"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_pre <- data_logistic %>% 
  as.data.frame() %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB","/")) %>% 
  mutate(NPC_stage2 = case_when(
    NPC_stage %in% c("I","II") ~ "Early stage",
    NPC_stage %in% c("III","IVA") ~ "Localregionally advanced stage",
    NPC_stage %in% c("IVB") ~ "Advanced stage",
    NPC_stage %in% c("/") ~ "/"
  ))
stages <- c("I","II","III","IVA","IVB")

indicators <- c("P85", "gam_P85_dual_cut")
result_list <- list()
result_statistical_list <- list()

for (i in 1:length(stages)) {
  current_NPC_stage <- stages[i]
  print(current_NPC_stage)
  
  subset_data <- data_pre[data_pre$NPC_stage ==  current_NPC_stage , ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  result_df
  
  result_df$Stage <- current_NPC_stage
  
  result_df <- result_df %>% 
    dplyr::select(Indicator,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/%s_pooled_analysis_Model_NPC_5_stages_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,I,II,III,IVA,IVB)
write.csv(merged.data,
          sprintf("%s/%s_pooled_analysis_Model_NPC_5_stages.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined ,
          sprintf("%s/%s_pooled_analysis_Model_NPC_5_stages_compare_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/%s_pooled_analysis_Model_NPC_5_stages_compare.csv",
                  save_path,saveDir),
          row.names = F)

#### for article
value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator")

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part1 <- merged2 %>%
  dplyr::select(Indicator, all_of(value_cols))


#### three-stage 
stages2 <- c("Early stage","Localregionally advanced stage","Advanced stage")
result_list <- list()
result_statistical_list <- list()

for (i in 1:length(stages2)) {
  current_NPC_stage <- stages2[i]
  print(current_NPC_stage)
  
  subset_data <- data_pre[data_pre$NPC_stage2 ==  current_NPC_stage , ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df
  
  result_df$Stage <- current_NPC_stage
  
  result_df <- result_df %>% 
    dplyr::select(Indicator,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/%s_pooled_analysis_Model_NPC_3_stages_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity)
write.csv(merged.data,
          sprintf("%s/%s_pooled_analysis_Model_NPC_3_stages.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_pooled_analysis_Model_NPC_3_stages_compare_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/%s_pooled_analysis_Model_NPC_3_stages_compare.csv",
                  save_path,saveDir),
          row.names = F)

#### for article
value_cols <- setdiff(colnames(merged.data), "Indicator")

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = "Indicator")

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part2 <- merged2 %>%
  dplyr::select(Indicator, all_of(value_cols))

merged.final <- left_join(merged.final.part1,merged.final.part2) %>% 
  dplyr::select(Indicator,I:II,`Early stage`,III:IVA,`Localregionally advanced stage`,IVB,`Advanced stage`)

write.csv(merged.final, 
          sprintf("%s/%s_Article.csv",
                  save_path,saveDir),
          row.names = F)

## eTable 10-13 ----
saveDir <- "eTable10_13"
save_path <- file.path(save.path,saveDir)
ensure_dir(save_path)
load("./my_color.rdata")

data_pre <- data %>% 
  dplyr::filter(center %in% c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))

### Overall ----
data_combine <- data_pre

centers <- unique(data_combine$center)
indicators <- c("P85", "EBV_DNA0") 

result_list <- list()
result_statistical_list <- list()
AUC_list <- list()
delong_list <- list()

for (i in 1:length(centers)) {
  center_name <- centers[i]
  print(center_name)
  
  data_pre <- data_combine %>% 
    filter(.$center == center_name) %>%
    as.data.frame()
  
  if(nrow(data_pre) == 0) next 
  
  result_df <- calcuate_diagnosis_index2(data=data_pre,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")
  
  result_df$Center <- center_name
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=data_pre,
                             indicator1="P85",
                             compare_indicator=indicators[-1],
                             outcome="group_ref")
  result_compare$Center <- center_name
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$spec.p.adjust.raw <- p.adjust(result_compare$spe.p.raw, method = correct_method)
  result_compare$ppv.p.adjust.raw <- p.adjust(result_compare$ppv.p.raw, method = correct_method)
  result_compare$npv.p.adjust.raw <- p.adjust(result_compare$npv.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare$spec.p.adjust <- format_p_special(result_compare$spec.p.adjust.raw)
  result_compare$ppv.p.adjust <- format_p_special(result_compare$ppv.p.adjust.raw)
  result_compare$npv.p.adjust <- format_p_special(result_compare$npv.p.adjust.raw)
  result_statistical_list[[i]] <- result_compare
  
  #### ROC
  if(T){
    roc_p85 <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$P85_SCO),direction="<",
                         levels = c(0, 1))
    value_p85 <- print_auc_and_ci(roc_p85)
    
    roc_ebv <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$log_EBV_DNA),
                         direction="<",
                         levels=c(0, 1))
    value_ebv <- print_auc_and_ci(roc_ebv)
    
    compare4 <- roc.test(roc_p85,roc_ebv)
    compare4_pvalue.raw <- compare4$p.value
    compare4_pvalue <- format_p_special(compare4$p.value)
    
    delong_df <- data.frame(
      center=center_name,
      Indicator = c("P85 vs EBV-DNA"),
      P_value.raw = c(compare4_pvalue.raw),
      P_value = c(as.character(compare4_pvalue)))
    
    delong_df$p.adjust.raw <- p.adjust(delong_df$P_value.raw, method = correct_method)
    delong_df$p.adjust <- format_p_special(delong_df$p.adjust.raw)
    delong_list[[i]] <- delong_df
    
    AUC_df <- data.frame(
      center=center_name,
      Indicator = c("P85","EBV-DNA"), # "EBV-DNA"
      AUC = c(value_p85,value_ebv) 
    )
    
    AUC_list[[i]] <- AUC_df 
  }
}

AUC.merge <- map_dfr(AUC_list,bind_rows)
write.csv(AUC.merge,
          sprintf("%s/%s_Data_AUC.merge.csv",save_path,saveDir),
          row.names = F)

delong.merge <- map_dfr(delong_list,bind_rows)
write.csv(delong.merge,
          sprintf("%s/%s_Data_delong.merge.csv",save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  dplyr::select(Indicator,Center,everything())

write.csv(merged.data,
          sprintf("%s/%s_P85_vs_EBVDNA_Raw_Data.csv",
                  save_path,saveDir),
          row.names = F)

merged.data.statistics.combined <- map_dfr(result_statistical_list,bind_rows) %>% 
  dplyr::select(compare,Center,everything())
merged.data.statistics.combined <- merged.data.statistics.combined[order(merged.data.statistics.combined$Center),]
merged.data.statistics.combined

write.csv(merged.data.statistics.combined,
          sprintf("%s/%s_P85_vs_EBVDNA_Raw_Data_statistics.csv",
                  save_path,saveDir),
          row.names = F)

merged.data.statistics <- map_dfr(result_statistical_list,bind_rows) %>% 
  dplyr::select(compare,Center,everything()) %>% 
  dplyr::mutate(sens_p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b")),
    spec_p_value=case_when(
      spe_method=="McNemar"~spec.p.adjust,
      spe_method=="McNemar mid-P"~str_c(spec.p.adjust,"b"),
    )) 
merged.data.statistics <- merged.data.statistics[order(merged.data.statistics$Center),]
merged.data.statistics
write.csv(merged.data.statistics,
          sprintf("%s/%s_P85_vs_EBVDNA_Raw_Data_statistics.csv",
                  save_path,saveDir),
          row.names = F)

### By Stages ----
indicators <- c("P85", "EBV_DNA0") 

data_pre <- data %>% 
  as.data.frame() %>%
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  filter(group_ref == 1) %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB")) 

centers <- unique(data_pre$center)

unique_combinations <- unique(data_pre[, c("NPC_stage", "center")])

result_list <- list()
result_statistical_list <- list()

for (i in 1:nrow(unique_combinations)) {
  current_NPC_stage <- unique_combinations[i, "NPC_stage"]
  current_center <- unique_combinations[i, "center"]
  print(str_c(current_NPC_stage,"_",current_center))
  
  subset_data <- data_pre[data_pre$NPC_stage == current_NPC_stage & 
                            data_pre$center == current_center, ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")  
  result_df$Stage <- current_NPC_stage
  result_df$Center <- current_center
  
  result_df <- result_df %>% 
    select(Indicator,Center,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=c("EBV_DNA0"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$Center <- current_center
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_5_stages_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,I,II,III,IVA,IVB)
write.csv(merged.data,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_5_stages.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_5_stages_compare_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols =-c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare,Center,I,II,III,IVA,IVB)

write.csv(merged.statistical.data,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_5_stages_compare.csv",
                  save_path,saveDir),
          row.names = F)

#### for article
value_cols <- setdiff(colnames(merged.data), c("Indicator","Center"))

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = c("Indicator","Center"))

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part1 <- merged2 %>%
  dplyr::select(Indicator,Center, all_of(value_cols))

#### 3 centres-3 stages 
data_pre <- data %>% 
  as.data.frame() %>%
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  filter(group_ref == 1) %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB")) %>% 
  mutate(NPC_stage = case_when(
    NPC_stage %in% c("I","II") ~ "Early stage",
    NPC_stage %in% c("III","IVA") ~ "Localregionally advanced stage",
    NPC_stage %in% c("IVB") ~ "Advanced stage"
  ))

indicators <- c("P85", "EBV_DNA0") 
centers <- unique(data_pre$center)

unique_combinations <- unique(data_pre[, c("NPC_stage", "center")])

result_list <- list()
result_statistical_list <- list()

for (i in 1:nrow(unique_combinations)) {
  current_NPC_stage <- unique_combinations[i, "NPC_stage"]
  current_center <- unique_combinations[i, "center"]
  print(str_c(current_NPC_stage,"_",current_center))
  
  subset_data <- data_pre[data_pre$NPC_stage == current_NPC_stage & 
                            data_pre$center == current_center, ]
  
  result_df <- calcuate_diagnosis_index2(data=subset_data,
                                         indicators = indicators,
                                         binomMethod_use = "wilson",
                                         outcome = "group_ref")  
  result_df$Stage <- current_NPC_stage
  result_df$Center <- current_center
  
  result_df <- result_df %>% 
    select(Indicator,Center,Stage,Sensitivity,TP,FN)
  
  result_list[[i]] <- result_df
  
  result_compare <- DTC_pair(data=subset_data,
                             indicator1="P85",
                             compare_indicator=c("EBV_DNA0"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$Center <- current_center
  result_compare$sens.p.adjust.raw <- p.adjust(result_compare$sen.p.raw, method = correct_method)
  result_compare$sens.p.adjust <- format_p_special(result_compare$sens.p.adjust.raw)
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value,sens.p.adjust.raw:sens.p.adjust)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows) 
write.csv(merged.data.combined,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_3_stages_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.data,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_3_stages.csv",
                  save_path,saveDir),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
write.csv(merged.statistical.data.combined,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_3_stages_compare_combined.csv",
                  save_path,saveDir),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sens.p.adjust,
    sen_method=="McNemar mid-P"~str_c(sens.p.adjust,"b"),
  )) %>% 
  pivot_wider(id_cols =-c(sen1:sens.p.adjust),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare ,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.statistical.data,
          sprintf("%s/%s_Individual_centers_Three_centers_NPC_3_stages_compare.csv",
                  save_path,saveDir),
          row.names = F)

#### for article
value_cols <- setdiff(colnames(merged.data), c("Indicator","Center"))

stat2 <- merged.statistical.data %>%
  separate(compare, into = c("refGroup", "Indicator"), sep = " vs ")

merged2 <- merged.data %>%
  dplyr::left_join(stat2, by = c("Indicator","Center"))

merged2 <- merged2 %>%
  mutate(across(all_of(paste0(value_cols, ".y")),
                ~ ifelse(Indicator == "P85", "ref", .)))

for (col in value_cols) {
  
  value_col <- paste0(col, ".x")
  p_col     <- paste0(col, ".y")
  
  merged2[[col]] <- ifelse(
    is.na(merged2[[p_col]]),
    merged2[[value_col]],
    str_c(merged2[[value_col]], "\n", merged2[[p_col]])
  )
}

merged.final.part2 <- merged2 %>%
  dplyr::select(Indicator,Center, all_of(value_cols))

merged.final <- left_join(merged.final.part1,merged.final.part2) %>% 
  dplyr::select(Indicator,Center,I:II,`Early stage`,III:IVA,`Localregionally advanced stage`,IVB,`Advanced stage`)

write.csv(merged.final, 
          sprintf("%s/%s_Article.csv",save_path,saveDir),
          row.names = F)

