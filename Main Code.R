#### R code for BNLF2b cohort study
rm(list=ls())

#################### Environmental preparation #################################
options("repos" = c(CRAN="https://mirrors.ustc.edu.cn/CRAN/")) 
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")

if (!require("pak", quietly = TRUE)) install.packages("pak")

CRAN_packages_install <- c("remotes","devtools",
                           "readxl","writexl",
                           "tableone","compareGroups", # baseline table
                           "tidyverse","tidyplots","pROC","patchwork",# plot
                           "DTComPair","ocbe-uio/contingencytables","rstatix","statpsych",
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
save.path <- sprintf("%s/results/",current_path)
fig.path <- sprintf("%s/results/Figures",current_path)
table.path <- sprintf("%s/results/Tables",current_path)

ensure_dir(save.path)
ensure_dir(fig.path)
ensure_dir(table.path)

load("./my_color.rdata")

data <- readxl::read_xlsx("./data/原始数据整理/前瞻性数据汇总/combined_data.xlsx") %>%
  dplyr::filter(enroll==1) %>% 
  dplyr::filter(time=="pre-treatment")

data$age <- as.numeric(data$age)
data$gender <- factor(ifelse(data$gender=="男","Male","Female"),levels = c("Male","Female"))
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
writexl::write_xlsx(data,sprintf("%s/combined_data.xlsx",save.path))


##### Figure2 ##################################################################

##### * Fig2A * ######
data_pre <- data  %>% as.data.frame()

##### ROC 
indicators <- c("P85_SCO", "VCA_IgA_SCO","EA_IgA_SCO","NA1_IgA_SCO") 

result_roc_df <- data.frame(
  Indicator = character(0),
  AUC = character(0),
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
  AUC_value_with_3_digits <- sprintf("%0.3f",AUC_value_with_DeLongCI[2])
  AUC_CI_lower <- sprintf("%0.3f",AUC_value_with_DeLongCI[1])
  AUC_CI_upper <- sprintf("%0.3f",AUC_value_with_DeLongCI[3])
  AUC_formatted <- str_c(AUC_value_with_3_digits, 
                         "\n(", AUC_CI_lower, "-", AUC_CI_upper, ")")
  
  if (indicator != "P85_SCO") {
    roc_test_result <- roc.test(roc_p85, roc_object)
    P_value <- ifelse(roc_test_result$p.value < 0.0001, "p < 0.0001", 
                      signif(roc_test_result$p.value,2))
  } else {
    P_value <- "ref"
  }
  
  result_roc_df <- rbind(result_roc_df, c(indicator, AUC_formatted, P_value))
}

colnames(result_roc_df) <- c("Indicator", "AUC", "P_value")
result_roc_df

save_path <- sprintf("%s/Fig2",fig.path)
ensure_dir(save_path)
write.csv(result_roc_df,
          sprintf("%s/Fig2A_supplementary_table_Delong_P_value.csv",save_path),
          row.names = F)


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

pdf(sprintf("%s/Fig2A.pdf",save_path),width = 6,height = 6)
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
            paste0("EBNA1-IgA\n",value_na1) )
  
  legend(0.8, 0.5, 
         bty = "n", 
         legend=text, 
         text.col = c(my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")]),
         col= c(my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")]),
         lwd=3) 
}
dev.off()

#### delong
roc.test(roc_p85,roc_vca)
roc.test(roc_p85,roc_ea)
roc.test(roc_p85,roc_na1)



#####  * Fig2B *######
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

ggVennDiagram(sets,set_color = c("#F94141","#F3B169" ,"#37AB78", "#1685a9"),label="count")+
  scale_fill_gradient(low="grey90",high = "red")

save_path <- sprintf("%s/Fig2",fig.path)
ensure_dir(save_path)
ggsave(sprintf("%s/Fig2B_Veen_NPC.pdf",save_path),width = 5,height = 5)

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

ggVennDiagram(sets,set_color = c("#F94141","#F3B169" ,"#37AB78", "#1685a9"),label="count")+
  scale_fill_gradient(low="grey90",high = "red")

ggsave(sprintf("%s/Fig2B_Veen_nonNPC.pdf",save_path),width = 5,height = 5)

##### * Fig2C-F & FigS3 * ######
data_combine <- data 

# 获取独立的中心列表
centers <- unique(data_combine$center)
centers
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
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare$Center <- center_name
  result_statistical_list[[i]] <- result_compare
  
  #### ROC 
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
    
    save_path <- sprintf("%s/FigS3",fig.path)
    ensure_dir(save_path)
    
    pdf(sprintf("%s/ROC_curve_%s.pdf",save_path,center_name),width = 6,height = 6)
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
              paste0("EBNA1-IgA\n",value_na1)#
      )
      
      legend(0.8, 0.5, 
             bty = "n", 
             legend=text,  
             text.col = my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")],
             col= my_color$single[c("P85-Ab","VCA-IgA","EA-IgA","EBNA1-IgA")], 
             lwd=3) 
    }
    dev.off()
    
    compare1 <- roc.test(roc_p85,roc_vca)
    compare1_pvalue <- ifelse(compare1$p.value < 0.0001, "p < 0.0001", 
                              signif(compare1$p.value,2))
    compare2 <- roc.test(roc_p85,roc_ea)
    compare2_pvalue <- ifelse(compare2$p.value < 0.0001, "p < 0.0001", 
                              signif(compare2$p.value,2))
    
    compare3 <- roc.test(roc_p85,roc_na1)
    compare3_pvalue <- ifelse(compare3$p.value < 0.0001, "p < 0.0001", 
                              signif(compare3$p.value,2))
    
    delong_df <- data.frame(
      center=center_name,
      Indicator = c("P85 vs VCA","P85 vs EA","P85 vs NA1"),
      P_value = c(compare1_pvalue,compare2_pvalue,compare3_pvalue))
    
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
save_path <- sprintf("%s/FigS3",fig.path)
ensure_dir(save_path)
write.csv(AUC.merge,
          sprintf("%s/AUC.merge.csv",save_path),
          row.names = F)

delong.merge <- map_dfr(delong_list,bind_rows)
write.csv(delong.merge,
          sprintf("%s/delong.merge.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  dplyr::select(Indicator,Center,everything())

save_path <- sprintf("%s/Fig2",fig.path)
ensure_dir(save_path)

write.csv(merged.data,
          sprintf("%s/02-Center_allmerge_with_predefined_cut-off_original.csv",save_path),
          row.names = F)

merged.data.statistics <- map_dfr(result_statistical_list,bind_rows) %>% 
  dplyr::select(compare,Center,everything())

write.csv(merged.data.statistics,
          sprintf("%s/02-Center_allmerge_with_predefined_cut-off_original_statistics.csv",
                  save_path),
          row.names = F)

data_forest_center <- read.csv(sprintf("%s/Fig2/02-Center_allmerge_with_predefined_cut-off_original.csv",
                                       fig.path)) %>% 
  as.data.frame()

data_forest <- data_forest_center
data_forest$Indicator %>% unique()
indicators_used <- c("P85","VCA_IgA","EA_IgA","NA1_IgA")

data_forest$Center <- factor(data_forest$Center,
                             levels = c("SYSUCC", "ZSCPH", "WZRCH", "TJH-HUST", "FJCH"))

matrics <- c("Sensitivity","Specificity","PPV","NPV")
save_path <- sprintf("%s/Fig2",fig.path)
ensure_dir(save_path)

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
      strip.background = element_blank(),  # 去掉背景
      strip.placement = "outside",
      panel.spacing.y = unit(0.5, "lines"),
      axis.text = element_text(color="black")
    ) +
    scale_color_manual(values = my_color$single)+
    xlab(str_c(used_matrics," (%)"))+
    ylab("")
  
  plot_list[[i]] <- p
  # ggsave(sprintf("%s/Fig2C-F_%s_across_centers.pdf",save_path,used_matrics),plot = p,
  #        width = 6,height = 8)
}

library(patchwork)

(plot_list[[1]]+plot_list[[2]])/(plot_list[[3]]+plot_list[[4]])
ggsave(sprintf("%s/Fig2C-F_combined_across_centers.pdf",save_path),width =12,height = 10)

##### Figure3 ##################################################################
data_analysis <- data  %>% as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))

##### * Fig3A-D & FigS4 * ######
data_combine <- data_analysis

centers <- unique(data_combine$center)
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA","EBV_DNA0") 

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
                             compare_indicator=c("EBV_DNA0"),
                             outcome="group_ref")
  result_compare$Center <- center_name
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
    
    save_path <- sprintf("%s/FigS4",fig.path)
    ensure_dir(save_path)
    
    pdf(sprintf("%s/Three_markers_ROC_curve_%s.pdf",save_path,center_name),
        width = 6,height = 6)
    {
      plot(roc_p85, 
           col=my_color$single["P85-Ab"],
           legacy.axes=TRUE)
      
      plot.roc(roc_ebv,
               add=TRUE, 
               col=my_color$single["EBV-DNA"], 
               print.auc.x=0.35,
               print.auc.y=0.50)
      
      text<-c(paste0("P85-Ab\n",value_p85),
              paste0("EBV-DNA\n",value_ebv))
      
      legend(0.8, 0.5, 
             bty = "n",  
             legend=text, 
             text.col =my_color$single[c("P85-Ab","EBV-DNA")],
             col=my_color$single[c("P85-Ab","EBV-DNA")],
             lwd=3) 
    }
    dev.off()
    
    compare4 <- roc.test(roc_p85,roc_ebv)
    compare4_pvalue <- ifelse(compare4$p.value < 0.0001, "p < 0.0001", 
                              signif(compare4$p.value,2))
    delong_df <- data.frame(
      center=center_name,
      Indicator = c("P85 vs EBV-DNA"),
      P_value = c(as.character(compare4_pvalue)))
    
    delong_list[[i]] <- delong_df
    
    # write.csv(delong_df,
    #           sprintf("%s/Three_markers_Center_%s_Delong_test.csv",
    #                   save_path,center_name),
    #           row.names = F)
    
    AUC_df <- data.frame(
      center=center_name,
      Indicator = c("P85","EBV-DNA"), # "EBV-DNA"
      AUC = c(value_p85,value_ebv) 
    )
    
    AUC_list[[i]] <- AUC_df 
    # write.csv(AUC_df,
    #           sprintf("%s/Three_markers_Center_%s_AUC.csv",save_path,center_name),row.names = F)
  }
}

AUC.merge <- map_dfr(AUC_list,bind_rows)
save_path <- sprintf("%s/FigS4",fig.path)
ensure_dir(save_path)
write.csv(AUC.merge,
          sprintf("%s/AUC.merge.csv",save_path),
          row.names = F)

delong.merge <- map_dfr(delong_list,bind_rows)
write.csv(delong.merge,
          sprintf("%s/delong.merge.csv",save_path),
          row.names = F)


save_path <- sprintf("%s/Fig3",fig.path)
ensure_dir(save_path)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  dplyr::select(Indicator,Center,everything())

write.csv(merged.data,
          sprintf("%s/02-Center_allmerge_with_predefined_cut-off_original_with_EBVDNA.csv",
                  save_path),
          row.names = F)

merged.data.statistics <- map_dfr(result_statistical_list,bind_rows) %>% 
  dplyr::select(compare,Center,everything())

write.csv(merged.data.statistics,
          sprintf("%s/02-Center_allmerge_with_predefined_cut-off_original_statistics_with_EBVDNA.csv",
                  save_path),
          row.names = F)

#### plot
data_forest <- read.csv(sprintf("%s/Fig3/02-Center_allmerge_with_predefined_cut-off_original_with_EBVDNA.csv",
                                fig.path)) %>% as.data.frame()


data_forest$Indicator %>% unique()
indicators_used <- c("P85","EBV_DNA0")
data_forest$Center <- factor(data_forest$Center,
                             levels = c("SYSUCC", "WZRCH", "FJCH"))

#### 
matrics <- c("Sensitivity","Specificity","PPV","NPV")

save_path <- sprintf("%s/Fig3",fig.path)
ensure_dir(save_path)

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
    dplyr::filter(Indicator%in%c("P85","EBV_DNA0")) %>%  
    mutate(
      est = as.numeric(est), 
      low_ci = as.numeric(low_ci), 
      high_ci = as.numeric(high_ci) 
    ) %>% 
    mutate(Indicator=factor(.$Indicator,levels = rev(c("P85","EBV_DNA0")),
                            labels = rev(c("P85-Ab","EBV-DNA")))) %>% 
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
      strip.background = element_blank(),  # 去掉背景
      strip.placement = "outside",
      panel.spacing.y = unit(0.5, "lines"),
      axis.text = element_text(color="black")
    ) +
    scale_color_manual(values = my_color$single)+
    xlab(str_c(used_matrics," (%)"))+
    ylab("") 
  
  plot_list[[i]] <- p
  # ggsave(sprintf("%s/Fig3B_%s_across_centers.pdf",save_path,used_matrics),plot = p,
  #        width = 6,height = 8)
  
}

library(patchwork)
(plot_list[[1]]+plot_list[[2]])/(plot_list[[3]]+plot_list[[4]])

ggsave(sprintf("%s/Fig3A_D_combined_across_centers.pdf",save_path),width =8,height = 5)

##### * Fig3E * ######
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA","EBV_DNA0") 

data_pre <- data_analysis

save_path <- sprintf("%s/Fig3",fig.path)
ensure_dir(save_path)

##### ROC 
indicators <- c("P85_SCO","log_EBV_DNA") 

result_roc_df <- data.frame(
  Indicator = character(0),
  AUC = character(0),
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
  AUC_value_with_3_digits <- sprintf("%0.5f",AUC_value_with_DeLongCI[2])
  AUC_CI_lower <- sprintf("%0.5f",AUC_value_with_DeLongCI[1])
  AUC_CI_upper <- sprintf("%0.5f",AUC_value_with_DeLongCI[3])
  AUC_formatted <- str_c(AUC_value_with_3_digits, 
                         "\n(", AUC_CI_lower, "-", AUC_CI_upper, ")")
  
  if (indicator != "P85_SCO") {
    roc_test_result <- roc.test(roc_p85, roc_object)
    P_value <- ifelse(roc_test_result$p.value < 0.0001, "p < 0.0001", 
                      signif(roc_test_result$p.value,2))
  } else {
    P_value <- "ref"
  }
  
  result_roc_df <- rbind(result_roc_df, c(indicator, AUC_formatted, P_value))
}

colnames(result_roc_df) <- c("Indicator", "AUC", "P_value")

result_roc_df

write.csv(result_roc_df,sprintf("%s/Supplementary_table_Delong_P_value.csv",save_path),
          row.names = F)

roc_p85 <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$P85_SCO),direction="<",
                     levels = c(0, 1))
value_p85 <- print_auc_and_ci(roc_p85)

roc_EBV <- pROC::roc(data_pre$group_ref, as.numeric(data_pre$log_EBV_DNA),direction="<",
                     levels=c(0, 1))
value_EBV <- print_auc_and_ci(roc_EBV)


pdf(sprintf("%s/Fig3E.pdf",save_path),width = 6,height = 6)
{
  plot(roc_p85, 
       col=my_color$single["P85-Ab"],
       legacy.axes=TRUE)

  plot.roc(roc_EBV,
           add=TRUE, 
           col=my_color$single["EBV-DNA"],
           print.auc.x=0.35,
           print.auc.y=0.50)
  
  text<-c(paste0("P85-Ab\n",value_p85),
          paste0("EBV-DNA\n",value_EBV)
  )
  
  legend(0.8, 0.5,
         bty = "n", 
         legend=text, 
         text.col = my_color$single[c("P85-Ab","EBV-DNA")],
         col=my_color$single[c("P85-Ab","EBV-DNA")],
         lwd=3) 
}
dev.off()

roc.test(roc_p85,roc_EBV)

##### * Fig3F * ######
data_plot <- data %>% as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1)) %>% 
  dplyr::filter(group=="NPC") %>% 
  mutate(`P85-Ab`=as.character(.$P85),
         `EBV-DNA`=as.character(.$EBV_DNA0))

sets <- list(
  `P85-Ab`= which(data_plot$`P85-Ab` == 1),
  `EBV-DNA`= which(data_plot$EBV_DNA0 == 1))

ggVennDiagram(sets,set_color = c("#F94141","#3B54A3"),label="count")+
  scale_fill_gradient(low="grey90",high = "red")

save_path <- sprintf("%s/Fig3",fig.path)
ensure_dir(save_path)
ggsave(sprintf("%s/Fig3F_Veen_NPC.pdf",save_path),width = 5,height = 5)

data_plot <- data %>% as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1)) %>% 
  dplyr::filter(group=="nonNPC") %>% 
  mutate(`P85-Ab`=as.character(.$P85),
         `EBV-DNA`=as.character(.$EBV_DNA0))

sets <- list(
  `P85-Ab`= which(data_plot$`P85-Ab` == 0),
  `EBV-DNA`= which(data_plot$EBV_DNA0 == 0))

ggVennDiagram(sets,set_color = c("#F94141","#3B54A3"),label="count")+
  scale_fill_gradient(low="grey90",high = "red")

ggsave(sprintf("%s/Fig3F_Veen_nonNPC.pdf",save_path),width = 5,height = 5)



##### FigS1 ####################################################################
data_kappa <- data %>% 
  as.data.frame() %>% 
  tidyr::drop_na(P85_plasma)

library(ggpubr)
ggscatter(data_kappa, x = "P85_SCO", 
          y = "P85_plasma_SCO",
          size = 1.5,
          color = "#0072B2",
          add = "reg.line",  # 添加回归线
          add.params = list(color = "black", fill = "grey50", size = 1),  
          conf.int = TRUE) +
  stat_cor(method = "pearson", label.x = 300, label.y = 100, label.sep = "\n", size=3) +
  xlab("Serum P85-Ab COI") +
  ylab("Plasam P85-Ab COI")+
  theme_classic()+theme(aspect.ratio = 1)

save_path <- sprintf("%s/FigS1",fig.path)
ensure_dir(save_path)
ggsave(sprintf("%s/FigS1_cor_plasma_serum.pdf",save_path),width = 6,height = 6)

suppressMessages(library(irr))
suppressMessages(library(statpsych))

data_kappa$P85_plasma <- factor(data_kappa$P85_plasma,levels = c("(+)","(-)"),
                                labels = c(1,0))
#### with CI version 
indicators <- c("P85", "P85_plasma")

kappa_matrix <- data.frame(matrix(ncol = length(indicators), 
                                  nrow = length(indicators)))
colnames(kappa_matrix) <- indicators
rownames(kappa_matrix) <- indicators

rater1 <- data_kappa[[indicators[1]]]
rater2 <- data_kappa[[indicators[2]]]

f00 <- sum(rater1 == 0 & rater2 == 0)
f01 <- sum(rater1 == 0 & rater2 == 1)
f10 <- sum(rater1 == 1 & rater2 == 0)
f11 <- sum(rater1 == 1 & rater2 == 1)

ci_result <- statpsych::ci.kappa(0.05, f00, f01, f10, f11)
ci_result

kappa_with_ci <- sprintf("%0.3f\n(%.3f - %.3f)", 
                         ci_result[2, 1], 
                         ci_result[2, 3], 
                         ci_result[2, 4])
kappa_with_ci
write.csv(as.data.frame(kappa_with_ci),
          file = sprintf("%s/Kappa_value.csv",save_path))

##### FigS6-S8 #################################################################
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 
data_combine <- data

###### ** five centers  #########
data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(NPC_stage %in% c("I","II","III","IVA","IVB"))

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
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$Center <- current_center
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

save_path <- sprintf("%s/Figures/FigS6_S8",save.path)
ensure_dir(save_path)

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/Individual_centers_NPC_5_stages_combined.csv",save_path),row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,I,II,III,IVA,IVB)

write.csv(merged.data,
          sprintf("%s/Individual_centers_NPC_5_stages.csv",save_path),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/Individual_centers_NPC_5_stages_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare,Center,I,II,III,IVA,IVB)

write.csv(merged.statistical.data,
          sprintf("%s/Individual_centers_NPC_5_stages_compare.csv",save_path),
          row.names = F)


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
                             compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare$Center <- current_center
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/Individual_centers_NPC_3_stages_combined.csv",save_path),row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.data,
          sprintf("%s/Individual_centers_NPC_3_stages.csv",save_path),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/Individual_centers_NPC_3_stages_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare ,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.statistical.data,
          sprintf("%s/Individual_centers_NPC_3_stages_compare.csv",save_path),
          row.names = F)

###### ** three centers #######
data_pre <- data %>% 
  as.data.frame() %>%
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  filter(group_ref == 1) %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB")) 

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
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/Individual_centers_Three_centers_NPC_5_stages_combined.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,I,II,III,IVA,IVB)
write.csv(merged.data,
          sprintf("%s/Individual_centers_Three_centers_NPC_5_stages.csv",save_path),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/Individual_centers_Three_centers_NPC_5_stages_compare_combined.csv",
                  save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare,Center,I,II,III,IVA,IVB)

write.csv(merged.statistical.data,
          sprintf("%s/Individual_centers_Three_centers_NPC_5_stages_compare.csv",save_path),
          row.names = F)

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
  result_compare <- result_compare %>% 
    dplyr::select(compare,Center,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows) 
write.csv(merged.data.combined,
          sprintf("%s/Individual_centers_Three_centers_NPC_3_stages_combined.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.data,
          sprintf("%s/Individual_centers_Three_centers_NPC_3_stages.csv",save_path),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
write.csv(merged.statistical.data.combined,
          sprintf("%s/Individual_centers_Three_centers_NPC_3_stages_compare_combined.csv",
                  save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value) %>% 
  dplyr::select(compare ,Center,`Early stage`,`Localregionally advanced stage`,`Advanced stage`)

write.csv(merged.statistical.data,
          sprintf("%s/Individual_centers_Three_centers_NPC_3_stages_compare.csv",save_path),
          row.names = F)

###### ** Plot #######
data_forest_original <- read.csv(sprintf("%s/FigS6_S8/Individual_centers_NPC_3_stages.csv",
                                         fig.path)) %>% 
  as.data.frame()

data_forest_stat <- read.csv(sprintf("%s/FigS6_S8/Individual_centers_NPC_3_stages_compare.csv",
                                     fig.path)) %>% 
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
    # add_data_labels(label = est) %>% 
    adjust_y_axis_title("Sensitivity (%)") %>% 
    remove_x_axis_labels() %>% 
    split_plot(by = Center, ncol = 2, nrow = 3)
  
  ggsave(sprintf("%s/Sensitivity_%s_across_centers_between_among_antibody.pdf",
                 save_path,three_stages[i]),
         width = 4,height = 12)
  
}

data_forest_original <- read.csv(sprintf("%s/FigS6_S8/Individual_centers_Three_centers_NPC_3_stages.csv",
                                         fig.path)) %>% 
  as.data.frame()

data_forest_stat <- read.csv(sprintf("%s/FigS6_S8/Individual_centers_Three_centers_NPC_3_stages_compare.csv",
                                     fig.path)) %>% 
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
                             levels = c("SYSUCC", "WZRCH", "FJCH"))

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
    dplyr::filter(Indicator%in%c("P85","EBV_DNA0")) %>%  
    mutate(
      est = as.numeric(est), 
      low_ci = as.numeric(low_ci), 
      high_ci = as.numeric(high_ci) 
    )
  
  data_sen %>% 
    mutate(Indicator=factor(.$Indicator,levels = c("P85","EBV_DNA0"),
                            labels = c("P85-Ab","EBV-DNA"))) %>% 
    mutate(label=str_c(ifelse(est==100,"100",sprintf("%0.1f",est)),
                       "\n",pvalue)) %>% 
    tidyplot(x = Center, y=est,color = Indicator,dodge_width = 0.6) %>%
    add_mean_dot(size = 1) %>% 
    add(ggplot2::geom_errorbar(aes(ymin =low_ci, ymax = high_ci,
                                   group=Indicator),width=0.5,
                               position = position_dodge(0.6))) %>% 
    remove_x_axis_title() %>% 
    add(ggplot2::theme(axis.text.x = element_text(angle=60,hjust = 1,vjust = 1))) %>% 
    adjust_colors(new_colors = my_color$single[c("P85-Ab","EBV-DNA")]) %>%
    adjust_y_axis(limits=c(0,100)) %>% 
    add(ggplot2::geom_label(aes(y=low_ci-15,label = label),fill="transparent",
                            label.size=NA,size=3,
                            position = position_dodge(0.6))) %>%
    add(ggplot2::theme(aspect.ratio = 0.8)) %>%
    adjust_y_axis_title("Sensitivity (%)") 
  
  ggsave(sprintf("%s/Sensitivity_%s_across_centers_between_P85_EBVDNA.pdf",
                 save_path,three_stages[i]),
         width = 5,height = 5)
  
}

##### FigS9-S10 ################################################################
###### * Sex ###########
data_combine <- data
sex_class <- unique(data_combine$gender)
sex_class
table(data_combine$gender)
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
  
  # write.csv(result_df,
  #           sprintf("./results/10-subgroup_%s.csv",sex_name),
  #           row.names = F)
  
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
    AUC = c(sprintf("%0.3f",pROC::ci.auc(roc_p85)[2]),
            sprintf("%0.3f",pROC::ci.auc(roc_vca)[2]),
            sprintf("%0.3f",pROC::ci.auc(roc_ea)[2]),
            sprintf("%0.3f",pROC::ci.auc(roc_na1)[2])),
    low = c(sprintf("%0.3f",pROC::ci.auc(roc_p85)[1]),
            sprintf("%0.3f",pROC::ci.auc(roc_vca)[1]),
            sprintf("%0.3f",pROC::ci.auc(roc_ea)[1]),
            sprintf("%0.3f",pROC::ci.auc(roc_na1)[1])),
    hi = c(sprintf("%0.3f",pROC::ci.auc(roc_p85)[3]),
           sprintf("%0.3f",pROC::ci.auc(roc_vca)[3]),
           sprintf("%0.3f",pROC::ci.auc(roc_ea)[3]),
           sprintf("%0.3f",pROC::ci.auc(roc_na1)[3])),
    format= c(value_p85,value_vca,value_ea,value_na1),
    group=sex_name)
  
  # write.csv(auc_df,
  #           sprintf("./results/10-subgroup_auc_%s.csv",sex_name),
  #           row.names = F)
  auc_list[[i]] <- auc_df
}

save_path <- sprintf("%s/FigS9_S10",fig.path)
ensure_dir(save_path)

merged.data <- map_dfr(result_list,bind_rows) 
write.csv(merged.data,
          sprintf("%s/subgroup_sex_diagnostic_performance.csv",save_path),row.names = F)

merged.auc_list <- map_dfr(auc_list,bind_rows) 
write.csv(merged.auc_list,
          sprintf("%s/subgroup_sex_AUC.csv",save_path),row.names = F)


data_analysis <- data %>% 
  dplyr::filter(center %in% c("SYSUCC","WZRCH","FJCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))

data_combine <- data_analysis

sex_class <- unique(data_combine$gender)
sex_class
table(data_combine$center)
indicators <- c("EBV_DNA0")

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
  
  # write.csv(result_df,
  #           sprintf("./results/10-subgroup_%s.csv",sex_name),
  #           row.names = F)
  
  result_list[[i]] <- result_df
  
  roc_ebv <- pROC::roc(data_used$group_ref, as.numeric(data_used$log_EBV_DNA),direction="<",
                       levels = c(0, 1))
  value_ebv <- print_auc_and_ci(roc_ebv)
  
  auc_df <- data.frame(
    Indicator = c("EBV_DNA"),
    AUC = c(sprintf("%0.3f",pROC::ci.auc(roc_ebv)[2])),
    low = c(sprintf("%0.3f",pROC::ci.auc(roc_ebv)[1])),
    hi = c(sprintf("%0.3f",pROC::ci.auc(roc_ebv)[3])),
    format= c(value_ebv),
    group=sex_name)
  
  # write.csv(auc_df,
  #           sprintf("./results/10-subgroup_auc_%s.csv",sex_name),
  #           row.names = F)
  auc_list[[i]] <- auc_df
}

merged.data <- map_dfr(result_list,bind_rows) 

write.csv(merged.data,
          sprintf("%s/subgroup_sex_diagnostic_performance_EBVDNA.csv",save_path),row.names = F)

merged.auc_list <- map_dfr(auc_list,bind_rows) 
write.csv(merged.auc_list,
          sprintf("%s/subgroup_sex_AUC_EBVDNA.csv",save_path),row.names = F)


###### * Age ###########
data_combine <- data
data_age_calculate <- data_combine %>% 
  mutate(age_cut = case_when(
    age < 45 ~ "less than 45",
    age >= 45 ~ "more tham 45"
  ))

age_class <- unique(data_age_calculate$age_cut)
table(data_age_calculate$age_cut)

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
  
  # write.csv(result_df,
  #           sprintf("./results/10-subgroup_%s_45cut.csv",age_name),
  #           row.names = F)
  
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
    AUC = c(sprintf("%0.3f",pROC::ci.auc(roc_p85)[2]),
            sprintf("%0.3f",pROC::ci.auc(roc_vca)[2]),
            sprintf("%0.3f",pROC::ci.auc(roc_ea)[2]),
            sprintf("%0.3f",pROC::ci.auc(roc_na1)[2])),
    low = c(sprintf("%0.3f",pROC::ci.auc(roc_p85)[1]),
            sprintf("%0.3f",pROC::ci.auc(roc_vca)[1]),
            sprintf("%0.3f",pROC::ci.auc(roc_ea)[1]),
            sprintf("%0.3f",pROC::ci.auc(roc_na1)[1])),
    hi = c(sprintf("%0.3f",pROC::ci.auc(roc_p85)[3]),
           sprintf("%0.3f",pROC::ci.auc(roc_vca)[3]),
           sprintf("%0.3f",pROC::ci.auc(roc_ea)[3]),
           sprintf("%0.3f",pROC::ci.auc(roc_na1)[3])),
    format= c(value_p85,value_vca,value_ea,value_na1),
    group=age_name)
  
  # write.csv(auc_df,
  #           sprintf("./results/10-subgroup_auc_%s_45cut.csv",age_name),
  #           row.names = F)
  
  auc_list[[i]] <- auc_df
}

merged.data <- map_dfr(result_list,bind_rows) 
write.csv(merged.data,
          sprintf("%s/subgroup_age_diagnostic_performance.csv",save_path),row.names = F)

merged.auc_list <- map_dfr(auc_list,bind_rows) 
write.csv(merged.auc_list,
          sprintf("%s/subgroup_age_AUC.csv",save_path),row.names = F)

data_analysis <- data %>% 
  dplyr::filter(center %in% c("SYSUCC","WZRCH","FJCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))

data_combine <- data_analysis

data_age_calculate <- data_combine %>% 
  mutate(age_cut = case_when(
    age < 45 ~ "less than 45",
    age >= 45 ~ "more tham 45"
  ))

age_class <- unique(data_age_calculate$age_cut)
table(data_age_calculate$age_cut)

indicators <- c("EBV_DNA0") 

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
  
  # write.csv(result_df,
  #           sprintf("./results/10-subgroup_%s_45cut.csv",age_name),
  #           row.names = F)
  
  result_list[[i]] <- result_df
  
  roc_ebv <- pROC::roc(data_used$group_ref, as.numeric(data_used$log_EBV_DNA),direction="<",
                       levels = c(0, 1))
  value_ebv <- print_auc_and_ci(roc_ebv)
  
  auc_df <- data.frame(
    Indicator = c("EBV_DNA"),
    AUC = c(sprintf("%0.3f",pROC::ci.auc(roc_ebv)[2])),
    low = c(sprintf("%0.3f",pROC::ci.auc(roc_ebv)[1])),
    hi = c(sprintf("%0.3f",pROC::ci.auc(roc_ebv)[3])),
    format= c(value_ebv),
    group=age_name)
  
  # write.csv(auc_df,
  #           sprintf("./results/10-subgroup_auc_%s_45cut.csv",age_name),
  #           row.names = F)
  
  auc_list[[i]] <- auc_df
}

merged.data <- map_dfr(result_list,bind_rows) 
write.csv(merged.data,
          sprintf("%s/subgroup_age_diagnostic_performance_EBVDNA.csv",save_path),row.names = F)

merged.auc_list <- map_dfr(auc_list,bind_rows) 
write.csv(merged.auc_list,
          sprintf("%s/subgroup_age_AUC_EBVDNA.csv",save_path),row.names = F)

###### * Plot ###########
data_age <- read.csv(sprintf("%s/FigS9_S10/subgroup_age_diagnostic_performance.csv",fig.path))
data_sex <- read.csv(sprintf("%s/FigS9_S10/subgroup_sex_diagnostic_performance.csv",fig.path))
data_age_EBVDNA <- read.csv(sprintf("%s/FigS9_S10/subgroup_age_diagnostic_performance_EBVDNA.csv",fig.path))
data_sex_EBVDNA <- read.csv(sprintf("%s/FigS9_S10/subgroup_sex_diagnostic_performance_EBVDNA.csv",fig.path))
data_subgroup <- rbind(data_age,data_age_EBVDNA,data_sex,data_sex_EBVDNA) 

suppressMessages(library(forestploter))

data_forest <- data_subgroup %>% 
  dplyr::select(Indicator,group,Total,TP:TN,sen_est1:spe_hi2)

custom_order <- c("Male", "Female", "less than 45", "more tham 45")
ordered_indices <- order(match(data_forest$group, custom_order))
sorted_data <- data_forest[ordered_indices, ]
data_forest <- sorted_data %>% as_tibble() %>% 
  add_row(group="Sex",.before = 1) %>% 
  add_row(group="Age,yr",.before = 12) 

writexl::write_xlsx(data_forest,
                    sprintf("%s/subgroup_combine_for_forest.xlsx",save_path))

indictors_all <- c("P85","VCA_IgA","EA_IgA","NA1_IgA","EBV_DNA0")
for (i in 1:length(indictors_all)){
  indicators_used <- indictors_all[i]
  print(indictors_all[i])
  data_forest_plot <- data_forest %>% 
    dplyr::filter(Indicator %in% indictors_all[i] | is.na(Indicator)) 
  data_forest_plot$group <- ifelse(is.na(data_forest_plot$Total), 
                                   data_forest_plot$group,
                                   paste0("   ", data_forest_plot$group))
  data_forest_plot <- data_forest_plot %>%
    mutate(across(c(Total, TP, FP, FN, TN), ~ifelse(is.na(.), "", .)))
  
  data_forest_plot$Sensitivity <- paste(rep(" ", 10), collapse = " ")
  data_forest_plot$Specificity <- paste(rep(" ", 10), collapse = " ")
  
  data_forest_plot$`Sensitivity (95%CI)` <- ifelse(is.na(data_forest_plot$sen_est1),"",
                                                   sprintf("%.3f (%.3f to %.3f)",
                                                           data_forest_plot$sen_est1,
                                                           data_forest_plot$sen_low1,
                                                           data_forest_plot$sen_hi1))
  
  data_forest_plot$`Specificity (95%CI)` <- ifelse(is.na(data_forest_plot$spe_est2),"",
                                                   sprintf("%.3f (%.3f to %.3f)",
                                                           data_forest_plot$spe_est2,
                                                           data_forest_plot$spe_low2,
                                                           data_forest_plot$spe_hi2))
  data_forest_plot <- data_forest_plot %>% 
    dplyr::select(group:TN,Sensitivity,`Sensitivity (95%CI)`,Specificity,`Specificity (95%CI)`,sen_est1:spe_hi2)
  
  tm <- forestploter::forest_theme(base_size = 10,  # 设置基本字体大小
                                   refline_gp = grid::gpar(lwd = 1, lty = "dashed", col = "red"),
                                   arrow_type = "closed",  # 设置箭头类型为闭合箭头
                                   footnote_gp = grid::gpar(col = "skyblue"))
  
  p <- forestploter::forest(data_forest_plot[,c(1:6,7:10)],    #  选择需要绘制的列
                            est=list(data_forest_plot$sen_est1, 
                                     data_forest_plot$spe_est2),
                            lower=list(data_forest_plot$sen_low1,
                                       data_forest_plot$spe_low2), 
                            upper=list(data_forest_plot$sen_hi1,
                                       data_forest_plot$spe_hi2),
                            sizes= list(as.numeric(data_forest_plot$Total)/3777,
                                        as.numeric(data_forest_plot$Total)/3777),
                            ci_column  = c(7,9),    #  置信区间的列索引
                            ref_line  = c(0.9,0.95),#  参考线位置
                            xlim = list(c(0.5, 1.0),
                                        c(0.6, 1.0)),
                            ticks_at=list(c(0.5, 0.8,1.0),
                                          c(0.6, 0.8,0.9,1.0)),
                            theme  =  tm) ;p
  
  pdf(sprintf("%s/forest_plot_subgroup_%s.pdf",save_path,indicators_used))
  print(p)
  dev.off()
}

### forest plot for AUC
data_age <- read.csv(sprintf("%s/FigS9_S10/subgroup_age_AUC.csv",fig.path))
data_sex <- read.csv(sprintf("%s/FigS9_S10/subgroup_sex_AUC.csv",fig.path))
data_age_EBVDNA <- read.csv(sprintf("%s/FigS9_S10/subgroup_age_AUC_EBVDNA.csv",fig.path))
data_sex_EBVDNA <- read.csv(sprintf("%s/FigS9_S10/subgroup_sex_AUC_EBVDNA.csv",fig.path))
data_subgroup <- rbind(data_age,data_age_EBVDNA,data_sex,data_sex_EBVDNA)

suppressMessages(library(forestploter))

data_forest <- data_subgroup %>% 
  select(Indicator,group,AUC:format)

custom_order <- c("Male", "Female", "less than 45", "more tham 45")
ordered_indices <- order(match(data_forest$group, custom_order))
sorted_data <- data_forest[ordered_indices, ]
data_forest <- sorted_data %>% as_tibble() %>% 
  add_row(group="Sex",.before = 1) %>% 
  add_row(group="Age,yr",.before = 12) %>% 
  dplyr::rename("AUC (95%CI)"="format")

writexl::write_xlsx(data_forest,
                    sprintf("%s/subgroup_combine_for_forest_AUC.xlsx",save_path))

totol_num <- readxl::read_xlsx(sprintf("%s/subgroup_combine_for_forest.xlsx",save_path)) %>% 
  tidyr::drop_na() %>% 
  distinct() %>% 
  mutate(Indicator=case_when(
    Indicator=="P85"~"P85-Ab",
    Indicator=="VCA_IgA"~"VCA-IgA",
    Indicator=="EA_IgA"~"EA-IgA",
    Indicator=="NA1_IgA"~"EBNA1-IgA",
    Indicator=="db_method2"~"db_method2",
    Indicator=="EBV_DNA0"~"EBV_DNA"
  )) %>% 
  mutate(connect_key=str_c(Indicator,group,sep = "_")) %>% 
  dplyr::select(connect_key,Total)

data_forest <- data_forest %>% 
  mutate(connect_key=str_c(Indicator,group,sep = "_")) %>% 
  left_join(totol_num,by="connect_key") %>% 
  dplyr::select(-connect_key)

data_forest$group <- ifelse(is.na(data_forest$Total), 
                            data_forest$group,
                            paste0("   ", data_forest$group))

data_forest <- data_forest %>%
  mutate(across(c(Total,`AUC (95%CI)`), ~ifelse(is.na(.), "", .)))

data_forest$AUC_ci  <- paste(rep(" ", 10), collapse = " ")

tm <- forestploter::forest_theme(base_size = 10,  # 设置基本字体大小
                                 refline_gp = grid::gpar(lwd = 1, lty = "dashed", col = "red"),
                                 arrow_type = "closed",  # 设置箭头类型为闭合箭头
                                 footnote_gp = grid::gpar(col = "skyblue"))

p <- forest(data_forest[,c(2,7,8,6)],    #  选择需要绘制的列
            est=data_forest$AUC,
            lower=data_forest$low,  #  下限
            upper=data_forest$hi,    #  上限
            ci_column  = 3,    #  置信区间的列索引
            ref_line  = 0.95,#  参考线位置
            xlim = c(0.70, 1.0),
            ticks_at=c(0.70,0.80,0.9,1.0),
            theme  =  tm) ;p
pdf(sprintf("%s/forest_plot_subgroup_AUC.pdf",save_path))
print(p)
dev.off()

##### Table2 ###################################################################
###### ** five centers ########
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA") 
data_pre <- data %>% as.data.frame()

result_df <- calcuate_diagnosis_index2(data=data_pre,
                                       indicators = indicators,
                                       binomMethod_use = "wilson",
                                       outcome = "group_ref")

result_df
save_path <- sprintf("%s/Table2",table.path)
ensure_dir(save_path)
write.csv(result_df,sprintf("%s/Table2_Pooled_analysis.csv",save_path),row.names = F)

result_compare <- DTC_pair(data=data_pre,
                           indicator1="P85",
                           compare_indicator=c("VCA_IgA","EA_IgA", "NA1_IgA"),
                           outcome="group_ref")
result_compare
write.csv(result_compare,
          sprintf("%s/Table2_Pooled_analysis_compare.csv",save_path),row.names = F)

###### ** three centers ########
indicators <- c("P85", "EBV_DNA0") 
data_analysis <- data  %>% as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))
table(data_analysis$center)

data_pre <- data_analysis

result_df <- calcuate_diagnosis_index2(data=data_pre,
                                       indicators = indicators,
                                       binomMethod_use = "wilson",
                                       outcome = "group_ref")
save_path <- sprintf("%s/Table2",table.path)
ensure_dir(save_path)

write.csv(result_df,
          sprintf("%s/Table2_Pooled_analysis_Three_centers.csv",save_path),row.names = F)

result_compare <- DTC_pair(data=data_pre,
                           indicator1="P85",
                           compare_indicator=c("EBV_DNA0"),
                           outcome="group_ref")
result_compare

write.csv(result_compare,
          sprintf("%s/Table2_Pooled_analysis_compare_three_centers.csv",save_path),row.names = F)


#### Table3 ####################################################################
###### ** five centers ######
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA")

data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(NPC_stage %in% c("I","II","III","IVA","IVB","/"))

###### 5 stages 
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
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=TP+FN) %>% 
  pivot_wider(id_cols =-c(sum,TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>%
  dplyr::select(Indicator,I,II,III,IVA,IVB)

save_path <- sprintf("%s/Table3",table.path)
ensure_dir(save_path)

write.csv(merged.data,
          sprintf("%s/pooled_analysis_NPC_5_stages.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/pooled_analysis_NPC_5_stages_compare.csv",save_path),
          row.names = F)


###### 3 stages 
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
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity)

write.csv(merged.data,
          sprintf("%s/pooled_analysis_NPC_3_stages.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/pooled_analysis_NPC_3_stages_compare.csv",save_path),
          row.names = F)

###### ** three centers #####
indicators <- c("P85", "EBV_DNA0")

data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  filter(NPC_stage %in% c("I","II","III","IVA","IVB","/"))

##### 5 stages 
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
                             compare_indicator=c("EBV_DNA0"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=TP+FN) %>% 
  pivot_wider(id_cols =-c(sum,TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,I,II,III,IVA,IVB)

write.csv(merged.data,
          sprintf("%s/Three_centers_NPC_5_stages.csv",save_path),row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/Three_centers_NPC_5_stages_compare.csv",save_path),
          row.names = F)


##### 3 stages 
data_pre <- data %>% 
  as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
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
                             compare_indicator=c("EBV_DNA0"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity)
write.csv(merged.data,sprintf("%s/Three_centers_NPC_3_stages.csv",save_path),row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/Three_centers_NPC_3_stages_compare.csv",save_path),row.names = F)



#### Table4 ####################################################################
###### ** five centers ########
indicators <- c("P85", "VCA_IgA", "EA_IgA", "NA1_IgA")

data_pre <- data %>% 
  as.data.frame() %>% 
  filter(group_ref == 0)

stages <- unique(data_pre$category_green)
# [1] "Lymphoma"                     "Other Benign Disease"         "EBV-related Benign Disease"   "Non-NPC Head and Neck Cancer"
# [5] "Other Malignancy" 

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
  
  result_statistical_list[[i]] <- result_compare
}


save_path <- sprintf("%s/Table4",table.path)
ensure_dir(save_path)

merged.data.combined <- map_dfr(result_list,bind_rows) 
write.csv(merged.data,
          sprintf("%s/pooled_analysis_Control_combined.csv",save_path),row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=FP+TN) %>% 
  pivot_wider(id_cols =-c(sum,FP:TN),
              names_from = Class,
              values_from = Specificity) %>% 
  dplyr::select(Indicator,Lymphoma,`EBV-related Benign Disease`,
                `Non-NPC Head and Neck Cancer`,`Other Malignancy`,`Other Benign Disease`)
write.csv(merged.data,
          sprintf("%s/pooled_analysis_Control.csv",save_path),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/pooled_analysis_Control_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    spe_method=="McNemar"~spe.p.value,
    spe_method=="McNemar mid-P"~str_c(spe.p.value,"a"),
  )) %>%
  pivot_wider(id_cols = -c(spe1:spe.p.value),
              names_from = Class,
              values_from = p_value) %>% 
  dplyr::select(compare,Lymphoma,`EBV-related Benign Disease`,
                `Non-NPC Head and Neck Cancer`,`Other Malignancy`,`Other Benign Disease`)

write.csv(merged.statistical.data,
          sprintf("%s/pooled_analysis_Control_compare.csv",save_path),
          row.names = F)

###### ** three centers #####
data_analysis <- data %>% as.data.frame() %>% 
  dplyr::filter(center %in% c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))

##### Diagnosis Performance 
indicators <- c("P85","EBV_DNA0")

data_pre <- data_analysis %>% 
  filter(group_ref == 0)

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
                             compare_indicator=c("EBV_DNA0"),
                             outcome="group_ref")
  result_compare <- result_compare %>% dplyr::select(compare,spe1:spe.p.value)
  
  result_compare$Class <- current_NPC_stage
  
  result_statistical_list[[i]] <- result_compare
}


merged.data.combined <- map_dfr(result_list,bind_rows) 
write.csv(merged.data.combined,
          sprintf("%s/Three_centers_Control_combined.csv",save_path),row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  mutate(sum=FP+TN) %>% 
  pivot_wider(id_cols =-c(sum,FP:TN),
              names_from = Class,
              values_from = Specificity) %>% 
  dplyr::select(Indicator,Lymphoma,`EBV-related Benign Disease`,
                `Non-NPC Head and Neck Cancer`,`Other Malignancy`,`Other Benign Disease`)
write.csv(merged.data,
          sprintf("%s/Three_centers_Control.csv",save_path),row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/Three_centers_Control_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    spe_method=="McNemar"~spe.p.value,
    spe_method=="McNemar mid-P"~str_c(spe.p.value,"a"),
  )) %>%
  pivot_wider(id_cols = -c(spe1:spe.p.value),
              names_from = Class,
              values_from = p_value) %>% 
  dplyr::select(compare,Lymphoma,`EBV-related Benign Disease`,
                `Non-NPC Head and Neck Cancer`,`Other Malignancy`,`Other Benign Disease`)
write.csv(merged.statistical.data,
          sprintf("%s/Three_centers_Control_compare.csv",save_path),
          row.names = F)

##### Model development and validation #########################################
suppressWarnings(library(mgcv))
set.seed(100)

data_logistic <- data
data_logistic$group_ref <- ifelse(data_logistic$group=="NPC",1,0)

###### * development ############
###### ** five centers #####
###### *** VCA_IgA+EBNA1-IgA ------
if(T){
  fit <- glm(group_ref~VCA_IgA_SCO+NA1_IgA_SCO,family=binomial(link = "logit"),
             data = data_logistic)
  summary(fit)
  
  data_logistic$db_method_new_probability <- stats::predict(newdata=data_logistic,fit,"response")
  roc_db <- pROC::roc(data_logistic$group_ref, data_logistic$db_method_new_probability)
  plot(roc_db,
       legacy.axes = TRUE,
       thresholds="best",
       print.thres="best") 
  roc_result <- coords(roc_db, "best");roc_result
  thresthold_db <- roc_result$threshold;thresthold_db
  data_logistic$db_method3 <- ifelse(data_logistic$db_method_new_probability>=thresthold_db,1,0)
}

###### *** P85-Ab + VCA-IgA ----------
if(T){
  fit.gam_P85_VCA <- mgcv::gam(group_ref==1~s(P85_SCO)+s(VCA_IgA_SCO), 
                               data = data_logistic, 
                               family=binomial)
  
  predicted.gam_P85_VCA <- mgcv::predict.gam(object = fit.gam_P85_VCA,newdata = data_logistic)
  apparent.auc.gam_P85_VCA <- pROC::auc(data_logistic$group_ref~predicted.gam_P85_VCA)
  apparent.auc.gam_P85_VCA
  data_logistic$gam_P85_VCA <- predicted.gam_P85_VCA
  
  roc_obj <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$gam_P85_VCA));roc_obj
  plot(roc_obj,
       legacy.axes = TRUE,
       thresholds="best",
       print.thres="best") 
  roc_result <- coords(roc_obj, "best");roc_result
  thresthold_gam_P85_VCA <- roc_result$threshold;thresthold_gam_P85_VCA
  data_logistic$P85_VCA_cut <- ifelse(data_logistic$gam_P85_VCA>=thresthold_gam_P85_VCA,1,0)
}


###### *** P85-Ab + EBNA1-IgA -------------
if(T){
  fit.gam_P85_EBNA1 <- mgcv::gam(group_ref==1~s(P85_SCO)+s(NA1_IgA_SCO), 
                                 data = data_logistic, 
                                 family=binomial)
  
  predicted.gam_P85_EBNA1 <- mgcv::predict.gam(object = fit.gam_P85_EBNA1,newdata = data_logistic)
  apparent.auc.gam_P85_EBNA1 <- pROC::auc(data_logistic$group_ref~predicted.gam_P85_EBNA1)
  apparent.auc.gam_P85_EBNA1
  data_logistic$gam_P85_EBNA1 <- predicted.gam_P85_EBNA1
  
  roc_obj <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$gam_P85_EBNA1),
                       direction="<",
                       levels = c(0, 1));roc_obj
  plot(roc_obj,
       legacy.axes = TRUE,
       thresholds="best",
       print.thres="best") 
  roc_result <- coords(roc_obj, "best");roc_result
  thresthold_gam_P85_EBNA1 <- roc_result$threshold;thresthold_gam_P85_EBNA1
  data_logistic$P85_EBNA1_cut <- ifelse(data_logistic$gam_P85_EBNA1>=thresthold_gam_P85_EBNA1,1,0)
}


###### *** P85-Ab + VCA-IgA + EBNA1-IgA --------------
if(T){
  fit.gam_P85_dual <- mgcv::gam(group_ref~s(P85_SCO)+s(NA1_IgA_SCO)+s(VCA_IgA_SCO),
                                data = data_logistic, 
                                family=binomial)
  
  predicted.gam_P85_dual <- mgcv::predict.gam(object = fit.gam_P85_dual,newdata = data_logistic)
  apparent.auc.gam_P85_dual <- pROC::auc(data_logistic$group_ref~predicted.gam_P85_dual)
  apparent.auc.gam_P85_dual
  data_logistic$gam_P85_dual <- predicted.gam_P85_dual
  
  roc_obj <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$gam_P85_dual),
                       direction="<",
                       levels = c(0, 1));roc_obj
  plot(roc_obj,
       legacy.axes = TRUE,
       thresholds="best", 
       print.thres="best") 
  roc_result <- coords(roc_obj, "best");roc_result
  thresthold_gam_P85_dual <- roc_result$threshold;thresthold_gam_P85_dual
  data_logistic$P85_dual_cut <- ifelse(data_logistic$gam_P85_dual>=thresthold_gam_P85_dual,1,0)
}

###### ** three centers ########
data_analysis <- data %>% as.data.frame() %>% 
  dplyr::filter(center%in%c("SYSUCC","FJCH","WZRCH")) %>% 
  dplyr::filter(!is.na(EBV_DNA_absolute)) %>% 
  dplyr::mutate(log_EBV_DNA = log(EBV_DNA_absolute + 1))
data_analysis$group_ref <- ifelse(data_analysis$group=="NPC",1,0)

####### *** P85-Ab + log EBV-DNA -------------
if(T){
  fit.gam_P85_EBVDNA <- mgcv::gam(group_ref==1~s(log_EBV_DNA)+s(P85_SCO),
                                  data = data_analysis,
                                  family=binomial)
  
  predicted.gam_P85_EBVDNA <- mgcv::predict.gam(object = fit.gam_P85_EBVDNA,newdata = data_analysis)
  apparent.auc.gam_P85_EBVDNA <- pROC::auc(data_analysis$group_ref~predicted.gam_P85_EBVDNA)
  apparent.auc.gam_P85_EBVDNA
  data_analysis$gam_P85_EBVDNA <- predicted.gam_P85_EBVDNA
  
  roc_obj <- pROC::roc(data_analysis$group_ref, as.numeric(data_analysis$gam_P85_EBVDNA));roc_obj
  plot(roc_obj,
       legacy.axes = TRUE,
       thresholds="best",
       print.thres="best") 
  roc_result <- coords(roc_obj, "best");roc_result
  thresthold_gam_P85_EBVDNA <- roc_result$threshold;thresthold_gam_P85_EBVDNA
  data_analysis$P85_EBVDNA_cut <- ifelse(data_analysis$gam_P85_EBVDNA>=thresthold_gam_P85_EBVDNA,1,0)
}


####### *** P85-Ab+dual+EBVDNA ---------------------------------
if(T){
  fit.gam_P85_dual_EBVDNA <- mgcv::gam(group_ref==1~s(log_EBV_DNA)+s(P85_SCO)+s(NA1_IgA_SCO)+s(VCA_IgA_SCO),
                                       data = data_analysis,family=binomial)
  
  predicted.gam_P85_dual_EBVDNA <- mgcv::predict.gam(object = fit.gam_P85_dual_EBVDNA,newdata = data_analysis)
  apparent.auc.gam_P85_dual_EBVDNA <- pROC::auc(data_analysis$group_ref~predicted.gam_P85_dual_EBVDNA)
  apparent.auc.gam_P85_dual_EBVDNA
  data_analysis$gam_P85_dual_EBVDNA <- predicted.gam_P85_dual_EBVDNA
  
  roc_obj <- pROC::roc(data_analysis$group_ref, as.numeric(data_analysis$gam_P85_dual_EBVDNA));roc_obj
  plot(roc_obj,
       legacy.axes = TRUE,
       thresholds="best", 
       print.thres="best") 
  roc_result <- coords(roc_obj, "best");roc_result
  thresthold_gam_P85_dual_EBVDNA <- roc_result$threshold;thresthold_gam_P85_dual_EBVDNA
  data_analysis$P85_dual_EBVDNA_cut <- ifelse(data_analysis$gam_P85_dual_EBVDNA>=thresthold_gam_P85_dual_EBVDNA,
                                              1,0)
}


####### *** p85+VCA+EBNA1 -----------------------------------
if(T){
  P85_dual_df <- data_logistic %>%
    dplyr::select(center,patientName,adminsionID,gam_P85_dual,P85_dual_cut)
  data_analysis <- data_analysis %>% left_join(P85_dual_df)
}


data_analysis$P85_EBVDNA_cut <- factor(data_analysis$P85_EBVDNA_cut,levels = c("1","0"))
data_analysis$P85_dual_EBVDNA_cut <- factor(data_analysis$P85_dual_EBVDNA_cut,levels = c("1","0"))
data_analysis$P85_dual_cut <- factor(data_analysis$P85_dual_cut,levels = c("1","0"))

#### IECV & TableS3 ###################################################################
###### ** five centers #########
data_IECV <- data_logistic
clusters <- unique(data_IECV$center)
N.clust <- length(clusters) # # 5 clusters
data.in <- data.leftout <- list()
results <- list()

# create the datasets
for(i in 1:N.clust){
  data.in[[i]] <- data_IECV[data_IECV$center!=clusters[i],]
  data.leftout[[i]] <- data_IECV[data_IECV$center==clusters[i],]
}

####### *** dual-antibody strategy -----------------------
model <- "dual-antibody strategy"
model.fit <- fit
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
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.logistic[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
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
                                        leftout.prediction.logistic.agg)
round(IECV.logistic,3)

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

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.3f (%0.3f to %0.3f)",auc,LCI,UCI)) %>% 
  left_join(leftout.performance.logistic.agg %>% 
              mutate(calibration.intercept=sprintf("%0.3f (%0.3f to %0.3f)",
                                                   calibration.intercept,
                                                   calibration.intercept_LCI,
                                                   calibration.intercept_UCI),
                     calibration.slope=sprintf("%0.3f (%0.3f to %0.3f)",
                                               calibration.slope,
                                               calibration.slope_LCI,
                                               calibration.slope_UCI)) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% 
  dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

results[[model]] <- data.save

####### *** P85-Ab+VCA-IgA -------------------------------
model <- "P85-Ab+VCA-IgA"
model.fit <- fit.gam_P85_VCA
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,
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
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
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
round(IECV.gam,3)

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

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.3f (%0.3f to %0.3f)",auc,LCI,UCI)) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.3f (%0.3f to %0.3f)",
                                                   calibration.intercept,
                                                   calibration.intercept_LCI,
                                                   calibration.intercept_UCI),
                     calibration.slope=sprintf("%0.3f (%0.3f to %0.3f)",
                                               calibration.slope,
                                               calibration.slope_LCI,
                                               calibration.slope_UCI)) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% 
  dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

results[[model]] <- data.save

####### *** P85-Ab+EBNA1-IgA -----------------------------
model <- "P85-Ab+EBNA1-IgA"
model.fit <- fit.gam_P85_EBNA1
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,
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
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
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
round(IECV.gam,3)

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

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.3f (%0.3f to %0.3f)",auc,LCI,UCI)) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.3f (%0.3f to %0.3f)",
                                                   calibration.intercept,
                                                   calibration.intercept_LCI,
                                                   calibration.intercept_UCI),
                     calibration.slope=sprintf("%0.3f (%0.3f to %0.3f)",
                                               calibration.slope,
                                               calibration.slope_LCI,
                                               calibration.slope_UCI)) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% 
  dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

results[[model]] <- data.save

####### *** triple-antibody strategy ---------------------
model <- "triple-antibody strategy"
model.fit <- fit.gam_P85_dual
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,
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
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
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
round(IECV.gam,3)

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

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.3f (%0.3f to %0.3f)",auc,LCI,UCI)) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.3f (%0.3f to %0.3f)",
                                                   calibration.intercept,
                                                   calibration.intercept_LCI,
                                                   calibration.intercept_UCI),
                     calibration.slope=sprintf("%0.3f (%0.3f to %0.3f)",
                                               calibration.slope,
                                               calibration.slope_LCI,
                                               calibration.slope_UCI)) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% 
  dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

results[[model]] <- data.save

####### *** combine --------------------------------
results_combine <- map_dfr(results,bind_rows)
save_path <- sprintf("%s/TableS3",table.path)
ensure_dir(save_path)
write.csv(results_combine,sprintf("%s/Five_centres_model_combined.csv",save_path))

###### ** three centers #########
data_IECV <- data_analysis
clusters <- unique(data_IECV$center)
N.clust <- length(clusters) # 3 clusters
data.in <- data.leftout <- list()
results <- list()

# create the datasets
for(i in 1:N.clust){
  data.in[[i]] <- data_IECV[data_IECV$center!=clusters[i],]
  data.leftout[[i]] <- data_IECV[data_IECV$center==clusters[i],]
}

####### *** P85-EBVDNA -----------------------
model <- "P85-EBVDNA strategy"
model.fit <- fit.gam_P85_EBVDNA
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,
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
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
    calibration.slope_LCI=calibration.slope.LCI,
    calibration.slope_UCI=calibration.slope.UCI,
    cluster=unique(data.leftout[[i]]$center)
  )
  
}

desired_order <- c("SYSUCC","WZRCH", "FJCH")
leftout.performance.gam.agg <- map_dfr(leftout.performance.gam,bind_rows) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)

leftout.prediction.gam.agg <- do.call(c, leftout.prediction.gam)
IECV.observed <- do.call(rbind, data.leftout)$group_ref
IECV.cluster <- do.call(rbind, data.leftout)$center
IECV.gam <- calculate_performance2(observed = IECV.observed,
                                   predicted = leftout.prediction.gam.agg)
round(IECV.gam,3)

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

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.3f (%0.3f to %0.3f)",auc,LCI,UCI)) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.3f (%0.3f to %0.3f)",
                                                   calibration.intercept,
                                                   calibration.intercept_LCI,
                                                   calibration.intercept_UCI),
                     calibration.slope=sprintf("%0.3f (%0.3f to %0.3f)",
                                               calibration.slope,
                                               calibration.slope_LCI,
                                               calibration.slope_UCI)) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% 
  dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

results[[model]] <- data.save

####### *** Combi-4 --------------------
model <- "Combi-4 strategy"
model.fit <- fit.gam_P85_dual_EBVDNA
leftout.prediction.gam <- list()
leftout.performance.gam <- list()
model.CV <- list()

for (i in 1:N.clust){
  model.CV[[i]] <- gam(model.fit$formula,
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
  calibration.intercept.LCI <- glm.fit$coef[1,1]-1.96*glm.fit$coef[1,2]
  calibration.intercept.UCI <- glm.fit$coef[1,1]+1.96*glm.fit$coef[1,2]
  
  calibration.slope <- glm.fit$coef[2,1] 
  calibration.slope.LCI <- glm.fit$coef[2,1]-1.96*glm.fit$coef[2,2]
  calibration.slope.UCI <- glm.fit$coef[2,1]+1.96*glm.fit$coef[2,2]
  
  leftout.performance.gam[[i]] <- tibble(
    AUC=AUC_value,
    AUC_LCI=LCI_value,
    AUC_UCI=UCI_value,
    calibration.intercept=calibration.intercept,
    calibration.intercept_LCI=calibration.intercept.LCI,
    calibration.intercept_UCI=calibration.intercept.UCI,
    calibration.slope=calibration.slope,
    calibration.slope_LCI=calibration.slope.LCI,
    calibration.slope_UCI=calibration.slope.UCI,
    cluster=unique(data.leftout[[i]]$center)
  )
  
}

desired_order <- c("SYSUCC","WZRCH", "FJCH")
leftout.performance.gam.agg <- map_dfr(leftout.performance.gam,bind_rows) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)

leftout.prediction.gam.agg <- do.call(c, leftout.prediction.gam)
IECV.observed <- do.call(rbind, data.leftout)$group_ref
IECV.cluster <- do.call(rbind, data.leftout)$center
IECV.gam <- calculate_performance2(observed = IECV.observed,
                                   predicted = leftout.prediction.gam.agg)
round(IECV.gam,3)

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

auc.save <- auc.clusters %>% 
  dplyr::select(-c(SE)) %>% 
  mutate(cluster = factor(cluster, levels = desired_order)) %>%
  arrange(cluster)
auc.save[, 1:3] <- round(auc.save[, 1:3], 3)
auc.save$model <- model

data.save <- auc.save %>% 
  mutate(AUC=sprintf("%0.3f (%0.3f to %0.3f)",auc,LCI,UCI)) %>% 
  left_join(leftout.performance.gam.agg %>% 
              mutate(calibration.intercept=sprintf("%0.3f (%0.3f to %0.3f)",
                                                   calibration.intercept,
                                                   calibration.intercept_LCI,
                                                   calibration.intercept_UCI),
                     calibration.slope=sprintf("%0.3f (%0.3f to %0.3f)",
                                               calibration.slope,
                                               calibration.slope_LCI,
                                               calibration.slope_UCI)) %>% 
              dplyr::select(cluster,calibration.intercept,calibration.slope)
  ) %>% 
  dplyr::select(model,cluster,AUC,calibration.intercept,calibration.slope)

results[[model]] <- data.save

####### *** combine -------------
results_combine <- map_dfr(results,bind_rows)
write.csv(results_combine,sprintf("%s/Three_centres_model_combined.csv",save_path))


#### TableS2 ###################################################################
indicators <- c("P85", "db_method3", "P85_VCA_cut","P85_EBNA1_cut","P85_dual_cut")
data_logistic$db_method3 <- factor(data_logistic$db_method3,levels = c("1","0"))
data_logistic$P85_VCA_cut <- factor(data_logistic$P85_VCA_cut,levels = c("1","0"))
data_logistic$P85_EBNA1_cut <- factor(data_logistic$P85_EBNA1_cut,levels = c("1","0"))
data_logistic$P85_dual_cut <- factor(data_logistic$P85_dual_cut,levels = c("1","0"))

data_pre <- data_logistic

result_df <- calcuate_diagnosis_index2(data=data_pre,
                                       indicators = indicators,
                                       binomMethod_use = "wilson",
                                       outcome = "group_ref")
result_df
save_path <- sprintf("%s/TableS2",table.path)
ensure_dir(save_path)
write.csv(result_df,sprintf("%s/TableS2_Pooled_analysis_model.csv",save_path),
          row.names = F)

result_compare <- DTC_pair(data=data_pre,
                           indicator1="P85",
                           compare_indicator=c("db_method3", "P85_VCA_cut",
                                               "P85_EBNA1_cut","P85_dual_cut"),
                           outcome="group_ref")
result_compare
write.csv(result_compare,
          sprintf("%s/TableS2_Pooled_analysis_model_compare.csv",save_path),
          row.names = F)

#### TableS4 ###################################################################
##### * five centers ############
#### five stages
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

indicators <- c("P85","db_method3", "P85_dual_cut") 
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
                             compare_indicator=c("db_method3", "P85_VCA_cut",
                                                 "P85_EBNA1_cut","P85_dual_cut"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

save_path <- sprintf("%s/TableS4",table.path)
ensure_dir(save_path)

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,sprintf("%s/pooled_analysis_Model_NPC_5_stages_combined.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,I,II,III,IVA,IVB)

write.csv(merged.data,sprintf("%s/pooled_analysis_Model_NPC_5_stages.csv",save_path),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined ,
          sprintf("%s/pooled_analysis_Model_NPC_5_stages_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/pooled_analysis_Model_NPC_5_stages_compare.csv",save_path),
          row.names = F)

#### three stages 
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
                             compare_indicator=c("db_method3","P85_dual_cut"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,sprintf("%s/pooled_analysis_Model_NPC_3_stages_combined.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity)
write.csv(merged.data,sprintf("%s/pooled_analysis_Model_NPC_3_stages.csv",save_path),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
write.csv(merged.statistical.data.combined,
          sprintf("%s/pooled_analysis_Model_NPC_3_stages_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/pooled_analysis_Model_NPC_3_stages_compare.csv",save_path),
          row.names = F)

##### * three centers ############
data_pre <- data_analysis %>%
  as.data.frame() %>%
  filter(NPC_stage %in% c("I","II","III","IVA","IVB","/")) %>%
  mutate(NPC_stage2 = case_when(
    NPC_stage %in% c("I","II") ~ "Early stage",
    NPC_stage %in% c("III","IVA") ~ "Localregionally advanced stage",
    NPC_stage %in% c("IVB") ~ "Advanced stage",
    NPC_stage %in% c("/") ~ "/"
  ))

stages <- c("I","II","III","IVA","IVB")
indicators <- c("P85", "P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut") 
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
                             compare_indicator=c("P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,
          sprintf("%s/Three_centers_Model_NPC_5_stages_combined.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity) %>% 
  dplyr::select(Indicator,I,II,III,IVA,IVB)
write.csv(merged.data,
          sprintf("%s/Three_centers_Model_NPC_5_stages.csv",save_path),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows)
write.csv(merged.statistical.data.combined,
          sprintf("%s/Three_centers_Model_NPC_5_stages_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/Three_centers_Model_NPC_5_stages_compare.csv",save_path),
          row.names = F)


#### three stages 
stages2 <- c("Early stage","Localregionally advanced stage","Advanced stage")
indicators <- c("P85", "P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut") 
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
                             compare_indicator=c("P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut"),
                             outcome="group_ref")
  result_compare$Stage <- current_NPC_stage
  result_compare <- result_compare %>% 
    dplyr::select(compare,Stage,sen1:sen.p.value)
  result_statistical_list[[i]] <- result_compare
  
}

merged.data.combined <- map_dfr(result_list,bind_rows)
write.csv(merged.data.combined,sprintf("%s/Three_centers_Model_NPC_3_stages_combined.csv",save_path),
          row.names = F)

merged.data <- map_dfr(result_list,bind_rows) %>% 
  pivot_wider(id_cols =-c(TP:FN),
              names_from = Stage,
              values_from = Sensitivity)
write.csv(merged.data,sprintf("%s/Three_centers_Model_NPC_3_stages.csv",save_path),
          row.names = F)

merged.statistical.data.combined <- map_dfr(result_statistical_list,bind_rows) 
write.csv(merged.statistical.data.combined,
          sprintf("%s/Three_centers_Model_NPC_3_stages_compare_combined.csv",save_path),
          row.names = F)

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) %>% 
  mutate(p_value=case_when(
    sen_method=="McNemar"~sen.p.value,
    sen_method=="McNemar mid-P"~str_c(sen.p.value,"a"),
  )) %>% 
  pivot_wider(id_cols = -c(sen1:sen.p.value),
              names_from = Stage,
              values_from = p_value)

write.csv(merged.statistical.data,
          sprintf("%s/Three_centers_Model_NPC_3_stages_compare.csv",save_path),
          row.names = F)



#### TableS5 ###################################################################
data_pre <- data_analysis
indicators <- c("P85", "P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut") 
result_df <- calcuate_diagnosis_index2(data=data_pre,
                                       indicators = indicators,
                                       binomMethod_use = "wilson",
                                       outcome = "group_ref")
result_df
save_path <- sprintf("%s/TableS5",table.path)
if(!dir.exists(save_path))dir.create(save_path,recursive = T)

write.csv(result_df,sprintf("%s/TableS5_three_centers_model.csv",save_path),
          row.names = F)

result_compare <- DTC_pair(data=data_pre,
                           indicator1="P85",
                           compare_indicator=c("P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut"),
                           outcome="group_ref")
result_compare
write.csv(result_compare,sprintf("%s/TableS5_three_centers_model_compare.csv",save_path),
          row.names = F)

result_compare1 <- DTC_pair(data=data_pre,
                            indicator1="P85_EBVDNA_cut",
                            compare_indicator=c("P85_dual_cut","P85_dual_EBVDNA_cut","P85"),
                            outcome="group_ref")
result_compare1
write.csv(result_compare1,sprintf("%s/TableS5_three_centers_PVE_vs_PVED.csv",save_path),
          row.names = F)

#### Table4 & TableS6 ##########################################################
symptom_data <- data_analysis %>% 
  dplyr::select(`Asymptomatic`,`NPC-specific Symptoms Cluster`,`NPC-nonspecific Symptoms Cluster`,
                `Neck mass`,`Nasal Symptoms`,`Aural Symptoms`,`Cranial Symptoms`,`Ophthalmic Symptoms`,
                `Neurological Deficits`,
                group,group_ref,P85,
                P85_dual_cut,P85_EBVDNA_cut,P85_dual_EBVDNA_cut)

# diagnostic performance
indicators <- c("P85", "P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut") 
candidate_symptom <- c('Asymptomatic','NPC-specific Symptoms Cluster','NPC-nonspecific Symptoms Cluster',
                       "Neck mass","Nasal Symptoms","Aural Symptoms",
                       "Cranial Symptoms","Ophthalmic Symptoms","Neurological Deficits")

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
                             compare_indicator=c("P85_dual_cut","P85_EBVDNA_cut","P85_dual_EBVDNA_cut"),
                             outcome="group_ref")
  result_compare$symptom <-  used_symptom
  result_statistical_list[[i]] <- result_compare
}

save_path <- sprintf("%s/Table4/",table.path)
ensure_dir(save_path)

merged.data <- map_dfr(result_list,bind_rows)
writexl::write_xlsx(merged.data,
                    sprintf("%s/Three_centers_Model_symptom_combined.xlsx",save_path))

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) 
writexl::write_xlsx(merged.statistical.data,
                    sprintf("%s/Three_centers_Model_symptom_compare_combined.xlsx",save_path))


#### TableS7 ###################################################################
symptom_data <- data_logistic %>% 
  dplyr::select(`Asymptomatic`,`NPC-specific Symptoms Cluster`,`NPC-nonspecific Symptoms Cluster`,
                group,group_ref,P85,P85_dual_cut)

# diagnostic performance
indicators <- c("P85","P85_dual_cut")
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
                             compare_indicator=c("P85_dual_cut"),
                             outcome="group_ref")
  result_compare$symptom <-  used_symptom
  result_statistical_list[[i]] <- result_compare
}

save_path <- sprintf("%s/TableS7/",table.path)
ensure_dir(save_path)

merged.data <- map_dfr(result_list,bind_rows)
writexl::write_xlsx(merged.data,
                    sprintf("%s/Five_centers_Model_symptom_combined.xlsx",save_path))

merged.statistical.data <- map_dfr(result_statistical_list,bind_rows) 
writexl::write_xlsx(merged.statistical.data,
                    sprintf("%s/Five_centers_Model_symptom_compare_combined.xlsx",save_path))


#### FigS12 ###################################################################
colnames(data_logistic)
roc_p85 <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$P85_SCO),
                     direction="<",
                     levels = c(0, 1))
value_p85 <- print_auc_and_ci(roc_p85);value_p85

## P85-Ab+vca-IgA(five centers)
roc_vca <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$gam_P85_VCA),
                     direction="<",
                     levels = c(0, 1))
value_vca <- print_auc_and_ci(roc_vca);value_vca

## P85-Ab+EBNA1-IgA (five centers)
roc_na1 <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$gam_P85_EBNA1),
                     direction="<",
                     levels = c(0, 1))
value_na1 <- print_auc_and_ci(roc_na1);value_na1

# VCA-IgA+EBNA1-IgA (five centers)
roc_db <- pROC::roc(data_logistic$group_ref, 
                    data_logistic$db_method_new_probability,
                    direction="<",levels=c(0, 1))
value_db <- print_auc_and_ci(roc_db);value_db

## P85-Ab+VCA-IgA+EBNA1-IgA (five centers)
roc_P85_db <- pROC::roc(data_logistic$group_ref, as.numeric(data_logistic$gam_P85_dual),
                        direction="<",levels = c(0, 1))
value_P85_db <- print_auc_and_ci(roc_P85_db);value_P85_db


## P85-Ab three centers
roc_P85_v2 <- pROC::roc(data_analysis$group_ref, 
                        as.numeric(data_analysis$P85_SCO),
                        direction="<",levels = c(0, 1))
value_P85_v2 <- print_auc_and_ci(roc_P85_v2);value_P85_v2

# P85-Ab+db (three centers)
roc_P85_db_v2 <- pROC::roc(data_analysis$group_ref,
                           as.numeric(data_analysis$gam_P85_dual),
                           direction="<",levels = c(0, 1))
value_P85_db_v2 <- print_auc_and_ci(roc_P85_db_v2);value_P85_db_v2

# P85-Ab+log EBV-DNA (three centers)
roc_ebvdna <- pROC::roc(data_analysis$group_ref,
                        as.numeric(data_analysis$gam_P85_EBVDNA),
                        direction="<",levels = c(0, 1))
value_ebvdna <- print_auc_and_ci(roc_ebvdna);value_ebvdna

## P85-Ab+logEBV-DNA+dual-antibody (three centers)
roc_dual_ebv <- pROC::roc(data_analysis$group_ref, 
                          as.numeric(data_analysis$gam_P85_dual_EBVDNA),
                          direction="<",levels = c(0, 1))
value_dual_ebv <- print_auc_and_ci(roc_dual_ebv);value_dual_ebv

save_path <- sprintf("%s/FigS12/",fig.path)
ensure_dir(save_path)

pdf(sprintf("%s/Fig4_Model_combination.pdf",save_path),
    width = 6,height = 6)
{
  plot(roc_p85, 
       col=my_color$single["P85-Ab"],
       legacy.axes=TRUE)
  
  plot.roc(roc_db,
           add=TRUE,
           col="#8cc8c4",
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  plot.roc(roc_vca,
           add=TRUE,
           col=my_color$single["VCA-IgA"],
           print.auc.x=0.35,
           print.auc.y=0.55) 
  
  plot.roc(roc_na1,
           add=TRUE, # 增加曲线
           col=my_color$single["EBNA1-IgA"],
           print.auc.x=0.35,
           print.auc.y=0.55) 
  
  
  plot.roc(roc_P85_db,
           add=TRUE, # 增加曲线
           col=my_color$single["EA-IgA"],
           print.auc.x=0.35,
           print.auc.y=0.55) 
  
  text<-c(paste0("P85-Ab\n",value_p85),
          paste0("VCA-IgA+EBNA1-IgA\n",value_db),
          paste0("P85-Ab+VCA-IgA\n",value_vca),
          paste0("P85-Ab+EBNA1-IgA\n",value_na1),
          paste0("P85-Ab+VCA-IgA+EBNA1-IgA\n",value_P85_db))
  
  legend(0.8, 0.5, 
         bty = "n", 
         legend=text, 
         col=my_color$model,
         lwd=3) 
}
dev.off()


roc.test(roc_p85, roc_vca)
roc.test(roc_p85, roc_na1)
roc.test(roc_p85, roc_db)
roc.test(roc_p85, roc_P85_db)

roc.test(roc_db, roc_vca)
roc.test(roc_db, roc_na1)
roc.test(roc_db, roc_P85_db)

roc.test(roc_vca, roc_na1)
roc.test(roc_vca, roc_P85_db)

roc.test(roc_na1,roc_P85_db)

pdf(sprintf("%s/Fig4_Model_combination_2.pdf",save_path),
    width = 6,height = 6)
{
  plot(roc_P85_v2, 
       col=my_color$single["P85-Ab"],
       legacy.axes=TRUE)
  
  plot.roc(roc_P85_db_v2,
           add=TRUE,
           col=my_color$single["EA-IgA"],
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  plot.roc(roc_ebvdna,
           add=TRUE, 
           col="#e29c45", 
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  plot.roc(roc_dual_ebv,
           add=TRUE,
           col="#6f2488",
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  text<-c(paste0("P85-Ab\n",value_P85_v2),
          paste0("P85-Ab+VCA-IgA+EBNA1-IgA\n",value_P85_db_v2),
          paste0("P85-Ab+EBV-DNA\n",value_ebvdna),
          paste0("P85-Ab+EBV-DNA+VCA-IgA+EBNA1-IgA\n",value_dual_ebv)
  )
  
  legend(0.8, 0.5,
         bty = "n",
         legend=text,
         col=c(my_color$single["P85-Ab"],my_color$single["EA-IgA"],"#e29c45","#6f2488"),
         lwd=3)
}
dev.off()

roc.test(roc_P85_v2, roc_P85_db_v2)
roc.test(roc_P85_v2, roc_ebvdna)
roc.test(roc_P85_v2, roc_dual_ebv)

roc.test(roc_P85_db_v2, roc_ebvdna)
roc.test(roc_P85_db_v2, roc_dual_ebv)

roc.test(roc_ebvdna,roc_dual_ebv)


symptom_data <- data_analysis %>% 
  dplyr::select(`NPC-specific Symptoms Cluster`,
                group,group_ref,P85_SCO,P85,
                gam_P85_dual,P85_dual_cut,
                gam_P85_EBVDNA,P85_EBVDNA_cut,
                gam_P85_dual_EBVDNA,P85_dual_EBVDNA_cut)

data_NPC_specific <- symptom_data %>% dplyr::filter(`NPC-specific Symptoms Cluster`=="1")

data_plot <- data_NPC_specific
roc_P85_v2 <- pROC::roc(data_plot$group_ref, 
                        as.numeric(data_plot$P85_SCO),
                        direction="<",levels = c(0, 1))
value_P85_v2 <- print_auc_and_ci(roc_P85_v2);value_P85_v2

# P85-Ab+db (three centers)
roc_P85_db_v2 <- pROC::roc(data_plot$group_ref,
                           as.numeric(data_plot$gam_P85_dual),
                           direction="<",levels = c(0, 1))
value_P85_db_v2 <- print_auc_and_ci(roc_P85_db_v2);value_P85_db_v2

# P85-Ab+log EBV-DNA (three centers)
roc_ebvdna <- pROC::roc(data_plot$group_ref,
                        as.numeric(data_plot$gam_P85_EBVDNA),
                        direction="<",levels = c(0, 1))
value_ebvdna <- print_auc_and_ci(roc_ebvdna);value_ebvdna

## P85-Ab+logEBV-DNA+dual-antibody (three centers)
roc_dual_ebv <- pROC::roc(data_plot$group_ref, 
                          as.numeric(data_plot$gam_P85_dual_EBVDNA),
                          direction="<",levels = c(0, 1))
value_dual_ebv <- print_auc_and_ci(roc_dual_ebv);value_dual_ebv

pdf(sprintf("%s/NPC-specific_roc.pdf",save_path),
    width = 6,height = 6)
{
  plot(roc_P85_v2, 
       col=my_color$single["P85-Ab"],
       legacy.axes=TRUE)
  
  plot.roc(roc_P85_db_v2,
           add=TRUE,
           col=my_color$single["EA-IgA"], 
           print.auc.x=0.35,
           print.auc.y=0.55) 
  
  plot.roc(roc_ebvdna,
           add=TRUE,
           col="#e29c45",
           print.auc.x=0.35,
           print.auc.y=0.55) 
  
  plot.roc(roc_dual_ebv,
           add=TRUE,
           col="#6f2488", 
           print.auc.x=0.35,
           print.auc.y=0.55)
  
  text<-c(paste0("P85-Ab\n",value_P85_v2),
          paste0("P85-Ab+VCA-IgA+EBNA1-IgA\n",value_P85_db_v2),
          paste0("P85-Ab+EBV-DNA\n",value_ebvdna),
          paste0("P85-Ab+EBV-DNA+VCA-IgA+EBNA1-IgA\n",value_dual_ebv)
  )
  
  legend(0.8, 0.5, 
         bty = "n", 
         legend=text, 
         col=c(my_color$single["P85-Ab"],my_color$single["EA-IgA"],"#e29c45","#6f2488"),
         lwd=3) 
}
dev.off()

roc.test(roc_P85_v2, roc_ebvdna)
roc.test(roc_P85_v2, roc_dual_ebv)
roc.test(roc_P85_v2, roc_P85_db_v2)

roc.test(roc_P85_db_v2, roc_ebvdna)
roc.test(roc_P85_db_v2, roc_dual_ebv)

roc.test(roc_ebvdna,roc_dual_ebv)