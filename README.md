# What is the study?
Anti-Epstein–Barr virus (EBV) BNLF2b (P85-Ab) has shown promise in nasopharyngeal carcinoma (NPC) screening. However, its diagnostic performance in outpatient settings remains unclear. In this prospective, multicenter, head-to-head comparison cohort study, we compared the diagnostic performance of P85-Ab with that of VCA-IgA, EA-IgA, EBNA1-IgA, and plasma EBV-DNA in outpatient clinics.

# Key findings of the study
- P85-Ab demonstrated the best diagnostic performance among the EBV-related biomarkers assessed for NPC detection.
- P85-Ab exhibited the highest sensitivity in detecting early-stage NPC.
- P85-Ab maintained robust specificity across different disease categories, particularly in the differential diagnosis of other EBV-related diseases.
- We explored optimal combination strategy and proposed a clinical application strategy using P85-Ab, VCA-IgA, EA-IgA, EBNA1-IgA, and plasma EBV-DNA.

# What is this repository?
This repository contains the R code used to generate our results and the trained models.

## Example: Using the Trained Model
To use the trained model, it is recommended to prepare an Excel file. The column names should include the following variables: `P85_SCO`, `VCA_IgA_SCO`, `NA1_IgA_SCO`, and `log_EBV_DNA`. Each row should represent a single participant.
**Note:** The required variables may vary depending on the specific model used—for example, some models only require two variables.
The detailed descriptions of these variables are as follows:
- P85_SCO: Cut-off index (COI) for the P85 antibody. This is a continuous variable.
- VCA_IgA_SCO: ratio of optical density to reference control (rOD) of VCA-IgA. This is a continuous variable.
- NA1_IgA_SCO: rOD of EBNA1-IgA. This is a continuous variable.
- log_EBV_DNA: Log-transformed value of plasma EBV-DNA load (log₁₀ of EBV-DNA + 1). This is a continuous variable.

```
library(mgcv)

# Load the pretrained model
model <- readRDS("./Model.rds")

# Check available strategies in the model object
names(model)
# [1] "P85-EBVDNA Strategy"      "Triple-Antibody Strategy" "Combi-4 Strategy"

# Prepare test data
test_data <- data.frame(
  P85_SCO = c(10, 0),
  VCA_IgA_SCO = c(5, 0),
  NA1_IgA_SCO = c(5, 0),
  log_EBV_DNA = c(2, 0)
)

# Alternatively, read from an Excel file
# Make sure the column names match exactly as required
# Example file: "example_input.xlsx" with variables in the first sheet
# test_data <- readxl::read_excel("example_input.xlsx")

# Display the test data
print(test_data)
#   P85_SCO VCA_IgA_SCO NA1_IgA_SCO log_EBV_DNA
# 1      10           5           5           2
# 2       0           0           0           0

# Run prediction using the Combi-4 Strategy
results <- mgcv::predict.gam(
  object = model$`Combi-4 Strategy`$model,
  newdata = test_data
)
print(results)
#         1         2 
#  5.312172 -5.397576 

# Classify based on the cutoff
group_infer <- ifelse(
  results > model$`Combi-4 Strategy`$cutoff,
  "NPC", "nonNPC"
)
print(group_infer)
#   1       2 
# "NPC" "nonNPC" 
```




# Contact us
- Su-Chen Li, lisc@sysucc.org.cn
- Lin-Quan Tang, tanglq@sysucc.org.cn
