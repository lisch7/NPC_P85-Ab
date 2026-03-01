# What is the study?
Anti-Epstein–Barr virus (EBV) BNLF2b total antibody (P85-Ab) has shown promise in nasopharyngeal carcinoma (NPC) screening. However, its diagnostic performance in outpatient settings remains unclear. In this prospective, multicenter, head-to-head comparison cohort study, we compared the diagnostic performance of P85-Ab with that of VCA-IgA, EA-IgA, and EBNA1-IgA in outpatient clinics. Also, we performed an exploratory analysis comparing the diagnostic performance of P85-Ab with plasma EBV DNA.

# Key findings of the study
- P85-Ab demonstrated the best diagnostic performance among the EBV-related biomarkers assessed for NPC detection.
- P85-Ab exhibited the highest sensitivity in detecting early-stage NPC.
- P85-Ab maintained robust specificity across different disease categories, particularly in the differential diagnosis of other EBV-related diseases.
- P85-Ab alone is suggested for populations with asymptomatic or NPC-nonspecific symptoms, and combining P85-Ab with VCA-IgA and EBNA1-IgA is suggested for populations with NPC-specific symptoms.

# What is this repository?
This repository contains the R code used to generate our results and the trained models.

## Example: Using the Trained Model
To use the trained model, it is recommended to prepare an Excel file. The column names should include the following variables: `P85_SCO`, `VCA_IgA_SCO`, and`NA1_IgA_SCO`. Each row should represent a single participant.

- P85_SCO: Cut-off index (COI) for the P85 antibody. This is a continuous variable.
- VCA_IgA_SCO: ratio of optical density to reference control (rOD) of VCA-IgA. This is a continuous variable.
- NA1_IgA_SCO: rOD of EBNA1-IgA. This is a continuous variable.

```R
library(mgcv)

# Load the pretrained model
model <- readRDS("./Model.rds")

# Check available strategies in the model object
names(model)
# [1] "P85-Ab+VCA-IgA"  "P85-Ab+EBNA1-IgA"  "Triplet-Antibody Strategy"

# Prepare test data
test_data <- data.frame(
  P85_SCO = c(10, 0),
  VCA_IgA_SCO = c(5, 0),
  NA1_IgA_SCO = c(5, 0)
)

# Alternatively, read from an Excel file
# Make sure the column names match exactly as required
# Example file: "example_input.xlsx" with variables in the first sheet
# test_data <- readxl::read_excel("example_input.xlsx")

# Display the test data
print(test_data)
# P85_SCO VCA_IgA_SCO NA1_IgA_SCO
# 1      10           5           5
# 2       0           0           0

# Run prediction using the Triplet-Antibody Strategy
results <- mgcv::predict.gam(
  object = model$`Triplet-Antibody Strategy`$model,
  newdata = test_data,
  type="response"
)
print(results)
# 1           2 
# 0.975855667 0.005126629 

# Classify based on the cutoff
group_infer <- ifelse(
  results > model$`Triplet-Antibody Strategy`$cutoff,
  "NPC", "nonNPC"
)
print(group_infer)
#     1        2 
# "NPC" "nonNPC" 
```




# Contact us
- Su-Chen Li, lisc@sysucc.org.cn
- Lin-Quan Tang, tanglq@sysucc.org.cn
