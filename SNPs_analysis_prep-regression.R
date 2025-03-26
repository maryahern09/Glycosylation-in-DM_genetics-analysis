#Code to Examine Glycogenes SNPs in DM Pt's
#created by Mary 3/18/25

#install packages
install.packages("data.table")
library(data.table)
library(readxl)
library(tidyverse)

#load data
data_SNPs = fread("GlycogenesSNPsOUT.raw")
data_phenotype = read_excel("Data Phenotype for Transfer 1.15.25-Sri (1).xlsx")

#examine data SNPs
table(data_SNPs$SEX)
table(data_SNPs$PHENOTYPE)
table(data_SNPs$MAT)
table(data_SNPs$PAT)
# looks like there is basically nothing in these variables, they are all the same

#merge datasets by IID
data_phenotype <- data_phenotype %>%
  rename(IID = biobank_id)

data_merged <- merge(data_phenotype, data_SNPs, by = "IID")

any(is.na(data_merged[, 18:ncol(data_merged)]))

table(data_merged$rs1405919033_C)
data_SNPs$rs1405919033_C

### Quality control- remove individuals with more than 2% missing and SNPs with more than 2% missing
#%missing for observations
for(i in 1:1026){
  data_merged$percentmissing[i] <-sum(is.na(data_merged[i,])/294442)*100
  print(i)
}

summary(data_merged$percentmissing)#none were above 2% so the dataset is good

#MAF
for(i in 1:5284) {
  data_map2$MAF[i] <- (mean(data3[, i + 12], na.rm = TRUE)) / 2
  print(i)
}

#prep dataset
data_merged <- data_merged %>%
  select(-c("FID", "PAT", "MAT", "SEX", "PHENOTYPE"))

#create DM variable
data_final <- data_merged %>%
  mutate(DM = case_when(
    hgb_a1c_glyco <5.7 ~ "Normal",
    hgb_a1c_glyco >=5.7 & hgb_a1c_glyco <=6.4 ~ "Pre-DM",
    hgb_a1c_glyco >6.5 ~ "Diabetes",
    TRUE ~ NA_character_
  ))
  
#create binary DM variable
data_final <- data_final %>%
  mutate(DM_binary = case_when(
    DM %in% c("Normal","Pre-DM") ~ 0,
    DM == "Diabetes" ~ 1,
    TRUE ~ NA_real_
  ))


#what are the allele frequencies for DM, preDM, and Normal
#rename gender category
data_final <- data_final %>%
  rename(gender_2 = "gender_2 (Female is 0)")

#remove all SNP columns that have all of the same values
# Identify columns (SNPs) with only one unique value
columns_to_remove <- sapply(data_final[, 19:294456], function(x) length(unique(x)) < 2)

# Remove those SNPs from the dataset
data_final_removed <- data_final[, !columns_to_remove]

# Check how many SNPs were removed
sum(columns_to_remove)

#remove SNPs where results are all NA
data_final_removed <- data_final_removed %>%
  select(where(~ !all(is.na(.))))

#logistic regression
#DM status
# Create blank OUT file
OUT <- as.data.frame(matrix(NA, nrow = 1026, ncol = 5))
colnames(OUT) <- c("SNP", "Coef", "SE", "OR", "PValue")

# Loop through all SNPs
for (i in 19:281980) {  # Assuming first column is ID
  SNP <- data_final_removed[, i]
  model <- glm(DM_binary ~ SNP + gender_2 + age_2, data = data_final_removed, family = binomial(link = "logit"))
  coef_vals <- summary(model)$coef
  OUT[i, ] <- c(colnames(data_final)[i],  # SNP name
                   coef_vals[2, 1],                      # Coefficient
                   coef_vals[2, 2],                      # Standard Error
                   exp(coef_vals[2, 1]),                 # Odds Ratio
                   coef_vals[2, 4])                      # P-value
}

# Save OUT file
write.table(OUT, "Diabetes_SNP_Logistic Regression_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")

#linear regression
#A1c continuous
# Create blank OUT file
OUT_linear <- as.data.frame(matrix(NA, nrow = 1026, ncol = 4))
colnames(OUT_linear) <- c("SNP", "Coef", "SE", "OR", "PValue")

# Loop through all SNPs
for (i in 19:294456) {  # Assuming first column is ID
  SNP <- data_final_removed[, i]
  model <- glm(hgb_a1c_glyco ~ SNP + gender_2 + age_2, data = data_final, family = gaussian)
  coef_vals <- summary(model)$coef
  OUT_linear[i, ] <- c(colnames(data_final)[i],  # SNP name
                coef_vals[2, 1],                      # Coefficient
                coef_vals[2, 2],                      # Standard Error
                coef_vals[2, 4])                      # P-value
}

# Save OUT file
write.table(OUT_linear, "Diabetes_SNP_Linear Regression_results.txt", quote = FALSE, row.names = FALSE, sep = "\t")
