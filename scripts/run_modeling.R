# scripts/run_modeling.R
# This script performs statistical modeling (linear regression and mixed models)
# to analyze the association between the Polygenic Risk Score (PRS) and health metrics.

message("Starting statistical modeling in R...")

# Define paths
health_path <- "Dataset/health_cleaned.csv"
prs_path <- "Dataset/prs_out.profile"
model_output_path <- "Dataset/model_results.RData"
summary_txt_path <- "Dataset/model_summary.txt"

# Load libraries
library(nlme)
library(dplyr)

# Read inputs
if (!file.exists(health_path)) {
  stop(paste("Cleaned health data not found at", health_path))
}
if (!file.exists(prs_path)) {
  stop(paste("PLINK PRS profile not found at", prs_path))
}

health_df <- read.csv(health_path, stringsAsFactors = FALSE)
# PLINK profile output is space-separated, we read it using read.table
prs_df <- read.table(prs_path, header = TRUE, stringsAsFactors = FALSE)

message(paste("Loaded", nrow(health_df), "health records and", nrow(prs_df), "PRS profiles."))

# Merge on user_id / IID
# In VCF/PLINK, the individual ID is in the IID column. In health data, it is user_id.
merged_df <- inner_join(health_df, prs_df, by = c("user_id" = "IID"))
message(paste("Merged dataset contains", nrow(merged_df), "records."))

if (nrow(merged_df) == 0) {
  stop("Merged dataset is empty! Check that VCF sample IDs match health raw user_ids.")
}

# Standardize PRS score for easier interpretation (mean=0, sd=1)
merged_df$prs_standardized <- scale(merged_df$SCORE)[,1]

# 1. Linear Regression Model: Predict cholesterol using PRS, age, gender, and BMI
message("Fitting linear regression model (cholesterol ~ prs + age + gender + bmi)...")
lm_model <- lm(cholesterol ~ prs_standardized + age + gender + bmi, data = merged_df)
lm_summary <- summary(lm_model)
print(lm_summary)

# 2. Mixed-effects Model (using nlme):
# We group patients into age decades as a grouping factor (random intercept)
message("Fitting mixed-effects model with random intercepts by age group (decade)...")
merged_df$age_group <- paste0(floor(merged_df$age / 10) * 10, "s")

# Fit mixed model using nlme::lme
lme_model <- lme(
  fixed = cholesterol ~ prs_standardized + bmi + gender,
  random = ~ 1 | age_group,
  data = merged_df,
  method = "REML"
)
lme_summary <- summary(lme_model)
print(lme_summary)

# 3. Save modeling results
# We extract coefficients for the DB loading script
coef_lm <- as.data.frame(lm_summary$coefficients)
names(coef_lm) <- c("Estimate", "StdError", "tValue", "pValue")
coef_lm$Variable <- rownames(coef_lm)
rownames(coef_lm) <- NULL

# Extract mixed model fixed effects
fixed_effects <- summary(lme_model)$tTable
coef_lme <- as.data.frame(fixed_effects)
names(coef_lme) <- c("Estimate", "StdError", "DF", "tValue", "pValue")
coef_lme$Variable <- rownames(coef_lme)
rownames(coef_lme) <- NULL

# Save summaries to a text file for documentation
sink(summary_txt_path)
cat("========================================================\n")
cat("LINEAR REGRESSION MODEL: CHOLESTEROL PREDICTION\n")
cat("========================================================\n")
print(lm_summary)
cat("\n\n========================================================\n")
cat("MIXED EFFECTS MODEL: RANDOM INTERCEPT BY AGE DECADE\n")
cat("========================================================\n")
print(lme_summary)
sink()

# Save objects to RData for database ingestion and plumber serving
save(lm_model, lme_model, coef_lm, coef_lme, merged_df, file = model_output_path)
message(paste("Saved statistical model results to RData:", model_output_path))
message(paste("Saved model text summaries to:", summary_txt_path))
message("Statistical modeling completed successfully!")
