# scripts/clean_health_data.R
# This script cleans the raw health screening data to prepare it for statistical modeling.

message("Starting raw health data cleaning...")

# Define paths
input_path <- "Dataset/health_raw.csv"
output_path <- "Dataset/health_cleaned.csv"

# Load libraries
library(dplyr)

# Read raw data
if (!file.exists(input_path)) {
  stop(paste("Raw health data not found at", input_path))
}
df <- read.csv(input_path, stringsAsFactors = FALSE)
initial_rows <- nrow(df)
message(paste("Loaded raw health data containing", initial_rows, "records."))

# 1. Clean Gender
# Standardize gender: "male" / "M" -> "M", "female" / "F" -> "F", uppercase everything
df$gender <- toupper(df$gender)
df$gender[df$gender %in% c("MALE", "M")] <- "M"
df$gender[df$gender %in% c("FEMALE", "F")] <- "F"

# Drop rows with invalid gender (not M or F)
df <- df %>% filter(gender %in% c("M", "F"))
gender_dropped <- initial_rows - nrow(df)
if (gender_dropped > 0) {
  message(paste("Dropped", gender_dropped, "rows due to invalid gender."))
}

# 2. Clean Age
# Invalid ages: age > 120 or age <= 0
# For testing, we drop these rows
valid_age_df <- df %>% filter(age > 0 & age <= 120)
age_dropped <- nrow(df) - nrow(valid_age_df)
df <- valid_age_df
if (age_dropped > 0) {
  message(paste("Dropped", age_dropped, "rows due to invalid age (>120 or <=0)."))
}

# 3. Clean BMI
# Negative or 0 BMI values are invalid. Let's replace them with the median of valid BMI values.
valid_bmis <- df$bmi[df$bmi > 0 & df$bmi < 100]
median_bmi <- median(valid_bmis, na.rm = TRUE)
invalid_bmi_count <- sum(df$bmi <= 0 | df$bmi >= 100 | is.na(df$bmi))
df$bmi[df$bmi <= 0 | df$bmi >= 100 | is.na(df$bmi)] <- median_bmi
message(paste("Imputed", invalid_bmi_count, "invalid BMI values with median BMI:", median_bmi))

# 4. Clean & Impute Cholesterol
# Missing cholesterol values (NAs) will be imputed using median cholesterol
valid_chol <- df$cholesterol[!is.na(df$cholesterol)]
median_chol <- median(valid_chol, na.rm = TRUE)
missing_chol_count <- sum(is.na(df$cholesterol))
df$cholesterol[is.na(df$cholesterol)] <- median_chol
message(paste("Imputed", missing_chol_count, "missing cholesterol values with median cholesterol:", median_chol))

# Save cleaned data
write.csv(df, output_path, row.names = FALSE, quote = FALSE)
message(paste("Cleaned data saved to", output_path))
message(paste("Final dataset contains", nrow(df), "records (dropped", initial_rows - nrow(df), "rows total)."))
