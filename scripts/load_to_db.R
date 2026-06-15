# scripts/load_to_db.R
# This script loads all analyzed data, model results, and metadata into a standardized SQLite database.

message("Starting database loading...")

# Load libraries
library(DBI)
library(RSQLite)
library(dplyr)

# Define paths
db_path <- "database/bio_analysis.db"
model_results_path <- "Dataset/model_results.RData"
health_cleaned_path <- "Dataset/health_cleaned.csv"

# Ensure directory exists
dir.create(dirname(db_path), showWarnings = FALSE, recursive = TRUE)

# Connect to database
con <- dbConnect(SQLite(), db_path)
message(paste("Connected to SQLite database at", db_path))

# Enable foreign keys
dbExecute(con, "PRAGMA foreign_keys = ON;")

# 1. Create tables
message("Creating tables...")

# Table: pipeline_metadata
dbExecute(con, "
CREATE TABLE IF NOT EXISTS pipeline_metadata (
    run_id TEXT PRIMARY KEY,
    version TEXT NOT NULL,
    executed_at TEXT NOT NULL,
    status TEXT NOT NULL,
    comments TEXT
);
")

# Table: codebook
dbExecute(con, "
CREATE TABLE IF NOT EXISTS codebook (
    column_name TEXT PRIMARY KEY,
    description TEXT NOT NULL,
    data_type TEXT NOT NULL,
    category TEXT NOT NULL
);
")

# Table: health_records
dbExecute(con, "
CREATE TABLE IF NOT EXISTS health_records (
    user_id TEXT PRIMARY KEY,
    age INTEGER NOT NULL,
    gender TEXT NOT NULL,
    bmi REAL NOT NULL,
    sbp REAL NOT NULL,
    dbp REAL NOT NULL,
    cholesterol REAL NOT NULL,
    glucose REAL NOT NULL,
    smoking_status TEXT NOT NULL,
    exercise_hours REAL NOT NULL,
    clean_timestamp TEXT NOT NULL
);
")

# Table: prs_results
dbExecute(con, "
CREATE TABLE IF NOT EXISTS prs_results (
    user_id TEXT PRIMARY KEY,
    prs_score REAL NOT NULL,
    prs_standardized REAL NOT NULL,
    run_id TEXT NOT NULL,
    FOREIGN KEY(user_id) REFERENCES health_records(user_id),
    FOREIGN KEY(run_id) REFERENCES pipeline_metadata(run_id)
);
")

# Table: model_coefficients
dbExecute(con, "
CREATE TABLE IF NOT EXISTS model_coefficients (
    model_type TEXT NOT NULL,
    variable TEXT NOT NULL,
    estimate REAL NOT NULL,
    std_error REAL NOT NULL,
    t_value REAL NOT NULL,
    p_value REAL NOT NULL,
    PRIMARY KEY(model_type, variable)
);
")

# 2. Insert metadata & codebook
run_id <- paste0("run_", format(Sys.time(), "%Y%m%d_%H%M%S"))
version <- "v1.0.0"
executed_at <- format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")

dbExecute(con, "
INSERT OR REPLACE INTO pipeline_metadata (run_id, version, executed_at, status, comments)
VALUES (?, ?, ?, ?, ?);
", list(run_id, version, executed_at, "SUCCESS", "Full Spectrum Bio-Analysis pipeline run"))

# Codebook entries
codebook_data <- data.frame(
  column_name = c("user_id", "age", "gender", "bmi", "sbp", "dbp", "cholesterol", "glucose", "smoking_status", "exercise_hours", "prs_score", "prs_standardized"),
  description = c(
    "Unique identifier for the patient matching 1000 Genomes VCF samples",
    "Age of the patient in years",
    "Biological gender (M: Male, F: Female)",
    "Body Mass Index (kg/m^2)",
    "Systolic blood pressure (mmHg)",
    "Diastolic blood pressure (mmHg)",
    "Total cholesterol (mg/dL)",
    "Fasting blood glucose (mg/dL)",
    "Smoking status (Never, Former, Current)",
    "Weekly physical exercise in hours",
    "Raw Polygenic Risk Score calculated by PLINK",
    "Standardized Polygenic Risk Score (mean=0, sd=1)"
  ),
  data_type = c("TEXT", "INTEGER", "TEXT", "REAL", "REAL", "REAL", "REAL", "REAL", "TEXT", "REAL", "REAL", "REAL"),
  category = c("Demographics", "Demographics", "Demographics", "Clinical", "Clinical", "Clinical", "Clinical", "Clinical", "Lifestyle", "Lifestyle", "Genomic", "Genomic"),
  stringsAsFactors = FALSE
)

dbWriteTable(con, "codebook", codebook_data, append = TRUE, row.names = FALSE)

# 3. Load clean health data and model results
load(model_results_path) # Loads merged_df, coef_lm, coef_lme

# Insert health records
health_db_df <- merged_df %>%
  select(user_id, age, gender, bmi, sbp, dbp, cholesterol, glucose, smoking_status, exercise_hours) %>%
  mutate(clean_timestamp = executed_at)

dbWriteTable(con, "health_records", health_db_df, append = TRUE, row.names = FALSE)

# Insert PRS results
prs_db_df <- merged_df %>%
  select(user_id, SCORE, prs_standardized) %>%
  rename(prs_score = SCORE) %>%
  mutate(run_id = !!run_id)

dbWriteTable(con, "prs_results", prs_db_df, append = TRUE, row.names = FALSE)

# Insert model coefficients
coef_lm_db <- coef_lm %>%
  select(Variable, Estimate, StdError, tValue, pValue) %>%
  rename(variable = Variable, estimate = Estimate, std_error = StdError, t_value = tValue, p_value = pValue) %>%
  mutate(model_type = "linear_regression")

coef_lme_db <- coef_lme %>%
  select(Variable, Estimate, StdError, tValue, pValue) %>%
  rename(variable = Variable, estimate = Estimate, std_error = StdError, t_value = tValue, p_value = pValue) %>%
  mutate(model_type = "mixed_effects")

dbWriteTable(con, "model_coefficients", rbind(coef_lm_db, coef_lme_db), append = TRUE, row.names = FALSE)

# 4. Create Indexes for optimization
message("Creating indexes for query optimization...")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_health_user_id ON health_records(user_id);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_prs_user_id ON prs_results(user_id);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_prs_run_id ON prs_results(run_id);")

# 5. Verify database contents
message("Verifying DB contents...")
tot_patients <- dbGetQuery(con, "SELECT COUNT(*) as count FROM health_records")$count
tot_prs <- dbGetQuery(con, "SELECT COUNT(*) as count FROM prs_results")$count
tot_coefs <- dbGetQuery(con, "SELECT COUNT(*) as count FROM model_coefficients")$count

message(paste("Database load verified!"))
message(paste("  - health_records rows:", tot_patients))
message(paste("  - prs_results rows:", tot_prs))
message(paste("  - model_coefficients rows:", tot_coefs))
message(paste("  - Current Run ID:", run_id))

# Disconnect
dbDisconnect(con)
message("Database connection closed. DB load completed successfully!")
