# api/app.R
# R Plumber API serving PRS and health data, generating automated reports,
# and hosting the system monitoring dashboard on port 8085.

library(plumber)
library(DBI)
library(RSQLite)
library(dplyr)

# Determine absolute project root dynamically
if (basename(getwd()) == "api") {
  project_root <- normalizePath(file.path(getwd(), ".."))
} else {
  project_root <- normalizePath(getwd())
}

db_path <- file.path(project_root, "database/bio_analysis.db")
model_results_path <- file.path(project_root, "Dataset/model_results.RData")

# Check database and model files
message("--- STARTUP DIAGNOSTICS ---")
message(paste("getwd():", getwd()))
message(paste("project_root:", project_root))
message(paste("db_path:", db_path))
message(paste("db_path exists:", file.exists(db_path)))
message(paste("model_results_path:", model_results_path))
message(paste("model_results_path exists:", file.exists(model_results_path)))
message("----------------------------")

#* @apiTitle Full Spectrum Biomedical Data API
#* @apiDescription REST API for querying individual genomics data (PRS), clinical parameters, and generating visual reports.

#* Serve the System Monitoring Dashboard UI
#* @get /dashboard
#* @serializer contentType list(type="text/html; charset=utf-8")
function(res) {
  dashboard_path <- file.path(project_root, "api/dashboard.html")
  if (!file.exists(dashboard_path)) {
    res$status <- 404
    return("<html><body><h3>Error: Dashboard UI template not found.</h3></body></html>")
  }
  
  html_content <- readChar(dashboard_path, file.info(dashboard_path)$size)
  return(html_content)
}

#* Get list of queryable sample IDs in the database
#* @get /samples
#* @serializer json
function() {
  if (!file.exists(db_path)) {
    return(list(error = "Database not found. Run the pipeline first."))
  }
  
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  
  samples <- dbGetQuery(con, "SELECT user_id FROM health_records ORDER BY user_id")
  return(samples$user_id)
}

#* Get system status and file diagnostics
#* @get /system_status
#* @serializer json
function() {
  files <- list(
    vcf_gz = file.exists(file.path(project_root, "Dataset/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")),
    vcf_ann = file.exists(file.path(project_root, "Dataset/Annotated.vcf.gz")),
    gwas_tsv = file.exists(file.path(project_root, "Dataset/gwas_catalog_associations.tsv")),
    health_raw = file.exists(file.path(project_root, "Dataset/health_raw.csv")),
    health_clean = file.exists(file.path(project_root, "Dataset/health_cleaned.csv")),
    prs_profile = file.exists(file.path(project_root, "Dataset/prs_out.profile")),
    model_rdata = file.exists(model_results_path),
    database = file.exists(db_path)
  )
  
  db_stats <- list(
    health_records_count = 0,
    prs_results_count = 0,
    last_run_id = "N/A",
    last_run_time = "N/A",
    last_run_status = "N/A"
  )
  
  if (files$database) {
    con <- tryCatch({
      dbConnect(SQLite(), db_path)
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(con)) {
      db_stats$health_records_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM health_records")$count
      db_stats$prs_results_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM prs_results")$count
      
      meta <- dbGetQuery(con, "SELECT * FROM pipeline_metadata ORDER BY executed_at DESC LIMIT 1")
      if (nrow(meta) > 0) {
        db_stats$last_run_id <- meta$run_id
        db_stats$last_run_time <- meta$executed_at
        db_stats$last_run_status <- meta$status
      }
      dbDisconnect(con)
    }
  }
  
  return(list(
    system_time = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
    project_root = project_root,
    files = files,
    database_stats = db_stats
  ))
}

#* Run Snakemake pipeline
#* @post /run_pipeline
#* @serializer json
function() {
  message("Running Snakemake pipeline via POST request...")
  # Reruns snakemake
  # Using system2 to capture output
  log_output <- system("snakemake --cores 1 --forceall 2>&1", intern = TRUE)
  
  # Check if database file exists and is updated
  db_exists <- file.exists(db_path)
  
  return(list(
    status = ifelse(db_exists, "SUCCESS", "FAILED"),
    log = log_output
  ))
}

#* Run system integrity tests
#* @get /run_tests
#* @serializer json
function() {
  test_runner_path <- file.path(project_root, "scripts/run_tests.R")
  if (!file.exists(test_runner_path)) {
    return(list(error = "Test runner script not found."))
  }
  
  # Source and execute the test runner function
  source(test_runner_path, local = TRUE)
  test_results <- run_test_suite()
  return(test_results)
}

#* Get patient PRS and health parameters
#* @param user_id The ID of the patient (e.g. HG00096)
#* @get /get_prs
#* @serializer json
function(user_id = "") {
  if (user_id == "") {
    return(list(error = "Missing user_id parameter. Please specify user_id (e.g. ?user_id=HG00096)"))
  }
  
  # Connect to database
  con <- dbConnect(SQLite(), db_path)
  on.exit(dbDisconnect(con))
  
  # 1. Fetch patient clinical record
  user_query <- dbSendQuery(con, "SELECT * FROM health_records WHERE user_id = ?")
  dbBind(user_query, list(user_id))
  user_data <- dbFetch(user_query)
  dbClearResult(user_query)
  
  if (nrow(user_data) == 0) {
    return(list(error = paste("Patient ID", user_id, "not found in database.")))
  }
  
  # 2. Fetch patient PRS score
  prs_query <- dbSendQuery(con, "SELECT * FROM prs_results WHERE user_id = ?")
  dbBind(prs_query, list(user_id))
  prs_data <- dbFetch(prs_query)
  dbClearResult(prs_query)
  
  # 3. Calculate PRS percentile in the population
  all_prs <- dbGetQuery(con, "SELECT prs_score FROM prs_results")$prs_score
  user_score <- prs_data$prs_score
  percentile <- round(mean(all_prs <= user_score) * 100, 1)
  
  # 4. Load models to predict cholesterol
  load(model_results_path) # Loads lm_model
  
  # Form input data frame for prediction
  pred_input <- data.frame(
    prs_standardized = prs_data$prs_standardized,
    age = user_data$age,
    gender = user_data$gender,
    bmi = user_data$bmi,
    stringsAsFactors = FALSE
  )
  
  pred_cholesterol <- round(predict(lm_model, newdata = pred_input), 1)
  
  # Compile response JSON
  response <- list(
    user_id = user_id,
    demographics = list(
      age = user_data$age,
      gender = user_data$gender
    ),
    clinical_metrics = list(
      bmi = user_data$bmi,
      blood_pressure = paste0(user_data$sbp, "/", user_data$dbp),
      systolic = user_data$sbp,
      diastolic = user_data$dbp,
      cholesterol_observed = user_data$cholesterol,
      glucose = user_data$glucose
    ),
    lifestyle = list(
      smoking_status = user_data$smoking_status,
      exercise_hours_weekly = user_data$exercise_hours
    ),
    genomic_prs = list(
      prs_raw_score = user_score,
      prs_standardized_score = prs_data$prs_standardized,
      percentile_rank = percentile,
      risk_category = if (percentile >= 80) "High Risk" else if (percentile <= 20) "Low Risk" else "Average"
    ),
    predictions = list(
      model_type = "Linear Regression Model",
      predicted_cholesterol = pred_cholesterol,
      deviation = round(user_data$cholesterol - pred_cholesterol, 1),
      implication = if (user_data$cholesterol > pred_cholesterol) {
        "Observed cholesterol is higher than predicted by clinical & genomic model. Lifestyle factors may play a larger role."
      } else {
        "Observed cholesterol is aligned with or lower than predicted by model."
      }
    )
  )
  
  return(response)
}

#* Get automated health report (HTML / PDF)
#* @param user_id The ID of the patient
#* @get /report
#* @serializer contentType list(type="text/html; charset=utf-8")
function(user_id = "", res) {
  if (user_id == "") {
    res$status <- 400
    return("<html><body><h3>Error: Missing user_id parameter.</h3></body></html>")
  }
  
  # Connect to database to verify user
  con <- dbConnect(SQLite(), db_path)
  user_query <- dbSendQuery(con, "SELECT user_id FROM health_records WHERE user_id = ?")
  dbBind(user_query, list(user_id))
  user_data <- dbFetch(user_query)
  dbClearResult(user_query)
  dbDisconnect(con)
  
  if (nrow(user_data) == 0) {
    res$status <- 404
    return(paste("<html><body><h3>Error: Patient ID", user_id, "not found in database.</h3></body></html>"))
  }
  
  # Create a temporary file to hold the rendered report
  tmp_file <- tempfile(fileext = ".html")
  
  # Render the R Markdown report template
  rmarkdown::render(
    input = file.path(project_root, "reports/report_template.Rmd"),
    output_file = tmp_file,
    params = list(user_id = user_id),
    quiet = TRUE
  )
  
  # Read and return the HTML file contents
  html_content <- readChar(tmp_file, file.info(tmp_file)$size)
  unlink(tmp_file)
  
  return(html_content)
}

# Run Plumber server on port 8085
if (sys.nframe() == 0) {
  pr <- plumber::plumb(file.path(project_root, "api/app.R"))
  pr$run(host = "0.0.0.0", port = 8085)
}
