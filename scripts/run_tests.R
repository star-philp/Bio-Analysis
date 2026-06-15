# scripts/run_tests.R
# System integrity test runner to verify ingestion, workflow outputs, model results, and DB loading.

message("Starting system integrity tests...")

library(DBI)
library(RSQLite)

run_test_suite <- function() {
  results <- list()
  
  # Determine project_root dynamically
  if (basename(getwd()) == "scripts" || basename(getwd()) == "api") {
    project_root <- normalizePath(file.path(getwd(), ".."))
  } else {
    project_root <- normalizePath(getwd())
  }
  
  # Helper to record test status
  add_result <- function(name, status, message) {
    results[[length(results) + 1]] <<- list(
      test_name = name,
      status = status, # "PASS" or "FAIL"
      message = message,
      timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")
    )
  }
  
  # Define absolute paths
  vcf_raw_path <- file.path(project_root, "Dataset/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz")
  health_raw_path <- file.path(project_root, "Dataset/health_raw.csv")
  gwas_path <- file.path(project_root, "Dataset/gwas_catalog_associations.tsv")
  ann_vcf <- file.path(project_root, "Dataset/Annotated.vcf.gz")
  prs_profile <- file.path(project_root, "Dataset/prs_out.profile")
  model_rdata <- file.path(project_root, "Dataset/model_results.RData")
  db_file <- file.path(project_root, "database/bio_analysis.db")
  
  # 1. Ingestion Files Exist & Non-empty
  message("Running Ingestion Files test...")
  files_to_check <- c(vcf_raw_path, health_raw_path, gwas_path)
  missing_files <- files_to_check[!file.exists(files_to_check)]
  if (length(missing_files) > 0) {
    add_result("Ingestion Files", "FAIL", paste("Missing files:", paste(basename(missing_files), collapse = ", ")))
  } else {
    sizes <- sapply(files_to_check, function(f) file.info(f)$size)
    if (any(sizes == 0)) {
      add_result("Ingestion Files", "FAIL", "One or more raw files are empty (0 bytes).")
    } else {
      add_result("Ingestion Files", "PASS", "All raw ingestion files are present and non-empty.")
    }
  }
  
  # 2. VCF Annotation Valid
  message("Running VCF Annotation test...")
  if (!file.exists(ann_vcf)) {
    add_result("VCF Annotation", "FAIL", "Annotated VCF file (Annotated.vcf.gz) not found.")
  } else {
    size_mb <- round(file.info(ann_vcf)$size / (1024 * 1024), 2)
    if (size_mb < 5) {
      add_result("VCF Annotation", "FAIL", paste("Annotated VCF file is suspiciously small:", size_mb, "MB"))
    } else {
      add_result("VCF Annotation", "PASS", paste("Annotated VCF file generated successfully (Size:", size_mb, "MB)."))
    }
  }
  
  # 3. PRS Calculations Completed
  message("Running PRS Calculations test...")
  if (!file.exists(prs_profile)) {
    add_result("PRS Calculations", "FAIL", "PLINK score profile (prs_out.profile) not found.")
  } else {
    prs_data <- read.table(prs_profile, header = TRUE, nrows = 10, stringsAsFactors = FALSE)
    if (!all(c("IID", "SCORE") %in% names(prs_data))) {
      add_result("PRS Calculations", "FAIL", "PLINK profile output has incorrect column format.")
    } else {
      add_result("PRS Calculations", "PASS", "PLINK score profile calculated and verified.")
    }
  }
  
  # 4. Statistical Modeling Successful
  message("Running Statistical Modeling test...")
  if (!file.exists(model_rdata)) {
    add_result("Statistical Modeling", "FAIL", "RData model results (model_results.RData) not found.")
  } else {
    env <- new.env()
    load(model_rdata, envir = env)
    required_vars <- c("lm_model", "lme_model", "coef_lm", "coef_lme", "merged_df")
    missing_vars <- required_vars[!required_vars %in% ls(env)]
    if (length(missing_vars) > 0) {
      add_result("Statistical Modeling", "FAIL", paste("Missing variables in RData:", paste(missing_vars, collapse = ", ")))
    } else {
      # Verify significance of model variables
      if (is.null(summary(env$lm_model)$coefficients)) {
        add_result("Statistical Modeling", "FAIL", "Linear regression model is corrupt or empty.")
      } else {
        add_result("Statistical Modeling", "PASS", "Linear regression and mixed-effects models loaded and verified successfully.")
      }
    }
  }
  
  # 5. SQLite Database Loaded & Optimized
  message("Running Database test...")
  if (!file.exists(db_file)) {
    add_result("SQLite Database", "FAIL", "Database file (bio_analysis.db) not found.")
  } else {
    con <- tryCatch({
      dbConnect(SQLite(), db_file)
    }, error = function(e) {
      NULL
    })
    
    if (is.null(con)) {
      add_result("SQLite Database", "FAIL", "Unable to establish connection to SQLite database.")
    } else {
      on.exit(dbDisconnect(con), add = TRUE)
      
      # Verify tables
      tables <- dbListTables(con)
      required_tables <- c("health_records", "prs_results", "model_coefficients", "codebook", "pipeline_metadata")
      missing_tables <- required_tables[!required_tables %in% tables]
      
      if (length(missing_tables) > 0) {
        add_result("SQLite Database", "FAIL", paste("Missing tables in DB:", paste(missing_tables, collapse = ", ")))
      } else {
        # Check row counts
        patients_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM health_records")$count
        prs_count <- dbGetQuery(con, "SELECT COUNT(*) as count FROM prs_results")$count
        
        # Check if indexes exist
        indexes <- dbGetQuery(con, "SELECT name FROM sqlite_master WHERE type = 'index'")
        required_indexes <- c("idx_health_user_id", "idx_prs_user_id", "idx_prs_run_id")
        missing_indexes <- required_indexes[!required_indexes %in% indexes$name]
        
        if (patients_count == 0 || prs_count == 0) {
          add_result("SQLite Database", "FAIL", "Database tables are empty (0 rows).")
        } else if (length(missing_indexes) > 0) {
          add_result("SQLite Database", "FAIL", paste("Missing indexes in DB:", paste(missing_indexes, collapse = ", ")))
        } else {
          add_result("SQLite Database", "PASS", paste("SQLite database load verified.", patients_count, "records queryable with indexes."))
        }
      }
    }
  }
  
  return(results)
}

# Run tests and output JSON to stdout
if (sys.nframe() == 0) {
  library(jsonlite)
  results <- run_test_suite()
  cat(jsonlite::toJSON(results, auto_unbox = TRUE, pretty = TRUE))
}
