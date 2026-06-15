# scripts/generate_mock_data.R
# This script generates mock health data and mock GWAS weights for pipeline testing.

message("Starting mock data generation...")

# Define paths
vcf_path <- "Dataset/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
health_raw_path <- "Dataset/health_raw.csv"
gwas_path <- "Dataset/gwas_catalog_associations.tsv"

# 1. Fetch sample IDs from VCF using bcftools
message("Fetching VCF samples...")
samples <- system(paste("bcftools query -l", vcf_path), intern = TRUE)
num_samples <- length(samples)
message(paste("Found", num_samples, "samples in VCF."))

if (num_samples == 0) {
  stop("No samples found in VCF. Check if bcftools is installed and path is correct.")
}

# 2. Fetch some variants from VCF to construct a realistic GWAS weights file
message("Fetching variants for GWAS weights...")
# We query a region of chr 22 and take 50 variants
variants_raw <- system(
  paste("bcftools query -r 22:16050000-17000000 -f '%CHROM\\t%POS\\t%REF\\t%ALT\\n'", vcf_path, "| head -n 50"),
  intern = TRUE
)

if (length(variants_raw) == 0) {
  stop("No variants found in VCF region.")
}

variants_df <- do.call(rbind, lapply(strsplit(variants_raw, "\t"), function(x) {
  data.frame(CHROM = x[1], POS = as.integer(x[2]), REF = x[3], ALT = x[4], stringsAsFactors = FALSE)
}))

# Filter variants to keep only simple SNPs (REF and ALT length = 1)
variants_df <- variants_df[nchar(variants_df$REF) == 1 & nchar(variants_df$ALT) == 1, ]

if (nrow(variants_df) < 10) {
  stop("Too few single nucleotide variants found in VCF region.")
}

# Keep top 20 variants for the GWAS panel
gwas_variants <- head(variants_df, 20)
gwas_variants$panel_variant_id <- paste(gwas_variants$CHROM, gwas_variants$POS, gwas_variants$REF, gwas_variants$ALT, sep = ":")
gwas_variants$rsid <- paste0("rs_mock_", gwas_variants$POS)
gwas_variants$risk_allele <- gwas_variants$ALT
gwas_variants$non_risk_allele <- gwas_variants$REF

# Generate random effect sizes (beta) and p-values
set.seed(42)
gwas_variants$beta <- round(rnorm(nrow(gwas_variants), mean = 0.2, sd = 0.1), 3)
# Make sure some beta are negative, some positive
gwas_variants$beta[seq(2, nrow(gwas_variants), by=2)] <- -gwas_variants$beta[seq(2, nrow(gwas_variants), by=2)]
gwas_variants$p_value <- format(runif(nrow(gwas_variants), min = 1e-15, max = 5e-8), scientific = TRUE)

# Write GWAS Catalog mockup file (TSV)
# In PLINK, the file format will be: variant_id, risk_allele, beta
# We write the full TSV as mock GWAS catalog metadata, and we can extract the PLINK score file from it
write.table(gwas_variants, gwas_path, sep = "\t", row.names = FALSE, quote = FALSE)
message(paste("Generated GWAS catalog at", gwas_path))

# 3. Generate mock health data matching the sample IDs
# To test R data cleaning, we will inject some messy/dirty rows:
# - A few missing values (NA)
# - Some negative BMI values
# - Extreme values
# - Gender mismatch characters
message("Generating health raw data...")

ages <- sample(20:80, num_samples, replace = TRUE)
genders <- sample(c("M", "F"), num_samples, replace = TRUE)
exercise <- round(runif(num_samples, 0, 15), 1)

# Generate cholesterol that is correlated with age, exercise, and a genetic tendency
# We will simulate a genetic score from the VCF later, but here we generate base values
bmi <- round(rnorm(num_samples, mean = 24, sd = 4.5), 1)
sbp <- round(110 + 0.5 * ages + rnorm(num_samples, 0, 10), 0)
dbp <- round(70 + 0.3 * ages + rnorm(num_samples, 0, 8), 0)
glucose <- round(80 + 0.4 * ages + rnorm(num_samples, 0, 12), 0)

# Simulate cholesterol
# We add a correlation factor to age and BMI
cholesterol <- round(170 + 0.8 * ages + 1.2 * bmi - 1.5 * exercise + rnorm(num_samples, 0, 20), 0)

smoking <- sample(c("Never", "Former", "Current"), num_samples, replace = TRUE, prob = c(0.6, 0.25, 0.15))

health_df <- data.frame(
  user_id = samples,
  age = ages,
  gender = genders,
  bmi = bmi,
  sbp = sbp,
  dbp = dbp,
  cholesterol = cholesterol,
  glucose = glucose,
  smoking_status = smoking,
  exercise_hours = exercise,
  stringsAsFactors = FALSE
)

# Inject dirty/messy data to test cleaning:
# 1. 10 rows with negative or zero BMI
health_df$bmi[sample(1:num_samples, 10)] <- -9.0
# 2. 5 rows with extremely high/invalid age
health_df$age[sample(1:num_samples, 5)] <- 999
# 3. 10 rows with NAs in cholesterol
health_df$cholesterol[sample(1:num_samples, 10)] <- NA
# 4. Gender codes in lowercase or other abbreviations for 15 rows
health_df$gender[sample(1:num_samples, 15)] <- "male"
health_df$gender[sample(1:num_samples, 15)] <- "female"

# Save raw health data
write.csv(health_df, health_raw_path, row.names = FALSE, quote = FALSE)
message(paste("Generated mock raw health data at", health_raw_path))
message("Mock data generation completed successfully!")
