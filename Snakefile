# Snakefile
# Snakemake pipeline for Full Spectrum Bio-Analysis data processing.

rule all:
    input:
        "database/bio_analysis.db",
        "Dataset/model_summary.txt"

rule generate_mock_data:
    input:
        vcf = "Dataset/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    output:
        health_raw = "Dataset/health_raw.csv",
        gwas = "Dataset/gwas_catalog_associations.tsv"
    shell:
        "Rscript scripts/generate_mock_data.R"

rule annotate_vcf:
    input:
        vcf = "Dataset/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        gwas = "Dataset/gwas_catalog_associations.tsv"
    output:
        vcf_ann = "Dataset/Annotated.vcf.gz"
    shell:
        "bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -O z -o {output.vcf_ann} {input.vcf}"

rule calculate_prs:
    input:
        vcf_ann = "Dataset/Annotated.vcf.gz",
        gwas = "Dataset/gwas_catalog_associations.tsv",
        plink = "bin/plink"
    output:
        profile = "Dataset/prs_out.profile"
    shell:
        """
        # Extract columns for PLINK score file: variant_id, risk_allele, beta
        Rscript -e "df <- read.table('{input.gwas}', header=TRUE, sep='\\t'); write.table(df[, c('panel_variant_id', 'risk_allele', 'beta')], 'Dataset/gwas_weights.txt', col.names=FALSE, row.names=FALSE, quote=FALSE, sep='\\t')"
        
        # Run PLINK scoring
        ./bin/plink --vcf {input.vcf_ann} --score Dataset/gwas_weights.txt 1 2 3 --out Dataset/prs_out
        """

rule clean_health_data:
    input:
        raw = "Dataset/health_raw.csv"
    output:
        clean = "Dataset/health_cleaned.csv"
    shell:
        "Rscript scripts/clean_health_data.R"

rule run_modeling:
    input:
        clean = "Dataset/health_cleaned.csv",
        profile = "Dataset/prs_out.profile"
    output:
        rdata = "Dataset/model_results.RData",
        summary = "Dataset/model_summary.txt"
    shell:
        "Rscript scripts/run_modeling.R"

rule load_to_db:
    input:
        clean = "Dataset/health_cleaned.csv",
        rdata = "Dataset/model_results.RData"
    output:
        db = "database/bio_analysis.db"
    shell:
        "Rscript scripts/load_to_db.R"
