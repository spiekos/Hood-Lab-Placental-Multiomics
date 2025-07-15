# Data Cleaning Script Calls
Goes from raw files provided by data processing companies like Metabolon and Olink and puts them in a standardized format ready for downstream processing and analyses.

1. Clean patient clinical data for interomics analyses\
`python3 format_cohort_metadata.py ../../data/raw/Final_Cohort.csv ../../data/normalized_cleaned/metadata-final.csv`

2. Clean miRNA data for interomics analyses\
`python3 clean_mirna_data.py ../../data/raw/TranscriptSampleInformation.csv ../../data/raw/mini-DP3_miRNA_rds.csv ../../data/normalized_cleaned/mirna-final.csv`

3. Format placental histopathology features into two csv files\
`python3 format_placental_histopathology_report.py../../data/raw/Path\ report\ tabulation.xlsx   ../../data/raw/additional_histopathology_features_from_slides.csv  ../../data/normalized_cleaned/PlacentalHistopathologyReportsCleaned.csv ../../data/normalized_cleaned/PlacentalHistopathologySlideFeaturesCleaned.csv`

4. Clean metabolomics data for downstream analysis\
`python3 clean_metabolomics_data.py ../../data/raw/metabolomics_batch_normalized_data.csv ../../data/raw/metabolomics_patient_id_key.csv ../../data/raw/metabolomics_chemical_annotation.csv ../../data/normalized_cleaned/metabolites-final.csv`

5. Convert all red cell values (olink determined below the level of detection) and red text values (sample failed quality control for that olink panel) using VBA in excel for the `../../data/raw/Olink_full_OMICS.csv` proteomics data received from Olink. Save this as ../../data/raw/Olink_full_OMICS_Undetectable_Analysis_Missing_Values_NA.csv
   
6. For all proteins with more than 20% of the samples missing, evaluate the distribution of which samples are missing between the cohorts. Check if any protein is disproportionately missing in an obstetric condition following Bonferonni correction. Record analysis in a log file. Clean and format proteomics data for downstream analysis.\
`python3 clean_proteomics_data.py ../../data/raw/Olink_full_OMICS_Undetectable_Analysis_Missing_Values_NA.csv ../../data/normalized_cleaned/proteins-final.csv > ../../data/normalized_cleaned/logs/proteins_missingness_report.log`

7. For all transcripts with more than 20% of the samples missing, evaluate the distribution of which samples are missing between the cohorts. Check if any transcript is disproportionately missing in an obstetric condition following Bonferonni correction. Record analysis in a log file.\
`python3 clean_and_normalize_transcriptomics_data.py ../../data/raw/TranscriptSampleInformation.csv ../../data/raw/miniDP3_miRNA_rds.csv ../../data/normalized_cleaned/transcripts-min-imputed.csv > ../../data/normalized_cleaned/logs/transcripts_missingness_report.log`

8. Clean RNA-seq data for interomics analyses\
`python3 clean_rna-seq_data.py ../../data/raw/TranscriptSampleInformation.csv ../../data/raw/mini-DP3_RNA-seq.csv ../../data/raw/mini-DP3_RNA-seq_gProfiler.csv ../../data/normalized_cleaned/transcripts-processed.csv ../../data/normalized_cleaned/transcripts-meta-final.csv`

9. Filter transcriptomics data for thresholds\
`python3 process_transcriptomics_thresholds.py ../../data/normalized_cleaned/transcripts-processed.csv ../../data/raw/Final_Cohort.csv ../../data/normalized_cleaned/transcripts-final.csv`
