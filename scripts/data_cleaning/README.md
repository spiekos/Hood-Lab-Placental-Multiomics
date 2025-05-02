# Data Cleaning Script Calls
Goes from raw files provided by data processing companies like Metabolon and Olink and puts them in a standardized format ready for downstream processing and analyses.

1. Clean patient clinical data for interomics analyses
\n`python3 format_cohort_metadata.py ../data/Final_Cohort.csv ../data/metadata-final.csv`

2. Clean miRNA data for interomics analyses
`python3 clean_mirna_data.py ../data/TranscriptSampleInformation.csv ../data/mini-DP3_miRNA_rds.csv ../data/mirna-final.csv`

3. Clean and normalize protein data for general downstream analysis
`python3 clean_proteomics_data.py ../data/Olink_full_OMICS_Undetectable_Analysis_Missing_Values_NA.csv ../data/Olink_full_OMICS_min_imputed.csv > proteins_missingness_report.log`

4. Clean and normalize transcriptomics data for general downstream analysis
`python3 clean_and_normalize_transcriptomics_data.py ../data/TranscriptSampleInformation.csv ../data/miniDP3_miRNA_rds.csv ../data/transcripts-min-imputed.csv > transcripts_missingness_report.log`

5. Format placental histopathology features into two csv files
`python3 format_placental_histopathology_report.py ../data/Path\ report\ tabulation.xlsx   ../data/additional_histopathology_features_from_slides.csv  ../data/PlacentalHistopathologyReportsCleaned.csv ../data/PlacentalHistopathologySlideFeaturesCleaned.csv`

6. Formatted Olink_full_OMICS_Undetectable_Analysis_Missing_Values_NA.csv manually as proteins-final.csv for input into interomics analyses. This included dropping the conditions column and changing the column name from 'Sample' to 'Patient-ID'. Formatted For/ Sam.csv manually as metabolites-final.csv for input into interomics analyses. This includes dropping unneaded descriptor columns and rows and labeling the samples column as "Patient-ID"

7. Clean RNA-seq data for interomics analyses
`python3 clean_rna-seq_data.py ../data/TranscriptSampleInformation.csv ../data/mini-DP3_RNA-seq.csv ../data/mini-DP3_RNA-seq_gProfiler.csv ../data/transcripts-final.csv ../data/transcripts-meta-final.csv`

8. Clean metabolomics data for downstream analysis
`python3 clean_metabolomics_data_3.py ../data/metabolomics_batch_normalized_data.csv ../data/metabolomics_patient_id_key.csv ../data/metabolomics_chemical_annotation.csv ../data/metabolites-final.csv`

9. Filter transcriptomics data for thresholds
`python3.9 process_transcriptomics_thresholds.py ../data/transcripts-final.csv`
