`AdjustAnalytesForCovariates.ipynb`:
1. Takes in clinical, metabolite, protein, miRNA, and transcript (mRNA) data that is cleaned and deposited in [`data/normalized_cleaned`](../../data/normalized_cleaned).
2. Filters for placentas for which the fetal sex is recorded.
3. Filters for analytes that meet the following criteria: missingness < 20% of samples, unique value in > 50% of samples, miRNA minimum read count of 10, and transcript mRNA minimum read count of 500.
4. For each analyte for each fetal sex fit a generalized linear model: `analyte ~ weeks gestation at delivery + pregravid BMI2 + C(Labor Initiation) + C(Smoking Use) + C(Illicit Drug Use)`.
5. Perform benjamini-hochberg multiple hypothesis correction
6. For analytes that are significantly regulated by gestational week, adjust the analytes value to 40 weeks gestation for those analytes for those fetal sex. Adjusted analyte values are saved to [`data/normalized_cleaned_adjusted'](../../data/normalized_cleaned_adjusted)
7. Compare analytes significantly regulated by gestational weeks in female fetuses vs male fetuses.

`AdjustAnalytesForCovariates-Randomized-FDR.ipynb`:
1. Perform random simulation in which analytes are randomly assigned to group 1 vs 2 and then fit each analyte in each group with the generalized linear model: `analyte ~ weeks gestation at delivery + pregravid BMI2 + C(Labor Initiation) + C(Smoking Use) + C(Illicit Drug Use)`.
2. Adjust p-value for weeks of gestation variable with Benjamini-Hochberg multiple hypothesis correction. Record the number of significant analytes detected in each group.
3. Repeat 1000 times to build random distribution.
4. Evaluate number of anlaytes that are significantly regulated in male vs female and the overlap in significant analytes between the two compared to what we would expect by random chance.
