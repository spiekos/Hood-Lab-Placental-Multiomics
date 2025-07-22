`deconvolution_functions.py`: 

A small module of functions necessary for carrying out the statistical analysis of the deconvolution results. These include CLR functions for transforming the data, functions for running the statistical tests, and  for plotting and formatting figures. 

`deconvolution_analysis.py`: 

A python script whereby the statistical functions and plotting are run. This runs the functions in the necessary order: transform the data using the CLR functions, run the statistical analysis, and plot the results. 

`plot_significant_cell_types_boxplot.py`

Makes boxplots of cell type percentages of control vs each obstetric syndrome placentas. Each box displays the IQR from 25th percentile (Q1) to 75th percentile (Q3). The line in each box is the median. Whiskers extend to 1.5 times the IQR and the circles outside the whiskers are outliers. To run:

`python3 plot_significant_cell_types_boxplot.py ../../data/cellular_deconvolution/cibersortx_results_placenta.csv ../../data/cellular_deconvolution/cond_test_clr_imp_1_v_control.csv ../../data/normalized_cleaned/metadata-final.csv all_cell_types_control_vs significant_cell_types_control_vs`
