# set up environment
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal, sem
import sys


# declare universal variables
DICT_RENAME = {
	'FGR+hyp': 'FGR+HDP',
	'Severe PE': 'PE'
}
DROP_GROUPS = ['Control PTD']
ID = 'Patient-ID'
NORMALIZE_PHENOTYPE = 'Control'
PHENOTYPE = 'Condition'


def check_anlyte_in_df(l, analyte):
	# check to make sure analyte of interest is in input file
	l = list(l)
	if analyte not in l:
		print('Error! ', analyte, ' not in input file!')
		sys.exit(1)
	return


def replace_values(df, col, dict_replacement):
	for key, value in dict_replacement.items():
		df.loc[df[col] == key, col] = value
	return df


def create_df(file_input, file_input_metadata, analyte):
	# create pandas df of condition and anlayte of interest
	# from two input csv files
	df = pd.read_csv(file_input)
	df_metadata = pd.read_csv(file_input_metadata)
	# check that analyte is a column in the input file
	check_anlyte_in_df(df.columns, analyte)
	# Merge DataFrames on patient-id
	df_merged = pd.merge(df, df_metadata, on=ID)
	# filter for condition and analtye of interest columns only
	cols_of_interest = [PHENOTYPE, analyte]
	df_final = df_merged[cols_of_interest]
	# rename phenotypes
	df_final = replace_values(df_final, PHENOTYPE, DICT_RENAME)
	return df_final


def normalize_count(df, analyte):
	# normalize counts to the median of the control
	# filter for control rows
	control_values = df[df[PHENOTYPE] == NORMALIZE_PHENOTYPE]
	# calculate median of control values for the analyte
	median = control_values[analyte].median()
	# normalize every value for the analyte against the median
	df[analyte] = df[analyte] / median
	# impute missing values with half min value
	min_imputation_value = df[analyte].min() / 2
	df[analyte] = df[analyte].fillna(min_imputation_value)
	return df


def unique_values(df, col):
	# return list of unique values in specified column of df
	l = list(set(df[col]))
	l.sort()
	return l


def create_analyte_dict(df, analyte):
	# create dictionary with conditions and keys
	# and the a list of the measurements of the 
	# analyte of interest as the value
	d = {}
	# get a list of the unique phenotypes
	groups = unique_values(df, PHENOTYPE)
	# remove undesired groups from the list
	groups = [item for item in groups if item not in DROP_GROUPS]
	# for each phenotype get the list of measurements
	# for patients belonging to that group and
	# write them to dict
	for group in groups:
		df_group = df[df[PHENOTYPE] == group]
		measurements = list(df_group[analyte])
		d[group] = measurements
	return d


def perform_kruskal_wallis_test(d):
	# Extract measurements from the dictionary
	groups = [d[key] for key in d]

	# Perform Kruskal-Wallis test
	statistic, p_value = kruskal(*groups)

	# Print the statistic for p-value
	print(f'Statistic for the Kruskal-Wallis test: {statistic}')
	print(f'P-value for the Kruskal-Wallis test: {p_value}')
	return


def plot_values(d, analyte):
	# Calculate mean and standard error of the mean (SEM) for each group
	means = [np.mean(data) for data in d.values()]
	sem_values = [sem(data) for data in d.values()]

	# calcualte p-value and output to console
	perform_kruskal_wallis_test(d)

	# Plotting
	fig, ax = plt.subplots()
	bar_width = 0.35
	index = np.arange(len(d))

	# Plot the mean as a dot
	dot = ax.scatter(index, means, color='black', label='Mean')

	# Plot the 95% confidence interval as whiskers
	for i, (mean, sem_val) in enumerate(zip(means, sem_values)):
		ax.plot([i - bar_width / 4, i + bar_width / 4], [mean + 1.96 * sem_val] * 2,\
		 color='black', linewidth=2, linestyle='solid')
		ax.plot([i - bar_width / 4, i + bar_width / 4], [mean - 1.96 * sem_val] * 2,\
		 color='black', linewidth=2, linestyle='solid')
		ax.vlines(i, mean - 1.96 * sem_val, mean + 1.96 * sem_val, color='black', linewidth=2)

	# format labels
	ax.set_ylabel(analyte + '\nRelative Expression Level')
	ax.set_title('95% Confidence Interval\n' + analyte)
	ax.set_xticks(index)
	ax.set_xticklabels(d.keys())

	# Set x-axis to start at 0
	ax.set_ylim(bottom=0)

	# Save the plot to a PDF file
	filepath = 'output/' + analyte + '_95%CI_plot.pdf'
	plt.savefig(filepath)


def main():
	file_input = sys.argv[1]  # read in arguments
	file_input_metadata = sys.argv[2]
	analyte = str(sys.argv[3])
	# create pandas df of phenotypes and analyte measurements
	df = create_df(file_input, file_input_metadata, analyte)
	# normalize count of analyte
	df = normalize_count(df, analyte)
	# make a dictionary of analyte values by phenotype
	d = create_analyte_dict(df, analyte)
	# make plot of values
	plot_values(d, analyte)


if __name__ == '__main__':
	main()
