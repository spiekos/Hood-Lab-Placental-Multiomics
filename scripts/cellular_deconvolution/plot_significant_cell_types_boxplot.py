import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import sys


DICT_COLORS = {
	'Control': 'black',
	'FGR': 'cornflowerblue',
	'FGR+HDP': 'mediumorchid',
	'PE': 'pink',
	'PTD': 'gold'
}


DICT_RENAME = {
	'FGR': 'FGR',
	'FGR+hyp': 'FGR+HDP',
	'PTD': 'PTD',
	'Severe PE': 'PE'
}


def format_input(file_cibersortx, file_cohort):
	df_condition = pd.read_csv(file_cohort)[['Patient-ID', 'Condition']]
	# replace group names in df
	df_condition = rename_groups(df_condition)
	df_cibersortx = pd.read_csv(file_cibersortx).rename(columns={'Mixture': 'Patient-ID'})
	df_cibersortx = pd.merge(df_condition, df_cibersortx, on='Patient-ID', how='inner')
	return df_cibersortx


def rename_groups(df):
	# rename group names in df using universally declared dictionary
	for k, v in DICT_RENAME.items():
		df = df.replace(k, v)
	return df


def plot_cell_types_by_condition(df, file_cond_test, file_output_1, file_output_2):
	# for each condition make a box plot
	for condition in DICT_RENAME.values():
		# limit df to only those in the control or condition
		mask = df['Condition'].isin(['Control', condition])
		df_condition = df[mask]
		df_condition = df_condition
		df_condition = df_condition.drop(['P-value', 'Correlation', 'RMSE'], axis=1)
		# make boxplot of all cell types
		plot_data_box(df_condition, file_output_1, condition)
		# make boxplot of only significant cell types
		plot_signicant_cell_types(condition, df_condition, file_cond_test, file_output_2)
	return


def plot_data_box(df, file_output, condition):
    """
    Plots a box plot for the percentage of different cell types,
    comparing a specific condition against the overall distribution.

    Args:
        df (pd.DataFrame): DataFrame containing cell type percentages.
                           Expected columns: 'Patient-ID', 'Condition', and
                           columns for each cell type (e.g., 'Fetal Mesenchymal Stem Cells').
        file_output (str): Base filename for saving the plot.
        condition (str): The specific condition to highlight (e.g., 'FGR').
    """

    plt.figure(figsize=(12, 8))

    # Melt the DataFrame to long format for seaborn boxplot
    # We want to compare cell types, so 'cell_type' will be our x-axis.
    # The 'value' column will contain the percentages.
    # Exclude 'Patient-ID' and 'Condition' from the melt to get only cell type columns
    cell_type_columns = [col for col in df.columns if col not in ['Patient-ID', 'Condition']]
    df_melted = df.melt(id_vars=['Patient-ID', 'Condition'],
                        value_vars=cell_type_columns,
                        var_name='Cell Type',
                        value_name='Percentage of Cells')

    # Create a 'Group' column for easier distinction in the plot
    # Here, we're assuming the 'condition' you pass is the 'group' and
    # anything else is implicitly 'control' or just part of the overall distribution.
    # For a true "control vs. group" box plot, you'd ideally have a 'Control' condition
    # explicitly in your 'Condition' column.
    # For this example, we'll plot all data and highlight the specified condition.

    # Filter for the specified condition to potentially plot separately or highlight
    df_condition = df_melted[df_melted['Condition'] == condition]
    df_other = df_melted[df_melted['Condition'] != condition]

    # Use seaborn to create the box plot
    # We'll plot all data and differentiate by 'Condition' if there are multiple,
    # or just show the distribution if only one condition is present.
    sns.boxplot(data=df_melted, x='Percentage of Cells', y='Cell Type', hue='Condition', palette=DICT_COLORS)

    plt.title(f'Percentage of Cell Types by Condition: {condition}')
    plt.xlabel('Percentage of Cells in Placenta')
    plt.ylabel('Cell Type')
    plt.gca().invert_yaxis()
    plt.legend(title='Condition')
    plt.tight_layout() # Adjust layout to prevent labels from overlapping

    # Save plot
    plt.savefig(file_output + '_' + condition + '_boxplot.pdf')

    # Show plot
    plt.show()

    return


def plot_signicant_cell_types(condition, df_condition, file_cond_test, file_output_2):
	# read cond test csv in as pandas df
	df_sig_cond = rename_groups(pd.read_csv(file_cond_test))
	# limit to conditions and cell types that are significant
	df_sig_cond = df_sig_cond[df_sig_cond['pvalue_adj'] < 0.05]

	# limit condition-specific df to celltypes that are significantly different
	df_sig_condition_limited = df_sig_cond[df_sig_cond['group'] == condition]
	# if any cell-types are significant, plot a boxplot of just these
	if df_sig_condition_limited.shape[0] > 0:
		list_cell_types = list(df_sig_condition_limited['cell_type'])
		list_cell_types.append('Patient-ID')
		list_cell_types.append('Condition')
		df_condition_limited = df_condition[list_cell_types]
		plot_data_box(df_condition_limited, file_output_2, condition)
	return


def main():
	file_cibersortx = sys.argv[1]
	file_cond_test = sys.argv[2]
	file_cohort = sys.argv[3]
	file_output_1 = sys.argv[4]
	file_output_2 = sys.argv[5]
	df = format_input(file_cibersortx, file_cohort)
	plot_cell_types_by_condition(df, file_cond_test, file_output_1, file_output_2)


if __name__ == '__main__':
	main()