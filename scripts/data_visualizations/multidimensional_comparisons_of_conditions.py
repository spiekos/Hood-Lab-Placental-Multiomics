# import environment
import pandas as pd
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
import numpy as np
import sys

from matplotlib.patches import Ellipse
from scipy.spatial.distance import pdist, squareform
from sklearn.cross_decomposition import PLSRegression
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, euclidean_distances
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import QuantileTransformer, StandardScaler
from scipy.stats import chi2


# Declare Universal Variable
THRESH_MISSINGNESS = 0.8
THRESH_UNIQUENESS = 0.5
THRESH_MIRNA_READ_COUNTS = 10
THRESH_RNA_READ_COUNTS = 500


DICT_COLORS = {
	'Control': 'black',
	'FGR': 'cornflowerblue',
	'FGR+HDP': 'mediumorchid',
	'PE': 'pink',
	'PTD': 'gold'
}


DICT_CONDITIONS = {
	'Control': 'Control',
	'Control PTD': 'ControlPTD',
	'FGR': 'FGR',
	'FGR+hyp': 'FGR+HDP',
	'PTD': 'PTD',
	'Severe PE': 'PE'
}


def filter_data(df, data_type):
	df_filtered = df.shape[0] * THRESH_MISSINGNESS
	df.dropna(thresh=df_filtered, axis=1, inplace=True)
	minUnique = df.nunique()/len(df)
	col_drop = minUnique[minUnique < THRESH_UNIQUENESS].index
	df.drop(columns=col_drop, inplace=True)
	if data_type.lower() == 'mirna':
		col_sums = df.sum()
		col_sums.drop(col_sums.index[1], inplace=True)
		col_drop = col_sums[col_sums < THRESH_MIRNA_READ_COUNTS].index
		df.drop(columns=col_drop, inplace=True)
	if data_type.lower() == 'transcript':
		col_sums = df.sum()
		col_sums.drop(col_sums.index[1], inplace=True)
		col_drop = col_sums[col_sums < THRESH_RNA_READ_COUNTS].index
		df.drop(columns=col_drop, inplace=True)
	return df


def impute_missing_data(df):
	float_columns = [col for col in df.columns if col != 'Patient-ID']
	df[float_columns] = df[float_columns].astype(float)
	df.fillna(0, inplace=True)

	# Calculate the minimum values for each numeric column
	df_filtered= df[df[float_columns] != 0]
	min_values = df_filtered[float_columns].min()/2

	# Replace 0 values in the numeric columns with their respective minimum values
	df[float_columns] = df[float_columns].apply(lambda col: np.where(col == 0, min_values[col.name], col))

	return df


def divide_by_column_average(df, exclude_column='Patient-ID'):
    # Create a copy of the dataframe to avoid modifying the original one
    divided_df = df.copy()
    
    # Calculate the average value for each column, excluding the specified column
    column_means = df.loc[:, df.columns != exclude_column].mean()
    
    # Divide each cell by the average value of its column for all columns except the excluded one
    divided_df.loc[:, df.columns != exclude_column] = df.loc[:, df.columns != exclude_column].divide(column_means)
    
    return divided_df


def initiate_dict(l):
	d = {}
	for item in l:
		d[item] = []
	return d


def add_item_to_dict(d, k, v):
	if len(v) > 0:
		d[k].append(v)
	return d


def create_limited_analyte_dict(file_group):
	df = pd.read_table(file_group, sep=',').fillna('')
	d = initiate_dict(list(df.columns))
	for index, row in df.iterrows():
		d = add_item_to_dict(d, 'metabolites', row['metabolites'])
		d = add_item_to_dict(d, 'mirna', row['mirna'])
		d = add_item_to_dict(d, 'proteins', row['proteins'])
		d = add_item_to_dict(d, 'transcripts', row['transcripts'])
	return d


def read_analyte_df(file_input, dict_group, data_type):
	df = pd.read_table(file_input, sep=',', dtype={'Patient-ID':str})
	if len(dict_group[data_type]) > 0:
		df_temp = pd.DataFrame(df['Patient-ID'])
		df = df[dict_group[data_type]]
		df = filter_data(df, data_type)
		df = impute_missing_data(df)
		#df = df_temp.concat(df, axis=1)
		df = pd.concat([df_temp, df], axis=1)
	else:
		df = pd.DataFrame(df['Patient-ID'])
	return df


def format_limited_df(file_metabolites, file_mirna, file_proteins, file_transcripts, file_group):
	dict_group = create_limited_analyte_dict(file_group)
	df_metabolites = read_analyte_df(file_metabolites, dict_group, 'metabolites')
	df_mirna = read_analyte_df(file_mirna, dict_group, 'mirna')
	df_proteins = read_analyte_df(file_proteins, dict_group, 'proteins')
	df_transcripts = read_analyte_df(file_transcripts, dict_group, 'transcripts')
	df = pd.merge(df_metabolites, df_mirna, how='inner', on=['Patient-ID'])
	df = pd.merge(df, df_proteins, how='inner', on=['Patient-ID'])
	df = pd.merge(df, df_transcripts, how='inner', on=['Patient-ID'])
	return df


def rename_conditions(df):
	for k, v in DICT_CONDITIONS.items():
		df = df.replace(k, v)
	return df


def limit_conditions_to_two(df, group1, group2):
	df = rename_conditions(df)
	for i in DICT_CONDITIONS.values():
		if i != group1 and i != group2:
			df = df[df['Condition']!=i]
	return df


def format_data(file_metabolites, file_mirna, file_proteins, file_transcripts, file_group, file_clinical, group1, group2):
	df = format_limited_df(file_metabolites, file_mirna, file_proteins, file_transcripts, file_group)
	df = divide_by_column_average(df)
	df_clinical = pd.read_table(file_clinical, sep=',', dtype={'Patient-ID':str})
	df_clinical = limit_conditions_to_two(df_clinical, group1, group2)
	df_clinical = df_clinical[df_clinical['InfSex']>=0]
	df = df.merge(df_clinical[['Patient-ID', 'Condition']], on=['Patient-ID'], how='inner')
	d = {group1: 0, group2: 1}
	df['Condition'] = df['Condition'].map(d)
	return df


def plot_pca(df, data_type, group1, group2):
	# Normalize data
	n=df.shape[0]
	scaler = QuantileTransformer(n_quantiles=n)
	scaled_data = scaler.fit_transform(df.drop(['Patient-ID', 'Condition'], axis=1))

	# Perform PCA
	pca = PCA(n_components=2)
	principal_components = pca.fit_transform(scaled_data)

	# Access the principal components
	pc1 = principal_components[:, 0]
	pc2 = principal_components[:, 1]

	## Create a custom colormap for male (blue) and female (pink)
	condition = df['Condition']

	# Create a scatter plot with different colors for male and female
	fig, ax = plt.subplots(figsize=(8, 6))

	# Define your custom colormap (cmap)
	color1 = DICT_COLORS[group1]
	color2 = DICT_COLORS[group2]
	custom_colors = [color1,  color2]   # Customize with your desired colors
	custom_cmap = mcolors.ListedColormap(custom_colors)

	# Create the scatter plot with colors based on the 'InfSex' column
	scatter = plt.scatter(pc1, pc2, c=condition, cmap=custom_cmap)
	legend_handles = [mpatches.Patch(color=color, label=f'Class {i+1}') for i, color in enumerate(custom_colors)]
	ax.legend(handles=legend_handles, labels=[group1, group2], title='Conditions', loc='best')

	# add labels
	explained_variance_ratios = pca.explained_variance_ratio_
	pc1 = str(round(explained_variance_ratios[0]*100, 1))
	pc2 = str(round(explained_variance_ratios[1]*100, 1))
	ax.set_xlabel('PC 1 (' + pc1 + '%)')
	ax.set_ylabel('PC 2 (' + pc2 + '%)')
	title = data_type + ' PCA Plot'
	plt.title(title)

	# save and plot figure
	file_output = 'condition_limited_' + group1 + '_vs_' + group2 + '_' + data_type + '_pca.pdf'
	plt.savefig(file_output, bbox_inches='tight')
	plt.show()


def create_similarity_matrix(df):
	# Normalize data
	n=df.shape[0]
	scaler = QuantileTransformer(n_quantiles=n)
	df_scaled = scaler.fit_transform(df.drop(['Patient-ID', 'Condition'], axis=1))

	# Calculate Euclidean distance matrix
	distance_matrix = squareform(pdist(df_scaled))

	return distance_matrix


def calculate_stress_reduction(distance_matrix):
	# Compute the dissimilarity matrix for the original data
	dissimilarity_original = distance_matrix

	# Perform MDS on the original dissimilarity matrix
	mds = MDS(n_components=2, dissimilarity='precomputed', random_state=611, normalized_stress=False)
	embedding = mds.fit_transform(dissimilarity_original)

	# Compute the stress of the original data
	stress_initial = mds.stress_

	# Compute the pairwise Euclidean distances between points in the lower-dimensional space
	distance_embedding = euclidean_distances(embedding)

	# Compute the stress of the embedded data
	stress_final = np.sqrt(np.sum((distance_embedding - dissimilarity_original) ** 2))

	# Calculate the variance explained for the first two dimensions
	variance_explained = ((stress_initial - stress_final) / stress_initial) * 100

	return str(round(variance_explained, 2))


def plot_plsda(df, data_type, group1, group2):
    # Split the data into features (X) and labels (y)
    X = df.drop(['Condition', 'Patient-ID'], axis=1)
    y = df['Condition']

    # Split the data
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    # PLS-DA model
    n_components = 2
    plsda = PLSRegression(n_components=n_components)
    plsda.fit(X_train, y_train)

    # Project data
    X_test_plsda = plsda.transform(X_test)
    y_pred = plsda.predict(X_test)
    y_pred_labels = [1 if val >= 0.5 else 0 for val in y_pred]
    accuracy = str(round(accuracy_score(y_test, y_pred_labels)*100, 1))

    # Custom colors
    color1 = DICT_COLORS[group1]
    color2 = DICT_COLORS[group2]
    custom_colors = [color1, color2]
    custom_cmap = mcolors.ListedColormap(custom_colors)

    # Create scatter plot
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(X_test_plsda[:, 0], X_test_plsda[:, 1], c=y_pred_labels, cmap=custom_cmap)
    plt.legend(handles=[scatter.legend_elements()[0][0], scatter.legend_elements()[0][1]], labels=[group1, group2])

    # Confidence ellipse function
    def plot_confidence_ellipse(data, ax, n_std=2.0, facecolor='none', edgecolor='black', label=None):
        if data.shape[0] < 2:
            return
        cov = np.cov(data, rowvar=False)
        mean = np.mean(data, axis=0)
        lambda_, v = np.linalg.eig(cov)
        lambda_ = np.sqrt(lambda_)
        chi2_val = np.sqrt(chi2.ppf(0.95, df=2))  # 95% CI
        width, height = 2 * chi2_val * lambda_
        angle = np.degrees(np.arctan2(*v[:, 0][::-1]))
        ellipse = Ellipse(xy=mean, width=width, height=height, angle=angle, facecolor=facecolor,
                          edgecolor=edgecolor, label=label, linewidth=2, linestyle='--')
        ax.add_patch(ellipse)

    # Add confidence ellipses
    ax = plt.gca()
    for label, color in zip([0, 1], [color1, color2]):
        group_data = X_test_plsda[np.array(y_pred_labels) == label]
        plot_confidence_ellipse(group_data, ax, edgecolor=color)

    # Labels
    plt.xlabel('PLS-DA Component 1')
    plt.ylabel('PLS-DA Component 2')
    plt.title(f'{data_type} PLS-DA (Accuracy: {accuracy}%)')

    # Save and plot
    file_output = f'condition_limited_{group1}_vs_{group2}_{data_type}_pls-da.pdf'
    plt.savefig(file_output, bbox_inches='tight')
    plt.show()


def main():
	file_metabolites = sys.argv[1]
	file_mirna = sys.argv[2]
	file_proteins = sys.argv[3]
	file_transcripts = sys.argv[4]
	file_group = sys.argv[5]
	file_clinical = sys.argv[6]
	subset = sys.argv[7]
	group1 = sys.argv[8]
	group2 = sys.argv[9]
	df = format_data(file_metabolites, file_mirna, file_proteins, file_transcripts, file_group, file_clinical, group1, group2)
	plot_pca(df, subset, group1, group2)
	plot_plsda(df, subset, group1, group2)


if __name__ == '__main__':
	main()
