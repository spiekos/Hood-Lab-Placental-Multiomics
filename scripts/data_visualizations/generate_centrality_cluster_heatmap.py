import sys
import pandas as pd
import numpy as np
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt

# declare universal variables
GROUPS = [
			'control',
			'fgr',
		 	'fgr+hyp',
			'severePE',
			'ptd'
		]

DICT_GROUP_NAMES = {
			'control': 'Control',
			'fgr': 'FGR',
		 	'fgr+hyp': 'FGR+HDP',
			'severePE': 'PE',
			'ptd': 'PTD'
			}


def plot_heatmap(data, list_analytes, file_output):
	plt.rcParams['pdf.fonttype'] = 42
	cg = sns.clustermap(data, center=0, linewidths=.5, figsize=(14,12), cmap="viridis", xticklabels=True, yticklabels=True)
	plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=25, fontsize=16)
	plt.setp(cg.ax_heatmap.get_ymajorticklabels(), fontsize=8)
	cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize=8)
	# set col names
	col_names = []
	for i, ticklabel in enumerate(cg.ax_heatmap.xaxis.get_majorticklabels()):
		col_names.append(list(DICT_GROUP_NAMES.values())[int(ticklabel.get_text())])
	cg.ax_heatmap.set_xticklabels(col_names)
	# set row names
	row_names = []
	for i, ticklabel in enumerate(cg.ax_heatmap.yaxis.get_majorticklabels()):
		row_names.append(list_analytes[int(ticklabel.get_text())])
	cg.ax_heatmap.set_yticklabels(row_names)

	plt.savefig(file_output, bbox_inches='tight')
	plt.show()


def create_cluster_analyte_list(file_cluster):
	# opening the file in read mode 
	my_file = open(file_cluster, "r")  
	# reading the file 
	data = my_file.read() 
	# replacing end splitting the text when newline ('\n') is seen
	data_into_list = data.split("\n") 
	my_file.close()
	return data_into_list


def read_centrality_file_to_df(group):
	# read centrality scores for a given outcome to a csv
	filepath = 'output/' + group + '/' + group + '_centrality_communities.csv'
	df = pd.read_csv(filepath, index_col=0, header=0)
	# rename column to outcome
	df.rename(columns={df.columns[0]: DICT_GROUP_NAMES[group]}, inplace=True)
	return df


def identify_centrality_values(list_cluster):
	# create a df for centrality values for each analyte with each column as an outcome
	df = read_centrality_file_to_df(GROUPS[0])
	for g in GROUPS[1:]:
		df_temp = read_centrality_file_to_df(g)
		group_name = DICT_GROUP_NAMES[g]
		df[group_name] = df_temp[group_name]
		print(df.head())
	# limit to those analytes that are in the cluster of interest
	df = df.filter(items=list_cluster, axis=0)
	# fill na with 0
	df = df.fillna(0)
	return df


def generate_centrality_heatmap(file_cluster, file_output):
	# identify the subset of analytes that are of interest
	analytes = create_cluster_analyte_list(file_cluster)
	# get the centrality values for each group for those analytes
	df_centrality = identify_centrality_values(analytes)
	centrality_values = df_centrality.to_numpy()
	plot_heatmap(centrality_values, analytes, file_output)


def main():
	file_cluster = sys.argv[1]
	file_output = sys.argv[2]
	generate_centrality_heatmap(file_cluster, file_output)


if __name__ == "__main__":
    main()
