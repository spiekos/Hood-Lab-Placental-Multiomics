""" Generates a horizontal bar graph """
import sys
import numpy as np
import matplotlib.pyplot as plt

outputFile = sys.argv[1]

fig = plt.figure(figsize=(10.5, 6))

N = 5
index = np.arange(N)    # the x locations for the groups
width = 0.6    # the width of the bars: can also be len(x) sequence
p_values = (
	0.01424,
	0.01424,
	0.01424,
	0.01424,
	0.006896
	)
neg_log_p_values = list(-np.log10(p_values))
labels =  (
	'Transmembrane Receptor Protein\nTyrosine Kinase Signaling Pathway',
	'Positive Regulation Of Protein\nSerine/Threonine Kinase Activity',
	'Positive Regulation Of\nMAP Kinase Activity',
	'Positive Regulation Of\nCell Migration',
	'Phagocytosis'
	)

p1 = plt.barh(index, neg_log_p_values, width, align = 'center', color = 'lightpink', tick_label = labels)

plt.xlabel('-log10 p-value', fontsize = 20)
plt.title('Female Fetuses Only\nGene Ontology Terms', fontsize = 28)
plt.yticks(index)
plt.xlim(0, 3)
plt.yticks(fontsize = 20)
plt.xticks(fontsize = 20)

plt.tight_layout() # Adjust layout to prevent clipping of y-axis labels
plt.savefig(outputFile)
plt.close()
