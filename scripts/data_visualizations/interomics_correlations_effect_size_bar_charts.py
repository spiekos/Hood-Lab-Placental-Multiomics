from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns


GROUP_ORDER = ['CONTROL', 'PTD', 'FGR', 'FGR_HYPERTENSION', 'SEVERE_PE']
COLORS={
    'CONTROL': 'mediumblue', 
    'FGR': 'firebrick', 
    'SEVERE_PE': 'green', 
    'FGR_HYPERTENSION': 'darkviolet', 
    'PTD': 'orange'}
VALUES = {
    'CONTROL': [0, 0, 0, 0, 16532, 120527,11124],
    'FGR': [0, 0, 0, 0, 0, 26686, 10818],
    'SEVERE_PE': [0, 0, 0, 0, 108471, 24474, 14837],
    'FGR_HYPERTENSION': [0, 0, 0, 0, 0, 21415, 59197],
    'PTD': [0, 0, 0, 0, 25868, 25403, 7101]
}


def plot_bar_chart():
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams.update({'font.size': 18})
    N = 7

    ind = np.arange(N) 
    width = 0.125
    #width=0.15

    c = 1
    for group in GROUP_ORDER:
        plt.bar(ind + width*c, VALUES[group], width, color=COLORS[group])
        c+=1       
    plt.axhline(0, color='black')
    plt.axvline(x=4, color='black', linestyle='dashed')
    plt.ylabel('Count')

    plt.xticks(ind+width*3.5, ('Very Small', 'Small', 'Medium', 'Medium-Large', 'Large', 'Very Large', 'Huge'), rotation=90)
    plt.xlabel('Effect Size')
    plt.ylim((0, 130000))
    plt.savefig('../output/interomics_correlations_effect_size_bar_chart.pdf', bbox_inches='tight')
    plt.show()


def main():
    plot_bar_chart()


if __name__ == "__main__":
    main()
