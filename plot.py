import pandas as pd
import matplotlib.pyplot as plt


def show_plot(file_path, geneA, geneB):
    df = pd.read_csv(file_path, sep="\t", index_col=0)
    dt = df.T
    dt.plot(y=[geneA, geneB], style=".")
    plt.show()


def corr_genes(file_path, geneA, geneB):
    df = pd.read_csv(file_path, sep="\t", index_col=0)
    sd = df.T[[geneA, geneB]]   
    return sd.corr().values[0][1]
