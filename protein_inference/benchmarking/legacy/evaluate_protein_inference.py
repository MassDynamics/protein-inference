import matplotlib.pyplot as plt
from numpy import linspace


class EvaluateProteinInference():
    def __init__(self):

        return

    def plot_n_over_fdr(self, protein_table):

        fdr_values = linspace(0, 1, 100)
        n_identified = []
        for value in fdr_values:
            n_identified.append(sum(protein_table.FDR <= value))

        plt.plot(fdr_values, n_identified)
        plt.xlabel("FDR")
        plt.ylabel("Number of Proteins Inferred")
        plt.title("Inference per Discovery Rate")

        return plt.show()
