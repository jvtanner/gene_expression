
import matplotlib.pyplot as plt
import knn
import gsea
import functools as ft
import itertools
from collections import Counter

def q1():

    KNN = knn.KNN()
    KNN.load_data("GSE25628_filtered_expression.txt", "GSE25628_samples.txt")

    k = 3
    FN = [.05, .1, .25, .5, .75, .9, 1]

    vals = []
    xs = []
    ys = []
    for fn in FN:
        (s, sp) = KNN.calc_metrics(k, fn)
        xs.append(1-sp)
        ys.append(s)

    # print(xs, ys)

    plt.scatter(xs,ys)
    plt.title("ROC curve")
    plt.xlabel("1 - Specificity")
    plt.ylabel("Sensitivity")
    plt.show()


def q12():
    GSEA = gsea.GSEA()
    GSEA.load_data('GSE25628_filtered_expression.txt', 'GSE25628_samples.txt', 'c2.cp.kegg.v6.2.symbols.filtered.gmt')

    # all_genes = list(GSEA.all_gene_sets.values())
    #
    # unique_genes = set()
    # for gene_set in all_genes:
    #     for gene in gene_set:
    #         unique_genes.add(gene.strip())

    unique_genes = ft.reduce(set.union, GSEA.all_gene_sets.values(), set())
    unique_genes = set([s.strip() for s in unique_genes])

    print(len(unique_genes))
    return unique_genes

def q13():
    GSEA = gsea.GSEA()
    GSEA.load_data('GSE25628_filtered_expression.txt', 'GSE25628_samples.txt', 'c2.cp.kegg.v6.2.symbols.filtered.gmt')

    all_genes = [list(s) for s in list(GSEA.all_gene_sets.values())]
    big_group = []
    for gene_set in all_genes:
        big_group.extend(gene_set)
    big_group = [s.strip() for s in big_group]

    plt.hist(list(Counter(big_group).values()), bins=30)
    plt.title('Gene Frequency')
    plt.xlabel('Number of Genes')
    plt.ylabel('Occurrences')
    plt.show()



    unique_genes = q12()

    # for gene_set in all_genes:
    #     for gene in gene_set:
    #         unique_genes.append(gene.strip())
    #
    # for gene in unique_genes:
    #     for exp_gene in


q13()
