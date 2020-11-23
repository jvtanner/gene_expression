import pandas as pd
import math
import numpy as np
import sys
import time

group = 'c2.cp.kegg.v6.2.symbols.filtered.gmt'
expression = 'GSE25628_filtered_expression.txt'
sample = 'GSE25628_samples.txt'

def debug(input):
    if False:
        print(input)

def timer(start, end):
    print(end - start)

class GSEA(object):
    """
    asdfasdfa
    """

    def __init__(self):
        """
        Initialize GSEA object
        """

        # What needs to be returned at the very end
        self.list_sigs = []

    def load_data(self, expfile: str, sampfile: str, keggfile: str):
        """
        Take the file paths to the expression, sample and gene set data,
        read them in and store within the GSEA instance.
        :return:
        """
        self.expfile = expfile
        self.sampfile = sampfile
        self.keggfile = keggfile
        self.labels = {}

        # Parse the expression file to a Dataframe
        with open(self.expfile, 'r') as f:
            self.express = pd.read_csv(f, delimiter='\t', index_col=0)

        # Parse the sample file to a dictionary
        with open(self.sampfile, 'r') as g:
            for line in g:
                (key, value) = line.split()
                self.labels[key] = int(value)

        # Parse the kegg file to a dictionary of sets
        self.all_gene_sets = {}
        with open(self.keggfile, 'r') as h:
            lines = h.readlines()
            for line in lines:
                self.all_gene_sets[line.split('\t')[0]] = set(line.split('\t')[2:])

        # List of all the relevant gene groups, to be used in other functions
        self.gene_groups = list(self.all_gene_sets.keys())
        # List of genes we'll be examining from groups, to be used in other functions
        self.genes_in_groups = list(self.all_gene_sets.values())


    def get_gene_rank_order(self):
        """
        Return a list of all genes (as strings) ranked by their logFC between patient and control,
        with the gene with the highest logFC ordered first.
        :return: list of genes ranked by logFC difference in descending order
        """

        # Split expression dataframe into 2 parts: 1 for healthy, 1 for disease:
        # The DF of healthy people will be the original express DF - diseased columns

        start = time.time()

        healthy_expression = self.express.drop(columns=[disease for disease in self.disease_group])
        # The DF of diseased people will be the original express DF - healthy columns
        disease_expression = self.express.drop(columns=[healthy for healthy in self.healthy_group])

        # debug('splitting time: {}'.format(time.time() - start))

        # New Dataframe with averages for healthy and disease group's average gene expression.
        averages = pd.DataFrame(healthy_expression.mean(axis=1), columns=['healthy_mean'])
        averages['disease_mean'] = disease_expression.mean(axis=1)

        # Add column for logFC
        ### POSSIBLE THAT NO LOGS ARE NEEDED HERE ###

        start1 = time.time()

        averages['logFC'] = averages['disease_mean'] - averages['healthy_mean']

        print(averages.at['BMP4', 'healthy_mean'])

        # debug('make avg col: {}'.format(time.time() - start1))

        # Rank genes according to logFC value
        sorted_by_rank = averages.sort_values(by=['logFC'], ascending= False)

        # debug('sorted by rank: {}'.format(sorted_by_rank))
        # Debug
        # debug("done gene_rank_order")

        # Return list of genes in order of significance
        # debug('list of genes by ES: {}'.format(sorted_by_rank.index.values))
        return sorted_by_rank.index.values

    def get_enrichment_score(self, gene_group: str):
        """
        Return the enrichment score, a float correct to two decimal places for a given gene set,
        such as ‘KEGG_CITRATE_CYCLE_TCA_CYCLE’ (which is the string for the gene set name corresponding to the gene set).
        This method should run get_gene_rank_order at some point to initiate enrichment calculations.
        :return: enrichment score, float rounded to 2 decimal places
        """
        # List of all ranked genes

        start0 = time.time()

        ranked_genes = self.get_gene_rank_order()

        # debug('gene ranking: {}'.format(time.time() - start0))

        # list of all genes in expression data
        all_expression_genes = self.express.index.values

        # start = time.time()
        # Exclude the genes in the group that aren't in the expression data
        rel_group_genes = set()
        # For each gene within the group,
        for gene in self.all_gene_sets[gene_group]:
            # If the gene is also in the list of genes from the expression data, keep it
            if gene in all_expression_genes:
                rel_group_genes.add(gene)

        # Debug
        # print(time.time() - start)

        # Total number of genes in ranking
        NT = len(ranked_genes)
        # Total number of genes in the group we care about
        NG = len(rel_group_genes)
        # Establish penalty values
        penalty_hit = math.sqrt(float(NT-NG)/float(NG))
        penalty_miss = -math.sqrt(float(NG)/float(NT-NG))

        # Move down ranked list, add to the running sum the appropriate penalty
        # depending on whether the ranked gene is part of the gene group
        start = time.time()

        running_score = 0.0
        supremum = 0.0
        for ranked_gene in ranked_genes:
            if ranked_gene in rel_group_genes:
                running_score += penalty_hit
                if running_score > supremum:
                    supremum = running_score
            else:
                running_score += penalty_miss

        # debug('ES for gene group {}: {}'.format(gene_group, (time.time() - start)))

        return float(round(supremum, 2))

    def get_sig_sets(self, p):
        """
        Return the list of significant gene sets (as strings), at a corrected threshold of p, by name.
        If no gene sets are significant, return an empty list.
        This method should run get_gene_rank_order and/or get_enrichment_score at some point to initiate
        enrichment calculations and then identify significant gene sets.
        :return: significant gene groups, list of strings
        """
        highest_score = 0.0
        baddest_gs = ''

        # Store copy of correct self.express
        express = self.express

        # Corrected p-value, Bonferroni method
        bonferroni_p = p / len(self.gene_groups)

        self.healthy_group = set()
        self.disease_group = set()
        # Get the mean gene expression values for both healthy and diseased groups
        # Separate the healthy from the disease groups in sets
        for patient in self.labels:
            if self.labels[patient] == 0:
                self.healthy_group.add(patient)
            else:
                self.disease_group.add(patient)

        # Make a list of 100 permutations of self.express
        permutation_list = []

        for i in range(100):
            # Make a copy of the expression table
            permutation = self.express.copy()
            # Shuffle the column names
            permutation.columns = np.random.permutation(self.express.columns)
            # Append permutation version to main list
            permutation_list.append(permutation)

        # For each gene group,
        for gene_group in self.gene_groups:

            start = time.time()
            debug("Processing gene group: {}".format(gene_group))

            # Get the enrichment score
            ES_score = self.get_enrichment_score(gene_group)

            if ES_score > highest_score:
                baddest_gs = gene_group
                highest_score = ES_score

            # debug('Enrichment score complete for {}: {}'.format(gene_group, (time.time() - start)))

            start1 = time.time()
            false_score_list = []
            for perm in permutation_list:
                self.express = perm
                start = time.time()
                false_score = self.get_enrichment_score(gene_group)
                # debug('false score took: {}'.format(time.time() - start))
                false_score_list.append(false_score)
            # Reset self.express to be ready for next gene group
            self.express = express
            debug('Actual ES: {}   False ES list: {}'.format(ES_score, false_score_list))
            # debug('FALSE SCORE COMPLETE FOR {}: {}'.format(gene_group, (time.time() - start1)))

            p_count = 0
            # From the list of false scores, check how many are greater than the true Enrichment Score
            for score in false_score_list:
                if score >= ES_score:
                    p_count += 1
            p_value = p_count / (len(false_score_list) * len(self.gene_groups))
            # If the p-value is less than the corrected p-value, add the name of gene group to list sig_ESs
            debug('p-value: {},    Bonferroni corrected p-value: {}'.format(p_value, bonferroni_p))
            if p_value <= bonferroni_p:
                self.list_sigs.append(gene_group)
        debug(self.list_sigs)
        print('highest scoring gene set: ', baddest_gs)
        return self.list_sigs


def main():
    if (len(sys.argv) != 4):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    expfile = sys.argv[1]
    sampfile = sys.argv[2]
    groupfile = sys.argv[3]

    # Create instance of class
    gsea = GSEA()
    gsea.load_data(expfile, sampfile, groupfile)

    gsea.get_sig_sets(.05)

if __name__ == '__main__':
    main()