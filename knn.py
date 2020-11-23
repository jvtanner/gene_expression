import pandas as pd
import sys
from math import sqrt
import numpy as np


def euclidian(col1, col2) -> float:
    """
    return the Euclidian distance between two different variables
    :param col1: gene expression of person1
    :param col2: gene expression of person2
    :return: Euc distance
    """
    return sqrt(sum([(col1[i]-col2[i])**2 for i in range(len(col1)-1)]))


class KNN(object):
    """
    Object to run K Nearest Neighbors algorithm
    """

    def __init__(self):
        """
        Initialize knn object.
        """

    def load_data(self, expfile: str, sampfile: str):
        """
        Reads in the paths to the expression and sample files, then stores them.
        :return:
        """

        self.expfile = expfile
        self.sampfile = sampfile

        self.labels = {}

        # Gene data
        with open(self.expfile, 'r') as f:
            self.express = pd.read_csv(f, delimiter='\t', index_col=0)

        # Labels
        with open(self.sampfile, 'r') as g:
            for line in g:
                (key, val) = line.split()
                self.labels[key] = int(val)

    def get_assignments(self, k: int, fn: float):
        """
        Returns the class assignments for all samples for given values of k and fn as a list of integers 0's and 1's
        :return: list of ints
        """

        # Store the distance between each column and all its neighbors

        self.k = k
        self.fn = fn
        names = self.express.columns.values

        all_distances = []
        dict = {}

        # For each of the patient's gene expression values,
        for h, i in self.express.iteritems():
            distances = []
            # Loop through the other patients' values,
            for _, j in self.express.iteritems():
                distances.append(euclidian(i, j))
            dict[h] = distances
            # Append to master list of distances for each patient, in correct order
        # print(dict)
        for label in self.labels:
            all_distances.append(dict[label])

        # Store the indices of the smallest k distances
        neighbors_index = []
        for i in range(len(all_distances)):
            # For each column store the indices of the k lowest distance values
            neighbors_index.append(np.argpartition(all_distances[i], self.k+1)[:self.k+1])

        # Make prediction based on the labels of the nearest neighbors
        end_predictions = []
        # Loop through each patient
        for h in range(len(neighbors_index)):
            score = []
            patient_bucket = neighbors_index[h]
            # Store the k distance values
            for i in range(1, self.k+1):
                name = names[patient_bucket[i]]
                score.append(self.labels[name])

            # Take average value of the neighbor labels
            raw_prediction = np.mean(score)
            # Make a prediction based on the average value and fn cutoff
            if raw_prediction > self.fn:
                end_predictions.append(1)
            else:
                end_predictions.append(0)
        # print(end_predictions)
        return end_predictions


    def calc_metrics(self, k: int, fn: float):
        """
        Return a list of float values [sensitivity, specificity} of a knn classifier using the values of k and fn.
        :return: list of float values for specificity and sensitivity
        """

        predictions = self.get_assignments(k, fn)
        reality = list(self.labels.values())

        TP = 0
        TN = 0
        FP = 0
        FN = 0

        # Get values for True Positives, True Negatives, False Positives, and False Negatives
        for i in range(len(predictions)):
            if predictions[i] == 1 and reality[i] == 1:
                TP += 1
            elif predictions[i] == 1 and reality[i] == 0:
                FP += 1
            elif predictions[i] == 0 and reality[i] == 1:
                FN += 1
            else:
                TN += 1

        sensitivity = float(TP/(TP+FN))
        specificity = float(TN/(TN+FP))

        # print([sensitivity, specificity])
        return [sensitivity, specificity]

def main():
    """
    Run KNN and return a prediction for all of the samples in the dataset.
    """

    # Check that arguments are in present in the command line as expected
    if (len(sys.argv) != 3):
        print("Please specify an input file and an output file as args.")
        return

    # input variables
    expfile = sys.argv[1]
    sampfile = sys.argv[2]

    # Define constants
    k = 3
    fn = 0.8

    # create an knn object and run
    classifier = KNN()
    classifier.load_data(expfile, sampfile)
    classifier.calc_metrics(k, fn)


if __name__ == "__main__":
    main()
