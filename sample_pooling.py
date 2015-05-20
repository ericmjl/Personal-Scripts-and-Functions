"""
Author: Eric J. Ma 
Date: 20 May 2015
For: Runstadler Lab

Assumptions:
1. prevalence is binomial where p = 0.1, n = 300.
2. prevalence is on a per cryo vial basis.
3. pooled samples are binned in order of collection with some bin size m.

Questions to answer:
1. Varying bin size m, what is the total number of extractions needed? 

To get help on this script, in the terminal run:

	python sample_pooling.py -h

To run the script, in the terminal run:

	python sample_pooling.py n p minb maxb

The parameters that are needed are:
- n: the total number of samples
- p: a prevalence estimate
- minb: minimum bin size you're willing to consider
- maxb: maximum bin size you're willing to consider

The dependencies apart from the standard Python 2.7 library that are needed are:
- numpy
- scipy
- pandas 
- matplotlib
"""

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import sys

from scipy.stats import bernoulli


def main(n=300, p=0.1, minb=1, maxb=30):
	# Input Parameters
	bin_sizes = np.arange(minb,maxb+1,1)

	data = dict()

	for bin_size in bin_sizes:
		data[bin_size] = list()
		# Perform x number of draws draws
		for i in range(200):
			# Prevalence draws
			rvs = bernoulli.rvs(p, size=n)
			num_extractions = len(rvs) / bin_size

			# Count the number of times a positive shows up in a pool
			num_round_1_positives = 0
			for i in range(0, len(rvs),bin_size):
				if np.sum(rvs[i:i+bin_size]) > 0:
					num_round_1_positives += 1

			# Compute the number of samples needed to re-screen.
			num_round_2_extractions = num_round_1_positives * bin_size

			# Compute the total number of extractions needed
			total = num_extractions + num_round_2_extractions
			data[bin_size].append(total)

	# Compute mean number of draws for each bin_size.
	data_mean = dict()
	for k, v in data.items():
		data_mean[k] = np.mean(v)

	# Make a dataframe and print to screen
	data = pd.DataFrame(data_mean.items())
	data.columns = ['Bin Size', 'Expected Num. of Extractions']
	data.set_index('Bin Size', inplace=True)
	print(data)
	print('')

	best_bin_size = data[data['Expected Num. of Extractions'] == data['Expected Num. of Extractions'].min()]
	print('Bin size that minimizes number of extractions: \n{0}'.format(best_bin_size))


	plt.plot(data_mean.keys(), data_mean.values())
	plt.show()


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Model parameters')
	parser.add_argument('n', metavar='n', type=int, help='an integer for the number of samples')
	parser.add_argument('p', metavar='p', type=float, help='a probability between 0 and 1 of finding a positive sample')
	parser.add_argument('minb', metavar='minb', type=int, help='the minimum bin size you would like to consider.')
	parser.add_argument('maxb', metavar='maxb', type=int, help='the maximum bin size you would like to consider.')

	args = parser.parse_args()
	n, p, minb, maxb = args.n, args.p, args.minb, args.maxb

	main(n, p, minb, maxb)