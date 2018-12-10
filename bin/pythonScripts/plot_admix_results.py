import matplotlib.ticker as plticker
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import numpy as np
import pandas as pd
from scipy import stats

admixed_lst = []
pure_lst = []
r_squared_lst = []

def plot_admix_results(admixed, pure):
	estimates = pd.read_csv(
		'pops_data/admixture/CEU_YRI_admixed_{}admixed_{}pure.2.Q'.format(admixed, pure),
		header=None,
		delim_whitespace=True,
		names=['x','y'])

	true = pd.read_csv(
		'pops_data/admixed/CEU_YRI_admixed_{}admixed_{}pure_proportions.txt'.format(admixed, pure),
		header=None,
		delim_whitespace=True,
		names=['x','y'])

	_, _, r_value, _, _ = stats.linregress(estimates.x, true.x)
	r_squared = round(r_value**2, 2)

	admixed_lst.append(admixed)
	pure_lst.append(pure)
	r_squared_lst.append(r_squared)

	fig, ax = plt.subplots()
	ax.scatter(estimates.x, true.x)
	plt.xlabel('estimate')
	plt.ylabel('true')

	blank = Rectangle((0, 0), 0.1, 0.1, fc="w", fill=False, edgecolor='none', linewidth=0)
	ax.legend([blank], ['r^2: {}'.format(r_squared)], loc='lower center')

	plt.suptitle('admix results with {} test chromosomes and {} total pure chromosomes'.format(admixed, pure))
	# plt.show()
	plt.savefig('etc/admix/{}admixed_{}totalpure_MM.png'.format(admixed, pure))
	plt.clf()

def plot_r_squared():
	fig = plt.figure()
	ax = fig.add_subplot(111)

	plt.scatter(pure_lst, admixed_lst)
	for (adm, pur, rsq) in zip(admixed_lst, pure_lst, r_squared_lst):
		ax.annotate(str(rsq), xy=(pur, adm), textcoords='data') # label with r-squared values

	plt.xlabel('total pure chromosomes')
	plt.ylabel('total admixed chromosomes')
	plt.title('r^2 values for ADMIXTURE results\non test chromosome sets')

	loc = plticker.MultipleLocator(base=20) # this locator puts ticks at regular intervals
	ax.xaxis.set_major_locator(loc)

	plt.grid()

	filename = 'etc/admix/r_squared_results_MM.png'
	print('saving r^2 figure to {}'.format(filename))
	plt.savefig(filename)
	
	plt.show()

to_plot = (
	(10,300),
	(200,20)
	# (100, 20),
	# (100,30),
	# (100,40),
	# 
	# (200,20),
	# (200,30),
	# (200,40),
	# (200,80),
	# 
	# (250,20),
	# 
	# (300,20),
	# (300,30),
	# (300,40),
	# 
	# (350,20),
	# 
	# (400,20),
	# (400,30),
	# (400,40),
	# 
	# (500, 0),
	# (500, 20),
	# (500,40),
	# (500, 100),
	#(500, 200),
)

for admixed, pure in to_plot:
	plot_admix_results(admixed, pure)

plt.close('all')
plot_r_squared()