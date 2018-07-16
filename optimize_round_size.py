#!/usr/bin/env python3

import os
import subprocess as sp
import numpy as np
import scipy.stats as st
from sys import argv

def get_auto_cov(x, mean):
	n = len(x)
	xi_sum = np.sum([(x[i] - mean) * (x[i+1] - mean) for i in range(n-1)])
	return xi_sum/(n-2)

nm = []

if __name__ == "__main__":
	try:
		metric = argv[1]
	except:
		raise ValueError()

	devnull = open(os.devnull, 'w')

	# sp.run('make clean', shell=True, stdout=devnull)
	sp.run('make CPPFLAGS="-D PYTHON_SCRIPT"', shell=True, stdout=devnull)
	
	warmup = 500000
	rho = 0.7
	interrupt = 0

	def next_round_size():
		k = 50
		for i in range(500):
			yield int(k)
			k *= 1.5

	interrupt_str = "ON" if interrupt else "OFF"
	print('Parameter: ' + metric + f' (warmup={warmup}, rho1={rho}, interruption {interrupt_str})')
	print('{:<10} {:<10} {:<10} {:<14} {:<12} {:<10} {:^5}'.format('Rounds', 'Round size', 'Mean', 'Autocovariance', 'Variance', 'IC', 'Prec'))
	for round_size in next_round_size():
		acc_avg = dict()
		op = sp.run('./simulador {} {} {} {}'.format(warmup, round_size, rho, interrupt), stderr=sp.PIPE, stdout=devnull, shell=True, check=True)

		x = []
		for line in op.stderr.decode('utf8').strip().split('\n'):
			round_info = line.split()
			metrics = [round_info[i].replace("rm.", "") for i in range(0, len(round_info), 2)]
			round_means = [float(round_info[i]) for i in range(1, len(round_info), 2)]

			x.append(dict(zip(metrics, round_means)))

		xs = dict()
		mean = dict()
		autocov = dict()
		var = dict()
		for param in x[0].keys():
			xs[param] = [xi[param] for xi in x]
			mean[param] = np.mean(xs[param])
			autocov[param] = get_auto_cov(xs[param], mean[param])
			var[param] = np.var(xs[param], ddof=1)
		
		ic = st.t.interval(0.9, len(xs[metric])-1, scale=np.sqrt(var[metric]/len(xs[metric])))
		prec = ic[1]/mean[metric] * 100

		print('{:<10} {:<10} {:<10.6f} {:<14.6f} {:<12.6f} ±{:<10.6f}{:>5.2f}%'.format(len(xs[metric]), round_size, mean[metric], autocov[metric], var[metric], ic[1], prec))

		if any(var[param] / abs(autocov[param]) < 5 for param in xs.keys()):
			continue
		else:
			break

	# sp.run('make clean', shell=True, stdout=devnull)
