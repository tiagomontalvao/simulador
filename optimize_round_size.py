#!/usr/bin/env python3

import subprocess as sp
import numpy as np
import scipy.stats as st
from sys import argv

def get_auto_cov(x, mean):
	n = len(xs)
	xi_sum = np.sum([(x[i] - mean) * (x[i+1] - mean) for i in range(n-1)])
	return xi_sum/(n-2)

if __name__ == "__main__":
	try:
		metric = argv[1]
	except:
		raise ValueError()
	
	warmup = 8000
	rounds = 50
	rho = 0.1
	interrupt = 0

	def next_round_size():
		k = 500
		for i in range(500):
			yield int(k)
			k *= 1.7

	interrupt_str = "ON" if interrupt else "OFF"
	print('Parameter: ' + metric + f' (warmup={warmup}, #rounds={rounds}, rho1={rho}, interruption {interrupt_str})')
	print('{:<10} {:<14} {:<14} {:<12} {:<10} {:^5}'.format('Round size', 'Mean', 'Autocovariance', 'Variance', 'IC', 'Prec'))
	for round_size in next_round_size():
		acc_avg = dict()
		op = sp.run('/home/renato/Desktop/simulador/simulador {} {} {} {} {}'.format(warmup, rounds, round_size, rho, interrupt), stdout=sp.PIPE, shell=True, check=True)

		x = []
		for line in op.stdout.decode('utf8').strip().split('\n'):
			round_info = line.split()
			metrics = [round_info[i].replace("rm.", "") for i in range(0, len(round_info), 2)]
			round_means = [float(round_info[i]) for i in range(1, len(round_info), 2)]

			x.append(dict(zip(metrics, round_means)))

		xs = [xi[metric] for xi in x]
		mean = np.mean(xs)
		var = np.var(xs, ddof=1)

		ic = st.t.interval(0.9, len(xs)-1, scale=np.sqrt(var))
		prec = ic[1]/mean * 100

		print('{:<10} {:<14.6f} {:<14.6f} {:<12.6f} ±{:<10.6f}{:>5.2f}%'.format(round_size, mean, get_auto_cov(xs, mean), var, ic[1], prec))

		if (prec <= 5.0): break


