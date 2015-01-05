#!/usr/bin/python
"""
Author: Jian Dai
"""
import math
import numpy as np
import matplotlib.pyplot as plt

class PisemaEllipse:
	def getSigmaTilde(self, sigma):
		return (sigma - self.sigma11)/(self.sigma33 - self.sigma11)
	
	def getNuTilde(self, nu):
		return (2.0*nu/self.nuParallel + 1.0)/3.0

	def getPointTilde(self, p):
		sigma = p[0]
		nu = p[1]
		return [self.getSigmaTilde(sigma), self.getNuTilde(nu)]

	def __init__(self, s11, s22, s33, nuParallel, infilename):
		self.sigma11 = s11
		self.sigma22 = s22
		self.sigma33 = s33
		self.nuParallel = nuParallel

		betaNHRad = 17.0*math.pi/180.0
		self.sinBetaNHRad = math.sin(betaNHRad)
		self.cosBetaNHRad = math.cos(betaNHRad)
		self.sinBetaNHRad2 = math.sin(2.0*betaNHRad)
		self.cosBetaNHRad2 = math.cos(2.0*betaNHRad)

		sigma = np.linspace(0, 250, 100)
		self.sigmaTilde = self.getSigmaTilde(sigma)
		nu = np.linspace(-1.0*self.nuParallel, 2.0*self.nuParallel, 100)[:, np.newaxis]
		self.nuTilde = self.getNuTilde(nu)

		self.P = self.getPointTilde([self.sigma33, self.nuParallel*(3.0*self.cosBetaNHRad**2 - 1.0)/2.0])
		self.Q = self.getPointTilde([self.sigma22, -0.5*nuParallel])
		self.R = self.getPointTilde([self.sigma11, self.nuParallel*(3.0*self.sinBetaNHRad**2 - 1.0)/2.0])
		self.S = self.getPointTilde([self.sigma33 * self.sinBetaNHRad**2 + self.sigma11 * self.cosBetaNHRad**2, -0.5*self.nuParallel])
		self.T = self.getPointTilde([self.sigma33 * self.cosBetaNHRad**2 + self.sigma11 * self.sinBetaNHRad**2, self.nuParallel])

		self.expDataLst = []
		for line in open(infilename):
			resid = int(line.split()[0])
			cs = float(line.split()[1])
			dc = float(line.split()[2])
			oneRow = [resid, cs, dc]
			(self.expDataLst).append(oneRow)

	def plotEllipse(self):
		fig, ax = plt.subplots(figsize=(8,8))
		cs = ax.contour(self.sigmaTilde, (self.nuTilde).ravel(), \
				(self.nuTilde - self.cosBetaNHRad2 * self.sigmaTilde - self.sinBetaNHRad**2)**2 \
				- self.sigmaTilde * (1.0 - self.sigmaTilde) * self.sinBetaNHRad2**2, [0], colors='black')
		ps = cs.collections[0].get_paths()[0]
		vs = ps.vertices
		outfile = open('ellipse.dat', 'w')
		for v in vs:
			print>>outfile, '%8.3f%8.3f' % (v[0], v[1])
		outfile.close()
		for (point, color) in zip([self.P, self.Q, self.R, self.S, self.T], ['blue','red','green','yellow','orange']):
			ax.scatter(point[0], point[1], color=color, s=30)
		
		for expData in self.expDataLst:
			resid = int(expData[0])
			point = self.getPointTilde([expData[1], expData[2]])
			ax.scatter(point[0], point[1], color='magenta', marker='*', s=30)
			#if resid in [10,12,13,16,17,20,21,24]:
			#	ax.scatter(point[0], point[1], color='magenta', marker='*', s=30)
			#else:
			#	ax.scatter(point[0], point[1], color='gray', marker='*', s=30)

		plt.gca().invert_xaxis()
		plt.savefig('plot_ellipse.png', dpi=300)
		plt.show()
		
#pisemaEllipse = PisemaEllipse(31.0, 54.0, 202.0, 10.735, 'kdpf_pisema.dat')
pisemaEllipse = PisemaEllipse(57.3, 81.2, 227.8, 10.735, 'kdpf_pisema.dat')
pisemaEllipse.plotEllipse()
