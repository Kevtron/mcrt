#!/usr/bin/env python

import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from time import time
import sys


class Photon(object):
	def __init__(self,mubins,bins):
		self.escaped=False
		self.cost = np.sqrt(np.random.random_sample())
		self.phi = 2.0 * np.pi * np.random.random_sample()
		self.sint = np.sqrt(1.0 - self.cost * self.cost)
		self.cosp = np.cos(self.phi)
		self.sinp = np.sin(self.phi)
		self.xt = 0.0
		self.yt = 0.0
		self.zt = 0.0
		self.bins=bins

	def scatter(self):
		self.cost = 2.0 * np.random.random_sample() - 1.0
		self.sint = np.sqrt(1.0 - self.cost * self.cost)
		self.phi = 2.0 * np.pi * np.random.random_sample()
		self.cosp = np.cos(self.phi)
		self.sinp = np.sin(self.phi)

	def propogate(self,albedo,taumax):
		while self.zt <=1.0 and self.zt>=0.0:
			x=self.xt
			y=self.yt
			z=self.zt
			tau=-1.0*np.log(np.random.random_sample())
			s=tau/taumax
			self.xt=x+s*self.sint*self.cosp
			self.yt=y+s*self.sint*self.sinp
			self.zt=z+s*self.cost
			#print self.xt,self.yt,self.zt
			if np.random.random_sample() < albedo and self.zt>0.0 and self.zt<1.0:
				self.scatter()		
			if self.zt>=1.0:
				self.escaped=True
				#print self.cost,self.phi
			elif self.zt<=0.0:
			#	print "Re-entered"
				self.__init__(mubins,phot.bins)

	def bin(self,mubins):
		j=abs(int(self.cost*mubins))
		self.bins[j]=self.bins[j]+1

def binphotons(mu, mubins,n):#bin the photons uniformly from 0-90 degrees.
	dtheta=1.0/mubins
	width=0.5*dtheta #see ken wood's slab code.
	theta=np.zeros(mubins)
	bins=np.zeros(mubins)
	for i in range(mubins):
		theta[i]=np.arccos(i*dtheta+width)*(180.0/np.pi)
	for i in range(int(n)):
		#if mu[i]>=0:
		j=abs(int(mu[i]*mubins))
		bins[j]=bins[j]+1
	return bins, theta

if __name__ == "__main__":	
	nphotons=10000		
	i=0
	mubins=10
	phibins=10
	mulist=np.empty(nphotons)
	philist=np.empty(nphotons)
	dtheta=1.0/mubins
	width=0.5*dtheta #see ken wood's slab code.
	theta=np.zeros(mubins)
	bins=np.zeros(mubins)
	for i in range(mubins):
		theta[i]=np.arccos(i*dtheta+width)*(180.0/np.pi)
	i=0
	while i<nphotons:
		phot=Photon(mubins,bins)
		phot.propogate(1.0,71.0)
		if phot.escaped:
			mulist[i]=phot.cost
			philist[i]=phot.phi
			phot.bin(mubins)
			i+=1
			phot.__init__(mubins,phot.bins)
		if i%100==0:
			print i
	intensity = np.zeros(mubins)
	for j in range(mubins):
		intensity[j]=(muhist[j]*mubins)/(2*sum(muhist)*np.cos(theta[j]*np.pi/180.0))

	fig=plt.figure(figsize=(8,8))
	plt.ylabel(r'$\frac{I_\nu}{B_\nu}$')
	plt.xlabel(r'$\theta$')
	plt.plot(theta,intensity,'bo')
	f=open('intensity-oop.dat','w')
	for i in range(theta.size):
		f.write("%4.2f\t%6.4f\n" % (theta[i], intensity[i]))
	f.close()
	fig.savefig('intensity.png')
	fig.clf()
	plt.plot(theta,muhist/nphotons,'bo')
	fig.savefig('energy.png')
