import numpy as np
import matplotlib.pylab as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
from time import time
import sys

def newphoton():
	cost = np.sqrt(np.random.random_sample())
	phi = 2.0 * np.pi * np.random.random_sample()
	sint = np.sqrt(1-cost*cost)
	cosp = np.cos(phi)
	sinp = np.sin(phi)
	xt, yt, zt = 0.0,0.0,0.0
	return cost, sint, cosp, sinp, xt, yt, zt
	
def scatter():
	cost = 2.0*np.random.random_sample()-1.0
	sint = np.sqrt(1.0-cost*cost)
	phi = 2.0*np.pi*np.random.random_sample()
	cosp = np.cos(phi)
	sinp = np.sin(phi)
	return cost, sint, cosp, sinp, phi

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
	start = time()
	taumax = 71.0    
	nphotons = float(sys.argv[1])
	mubins = 10
	phibins = 10
	mu1=np.array(np.zeros(nphotons))
	phi1=np.array(np.zeros(nphotons))
	i=1
	albedo = 1
	while i <= nphotons:
		cost, sint, cosp, sinp, xt, yt, zt = newphoton()
		###get scattered on
		while zt>=0.0 and zt<=1.0:
			x=xt
			y=yt
			z=zt
			tau=-1.0*np.log(np.random.random_sample())
			s=tau/taumax
			xt=x+s*sint*cosp
			yt=y+s*sint*sinp
			zt=z+s*cost
			if np.random.random_sample() < albedo and zt>0.0 and zt<1.0:
				cost,sint,cosp,sinp,phi = scatter()
		if zt>=1.0:
			mu1[i-1]=cost
			phi1[i-1]=phi 
			i=i+1
	muhist,theta=binphotons(mu1,mubins,nphotons)
	intensity = np.zeros(mubins)
	for j in range(mubins):
		intensity[j]=(muhist[j]*mubins)/(2*sum(muhist)*np.cos(theta[j]*np.pi/180.0))
	print "%d scatterings completed in %f s" % (sum(muhist), (time()-start))
	fig=plt.figure(figsize=(8,8))
	plt.ylabel(r'$\frac{I_\nu}{B_\nu}$')
	plt.xlabel(r'$\theta$')
	plt.plot(theta,intensity,'bo')
	f=open('intensity.dat','w')
	for i in range(theta.size):
		f.write("%4.2f\t%6.4f\n" % (theta[i], intensity[i]))
	f.close()
	fig.savefig('intensity.png')
	fig.clf()
	plt.plot(theta,muhist/nphotons,'bo')
	fig.savefig('energy.png')
