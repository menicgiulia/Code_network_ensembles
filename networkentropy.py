from __future__ import division
from pylab import *
import networkx as nx
import time
from numexpr import evaluate as ev
from scipy.spatial.distance import cdist, pdist, squareform

#configuration ensemble
def entropy_conf(PP):

	N=PP.shape[0];
	connectivity = sum(PP,1); 
	avg_conn = mean(connectivity)
	#Lagrangians
	z = connectivity/(sqrt(avg_conn*N))
	old_z = z
	loops = 10000
	precision = 1e-5
	
	for idx in range(loops):
		zT = z[:,np.newaxis]  
		D = ev("(zT * z) + 1.")
		UD = ev("z/D")
		del D
		for i in xrange(N):  UD[i,i]=0.
		z = connectivity / sum(UD,1)
		rz= ev("abs(1.-z/old_z)")
		rz[np.isinf(rz)]=0.
		rz[np.isnan(rz)]=0.

		if max(rz)<precision:
			break
		old_z = z
	z2 = outer(z,z)
 	for i in xrange(N):  z2[i,i]=0.
	P = z2 / ( z2 + 1)
	Q=1.-P
	S = -ev("sum(log(P**P) + log(Q**Q) )")/2.
	print "number of loops ", idx+1
	return S







#spatial ensemble
def entropy_dist(PP, distance, Nbin):

	#linear binning
	mi,ma = np.min(distance[distance>0]),np.max(distance)
	limiti = linspace(mi,ma,Nbin+1)
	limiti[-1]+=0.1*ma
	limiti[0]-=0.1*mi


	b = searchsorted(limiti,distance)-1
	massimo = np.max(b)+1
	BC = [ sum(b[PP>0]==i)/2 for i in range(massimo) ]

	N=PP.shape[0];
	connectivity = sum(PP,1);
	avg_conn = mean(connectivity)

	#Lagrangians
	z = connectivity/(sqrt(avg_conn*N))
	w = BC / (avg_conn*N)

	old_z = z
	old_w = w

	loops = 10000
	precision = 1e-5


	for idx in range(loops):
		bigW = w.take(b)

		for i in xrange(N):  bigW[i,i]=0.
		
		U = ev("bigW * z")
		UT = U.T    
		D = ev("(UT * z) + 1.")
	
		UD = ev("U/D")
	
		del D,U,UT

		for i in xrange(N):  UD[i,i]=0.

		z = connectivity / sum(UD,1)
	
		zt = z[:].reshape(N,1)
		D2 = ev("(z*zt) / ( bigW *(z*zt) + 1.)")
	
		B2 = array([ ev("sum(where(b==i,D2,0.))") for i in range(len(w)) ])/2.
	
		w = ev("where( (BC!=0) & (B2!=0),BC/B2,0 )")
		rz= ev("abs(1.-z/old_z)")
		rw= ev("abs(1.-w/old_w)")
		rz[np.isinf(rz)]=0.
		rw[np.isinf(rw)]=0.
		rz[np.isnan(rz)]=0.
		rw[np.isnan(rw)]=0.
		
		if max(rz)<precision and max(rw)<precision:
			break
	
		old_z = z
		old_w = w


	bigW = w.take(b)
	for i in xrange(N):  bigW[i,i]=0.

	z2 = bigW * outer(z,z)
	P = z2 / ( z2 + 1)
	Q=1.-P
	S = -ev("sum(log(P**P) + log(Q**Q) )")/2.
	print "number of loops ", idx+1

	return S

  

#test function
if __name__ == "__main__":
    print"testing networkentropy.py"
    ba=nx.barabasi_albert_graph(100,2)
    barabasi=np.array(nx.to_numpy_matrix(ba))


    print "100 nodes, barabasi-albert"
    degree_sequence=sorted(nx.degree(ba).values(),reverse=True) # degree sequence

    # matrix with random positive weights
    vec_pos = np.random.rand(100)*5
    Y = pdist(vec_pos[:, np.newaxis], 'euclidean')
    Y[Y==0]+= 2*np.finfo(float).eps
    dmatrix=squareform(Y)
    
    
    #adjacency matrix
    print "Configuration ensemble"
    S_conf=entropy_conf(barabasi)
    print "Configuration entropy", S_conf
    
    
    print "Spatial ensemble"
    S_dist=entropy_dist(barabasi, dmatrix, 10)
    print "Spatial entropy", S_dist

