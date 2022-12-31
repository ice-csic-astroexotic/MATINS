import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.tri as tri

file = np.genfromtxt(fname="Bmap_surf.dat")
theta = file[:,0]
phi = file[:,1]
br = file[:,2]
N = len(theta)

# Try a dipole:
#for i in range(0,len(theta),1):
#  br[i] = np.cos(theta[i])*np.sin(phi[i])

theta = math.pi/2 - theta
phi = np.asarray(phi) - math.pi

p = []
t = []
v = []

for i in range(0,N,1):
  if (phi[i] == -math.pi):
    p.append(math.pi)
    t.append(theta[i])
    v.append(br[i])

theta = np.concatenate([theta,np.asarray(t)])
phi   = np.concatenate([phi,np.asarray(p)])
br    = np.concatenate([br,np.asarray(v)])


Ncont=100

fig = plt.figure() 



ax = fig.add_subplot(221, projection='aitoff')
#ax.scatter(phi,theta,c=br, cmap=plt.cm.Spectral_r)
ax.tricontour(phi, theta, br, Ncont, linewidths=0.5, cmap=plt.cm.Spectral_r)
#ax.tricontourf(phi, theta, br, Ncont, cmap=plt.cm.Spectral_r)
ax.set_title("aitoff")
ax.set_xticks([])
ax.set_yticks([])
#ax.grid(True)

ax = fig.add_subplot(222, projection='hammer') 
ax.tricontour(phi, theta, br, Ncont, linewidths=0.5, cmap=plt.cm.Spectral_r)
ax.set_title("hammer")
ax.set_xticks([])
ax.set_yticks([])
#ax.grid(True) 

ax = fig.add_subplot(223, projection='lambert') 
ax.tricontour(phi, theta, br, Ncont, linewidths=0.5, cmap=plt.cm.Spectral_r)
ax.set_title("lambert")
ax.set_xticks([])
ax.set_yticks([])
#ax.grid(True) 

ax = fig.add_subplot(224, projection='mollweide') 
ax.tricontour(phi, theta, br, Ncont, linewidths=0.5, cmap=plt.cm.Spectral_r)
ax.set_title("Mollweide")
ax.set_xticks([])
ax.set_yticks([])
#ax.grid(True) 

#plt.xlabel('R.A.') 
#plt.ylabel(r'Decl.') 
plt.savefig("proyecciones.pdf",  bbox_inches="tight")
plt.close()
