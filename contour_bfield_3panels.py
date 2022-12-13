"""
"""
import matplotlib 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap

def load_table_vertical(fname):  
  # Vertical slice
  lines = open(fname).readlines()

  #  Read the dimensions and the grid arrays (rows 0,1,2)
  nx, ny  = int(lines[0].split()[0]), int(lines[0].split()[1])

  r = np.array([float(x) for x in lines[1].split()]) - 12
  theta = np.array([float(x) for x in lines[2].split()])
  # read the B field components
  ntimes = int((len(lines)-3)/(nx*ny))
  ntimes = min(1000, ntimes)
  
  br = np.zeros([nx, ny, ntimes])
  bth = np.zeros([nx, ny, ntimes])
  bphi = np.zeros([nx, ny, ntimes])
  phi_sc = np.zeros([nx, ny, ntimes])
  psi_sc = np.zeros([nx, ny, ntimes])

  irow = 2
  for k in range(ntimes):
    for i in range(nx):
      for j in range(ny):
        #irow = ny*i + j 
        irow += 1
        br[i,j, k]=float(lines[irow].split()[0])
        bth[i,j, k]=float(lines[irow].split()[1])
        bphi[i,j, k]=float(lines[irow].split()[2]) 
        phi_sc[i,j, k]=float(lines[irow].split()[3]) 
        psi_sc[i,j, k]=float(lines[irow].split()[4]) 
  return nx, ny, ntimes, r, theta, br, bth, bphi, phi_sc, psi_sc

def load_table_equatorial(fname):  
  # Vertical slice
  lines = open(fname).readlines()
  #  Read the dimensions and the grid arrays (rows 0,1,2)
  nx, ny  = int(lines[0].split()[0]), int(lines[0].split()[1])
  r = np.array([float(x) for x in lines[1].split()]) - 12
  phi = np.array([float(x) for x in lines[2].split()])
  # read the B field components
  ntimes = int((len(lines)-3)/(nx*ny))
  ntimes = min(1000, ntimes)
  
  br = np.zeros([nx, ny, ntimes])
  bth = np.zeros([nx, ny, ntimes])
  bphi = np.zeros([nx, ny, ntimes])
  phi_sc = np.zeros([nx, ny, ntimes])
  psi_sc = np.zeros([nx, ny, ntimes])
  
  irow = 2
  for k in range(ntimes):
    for i in range(nx):
      for j in range(ny):
        #irow = ny*i + j 
        irow += 1
        br[i,j, k]=float(lines[irow].split()[0])
        bth[i,j, k]=float(lines[irow].split()[1])
        bphi[i,j, k]=float(lines[irow].split()[2]) 
        phi_sc[i,j, k]=float(lines[irow].split()[3]) 
        psi_sc[i,j, k]=float(lines[irow].split()[4]) 
  return nx, ny, ntimes, r, phi, br, bth, bphi, phi_sc, psi_sc

#def f(r,theta):
#    f = 1.e-2*((r-r[0])*(r[-1]-r))**2*np.sin(theta)*np.cos(theta)
#    return f

file_name = 'out/2D/b_merid_l0180_volume.dat' 
nx, ny, ntimes1, r, theta, br1, bth1, bphi1, phi1_sc, psi1_sc = load_table_vertical(file_name)
THETA, R1 = np.meshgrid(theta, r)
#
#  Change to cartesian coordinates 
#
X1=R1*np.sin(THETA)
Y1=R1*np.cos(THETA)

file_name = 'out/2D/b_merid_l90270_volume.dat'
nx, ny, ntimes2, r, theta, br2, bth2, bphi2, phi2_sc, psi2_sc = load_table_vertical(file_name)

file_name = 'out/2D/b_equator_volume.dat'
nx, ny, ntimes3, r, phi, br3, bth3, bphi3, phi3_sc, psi3_sc = load_table_equatorial(file_name)
PHI, R3 = np.meshgrid(phi, r)
#
#  Change to cartesian coordinates 
#
X3=R3*np.cos(PHI)
Y3=R3*np.sin(PHI)

# Take the minimum in case the files have different numbers of snapshots
ntimes = min(ntimes1, ntimes2, ntimes3)
#
# Plot figure
#
fig, (ax1, ax2, ax3) = plt.subplots(1,3,figsize=(17,8))
#
#   Choose the function you want to plot
#
func1 = psi1_sc
func2 = psi2_sc
func3 = psi3_sc
#print(np.max(func2),np.min(func2))
nlevels = 100

mylevels1 = np.linspace(np.min(func1), np.max(func1), nlevels) 
mylevels2 = np.linspace(np.min(func2), np.max(func2), nlevels) 
mylevels3 = np.linspace(np.min(func3), np.max(func3), nlevels) 

def animate(frame):
  ax1.cla()
  ax1.contourf(X1, Y1, func1[:,:,frame], levels=mylevels1, alpha=1, cmap=plt.cm.jet)
  C = ax1.contour(X1, Y1, func1[:,:,frame], mylevels1[0:nlevels:8] , colors='black')
  ax1.clabel(C, inline=2, fontsize=10)
  ax1.set_title("Meridional Cut (longitud. 0$^o$-180$^o$), Time [kyr] ={:.2f}".format(frame*1e0), fontsize = 10)
  ax1.set_xlabel(r"$X$ (km)", fontsize = 10)
  ax1.set_ylabel(r"$Z$ (km)", fontsize = 10)

  ax2.cla()
  ax2.contourf(X1, Y1, func2[:,:,frame], levels=mylevels2, alpha=1, cmap=plt.cm.jet)
  C = ax2.contour(X1, Y1, func2[:,:,frame], mylevels2[0:nlevels:8] , colors='black')
  ax2.clabel(C, inline=2, fontsize=10)
  ax2.set_title("Meridional Cut (longitud. 90$^o$-270$^o$), Time [kyr] ={:.2f}".format(frame*1e0), fontsize = 10)
  ax2.set_xlabel(r"$X$ (km)", fontsize = 10)
  ax2.set_ylabel(r"$Z$ (km)", fontsize = 10)
  
  ax3.cla()
  ax3.contourf(X3, Y3, func3[:,:,frame], levels=mylevels3, alpha=1, cmap=plt.cm.jet)
  C = ax3.contour(X3, Y3, func3[:,:,frame], mylevels3[0:nlevels:8] , colors='black')
  ax3.clabel(C, inline=2, fontsize=10)
  ax3.set_title("Equatorial Cut, Time [kyr] ={:.2f}".format(frame*1e0), fontsize = 10)
  ax3.set_xlabel(r"$X$ (km)", fontsize = 10)
  ax3.set_ylabel(r"$Y$ (km)", fontsize = 10)

  ###################
  
  # Colorbar
# Manual input of MIN and MAX
min_b1 = np.min(func1)
min_b2 = np.min(func2)
min_b3 = np.min(func3)
min_b = min(min_b1, min_b2, min_b3)
cbar_min = np.min(min_b)

max_b1 = np.max(func1)
max_b2 = np.max(func2)
max_b3 = np.max(func3)
max_b = max(max_b1, max_b2, max_b3)
cbar_max = np.max(max_b)

# Colormap, Logaritmic normalization and ScalarMappable
cbar_cmap = matplotlib.cm.jet
cbar_norm = matplotlib.colors.TwoSlopeNorm(vmin=cbar_min, vcenter=0.0, vmax=cbar_max)
mappable = matplotlib.cm.ScalarMappable(norm=cbar_norm, cmap=cbar_cmap)

# Creation of new Axes for the Colorbar
cbar_ax = fig.add_axes([0.911, 0.15, 0.025, 0.7])  # (Dimensions [left, bottom, width, height] of the colorbar Axes)

# Creation of the Colorbar
cbar = fig.colorbar(mappable, cax=cbar_ax)
cbar.set_label('b-field [G]', size=12) 

anim = FuncAnimation(fig, animate, 
            frames=ntimes, interval=1, blit=False, repeat=False)

##################

# set the spacing between subplots
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.4)

plt.show()
