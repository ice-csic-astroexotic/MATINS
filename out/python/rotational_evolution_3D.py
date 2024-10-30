######################################################## 
# Neutron stars period and inclination angle evolution #
########################################################

#import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import csv
from typing import Tuple

#import the rcParams subclass, that contains matplotlib core settings
from matplotlib import rcParams

rcParams['mathtext.fontset'] = 'stix'
rcParams['font.family'] = 'Liberation serif'
rcParams['font.size']='20'
rcParams['figure.figsize']='10.0, 10.0'
rcParams['figure.autolayout']=False

rcParams['axes.linewidth'] = '1.7' 

rcParams['xtick.direction'] = 'in'
rcParams['xtick.top'] = True
rcParams['xtick.major.size'] = '10.0'
rcParams['xtick.major.width'] = '1.7'
rcParams['xtick.minor.size'] = '5.0'
rcParams['xtick.minor.width'] = '1.7'
rcParams['xtick.labelsize'] = '20'

rcParams['ytick.direction'] = 'in'
rcParams['ytick.right'] = True
rcParams['ytick.major.size'] = '10.0'
rcParams['ytick.major.width'] = '1.7'
rcParams['ytick.minor.size'] = '5.0'
rcParams['ytick.minor.width'] = '1.7'
rcParams['ytick.labelsize'] = '20'

from matplotlib.ticker import FuncFormatter
formatter = FuncFormatter(lambda y, _: '{:.16g}'.format(y))

###########################################################################################

# Unit conversion
YR_TO_S = 3600 * 24 * 365  # Convert from [yr] to [s].

# Constants
Msun = 2.e33    # Solar mass [g]
c = 2.99792458e10   # Speed of light [cm/s]

# Neutron stars physical parameters
NS_mass = 1.4*Msun  # Neutron star mass [g]
NS_radius = 1.1e6 # Neutron star radius [cm] 

# Canonical neutron star moment of inertia in [g cm^2] assuming a perfect solid sphere.
NS_inertia = 1E+45 
# Auxiliary quantity beta as defined in eq. (72) of Pons & Vigano (2019).
beta = np.pi ** 2 * NS_radius ** 6 / (NS_inertia * c ** 3)
beta_vac = 6.489238670467112e-40 
print(beta)
###########################################################################################
# Functions

def period_derivative(B: float, chi: float, P: float, k_0: float, k_1: float) -> float:
    """
    This function determines the change in the rotation period of a pulsar. It is taken
    from eq. (70) of Pons & Vigano (2019). For more details see, e.g., Spitkovsky (2006)
    or Philippov et al. (2014), who determine the coefficients k_0, k_1, k_2 for a pulsar 
    embedded in a force-free and resistive magnetosphere from numerical simulations. 
    k0 ~ k1 ~ 1 for plasma filled magnetosphere solutions. 
    k0 = 0 and k1 = 2/3 for vacuum solution.

    Args:
        B (float): value of the dipolar component of the magnetic field at the
        magnetic pole for a simulated neutron star, measured in [G].
        chi (float): angle between the magnetic dipolar moment, i.e., the magnetic
        field axis, and the rotation axis for a simulated pulsar, measured in [rad].
        P (float): spin period of a simulated pulsar, measured in [s].
        k_0 (float): coefficient regulating the behaviour of the spin down evolution.
        k_1 (float): coefficient regulating the behaviour of the spin down evolution.

    Returns:
        (float): period derivative of a simulated pulsar in [s/yr].
    """

    # Period derivative.
    P_deriv = (
        beta
        * B ** 2
        / P
        * (k_0 + k_1 * np.sin(chi) ** 2)
    ) * YR_TO_S
    
    return P_deriv

def period_derivative_vac(B: float, P: float) -> float:
    """
    This function determines the change in the rotation period of a pulsar. It is taken
    from eq. (70) of Pons & Vigano (2019), considering k0 = 1 and k1 = 0 for vacuum 
    solution, beta = 1/(6.4e10)**2 and we don't take into account 
    inclination angle chi, chi=0. Therefore Pdot = B**2 * beta_vac / P. 

    Args:
        B (float): value of the dipolar component of the magnetic field at the
        magnetic pole for a simulated neutron star, measured in [G].
        P (float): spin period of a simulated pulsar, measured in [s].

    Returns:
        (float): period derivative of a simulated pulsar in [s/yr].
    """

    # Period derivative.
    P_deriv = (
        beta_vac
        * B ** 2
        / P) * YR_TO_S
    
    return P_deriv


def inclination_angle_derivative(B: float, chi: float, P: float, k_2: float) -> float:
    """
    This function determines the change in the misalignment angle, i.e., the angle between the
    magnetic dipolar moment and the rotation axis of a pulsar. It is taken from eq. (71) of
    Pons & Vigano (2019). For more details see, e.g., Spitkovsky (2006) or Philippov et al.
    (2014), who determine the coefficients k_0, k_1, k_2 (defined in the configuration file)
    for a pulsar embedded in a force-free and resistive magnetosphere from numerical simulations.
    k2 ~ 1 for plasma filled magnetosphere solutions;
    k2 = 2/3 for vacuum solution.
    
    Args:
        B (float): values of the dipolar component of the magnetic field at the
        magnetic pole for the sample of simulated neutron stars, measured in [G].
        chi (float): angles between the magnetic dipolar moment, i.e., the magnetic
        field axis, and the rotation axis for all simulated pulsars, measured in [rad].
        P (float): spin periods of simulated pulsars, measured in [s].
        k_2 (float): coefficient regulating the behaviour of the inclination angle evolution.

    Returns:
        (float): misalignment angle derivatives for all simulated pulsars in [rad/yr].
    """

    # Misalignment angle derivative.
    chi_deriv = (
        -k_2
        * beta
        * B ** 2
        / (P ** 2)
        * np.sin(chi)
        * np.cos(chi)
    ) * YR_TO_S

    return chi_deriv

    
def Pdot_tage(P:np.ndarray, tau_age:float) -> np.ndarray:
    '''
    Pdot as a function of period for varying characteristic age.
    Args:
        P (np.ndarray): Array of period values in [s].
        tau_age (float): characteristic age in [yr].
    Return:
        (np.ndarray): Values of period derivatives in [s/s].
    '''
    tau_age = tau_age * YR_TO_S
    Pdot = P / (2. * tau_age)
    
    return Pdot
    
    
def rotational_evolution(t_mt: np.ndarray, B_mt: np.ndarray, L_mt: np.ndarray, P_in: float, chi_in: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    '''
    Magneto-rotational evolution of a neutron star for two types of solutions:
    1) plasma filled magnetosphere with inclination angle evolution;
    2) vacuum magnetosphere with no inclination angle evolution.
    This function takes as input the time evolution of the dipolar magnetic field from the 
    magneto-thermal evolution code and starting from some initial conditions it evolves in time
    the spin period P and the inclination angle chi for the two soultions 1) and 2).
    
    Args:
        t_mt (np.ndarray): Time grid where the magneto-thermal evolution has been performed [yr].
        B_mt (np.ndarray): Dipolar magnetic field evolution from the magneto-thermal evolution output [G].
        L_mt (np.ndarray): Thermal-luminosity evolution from the magneto-thermal evolution output [G].
        P_in (float): Initial spin period [s].
        chi_in (float): Initial inclination angle [rad].
    Return:
        (np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray): 
        Tuple containg the new time grid used for the evolution, the magnetic field evolution interpolated 
        on the new time grid, the thermal luminosity evolution interpolated on the new time grid, the spin 
        period P evolution for solution 1) and 2), the spin period derivative P_dot evolution for solution 
        1) and 2), the inclination angle chi evolution for solution 1).
    '''
   
    # Define the new logarithmic time grid to use for the evolution.
    t_grid = np.logspace(-1, np.log10(max(t_mt)), 1000)
    
    # Interpolate the B field on the new time grid  
    B_mt_int = np.interp(t_grid, t_mt, B_mt)
    print(B_mt)
    
    # Interpolate the luminosity L on the new time grid  
    L_mt_int = np.interp(t_grid, t_mt, L_mt)
    
    # Initialize the arrays containing the evolution information for solution 1) and 2).
    P_1 = np.zeros(len(t_grid))
    Pdot_1 = np.zeros(len(t_grid))
    chi = np.zeros(len(t_grid))
    P_2 = np.zeros(len(t_grid))
    Pdot_2 = np.zeros(len(t_grid))
    P_3 = np.zeros(len(t_grid))
    Pdot_3 = np.zeros(len(t_grid))
    
    # set the initial conditions
    P_1[0] = P_in
    chi[0] = chi_in
    P_2[0] = P_in
    P_3[0] = P_in
    
    # Set the initial values for the derivatives for solution 1) and 2).
    Pdot_1[0] = period_derivative(B_mt_int[0], chi[0], P_1[0], k_0=1, k_1=1)   # period derivative filled magnetosphere solution
    chidot = inclination_angle_derivative(B_mt_int[0], chi[0], P_1[0], k_2=1)   # inclination angle derivative
    Pdot_2[0] = period_derivative(B_mt_int[0], chi[0], P_2[0], k_0=0, k_1=2./3.)   # period derivative vacuum solution
    Pdot_3[0] = period_derivative_vac(B_mt_int[0], P_3[0])   # period derivative vacuum solution (Observers Formula)


    with open("test_extended.csv", 'w') as f:
        writer = csv.writer(f, delimiter=',', lineterminator='\n',)
        header = ['t[yr]','B[G]','L[erg/s]','P_fill[s]','Pdot_fill[s/s]','chi[rad]','P_vac[s]','Pdot_vac[s/s]','P_vacObs[s]','Pdot_vacObs[s/s]']
        writer.writerow(header) 
        
        # Evolve all the physical quantities using finite differences method.
        for i in range(1,len(t_grid)):

            dt = t_grid[i]-t_grid[i-1]

            P_1[i] = P_1[i-1]+Pdot_1[i-1]*dt
            P_2[i] = P_2[i-1]+Pdot_2[i-1]*dt
            P_3[i] = P_3[i-1]+Pdot_3[i-1]*dt
            chi[i] = chi[i-1]+chidot*dt

            Pdot_1[i] = period_derivative(B_mt_int[i], chi[i], P_1[i], k_0=1, k_1=1)   
            chidot = inclination_angle_derivative(B_mt_int[i], chi[i], P_1[i], k_2=1) 
            Pdot_2[i] = period_derivative(B_mt_int[i], chi[0], P_2[i], k_0=0, k_1=2./3.)
            Pdot_3[i] = period_derivative_vac(B_mt_int[i], P_3[i])
            
            # write on the output.csv file
            row = [t_grid[i],B_mt_int[i],L_mt_int[i],P_1[i],Pdot_1[i] / YR_TO_S,chi[i],P_2[i],Pdot_2[i] / YR_TO_S,P_3[i],Pdot_3[i] / YR_TO_S] 
            writer.writerow(row) 
    
    Pdot_1 = Pdot_1 / YR_TO_S
    Pdot_2 = Pdot_2 / YR_TO_S
    Pdot_3 = Pdot_3 / YR_TO_S
  
    
    # Interpolate back on the time grid provided by the magneto-thermal code
    P_1_mt = np.interp(t_mt, t_grid, P_1)
    P_2_mt = np.interp(t_mt, t_grid, P_2)
    P_3_mt = np.interp(t_mt, t_grid, P_3)
    
    print(P_1)
        
    Pdot_1_mt = np.interp(t_mt, t_grid, Pdot_1)
    Pdot_2_mt = np.interp(t_mt, t_grid, Pdot_2)
    Pdot_3_mt = np.interp(t_mt, t_grid, Pdot_3)
    
    chi_mt = np.interp(t_mt, t_grid, chi)
    
    header = ['t[yr]','B[G]','L[erg/s]','P_fill[s]','Pdot_fill[s/s]','chi[rad]','P_vac[s]','Pdot_vac[s/s]','P_vacObs[s]','Pdot_vacObs[s/s]']

    df = pd.DataFrame(data=np.array([t_mt, B_mt, L_mt, P_1_mt, Pdot_1_mt, chi_mt, P_2_mt, Pdot_2_mt, P_3_mt, Pdot_3_mt]).T, columns=header)

   # adding column with constant value
    df['bpol'] = 3e14
    df['etor'] = 50
    df['btorAvg'] = 5.2e13

    # Save the data frame as .csv file.
    df.to_csv("cooling_curve.csv", index=False)
    return t_grid,B_mt_int,L_mt_int,P_1,P_2,Pdot_1,Pdot_2,chi
    

t_mt = np.loadtxt("cooling_curve.d", usecols=[1])
L_mt = np.loadtxt("cooling_curve.d", usecols=[2])
B_mt = np.loadtxt("cooling_curve.d", usecols=[7])

# Set inital period and inclination angle.
P_in = 1.e-3
chi_in = np.pi/3

# Evolve in time>
t,B,L,P_1,P_2,Pdot_1,Pdot_2,chi = rotational_evolution(t_mt, B_mt, L_mt, P_in, chi_in)


# P-PDOT plot
period_derivative_vect = np.vectorize(period_derivative)

# Setup the P-Pdot plot.
B_grid = 10**(1.*(np.arange(10,18,2)))    # array of dipolar magnetic field values [gauss]
tau_age_grid = 10**(1.*(np.arange(2,14,4)))    # array of characteristic age values [yr]
P_grid = np.logspace(np.log10(2.e-3), np.log10(100),10)    # periods array [s]

fig, ax = plt.subplots(figsize=(11,8))

# Plot the lines of constant magnetic field and characteristic age.
for i in range(len(B_grid)):
    ax.plot(
        P_grid, 
        period_derivative_vect(B=B_grid[i], chi=0, P=P_grid, k_0=1, k_1=1) / YR_TO_S, 
        linestyle='-', 
        color='black',
        linewidth=1,
        zorder=1,
        rasterized=True
    )
    
for i in range(len(tau_age_grid)):
    ax.plot(
        P_grid, 
        Pdot_tage(P_grid,tau_age_grid[i]), 
        linestyle='--',
        color='black',
        linewidth=1,
        zorder=1,
        rasterized=True
    )

ax.text(3.e-3, 8e-18, '$10^{10}$ G', fontsize = 20, rotation=-25, color = 'black')
ax.text(3.e-3, 8e-14, '$10^{12}$ G', fontsize = 20, rotation=-25, color = 'black')
ax.text(1.5e-2, 1.5e-10, '$10^{14}$ G', fontsize = 20, rotation=-25, color = 'black')

ax.text(3.e-3, 8e-13, '$10^{2}$ yr', fontsize = 20, rotation=25, color = 'black')
ax.text(3.e-3, 8e-17, '$10^{6}$ yr', fontsize = 20, rotation=25, color = 'black')
ax.text(3., 8e-18, '$10^{10}$ yr', fontsize = 20, rotation=25, color = 'black')

ax.set_xlabel(r"Spin period [s]")
ax.set_ylabel(r"Spin period derivative [s/s]")
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(2.e-3, 100)
ax.set_ylim(5.e-18, 1.e-9)
ax.xaxis.set_major_formatter(formatter)

ax.plot(P_1, Pdot_1, linestyle='-',linewidth=3,color='tab:blue', label='filled')
ax.plot(P_2, Pdot_2, linestyle='-',linewidth=3,color='tab:orange', label='vacuum')

plt.legend(frameon=False, loc=0)
plt.show(block=False)

