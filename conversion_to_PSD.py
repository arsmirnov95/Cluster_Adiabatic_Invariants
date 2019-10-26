import numpy as np
import pandas as pd
import scipy.io as io
from scipy import optimize
import matplotlib.pyplot as plt


""" Assign some variables """

# pitch angles-related:
pitch_angles = range(10,180,20) #PA's from 10 to 170
alphas = np.deg2rad(pitch_angles)
pitch_angles5 = range(10,100,20) # PA's from 10 to 90
pitch_angles

# energies-related.
# effective energies on available channels 1 and 3:
e13    = np.array([44.4927, 80.2213]) # energies of channels 1 and 3
e13_lg = np.log10(np.array([44.4927, 80.2213])) # energies of channels 1 and 3
# target energies in keV:
e_t_keV    = np.arange(40,160,10) # target energies in keV
e_t_keV_lg = np.log10(np.arange(40,160,10)) # lg of target energies in keV
# target energies in MeV
e_t_MeV    = np.arange(40,160,10)/1000. #40,160,10)/1000. # target energies in MeV
e_t_MeV_lg = np.log10(np.arange(40,160,10)/1000.) #np.arange(40,160,10)/1000.) # lg of target energies in MeV

""" Processing functions """

def alpha_k_fit(K_0, K_eachtime):
    alpha_k = np.interp(np.log10(K_0), np.log10(K[i,::-1]), pitch_angles5[::-1], right=np.nan, left=np.nan)
    return alpha_k

def retrieve_flux_from_PAD(p, alpha):
    fitted_j = np.log10(p[0]) + np.log10(np.sin(alpha)) + p[1]*np.log10(np.sin(alpha))
    return 10**fitted_j

def energy_fit_flux(fluxes_1_and_3):
    # equation of straight line retrieved by 2 points:
    f1 = np.log10(fluxes_1_and_3[0])
    f3 = np.log10(fluxes_1_and_3[1])
    fitted_fluxes_lg = f1 + (e_t_keV_lg - e13_lg[0])*(f3-f1)/(e13_lg[1] - e13_lg[0])
    return fitted_fluxes_lg
    
def PSD_target_lg(fitted_fluxes_lg, magnetic_field, alpha_K, mu_0):
    flux = 10**(fitted_fluxes_lg)
    
    PSD = flux/(e_t_MeV*(e_t_MeV+2*0.511))*1.66*10**(-10)*200.3
    
    mu  = e_t_MeV*(e_t_MeV+2*0.511)*(np.sin(alpha_K))**2*10**5/(magnetic_field*2*0.511)
    
    if (mu_0<np.nanmax(mu)) & (mu_0>np.nanmin(mu)):
        PSD_0 = np.interp(np.log10(mu_0),np.log10(mu), np.log10(PSD), right=np.nan, left=np.nan)
    else:
        PSD_0 = np.nan
    return PSD_0
    
""" Main function for PSD Conversion """

def convert_flux_to_PSD(K_0, mu_0, K, PAD_params1, PAD_params3, corr1, corr3, B, L, t)

    PSD = np.zeros(len(L))*np.nan
    fluxes_ = np.zeros((len(L),2))*np.nan
    pitch_angles = np.zeros(len(L))*np.nan
    
    for i in range(len(L)):
    
        if (K_0>np.nanmax(K[i,:])) or (K_0<np.nanmin(K[i,:])) or (corr1[i]<0.8) or (corr3[i]<0.8) or (PAD_params1[i,0]==0.1) or (PAD_params1[i,1]==0.1):
    
            PSD[i] = np.nan
    
        else:
    
            alpha_k = np.interp(np.log10(K_0), np.log10(K[i,::-1]), pitch_angles5[::-1], right=np.nan, left=np.nan)
            pitch_angles[i] = alpha_k
            alpha_K_rad = np.deg2rad(alpha_k)
    
            flux_chan1 = retrieve_flux_from_PAD(PAD_params1[i,:], alpha_K_rad) # [cm-2sr-1keV-1]
            flux_chan3 = retrieve_flux_from_PAD(PAD_params3[i,:], alpha_K_rad) # [cm-2sr-1keV-1]
            fluxes_[i,:]=[flux_chan1, flux_chan3]
            fitted_fluxes_lg = energy_fit_flux([flux_chan1, flux_chan3])
    
            PSD[i] = PSD_target_lg(fitted_fluxes_lg, B[i], alpha_K_rad, mu_0)
            
    PSD[np.log10(fluxes_[:,0])<0]=np.nan
    PSD[np.log10(fluxes_[:,1])<0]=np.nan
    df = pd.DataFrame({'time':t, 'Lstar':L, 'PSD':PSD})
    df.head()
    
    df=df.dropna()
    return df