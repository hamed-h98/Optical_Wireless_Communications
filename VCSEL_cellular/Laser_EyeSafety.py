# Eye-safety calculation for VCSELs based on: 
# Zeng, Zhihong, et al. "A VCSEL array transmission system with novel beam activation mechanisms." IEEE Transactions on Communications 70.3 (2021): 1886-1900.

import numpy as np 
import matplotlib.pyplot as plt


class Eye_Safety: 
    def __init__(self,lambda_nm, FWHM,t_exp):
        self.lambda_nm = lambda_nm
        self.FWHM = FWHM
        self.t_exp = t_exp
    
    def Pt(self):
        theta_FWHM = np.deg2rad(self.FWHM)
        theta_beam = theta_FWHM / np.sqrt(2*np.log(2))
        w0 = self.lambda_nm / (np.pi * theta_beam)
        z_mhp = 0.1
        w_mhp = w0*np.sqrt(1 + ((self.lambda_nm*z_mhp)/(np.pi*(w0**2)))**2)
        if self.t_exp < 10: 
            if self.lambda_nm == 850e-9:
                C4 = 10**(0.002*(self.lambda_nm*1e9 - 700))
                MPE = 18*(self.t_exp**0.75)*C4/self.t_exp
                d_pupil = 7 / 1000 
            elif self.lambda_nm == 1550e-9:
                MPE = (1e4) / self.t_exp 
                d_pupil = 1.5*(self.t_exp**(3/8)) / 1000 
        elif self.t_exp >= 10: 
            if self.lambda_nm == 850e-9: 
                C4 = 10**(0.002*(self.lambda_nm*1e9 - 700))
                C7 = 1
                MPE = 10*C4*C7 
                d_pupil = 7 / 1000 
            elif self.lambda_nm == 1550e-9:
                MPE = 1000 
                d_pupil = 3.5 / 1000 
        
        Pt = (np.pi*(d_pupil**2)*MPE)/(4*(1 - np.exp(-(d_pupil**2)/(2*w_mhp**2))))
        return Pt

    def plot_eyesafety(self,t,Pt,idx):
        colors = [
            (1, 0, 0),    # Red
            (0, 1, 0),    # Green
            (0, 0, 1),    # Blue
            (0.5, 0, 0.5), # Purple
            (0, 0.75, 0.75), # Teal
            (1, 0.5, 0)  # Orange
        ]
        plt.grid(True)
        plt.xlabel('Exposure Time [s]')
        plt.ylabel('Maximum Transmit Power [mW]')
        plt.title(f"{self.lambda_nm*1e9} nm")
        plt.plot(t,Pt*1000,linewidth = 2,color = colors[idx])
        if self.lambda_nm == 850e-9: 
            plt.ylim([0, 15])
        elif self.lambda_nm == 1550e-9:
            plt.ylim([0,1000])


FWHM = [0.4,4,6]
t_exp = np.arange(0.01,30,0.01) 
Pt_max_values = np.zeros((len(FWHM),len(t_exp)))
Wavelength = [850*1e-9, 1550e-9]
for wavelength in Wavelength:
    plt.figure()
    for idx, fwhm in enumerate(FWHM): 
        Pt_max = []
        for t in t_exp: 
            eye_safety = Eye_Safety(wavelength, fwhm, t)
            Pt_max.append(eye_safety.Pt())
        Pt_max_values[idx, :] = np.array(Pt_max) 
        eye_safety.plot_eyesafety(t_exp, Pt_max_values[idx,:],idx)
    plt.legend([rf"$\theta_{{FWHM}}= {fwhm}^\circ$" for fwhm in FWHM], fontsize = 12) 

plt.show()