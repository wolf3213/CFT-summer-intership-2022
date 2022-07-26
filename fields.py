import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import scipy.integrate as scp_int

import h5py as hdf
import h5_harm_script_parallel as hrms
import time
import os
import numpy as np
import sys
import os.path
import math
import cmath
import re


#author Adam Gonstal 09/07/2022
# read https://stackoverflow.com/questions/28926813/create-a-meshgrid-for-polar-coordinates
# read  https://stackoverflow.com/questions/57246146/logarithmic-colorbar


#################
name_of_file='early_dumps'



#################
class arr(): #this class is used to store information about the name
    def __init__(self, array, name):
        self.array = array
        self.name = name

if not os.path.exists('mems'):
    os.mkdir('mems')
#################
#importing data
nx,ny,nz,r,h,ph,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,x1,x2,x3,float_type=hrms.csfn_gread(name_of_file)
#r is radius values
#ph is phi values
#h is theta values
ratio=[]
ind_t=[]
ind_step=[]
mag_field=[]
energy_hot=[]
energy_thermal=[]
energy_cold=[]
#step=0
     
#######################
for step in np.arange(0,50,2):
	t,nstep,rho,ug,vu,B,pg,divb,uu,ud,bu,bd,bsq,ktot,rhor,beta,Lambda_sim,P_dump,T_dump,Hd,Qnu,Tau,Yee = hrms.csfn_dumpread("dump%03d"%step,gam,nx,ny,nz,float_type,name_of_file)      
	gd=hrms.gdet
########################
	# adjusting scales and units:
	GNEWT = 6.6742867e-8
	MSUN = 1.9892e33
	YEAR = 365.25 * 24. * 3600.
	KBOL = 1.3e-16
	CL = 2.99792458e10
	MP = 1.67262163783e-24
	ME = 9.10956e-28
	MN = 1.67492729e-24
	MHE = 6.6465e-24        #2p+2e+binding energy
	EE = 4.80320680e-10
	SIGMATH = 6.65248e-25
	SIGMASBOL = 5.67051e-5
	MAV = (MP + MN + 2. * ME + MHE) / 5.    # avrage particle mass
	# set mass of the black hole- sets the lenght scale, and the initial mass of the disk
	M_BLH = 1.0 * MSUN
	# fundamental units
	L_UNIT = GNEWT * M_BLH / (CL * CL)
	RHO_SCALE = 1.5e-5
	M_UNIT = RHO_SCALE * MSUN
	T_UNIT = L_UNIT / CL        
	USCALE = 1.0;
	# scaling ratio
	RHO_UNIT = M_UNIT / (L_UNIT * L_UNIT * L_UNIT)
	U_UNIT = RHO_UNIT * CL * CL * USCALE
	P_UNIT = U_UNIT
	SIM_UNIT = RHO_UNIT * CL * CL * CL / L_UNIT
	# converting units:
	#rho = rho*RHO_UNIT
	#Lambda_sim = Lambda_sim*SIM_UNIT
	#Qnu = Qnu*SIM_UNIT   
#######################
	#preparing data
	r=np.array(r)
	rho=np.reshape(rho,(nx,ny)) 
	bsq=np.reshape(bsq,(nx,ny)) 
	gd=np.reshape(gd,(nx,ny))
	ug=np.reshape(ug,(nx,ny))
	pg=np.reshape(pg,(nx,ny))
	e_hot=np.divide(ug,rho)
	k=0.0001
	e_cold_up=k*np.power(rho,gam)#this are arrays
	e_cold_down=(gam-1)*rho
	e_cold=np.divide(e_cold_up,e_cold_down)
	uth=rho*(e_hot-e_cold)
	ub=1/2*bsq
	ub_by_uth=np.divide(ub,uth) #result
	Gdet = gd[:,:] #this are slices of arrays
	Rho=rho[:,:]
	Ub=ub[:,:]
	Uth=uth[:,:]
	E_hot=e_hot[:,:]
	E_cold=e_cold[:,:]
	Ub_by_uth=ub_by_uth[:,:]
########## #averaging
	f_integrand_r = Rho*Gdet*2*np.pi
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	volume=vol_theta
##########
	f_integrand_r = Rho*Gdet*2*np.pi*E_cold
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	e_avr_cold=vol_theta/volume
##########
	f_integrand_r = Rho*Gdet*2*np.pi*E_hot
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	e_avr_hot=vol_theta/volume		
##########
	f_integrand_r = Ub_by_uth*Gdet*2*np.pi*Rho
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	 	
	avr_up = vol_theta #avaraged ratio
##########
	f_integrand_r = Rho*Gdet*2*np.pi*Ub
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	ub_avr=vol_theta/volume
##########
	f_integrand_r = Rho*Gdet*2*np.pi*Uth
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	uth_avr=vol_theta/volume	
##########
	result=avr_up/volume
	print(result)
	print(t)
	print(ub_avr)
	print(uth_avr)
	energy_cold.append(e_avr_cold) #this are avaraged arrays 
	energy_hot.append(e_avr_hot)
	mag_field.append(ub_avr)
	energy_thermal.append(uth_avr)
	ratio.append(result)
	ind_t.append(t)
	ind_step.append(step)
plt.plot(ind_step,energy_hot,'o')
plt.title('step vs hot energy')
plt.savefig('early_ step vs hot energy')
#plt.show()
plt.clf()
plt.plot(ind_step,energy_cold,'o')
plt.title('step vs cold energy')
plt.savefig('early_ step vs cold energy')
#plt.show()
plt.clf()
plt.plot(ind_step,energy_thermal,'o')
plt.title('step vs total thermal energy')
plt.savefig('early_step vs total thermal energy')
#plt.show()
plt.clf()
plt.plot(ind_step,mag_field,'o')
plt.title('step vs magfield energy')
plt.savefig('early_ step vs magfield')
#plt.show()
plt.clf()
#plt.plot(ind_t,ratio,'o')
#plt.show()
#plt.savefig('time vs energies ratio_early')
#plt.clf()
plt.plot(ind_step,ratio,'o')
plt.title('step vs ratio of ub by uth')
plt.savefig('early_ step vs energies ratio')
plt.clf
#plt.show()
