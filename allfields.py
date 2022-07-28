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
name_of_file='alldumps'



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
mag_field=[]#this is not specific 
energy_hot=[]
energy_thermal=[]
energy_cold=[]
u_thermal=[]#this is not specific 
mag_field_specific=[]
ratio_specific_array=[]
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
	e_th=e_hot-e_cold
	u_th=rho*e_th
	ub=1/2*bsq
	e_b=np.divide(ub,rho)
	ub_by_uth=np.divide(ub,u_th) 
	eb_by_eth=np.divide(e_b,e_th)
	Gdet = gd[:,:] #this are slices of arrays
	Rho=rho[:,:]
	Ub=ub[:,:]
	E_th=e_th[:,:]
	E_hot=e_hot[:,:]
	E_cold=e_cold[:,:]
	E_b=e_b[:,:]
	U_th=u_th[:,:]
	Ub_by_uth=ub_by_uth[:,:]
	Eb_by_eth=eb_by_eth[:,:]	
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
	f_integrand_r = Eb_by_eth*Gdet*2*np.pi*Rho
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	 	
	ratio_spefic_avr = vol_theta/volume #avaraged ratio specific
##########
	f_integrand_r = Rho*Gdet*2*np.pi*Ub
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	ub_avr=vol_theta/volume
##########
	f_integrand_r = Rho*Gdet*2*np.pi*E_th
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	E_th_avr=vol_theta/volume
##########
	#f_integrand_r = Rho*Gdet*2*np.pi*U_th
	#vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	#f_integrand_theta = vol_r
	#vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	#U_th_avr=vol_theta/volume
##########
	f_integrand_r = Rho*Gdet*2*np.pi*E_b
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	E_b_avr=vol_theta/volume		
##########
	result=avr_up/volume
	print(result)
	print(t)
	print(ub_avr)
	print(E_th_avr)
	#u_thermal.append(U_th_avr)
	mag_field_specific.append(E_b_avr)
	energy_cold.append(e_avr_cold) #this are avaraged arrays 
	energy_hot.append(e_avr_hot)
	mag_field.append(ub_avr)
	energy_thermal.append(E_th_avr)
	ratio.append(result)
	ind_t.append(t)
	ind_step.append(nstep)
	ratio_specific_array.append(ratio_spefic_avr)
#######################
for step in np.arange(60,410,10):
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
	e_th=e_hot-e_cold
	u_th=rho*e_th
	ub=1/2*bsq
	e_b=np.divide(ub,rho)
	ub_by_uth=np.divide(ub,u_th) 
	eb_by_eth=np.divide(e_b,e_th)
	Gdet = gd[:,:] #this are slices of arrays
	Rho=rho[:,:]
	Ub=ub[:,:]
	E_th=e_th[:,:]
	E_hot=e_hot[:,:]
	E_cold=e_cold[:,:]
	E_b=e_b[:,:]
	U_th=u_th[:,:]
	Ub_by_uth=ub_by_uth[:,:]
	Eb_by_eth=eb_by_eth[:,:]	
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
	f_integrand_r = Ub_by_uth*2*np.pi*Rho*Gdet
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	 	
	avr_up = vol_theta #avaraged ratio
##########
	f_integrand_r = Eb_by_eth*Gdet*2*np.pi*Rho
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	 	
	ratio_spefic_avr = vol_theta/volume #avaraged ratio specific
##########
	f_integrand_r = Rho*Gdet*2*np.pi*Ub
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	ub_avr=vol_theta/volume
##########
	f_integrand_r = Rho*Gdet*2*np.pi*E_th
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	E_th_avr=vol_theta/volume
##########
	#f_integrand_r = Rho*Gdet*2*np.pi*U_th
	#vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	#f_integrand_theta = vol_r
	#vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	#U_th_avr=vol_theta/volume
##########
	f_integrand_r = Rho*Gdet*2*np.pi*E_b
	vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)	        
	f_integrand_theta = vol_r
	vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)	
	E_b_avr=vol_theta/volume		
##########
	result=avr_up/volume
	print(result)
	print(t)
	print(ub_avr)
	print(E_th_avr)
	#u_thermal.append(U_th_avr)
	mag_field_specific.append(E_b_avr)
	energy_cold.append(e_avr_cold) #this are avaraged arrays 
	energy_hot.append(e_avr_hot)
	#mag_field.append(ub_avr)
	energy_thermal.append(E_th_avr)
	ratio.append(result)
	ind_t.append(t)
	ind_step.append(nstep)
	ratio_specific_array.append(ratio_spefic_avr)
ratio_after_avarage=np.divide(mag_field_specific,energy_thermal)
print(ind_t)
plt.figure(figsize=(12,8))
plt.plot(ind_t,energy_thermal,'b.')
plt.plot(ind_t,energy_thermal,'b')
plt.plot(ind_t,mag_field_specific,'r.')
plt.plot(ind_t,mag_field_specific,'r')
plt.title('time vs specific mag field && thermal energy')
plt.grid(True)
plt.xlabel('time')
plt.ylabel('total mag field && thermal energy')
plt.savefig('all time vs specific mag field and thermal energy')
#plt.show()
plt.clf()
#plt.plot(ind_t,u_thermal,'b.')
#plt.plot(ind_t,u_thermal,'b')
#plt.plot(ind_t,mag_field,'r.')
#plt.plot(ind_t,mag_field,'r')
#plt.title('time vs mag field && thermal energy')
#plt.xlabel('time')
#plt.ylabel('total mag field && thermal energy')
#plt.savefig('all time vs mag field and thermal energy')
#plt.show()
#plt.clf()
plt.plot(ind_t,energy_hot,'b.')
plt.plot(ind_t,energy_hot,'b')
plt.title('time vs hot energy')
plt.xlabel('time')
plt.ylabel('total internal energy')
plt.grid(True)
plt.savefig('all time vs hot energy')
#plt.show()
plt.clf()
##################
plt.plot(ind_t,energy_cold,'b.')
plt.plot(ind_t,energy_cold,'b')
plt.title('time vs cold energy')
plt.xlabel('time')
plt.ylabel('cold energy')
plt.grid(True)
plt.savefig('all time vs cold energy')
#plt.show()
plt.clf()
################
plt.plot(ind_t,np.log10(energy_thermal),'b.')
plt.plot(ind_t,np.log10(energy_cold),'g.')
plt.plot(ind_t,np.log10(energy_hot),'r.')
plt.plot(ind_t,np.log10(energy_thermal),'b',label='thermal energy')
plt.plot(ind_t,np.log10(energy_cold),'g',label='internal energy')
plt.plot(ind_t,np.log10(energy_hot),'r',label='cold energy')
plt.grid(True)
plt.legend()
plt.xlabel('time')
plt.ylabel('log10 of energy')
plt.title('time vs log of multiple energies wo mag')
plt.savefig('all time vs log_energies wo mag')
#plt.show()
plt.clf()
################
plt.plot(ind_t,np.log10(energy_thermal),'b.')
plt.plot(ind_t,np.log10(energy_cold),'g.')
plt.plot(ind_t,np.log10(energy_hot),'r.')
plt.plot(ind_t,np.log10(mag_field_specific),'k.')
plt.plot(ind_t,np.log10(energy_thermal),'b',label='thermal energy')
plt.plot(ind_t,np.log10(energy_cold),'g',label='cold energy')
plt.plot(ind_t,np.log10(energy_hot),'r',label='internal energy')
plt.plot(ind_t,np.log10(mag_field_specific),'k',label='magnetic energy')
plt.grid(True)
plt.legend()
plt.xlabel('time')
plt.ylabel('log10 of energy')
plt.title('time vs log of multiple energies')
plt.savefig('all time vs log_energies')
#plt.show()
plt.clf()
################
plt.plot(ind_t,energy_thermal,'b.')
plt.plot(ind_t,energy_cold,'g.')
plt.plot(ind_t,energy_hot,'r.')
plt.plot(ind_t,energy_thermal,'b',label='thermal energy')
plt.plot(ind_t,energy_cold,'g',label='cold energy')
plt.plot(ind_t,energy_hot,'r',label='internal energy')
plt.legend()
plt.grid(True)
plt.xlabel('time')
plt.ylabel('energy')
plt.title('time vs multiple energies energy')
plt.savefig('all time vs energies')
#plt.show()
plt.clf()
################
plt.plot(ind_t,energy_thermal,'b.')
plt.plot(ind_t,energy_thermal,'b')
plt.xlabel('time')
plt.ylabel('total thermal energy')
plt.title('time vs total thermal energy')
plt.grid(True)
plt.savefig('all time vs total thermal energy')
#plt.show()
plt.clf()
###################3
plt.plot(ind_t,np.log10(mag_field_specific),'b.')
plt.plot(ind_t,np.log10(mag_field_specific),'b')
plt.xlabel('time')
plt.ylabel('total log of specfic magnetic field energy')
plt.title('time vs log of spefic magfield energy')
plt.grid(True)
plt.savefig('all time vs log of magfield')
#plt.show()
plt.clf()
###################3
plt.plot(ind_t,mag_field_specific,'b.')
plt.plot(ind_t,mag_field_specific,'b')
plt.xlabel('time')
plt.ylabel('total specfic magnetic field energy')
plt.title('time vs specfic magfield energy')
plt.grid(True)
plt.savefig('all time vs magfield')
#plt.show()
plt.clf()
#plt.plot(ind_t,ratio,".")
#plt.show()
#plt.savefig('time vs energies ratio_early')
#plt.clf()

############
plt.plot(ind_t,ratio,'b.')
plt.plot(ind_t,ratio,'b')
plt.xlabel('time')
plt.ylabel('ub by uth ratio')
plt.title('time vs ratio of ub by uth')
plt.grid(True)
plt.savefig('all time vs energies ratio')
plt.clf()
############
plt.plot(ind_t,ratio_after_avarage,'b.')
plt.plot(ind_t,ratio_after_avarage,'b')
plt.xlabel('time')
plt.ylabel('<ub> by <uth> ratio')
plt.title('time vs ratio of <ub> by <uth> ')
plt.grid(True)
plt.savefig('all time vs  ratio of averaged energies')
plt.clf()
############
plt.plot(ind_t,np.log10(ratio_after_avarage),'b.')
plt.plot(ind_t,np.log10(ratio_after_avarage),'b')
plt.xlabel('time')
plt.ylabel('<ub> by <uth> ratio')
plt.title('time vs ratio of log( <ub> by <uth>) ')
plt.grid(True)
plt.savefig('all time vs  ratio of log averaged energies')
plt.clf()
################
plt.plot(ind_t,ratio_specific_array,'b.')
plt.plot(ind_t,ratio_specific_array,'b')
plt.xlabel('time')
plt.ylabel('eb by eth ratio')
plt.title('time vs ratio of eb by eth (specific energies ratio)')
plt.grid(True)
plt.savefig('all time vs specfic energies ratio')
plt.clf()
plt.close()
#plt.show()
