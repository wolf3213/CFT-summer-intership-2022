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
name_of_file='M1.0-a0.9-dumps'



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

step = 128

t,nstep,rho,ug,vu,B,pg,divb,uu,ud,bu,bd,bsq,ktot,rhor,beta,Lambda_sim,P_dump,T_dump,Hd,Qnu,Tau,Yee = hrms.csfn_dumpread("dump%03d"%step,gam,nx,ny,nz,float_type,name_of_file)  

##############################
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
rho = rho*RHO_UNIT
Lambda_sim = Lambda_sim*SIM_UNIT
Qnu = Qnu*SIM_UNIT        
#######################        
gd=hrms.gdet
print(nx) #so ive checked and this are dimensions of grid
print(ny)
print(nz)
print(np.size(ud))
print(np.size(ud,0))
print(np.size(ud,1))
print(np.size(ud,2))
#######################
#preparing data
r=np.array(r)
h=np.array(h)
rho=np.reshape(rho,(nx,ny)) #I dont understand why but without transporing doesnt work
rho=np.transpose(rho)
Yee=np.reshape(Yee,(nx,ny)) #I dont understand why but without transporing doesnt work
Yee=np.transpose(Yee)
bsq=np.reshape(bsq,(nx,ny)) #I dont understand why but without transporing doesnt work
bsq=np.transpose(bsq)
Tau=np.reshape(Tau,(nx,ny)) #I dont understand why but without transporing doesnt work
Tau=np.transpose(Tau)
h=np.reshape(h,(nx,ny)) #I dont understand why but without transporing doesnt work
h=np.transpose(h)
gd=np.reshape(gd,(nx,ny))
gd=np.transpose(gd)
ug=np.reshape(ug,(nx,ny))
#ud=np.reshape(ud,(4,nx,ny))
ug=np.reshape(ug,(nx,ny))
pg=np.reshape(pg,(nx,ny))

r0=[a[0,0] for a in r]
h0=[a[0] for a in h]
r0=np.array(r0)
h0=np.array(h0)

r0=r0
#print(h0)
#print(r0)


radius_matrix, theta_matrix = np.meshgrid(r0,h0)
#print(np.meshgrid(r0,h0)) #this makes grid in and radius and angle

X = radius_matrix * np.cos(theta_matrix) #this converts to cartesian cordinates
Y = radius_matrix * np.sin(theta_matrix)
fig = plt.figure()

#plt.plot(X,Y,marker="o",markersize='1')
#plt.show()
#############

#sizes of pictures
x_min=-300 #default are 25,50 and 100
x_max=300
y_min=0
y_max=600 #default are 50,100,200

##############
mass=0
i=0#radial cordinate
j=0#theta cordinate
for i in  np.arange(nx):
	for j in np.arange(ny):
		mass=mass+rho[j,i]*L_UNIT**3*_dx1*_dx2*gd[j,i] #mass=mass+rhov*L_UNIT**3*dr*dh*gd[i,j]
#rho and gd  here is 300 by 384
mass=mass*2*math.pi
masskg=mass/1000
masssolar=mass/MSUN
print('{} kg'.format(masskg))
print('{} MSUN'.format(masssolar))
##############
print("my code:")
mass=0
i=0#radial cordinate
j=0#theta cordinate
for i in  np.arange(nx):
	for j in np.arange(ny):
		if -ud[0][i,j,0]>1:
			mass=mass+rho[j,i]*L_UNIT**3*_dx1*_dx2*gd[j,i] #mass=mass+rhov*L_UNIT**3*dr*dh*gd[i,j]
mass=mass*2*math.pi
masskg=mass/1000
masssolar=mass/MSUN
print('{} kg'.format(masskg))
print('{} MSUN'.format(masssolar))
############
#Rho=np.transpose(rho)
Rho = rho[:,:]
Rho=np.transpose(Rho)
print(np.size(Rho,0))
#Gdet=np.transpose(gd)
Gdet = gd[:,:]
Gdet=np.transpose(Gdet)
print(np.size(Gdet,0))
Mask = np.zeros_like(Rho)
print(np.size(Mask,0))
print("simpson method")
i=0#radial cordinate
j=0#theta cordinate
for i in  np.arange(nx):
	for j in np.arange(ny):
		if -ud[0][i,j,0]>1:
			Mask[i][j]=1			



mass=0
f_integrand_r = Rho*Gdet*2*np.pi*Mask
vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)
        
f_integrand_theta = vol_r
vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)
 	
mass = vol_theta
mass = mass*L_UNIT*L_UNIT*L_UNIT
masskg=mass/1000
masssolar=mass/MSUN
print('{} kg'.format(masskg))
print('{} MSUN'.format(masssolar))


##############
#mass=0
#i=0#radial cordinate
#j=0#theta cordinate
#for a in rho:
#	for rhov in a:
#		if j==299:#why doesnt work with 300?
#			break
#		if ud[0,i,j]<-1 and r0[i]>80:
#			mass=mass+rhov*L_UNIT**3*_dx1*_dx2*gd[i,j] #mass=mass+rhov*L_UNIT**3*dr*dh*gd[i,j]
#		j=j+1
#	i=i+1
#	j=0
#	if i==384:	
#		break
#mass=mass*2*math.pi
#masskg=mass/1000
#masssolar=mass/MSUN
#print('{} kg'.format(masskg))
#print('{} MSUN'.format(masssolar))
############

mass=0
f_integrand_r = Rho*Gdet*2*np.pi
vol_r = scp_int.simps(f_integrand_r,dx=_dx1,axis=0)
        
f_integrand_theta = vol_r
vol_theta = scp_int.simps(f_integrand_theta,dx=_dx2,axis=0)
 	
mass = vol_theta
mass = mass*L_UNIT*L_UNIT*L_UNIT
masskg=mass/1000
masssolar=mass/MSUN
print('{} kg'.format(masskg))
print('{} MSUN'.format(masssolar))

##############

#rho_code=rho/RHO_UNIT
#mass=0
#i=0
#j=0
#for a in rho:
#	for rhov in a:
#		if j==299:#why doesnt work with 300?
#			break
#		dr=r0[i+1]-r0[i]
#		dh=h0[j+1]-h0[j]
#		if ud[0,i,j]*(1+ug[i,j]/rhov+pg[i,j]/rho_code[j,i])<-1 and r0[i]>80:
#			mass=mass+rhov*L_UNIT**3*_dx1*_dx2*gd[i,j] #mass=mass+rhov*L_UNIT**3*dr*dh*gd[i,j] should be there [j,i]????????/
#		j=j+1
#	i=i+1
#	j=0
#	if i==384:	
#		break
#mass=mass*2*math.pi
#masskg=mass/1000
#masssolar=mass/MSUN
#print('{} kg'.format(masskg))
#print('{} MSUN'.format(masssolar))


##############
#
#rho_code=rho/RHO_UNIT
#mass=0
#i=0
#j=0
#for a in rho:
#	for rhov in a:
#		if j==299:#why doesnt work with 300?
#			break
#		dr=r0[i+1]-r0[i]
#		dh=h0[j+1]-h0[j]
#		if ud[0,i,j]*(1+ug[i,j]/rhov+pg[i,j]/rho_code[j,i])*(0.9968 + 0.0085*Yee[j,i])<-1 and r0[i]>80:
#			mass=mass+rhov*L_UNIT**3*_dx1*_dx2*gd[i,j] #mass=mass+rhov*L_UNIT**3*dr*dh*gd[i,j] should be there [j,i]????????/
#		j=j+1
#	i=i+1
#	j=0
#	if i==384:	
#		break
#mass=mass*2*math.pi
#masskg=mass/1000
#masssolar=mass/MSUN
#print('{} kg'.format(masskg))
#print('{} MSUN'.format(masssolar))

##########################################################################################################################################################################


