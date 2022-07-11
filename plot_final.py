import numpy as np

import matplotlib
import matplotlib.pyplot as plt

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
#Units: 
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
M_BLH = 2.0 * MSUN

# fundamental units
L_UNIT = GNEWT * M_BLH / (CL * CL)
RHO_SCALE = 1.5e-5

M_UNIT = RHO_SCALE * MSUN
T_UNIT = L_UNIT / CL
         
#################
class arr(): #this class is used to store information about the name
    def __init__(self, array, name):
        self.array = array
        self.name = name

if not os.path.exists('mems'):
    os.mkdir('mems')
#################

nx,ny,nz,r,h,ph,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,x1,x2,x3,float_type=hrms.csfn_gread('dumps')
#r is radius values
#ph is phi values
#h is theta values

step = 0

t,nstep,rho,ug,vu,B,pg,divb,uu,ud,bu,bd,bsq,ktot,rhor,beta,Lambda_sim,P_dump,T_dump,Hd,Qnu,Tau,Yee = hrms.csfn_dumpread("dump%03d"%step,gam,nx,ny,nz,float_type,"./dumps")  


#print(nx) so ive checked and this are dimensions of grid
#print(ny)
#print(nz)


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
#############

#sizes of pictures
x_min=-50 #default are 25,50 and 100
x_max=50
y_min=0
y_max=100 #default are 50,100,200

################

#Z=Tau
#plt.pcolormesh(X,Y,Z) #makes figure

#plt.show()
#plt.colorbar()
#fig.savefig('all_grid_{}_{}.png'.format('Tau', step))
#plt.clf()
##################
#Z=rho
Z=arr(rho, "rho")


plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())

plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))

fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()

###################
Z=arr(Yee, "Yee")


plt.pcolormesh(X,Y,Z.array)
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()

fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
###################
Z=arr(bsq, "bsq")


plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()

#################
Z=arr(Tau, "Tau")


plt.contour(X,Y,Z.array)
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
#################
ZA=arr(Tau, "Tau")
ZB=arr(rho, "rho")

plt.contour(X,Y,ZA.array)
plt.pcolormesh(X,Y,ZB.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} vs {} for step {}'.format(ZA.name,ZB.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}_{}.png'.format(ZA.name,ZB.name, step))
plt.clf()
##############
mass=0
i=0
j=0
for a in rho:
	for rhov in a:
		if j==299:#why doesnt work with 300?
			break
		dr=r0[i+1]-r0[i]
		dh=h0[j+1]-h0[j]
		mass=mass+rhov*r0[i]**2*L_UNIT**3*math.sin(h0[j])*dr*dh
		j=j+1
	i=i+1
	j=0
	if i==384:	
		break
mass=mass*2*math.pi
masskg=mass/1000
masssolar=mass/MSUN
print('{} kg'.format(masskg))
print('{} MSUN'.format(masssolar))
##########################################################################################################################################################################

nx,ny,nz,r,h,ph,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,x1,x2,x3,float_type=hrms.csfn_gread('dumps')
step = 64

t,nstep,rho,ug,vu,B,pg,divb,uu,ud,bu,bd,bsq,ktot,rhor,beta,Lambda_sim,P_dump,T_dump,Hd,Qnu,Tau,Yee = hrms.csfn_dumpread("dump%03d"%step,gam,nx,ny,nz,float_type,"./dumps")  


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
#Z=Tau
#plt.pcolormesh(X,Y,Z) #makes figure

#plt.show()
#plt.colorbar()
#fig.savefig('all_grid_{}_{}.png'.format('Tau', step))
#plt.clf()
##################
#Z=rho
Z=arr(rho, "rho")


plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())

plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))

#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()

###################
Z=arr(Yee, "Yee")


plt.pcolormesh(X,Y,Z.array)
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()

fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
###################
Z=arr(bsq, "bsq")

plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()

#################
Z=arr(Tau, "Tau")

plt.contour(X,Y,Z.array)
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
#################
ZA=arr(Tau, "Tau")
ZB=arr(rho, "rho")

plt.contour(X,Y,ZA.array)
plt.pcolormesh(X,Y,ZB.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} vs {} for step {}'.format(ZA.name,ZB.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}_{}.png'.format(ZA.name,ZB.name, step))
plt.clf()
###############
mass=0
i=0
j=0
for a in rho:
	for rhov in a:
		if j==299:#why doesnt work with 300?
			break
		dr=r0[i+1]-r0[i]
		dh=h0[j+1]-h0[j]
		mass=mass+rhov*r0[i]**2*L_UNIT**3*math.sin(h0[j])*dr*dh
		j=j+1
	i=i+1
	j=0
	if i==384:	
		break
mass=mass*2*math.pi
masskg=mass/1000
masssolar=mass/MSUN
print('{} kg'.format(masskg))
print('{} MSUN'.format(masssolar))






