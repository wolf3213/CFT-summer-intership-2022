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

#################
class arr(): #this class is used to store information about the name
    def __init__(self, array, name):
        self.array = array
        self.name = name

if not os.path.exists('mems'):
    os.mkdir('mems')
#################

#author Adam Gonstal 09/07/2022
# read https://stackoverflow.com/questions/28926813/create-a-meshgrid-for-polar-coordinates
# read  https://stackoverflow.com/questions/57246146/logarithmic-colorbar
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
print(h0)
print(r0)




radius_matrix, theta_matrix = np.meshgrid(r0,h0)
print(np.meshgrid(r0,h0)) #this makes grid in and radius and angle

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

x_min=-25
x_max=25
y_min=0
y_max=50
plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())

plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))

#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()

###################
Z=arr(Yee, "Yee")

x_min=-25
x_max=25
y_min=0
y_max=50
plt.pcolormesh(X,Y,Z.array)
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()

fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
###################
Z=arr(bsq, "bsq")

x_min=-25
x_max=25
y_min=0
y_max=50
plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()

#################
Z=arr(Tau, "Tau")

x_min=-25
x_max=25
y_min=0
y_max=50
plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
##############
mass=0

for r in rho:
	for i in r:
		mass=mass+i
mass=mass/1000
print('{} kg'.format(mass))

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
print(h0)
print(r0)

x_min=-200
x_max=200
y_min=0
y_max=100


radius_matrix, theta_matrix = np.meshgrid(r0,h0)
print(np.meshgrid(r0,h0)) #this makes grid in and radius and angle

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

plt.pcolormesh(X,Y,Z.array,norm=matplotlib.colors.LogNorm())
plt.colorbar()
plt.axis([x_min, x_max, y_min, y_max]) #changes axis
plt.title('plot of {} for step {}'.format(Z.name, step))
#plt.show()
fig.savefig('part_of_grid_{}_{}.png'.format(Z.name, step))
plt.clf()
##############
mass=0

for r in rho:
	for i in r:
		mass=mass+i
mass=mass/1000
print('{} kg'.format(mass))






