import sys
from scipy.interpolate import interp1d
from numpy import sin, cos, tan, pi

from scipy import linalg


import pdb

import numpy as np
import glob
import os
from scipy.interpolate import griddata
from scipy.interpolate import interp1d
from numpy import ma
import h5py as hdf
import gc
import resource
import matplotlib.pyplot as plt #added 

def convert_to_single_file(startn=0,endn=-1,ln=10,whichi=0,whichn=1,**kwargs):
    which = kwargs.pop("which","convert_file")
    rg("gdump")
    flist1 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9]_0000") ) )
    flist2 = np.sort(glob.glob( os.path.join("dumps/", "dump[0-9][0-9][0-9][0-9]_0000") ) )
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist1,flist2))
    firsttime = 1
    for fldname in flist:
        #find the index of the file
        fldindex = np.int(fldname.split("_")[0].split("p")[-1])
        if fldindex < startn:
            continue
        if endn>=0 and fldindex >= endn:
            break
        if fldindex % whichn != whichi:
            #do every whichn'th snapshot starting with whichi'th snapshot
            continue
        print( "Reading " + fldname + " ..." )
        fname = "dump%03d" % fldindex
        if os.path.isfile( fname ):
            print("File %s exists, skipping..." % fname)
            continue
        if not os.path.isfile( fname ):
            rd(fname)
        

def ellk(a,r):
    ekval = ek(a,r)
    lkval = lk(a,r)
    return(lkval/ekval)

def ek(a,r):
    #-u_t, I suspect
    ek = (r**2-2*r+a*r**0.5)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return(ek)

def lk(a,r):
    udphi = r**0.5*(r**2-2*a*r**0.5+a**2)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return( udphi )

def Risco(ain):
    eps = np.finfo(np.float64).eps
    a = np.minimum(ain,1.)
    Z1 = 1 + (1. - a**2)**(1./3.) * ((1. + a)**(1./3.) + (1. - a)**(1./3.))
    Z2 = (3*a**2 + Z1**2)**(1./2.)
    risco = 3 + Z2 - np.sign(a)* ( (3 - Z1)*(3 + Z1 + 2*Z2) )**(1./2.)
    return(risco)

def Ebind(r,a):
    #1+u_t, I suspect
    Eb = 1 - (r**2-2*r+a*r**0.5)/(r*(r**2-3*r+2*a*r**0.5)**0.5)
    return( Eb )

def etaNT(a):
    return( Ebindisco(a) )

def Ebindisco(a):
    eps = np.finfo(np.float64).eps
    a0 = 0.99999 #1.-1e8*eps
    if a > a0: 
        a = a0
        Eb = Ebind( Risco(a), a )
        return((a-a0)/(1.-a0)*(1.-3.**(-0.5)) + (1.-a)/(1.-a0)*Eb)
    Eb = Ebind( Risco(a), a)
    #Eb = (1.-3.**(-0.5))*a**2
    return( Eb )



def convert_wrapper(**kwargs):
    if len(sys.argv[2:])==2 and sys.argv[2].isdigit() and sys.argv[3].isdigit():
        whichi = int(sys.argv[2])
        whichn = int(sys.argv[3])
    else:
        print( "Usage: %s %s <whichi> <whichn>" % (sys.argv[0], sys.argv[1]) )
        return
    convert_to_single_file(whichi = whichi, whichn = whichn, **kwargs)


#############

def Qmri(dir=2):
    """
    APPROXIMATELY Computes number of theta cells resolving one MRI wavelength
    """
    global bu,rho,uu,_dx2,_dx3
    #cvel()
    #corrected this expression to include both 2pi and dxdxp[3][3]
    #also corrected defition of va^2 to contain bsq+gam*ug term
    #need to figure out how to properly measure this in fluid frame
    vaudir = np.abs(bu[dir])/np.sqrt(rho+bsq+gam*ug)
    omega = dxpdx[3][3]*uu[3]/uu[0]+1e-15
    lambdamriudir = 2*np.pi * vaudir / omega
    if dir == 2:
        res=lambdamriudir/_dx2
    elif dir == 3:
        res=lambdamriudir/_dx3
    return(res)

def OmegaXX(dir=2):
    """
    APPROXIMATELY Computes number of theta cells resolving one MRI wavelength
    """
    global bu,rho,uu,_dx2,_dx3,omega
    #cvel()
    #corrected this expression to include both 2pi and dxdxp[3][3]
    #also corrected defition of va^2 to contain bsq+gam*ug term
    #need to figure out how to properly measure this in fluid frame
    vaudir = np.abs(bu[dir])/np.sqrt(rho+bsq+gam*ug)
    omega = dxpdx[3][3]*uu[3]/uu[0]+1e-15
    lambdamriudir = 2*np.pi * vaudir / omega
    if dir == 2:
        res=lambdamriudir/_dx2
    elif dir == 3:
        res=lambdamriudir/_dx3    
    return(omega)    

def goodlabs(fntsize=20):
    ax = plt.gca()
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontsize(fntsize)



def iofr(rval):
    rval = np.array(rval)
    if np.max(rval) < r[0,0,0]:
        return 0
    res = interp1d(r[:,0,0], ti[:,0,0], kind='linear', bounds_error = False, fill_value = 0)(rval)
    if len(res.shape)>0 and len(res)>0:
        res[rval<r[0,0,0]]*=0
        res[rval>r[nx-1,0,0]]=res[rval>r[nx-1,0,0]]*0+nx-1
    else:
        res = np.float64(res)
    return(np.floor(res+0.5).astype(int))








def mdot(a,b):
    """
    Computes a contraction of two tensors/vectors.  Assumes
    the following structure: tensor[m,n,i,j,k] OR vector[m,i,j,k], 
    where i,j,k are spatial indices and m,n are variable indices. 
    """
    if (a.ndim == 3 and b.ndim == 3) or (a.ndim == 4 and b.ndim == 4):
          c = (a*b).sum(0)
    elif a.ndim == 5 and b.ndim == 4:
          c = np.empty(np.maximum(a[:,0,:,:,:].shape,b.shape),dtype=b.dtype)
          for i in range(a.shape[0]):
                c[i,:,:,:] = (a[i,:,:,:,:]*b).sum(0)
    elif a.ndim == 4 and b.ndim == 5:
          c = np.empty(np.maximum(b[0,:,:,:,:].shape,a.shape),dtype=a.dtype)
          for i in range(b.shape[1]):
                c[i,:,:,:] = (a*b[:,i,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 5:
          c = np.empty((a.shape[0],b.shape[1],a.shape[2],a.shape[3],max(a.shape[4],b.shape[4])),dtype=a.dtype)
          for i in range(c.shape[0]):
                for j in range(c.shape[1]):
                      c[i,j,:,:,:] = (a[i,:,:,:,:]*b[:,j,:,:,:]).sum(0)
    elif a.ndim == 5 and b.ndim == 6:
          c = np.empty((a.shape[0],b.shape[1],b.shape[2],max(a.shape[2],b.shape[3]),max(a.shape[3],b.shape[4]),max(a.shape[4],b.shape[5])),dtype=a.dtype)
          for mu in range(c.shape[0]):
              for k in range(c.shape[1]):
                  for l in range(c.shape[2]):
                      c[mu,k,l,:,:,:] = (a[mu,:,:,:,:]*b[:,k,l,:,:,:]).sum(0)
    else:
           raise Exception('mdot', 'wrong dimensions')
    return c


def psicalc(B1=None):
    """
    Computes the field vector potential
    """
    global B
    if B1 is None: B1 = B[1]
    daphi = -(gdet*B1).mean(-1)*_dx2
    aphi=daphi[:,::-1].cumsum(axis=1)[:,::-1]
    aphi-=0.5*daphi #correction for half-cell shift between face and center in theta
    return(aphi)

def B_pcalc(B1=None, B2=None):
    global B
    if B1 is None: B1 = B[1]
    if B2 is None: B2 = B[2]
    B_poloidal = np.sqrt(B1*B1+B2*B2)
    return(B_poloidal)

def B_phicalc(B3=None):
    global B
    if B3 is None: B3 = B[3]
    B_toroidal = B3
    return(B_toroidal)


def myfloat(f,acc=1):
    """ acc=1 means np.float32, acc=2 means np.float64 """
    if acc==1:
        return( np.float32(f) )
    else:
        return( np.float64(f) )

def get_fracphi():
    fracphi = dxdxp[3,3,0,0,0]*_dx3*nz/(2*np.pi)
    return( fracphi )

def plco(myvar,**kwargs):
    global r,h,ph
    plt.clf()
    return plc(myvar,**kwargs)

def plc(myvar,**kwargs): #plc
    global r,h,ph
    #xcoord = kwargs.pop('x1', None)
    #ycoord = kwargs.pop('x2', None)
    if(np.min(myvar)==np.max(myvar)):
        print("The quantity you are trying to plot is a constant = %g." % np.min(myvar))
        return
    cb = kwargs.pop('cb', False)
    nc = kwargs.pop('nc', 15)
    k = kwargs.pop('k',0)
    mirrorx = kwargs.pop('mirrorx',0)
    mirrory = kwargs.pop('mirrory',0)
    symmx = kwargs.pop('symmx',0)
    #cmap = kwargs.pop('cmap',cm.jet)
    isfilled = kwargs.pop('isfilled',False)
    xy = kwargs.pop('xy',0)
    xcoord = kwargs.pop("xcoord",None)
    ycoord = kwargs.pop("ycoord",None)
    lin = kwargs.pop('lin',0)
    xmax = kwargs.pop('xmax',10)
    ymax = kwargs.pop('ymax',5)
    cbxlabel = kwargs.pop('cbxla',None)
    cbylabel = kwargs.pop('cbyla',None)
    fntsize = kwargs.pop("fntsize",20)
    cbgoodticks = kwargs.pop("cbgoodticks",1)
    xlabel = kwargs.pop("xla",None)
    ylabel = kwargs.pop("yla",None)
    dobh = kwargs.pop("dobh",1)
    pretty = kwargs.pop("pretty",0)
    ax = kwargs.pop("ax",None)
    cbticks = kwargs.pop("cbticks",None)
    domathify = kwargs.pop("domathify",0)
    if np.abs(xy)==1:
        if xcoord is None: xcoord = r * np.sin(h)
        if ycoord is None: ycoord = r * np.cos(h)
        if mirrory: ycoord *= -1
        if mirrorx: xcoord *= -1
    if xcoord is not None and ycoord is not None:
        xcoord = xcoord[:,:,None] if xcoord.ndim == 2 else xcoord[:,:,k:k+1]
        ycoord = ycoord[:,:,None] if ycoord.ndim == 2 else ycoord[:,:,k:k+1]
    if np.abs(xy)==1 and symmx:
        if myvar.ndim == 2:
            myvar = myvar[:,:,None] if myvar.ndim == 2 else myvar[:,:,k:k+1]
            myvar=np.concatenate((myvar[:,::-1],myvar),axis=1)
            xcoord=np.concatenate((-xcoord[:,::-1],xcoord),axis=1)
            ycoord=np.concatenate((ycoord[:,::-1],ycoord),axis=1)
        else:
            if myvar.shape[-1] > 1: 
                symmk = (k+nz/2)%nz 
            else: 
                symmk = k
            myvar=np.concatenate((myvar[:,ny-1:ny,k:k+1],myvar[:,::-1,symmk:symmk+1],myvar[:,:,k:k+1]),axis=1)
            xcoord=np.concatenate((xcoord[:,ny-1:ny,k:k+1],-xcoord[:,::-1],xcoord),axis=1)
            ycoord=np.concatenate((ycoord[:,ny-1:ny,k:k+1],ycoord[:,::-1],ycoord),axis=1)
    elif np.abs(xy) == 2 and symmx:
        #if fracphi == 0.5 done in a robust way
        if get_fracphi() < 0.75:
            r1 = np.concatenate((r,r,r[...,0:1]),axis=2)
            ph1 = np.concatenate((ph,ph+np.pi,ph[...,0:1]+2*np.pi),axis=2)
            myvar = np.concatenate((myvar,myvar,myvar[...,0:1]),axis=2)
        else:
            r1 = np.concatenate((r,r[...,0:1]),axis=2)
            ph1 = np.concatenate((ph,ph[...,0:1]+2*np.pi),axis=2)
            myvar = np.concatenate((myvar,myvar[...,0:1]),axis=2)
        xcoord=(r1*cos(ph1))[:,ny/2,:,None]
        ycoord=(r1*sin(ph1))[:,ny/2,:,None]
        myvar = myvar[:,ny/2,:,None]
    else:
        myvar = myvar[:,:,None] if myvar.ndim == 2 else myvar[:,:,k:k+1]
    if lin:
        xcoord = r
        ycoord = h
    if ax is None:
        ax = plt.gca()
    if  xcoord is None or ycoord is None:
        if isfilled:
            res = ax.contourf(myvar[:,:,0].transpose(),nc,**kwargs)
        else:
            res = ax.contour(myvar[:,:,0].transpose(),nc,**kwargs)
    else:
        if isfilled:
            res = ax.contourf(xcoord[:,:,0],ycoord[:,:,0],myvar[:,:,0],nc,**kwargs)
        else:
            res = ax.contour(xcoord[:,:,0],ycoord[:,:,0],myvar[:,:,0],nc,**kwargs)
    if xy>0 and not symmx:
        ax.set_xlim(0,xmax)
        ax.set_ylim(-ymax,ymax)
    if xy> 0 and symmx:
        ax.set_xlim(-xmax,xmax)
        ax.set_ylim(-ymax,ymax)
    if xlabel is not None:
        ax.set_xlabel(xlabel,fontsize=fntsize)
    if ylabel is not None:
        ax.set_ylabel(ylabel,fontsize=fntsize)
    if pretty:
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontsize(fntsize)
            if domathify: mathify_axes_ticks(ax,fontsize=fntsize)
    if cb: #use color bar
        cb = plt.colorbar(res,ax=ax)
        if pretty and cbgoodticks and cbticks is None:
            vmin = cb.vmin
            vmax = cb.vmax
            #this returns incorrect ticks! so ignore it
            #ticks = cb.ax.get_yticks()
            #nticks = len(ticks)
            #if not too many ticks, then pretty them up
            rvmin = np.round(vmin)
            rvmax = np.round(vmax)
            if rvmin == vmin and rvmax == vmax and vmax-vmin <= 10:
                ticks = np.arange(rvmin,rvmax+1)
                cb.set_ticks(ticks)
                mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=ticks)
            elif rvmin == vmin and rvmax == vmax and vmax-vmin <= 20:
                ticks = np.arange(rvmin,rvmax+1)[::2]
                cb.set_ticks(ticks)
                mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=ticks)
        if cbticks is not None:
            cb.set_ticks(cbticks)
            mathify_axes_ticks(cb.ax,fontsize=fntsize,yticks=cbticks)
        if cbxlabel is not None:
            cb.ax.set_xlabel(cbxlabel,fontsize=fntsize)
        if cbxlabel is not None:
            cb.ax.set_xlabel(cbxlabel,fontsize=fntsize)
        if cbylabel is not None:
            cb.ax.set_ylabel(cbylabel,fontsize=fntsize)
        if pretty:
            for label in cb.ax.get_yticklabels():
                label.set_fontsize(fntsize)
    if xy and dobh and "rhor" in globals(): 
        el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
        art=ax.add_artist(el)
        art.set_zorder(20)
    if cb:
        return res, cb
    else:
        return res

def faraday():
    global omegaf1, omegaf2
    if 'omegaf1' in globals():
        del omegaf1
    if 'omemaf2' in globals():
        del omegaf2
    omegaf1=fFdd(0,1)/fFdd(1,3)
    omegaf2=fFdd(0,2)/fFdd(2,3)

def Tcalcud():
    global Tud, TudEM, TudMA
    global mu, sigma
    global enth
    global unb, isunbound, ug
    pg = (gam-1)*ug
    w=rho+ug+pg
    eta=w+bsq
    if 'Tud' in globals():
        del Tud
    if 'TudMA' in globals():
        del TudMA
    if 'TudEM' in globals():
        del TudEM
    if 'mu' in globals():
        del mu
    if 'sigma' in globals():
        del sigma
    if 'unb' in globals():
        del unb
    if 'isunbound' in globals():
        del isunbound
    Tud = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    TudMA = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    TudEM = np.zeros((4,4,nx,ny,nz),dtype=np.float32,order='F')
    for kapa in np.arange(4):
        for nu in np.arange(4):
            if(kapa==nu): delta = 1
            else: delta = 0
            TudEM[kapa,nu] = bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta - bu[kapa]*bd[nu]
            TudMA[kapa,nu] = w*uu[kapa]*ud[nu]+pg*delta
            #Tud[kapa,nu] = eta*uu[kapa]*ud[nu]+(pg+0.5*bsq)*delta-bu[kapa]*bd[nu]
            Tud[kapa,nu] = TudEM[kapa,nu] + TudMA[kapa,nu]
    mu = -Tud[1,0]/(rho*uu[1])
    sigma = TudEM[1,0]/TudMA[1,0]
    enth=1+ug*gam/rho
    unb=enth*ud[0]
    isunbound=(-unb>1.0)


def aux():
    faraday()
    Tcalcud()

if __name__ == "__main__":
    if False:
        #1D plot example
        plt.clf()
        rg("gdump")
        rd("dump000")
        plt.plot(r[:,ny/2,0],rho[:,ny/2,0])
        plt.xscale("log")
        plt.yscale("log")
        plt.xlabel("r")
        plt.ylabel("rho")
    if True: #false->true
        #2D plot example
        plt.clf()
        #rg("gdump")
        #rd("dump000")
        #R-z plot of the logarithm of density distribution
        plc(r,np.log10(rho),cb=True,xy=1,xmax=100,ymax=50)


def bhole():

    ax = plt.gca()
    el = Ellipse((0,0), 2*rhor, 2*rhor, facecolor='k', alpha=1)
    art=ax.add_artist(el)
    art.set_zorder(20)
    plt.draw()



def testfail(fldname = "dump000"):
    try: 
        rd(fldname)
    except IOError as e:
        print("I/O error({0}): {1}".format(e.errno, e.strerror))


def get_sorted_file_list(prefix="dump"):
    flist0 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9]"%prefix) ) )
    flist1 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9][0-9]"%prefix) ) )
    flist2 = np.sort(glob.glob( os.path.join("dumps/", "%s[0-9][0-9][0-9][0-9][0-9]"%prefix) ) )
    flist0.sort()
    flist1.sort()
    flist2.sort()
    flist = np.concatenate((flist0,flist1,flist2))
    return flist
    


def fFdd(i,j):
    if i==0 and j==1:
        fdd =  gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # f_tr
    elif i==1 and j==0:
        fdd = -gdet*(uu[2]*bu[3]-uu[3]*bu[2]) # -f_tr
    elif i==0 and j==2:
        fdd =  gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # f_th
    elif i==2 and j==0:
        fdd = -gdet*(uu[3]*bu[1]-uu[1]*bu[3]) # -f_th
    elif i==0 and j==3:
        fdd =  gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # f_tp
    elif i==3 and j==0:
        fdd = -gdet*(uu[1]*bu[2]-uu[2]*bu[1]) # -f_tp
    elif i==1 and j==3:
        fdd =  gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # f_rp = gdet*B2
    elif i==3 and j==1:
        fdd = -gdet*(uu[2]*bu[0]-uu[0]*bu[2]) # -f_rp = gdet*B2
    elif i==2 and j==3:
        fdd =  gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # f_hp = gdet*B1
    elif i==3 and j==2:
        fdd = -gdet*(uu[0]*bu[1]-uu[1]*bu[0]) # -f_hp = gdet*B1
    elif i==1 and j==2:
        fdd =  gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # f_rh = gdet*B3
    elif i==2 and j==1:
        fdd = -gdet*(uu[0]*bu[3]-uu[3]*bu[0]) # -f_rh = gdet*B3
    else:
        fdd = np.zeros_like(uu[0])
    return fdd

delta = lambda kapa,nu: (kapa==nu)
fTudEM = lambda kapa,nu: bsq*uu[kapa]*ud[nu] + 0.5*bsq*delta(kapa,nu) - bu[kapa]*bd[nu]
fTudMA = lambda kapa,nu: (rho+gam*ug)*uu[kapa]*ud[nu]+(gam-1)*ug*delta(kapa,nu)
fTud = lambda kapa,nu: fTudEM(kapa,nu) + fTudMA(kapa,nu)
fRud = lambda kapa,nu: 4./3.*Erf*uradu[kapa]*uradd[nu]+1./3.*Erf*delta(kapa,nu)

def odot(a,b):
    """ Outer product of two vectors a^mu b_nu"""
    #the shape of the product is (4,4,nx,ny,max(a.nz,b.nz))
    outer_product = np.zeros(np.concatenate((np.array((4,4)),amax(a[0].shape,b[0].shape))),dtype=np.float32,order='F')
    for mu in np.arange(4):
        for nu in np.arange(4):
            outer_product[mu,nu] = a[mu]*b[nu]
    return(outer_product)


def amax(arg1,arg2):
    return(np.maximum(arg1,arg2))

def amin(arg1,arg2):
    return(np.minimum(arg1,arg2))


            

#############################
#
# Costas functions section
#
#############################


def csfn_l():
# calculates the specific angular momentum parameter using the definition l=-u_phiphi/u_tt
    global ud
    l_ang=np.empty((nx,ny));l_ang[:]=np.nan
        
    for i in range(0,nx):
          for j in range(0,ny):
               l_ang[i][j]=-ud[3][i][j][0]/ud[0][i][j][0]
    return l_ang

def csfn_Om():
# calculates the rotation velocity using the definition Om=u^phiphi/u^tt
    global uu
    Om=np.empty((nx,ny));Om[:]=np.nan
        
    for i in range(0,nx):
          for j in range(0,ny):
               Om[i][j]=uu[3][i][j][0]/uu[0][i][j][0]
    return Om

def csfn_lrz_slow():
# finds the lorentz factor using the covariant 4 velocity components uU^mu (another way is simply to divide 
# the spatial components with uU^0).This is a slow function, use the csfn_lrz() instead

    global uu

    lorentz=np.empty((nx,ny));lorentz[:]=np.nan
    for i in range(0,nx):
          for j in range(0,ny):
               lorentz[i][j]=np.sqrt(gdd[1][1][i][j][0]*np.square(uu[1][i][j][0])+gdd[2][2][i][j][0]*np.square(uu[2][i][j][0])+gdd[3][3][i][j][0]*np.square(uu[3][i][j][0])+1)
    return lorentz

def csfn_lrz():
# This is the normal function to calculate lorentz factor

    global uu,alpha

    lor=alpha*uu[0]
    lorentz=np.squeeze(lor)

    return lorentz

def csfn_hdf2harm_file(Quah_file,float_type,nx,ny,nz,n4=None,n5=None):
# This function just re-arrenges the h5 files read vectors to vectors used by the rest harm script

    if (n5==None):
        if (n4==None):
            print("not able to reshape object, not of rank 4 (vector) or rank 5 (2 order tensor)")
            sys.exit(42)
        else:
            Quahin=np.memmap(Quah_file, dtype=float_type, mode='r+', shape=(nx,ny,nz,n4))
            Varin=np.memmap(Quah_file+'.tmp', dtype=float_type, mode='w+', shape=(4,nx,ny,nz))
            if n4==3:
                for i in range(0,3):
                    Varin[i+1,:,:,:]=Quahin[:,:,:,i]
                Varin[0,:,:,:,]=0     
            else:
                Varin[:]=np.moveaxis(Quahin,3,0);

            del Quahin;Quahin=np.memmap(Quah_file, dtype=float_type, mode='w+', shape=(4,nx,ny,nz))
            Quahin[:]=Varin[:]; del Varin;

    else:
        Quahin=np.memmap(Quah_file, dtype=float_type, mode='r+', shape=(nx,ny,nz,n4,n5))
        Qouter=np.memmap("mems/Qu_outer_file.mymemmap", dtype=float_type, mode='w+', shape=(n4,n5,nx,ny,nz))
        #Quahin=np.moveaxis(Quahin,3,0);
        #Quahin=np.moveaxis(Quahin,4,1);
        for i in range(0,n4):
            for j in range(0,n5):
                Qouter[i,j,:,:,:]=Quahin[:,:,:,i,j]
    #Quahin=Qouter;del Qouter;
    del Quahin
   ##next line was commented out, due to err display in terminal 
   # os.system("mv ./mems/Qu_outer_file.mymemmap %s"%Quah_file);#os.system("rm ./mems/Qu_outer_file.mymemmap")

    return()
    #return("Qu_outer_file.mymemmap")

def csfn_hdf2harm(Quah5):
# This function just re-arrenges the h5 files read vectors to vectors used by the rest harm script
    ss=Quah5.shape;
    nx=ss[0];
    ny=ss[1];
    nz=ss[2];
    ncompa=ss[3];
    if (Quah5.ndim==5):
        ncompb=ss[4];
        Qouter=np.empty((ncompa,ncompb,nx,ny,nz));Qouter[:]=np.nan;
        for i in range(0,ncompa):
            for j in range(0,ncompb):
                Qouter[i,j,:,:,:]=Quah5[:,:,:,i,j]
            
    elif (Quah5.ndim==4):
        Qouter=np.empty((4,nx,ny,nz));Qouter[:]=np.nan;
        if ncompa==4:
            #Qouter=np.empty((ncompa,nx,ny,nz));Qouter[:]=np.nan;
            for i in range(0,ncompa):
                Qouter[i,:,:,:]=Quah5[:,:,:,i]
        else:
            Qouter[0,:,:,:]=0.
            for i in range(0,ncompa):
                Qouter[i+1,:,:,:]=Quah5[:,:,:,i]
    else:
        print("not able to reshape object, not of rank 4 (vector) or rank 5 (2 order tensor)")
    return(Qouter)

def csfn_gread(dump_directory=None):
    #global nx,ny,nz,_dx1,_dx2,_dx3,a,Rin,Rout,hslope,R0,x1,x2,x3,r,h,ph,gcov,gcon,alpha
    #global gdet,gam
    global gam,_dx2,nx,ny,nz,float_type,a,alpha

    '''
    The following global variables are missing from our version, 
    t,N1,N2,N3,ti,tj,tk,drdx,gn3,gv3,guu,gdd,dxdxp, games
    guu,gdd and gn3,gv3 are redudant since they are contained in the gcov,gcon.
    
    The following quantities are added
    alpha (lapse)
    
    execute the following command
    print(f2)
    for name in f2:
        print(name)
    to see what other variables might exist on the file
    '''
    
    if os.path.isfile(dump_directory+"/coords.h5"):
        print("Maximum memory use, before hdf.File ",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"KB")
        f2=hdf.File(dump_directory+"/coords.h5", 'r')
        print("Maximum memory use, after hdf.File ",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"KB")
        print(f2)
        #print("Under this database the following groups occur")
        #for name in f2:
        #    print(name)
        #print(" ")

        gcov_found=np.False_
        gcon_found=np.False_
        dxpdx_found=np.False_
        for name in f2:
            if(name=="gcov"):
                print("gcov found in %s"%(dump_directory+"/coords.h5"))
                gcov_found=np.True_
            if (name=="gcon"):
                print("gcon found in %s"%(dump_directory+"/coords.h5"))
                gcon_found=np.True_
            if (name=="dxpdx"):
                print("dxpdx found in %s"%(dump_directory+"/coords.h5"))
                dxpdx_found=np.True_



        nx=f2['N1'][()][0]#const
        ny=f2['N2'][()][0]#const
        nz=f2['N3'][()][0]#const
        if nz>1:
                nz = nz + 1

        R0=f2['R0'][()][0]#const
        float_type=type(R0).__name__

        _dx1=f2['dx1'][()];#const 
        _dx2=f2['dx2'][()];#const
        _dx3=f2['dx3'][()];#const
        gam=f2['gam'][()];#const
        hslope=f2['hslope'][()];#const

        Rin=f2['Rin'][()];#const
        Rout=f2['Rout'][()];#const
        a=f2['a'][()][0];#const

        r=np.memmap('mems/r.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        r[:]=f2['r'][()]; #del r

        h=np.memmap('mems/h.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        h[:]=f2['th'][()]; #del h

        ph=np.memmap('mems/ph.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        ph[:]=f2['ph'][()]; #del ph

        global gdet
        gdet=np.memmap('mems/gdet.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz));
        #ff2=f2['gdet'][()];gdet[:]=ff2[:];print("gdet in arg =0",np.argwhere(gdet==0));del gdet
        gdet[:]=f2['gdet'][()];#del gdet
        
        x1=f2['X1'][()];
        x2=f2['X2'][()];
        x3=f2['X3'][()];

        if (gcov_found):
            #print("Maximum memory use, after gdet",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"KB")
            global gcov
            gcov=np.memmap('mems/gcov.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz,4,4))
            gcov[:]=f2['gcov'][()];
            # here the reading takes lot of memory
            #print("Maximum memory use, before csfn_hdf2harm_file",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"KB")
            csfn_hdf2harm_file('mems/gcov.mymemmap',float_type,nx,ny,nz,n4=4,n5=4)

        
        if (gcon_found):
            gcon=np.memmap('mems/gcon.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz,4,4))
            gcon[:]=f2['gcon'][()];del gcon
            csfn_hdf2harm_file('mems/gcon.mymemmap',float_type,nx,ny,nz,n4=4,n5=4)
            gcon=np.memmap('mems/gcon.mymemmap', dtype=float_type, mode='r', shape=(4,4,nx,ny,nz))

        if (dxpdx_found):
            global dxpdx
            dxpdx=np.memmap('mems/dxpdx.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz,4,4))
            dxpdx[:]=f2['dxpdx'][()];
            csfn_hdf2harm_file('mems/dxpdx.mymemmap',float_type,nx,ny,nz,n4=4,n5=4)
            dxpdx=np.memmap('mems/dxpdx.mymemmap', dtype=float_type, mode='r', shape=(4,4,nx,ny,nz))
        
        alpha=np.memmap('mems/alpha.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz));alpha[:]=np.nan
        #alpha=np.empty((nx,ny,nz));alpha[:]=np.nan
        #rhor=1+(1-a**2)**0.5 #the horizon
        #print("Maximum memory use, after first reading ",resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,"KB")
        if (gcon_found):
            ''' That is a little insecure, but inside the horizon metric changes signature. Maybe we might examine 
            if the horizon is crossed properly'''
            with np.errstate(invalid='ignore'):
                alpha[:,:,:]=np.sqrt(-gcon[0,0,:,:,:])
            del gcon
        #del alpha
    else:
        print("ERROR: The files of the initial coords coords.h5 does not exist!")
    #return()
    #return(nx,ny,nz,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,x1,x2,x3,r,h,ph,gcov,gcon,alpha)
    f2.close()
    gc.collect()
    return(nx,ny,nz,r,h,ph,_dx1,_dx2,_dx3,a,gam,Rin,Rout,hslope,R0,x1,x2,x3,float_type)  

    
def csfn_dumpread(dump,ggam,nnx,nny,nnz,fl_type,dump_directory=None,mem_directory=None):
    global t,rho,ug,vu,B,pg,divb,uu,ud,bu,bd,bsq,ktot,rhor,Lambda_sim,P_dump,T_dump,Hd,Qnu,Tau,Yee
    #global gdet,gam
    global gam,nx,ny,nz,float_type,a
    gam=ggam;nx=nnx;ny=nny;nz=nnz;float_type=fl_type
    '''
    The following global variables are missing from our version,
    ti,tj,tk,cs2,Sden,U,gdetB,v1m,v1p,v2m,v2p 
    
    The following quantities are added (but not globalized)
    dt,dumps
    
    the following are redudant since they have already defined in the csfn_gread procedure
    nx,ny,nz,x1,x2,x3,r,h,ph,_dx1,_dx2,_dx3,hslope,a,R0,Rin,Rout,alpha
    
    The following are needed for the subsequent calculations
    gdet,gam
    
    notice that we added an extra component B[0]=0 and vu[0] to match Tschekhovskoy quantities
    and we devided by gdet to get the relevant values.
    
    execute the following command
    print(f1)
    for name in f1:
        print(name)
    to see what other variables might exist on the file
    '''
    if dump_directory==None:
        dump_directory="./dumps"
    if mem_directory==None:
        mem_directory="./mems/"
    
    #print(dump_directory + "/" + dump + ".h5")
    if os.path.isfile(dump_directory + "/" + dump + ".h5"):
        f1=hdf.File(dump_directory +  "/" + dump + ".h5", 'r')

        rho=np.memmap(mem_directory+'rho.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        rho[:]=f1['Rho'][()];

        dumps=f1['dump_cnt'][()]

        B1B2B3gdet=np.memmap(mem_directory +'B1B2B3gdet.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz,3))
        #B1B2B3gdet_out==np.memmap(mem_directory+'B1B2B3gdet_out.mymemmap', dtype=float_type, mode='w+', shape=(4,nx,ny,nz));
        B1B2B3gdet[:]=f1['B1B2B3gdet'][()]; 
        #ii=42;jj=96;kk=0;
        #del B1B2B3gdet; 
        #B1B2B3gdet_out[:]=csfn_hdf2harm(B1B2B3gdet);
        csfn_hdf2harm_file(mem_directory+'B1B2B3gdet.mymemmap',float_type,nx,ny,nz,n4=3)
        B1B2B3gdet=np.memmap(mem_directory+'B1B2B3gdet.mymemmap', dtype=float_type, mode='r', shape=(4,nx,ny,nz))
 

        B=np.memmap(mem_directory+'B.mymemmap', dtype=float_type, mode='w+', shape=(4,nx,ny,nz))
        gdet=np.memmap(mem_directory+'gdet.mymemmap', dtype=float_type, mode='r', shape=(nx,ny,nz))
        B[0,:,:,:]=0
        for i in range(0,3):
            B[i+1,:,:,:]=B1B2B3gdet[i+1,:,:,:]/gdet[:,:,:]
        #del B,B1B2B3gdet

        '''check if the U1U2U3gdet corresponds to vu'''
        U1U2U3gdet=np.memmap(mem_directory+'U1U2U3gdet.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz,3))
        U1U2U3gdet[:]=f1['U1U2U3gdet'][()]
        U1U2U3gdet=csfn_hdf2harm(U1U2U3gdet)
        vu=np.memmap(mem_directory+'vu.mymemmap', dtype=float_type, mode='w+', shape=(4,nx,ny,nz))
        #vu=np.empty((4,nx,ny,nz));vu[:]=np.nan;
        vu[0,:,:,:]=0
        for i in range(0,3):
            vu[i+1,:,:,:]=U1U2U3gdet[i+1,:,:,:]/gdet
        #del vu,U1U2U3gdet,gdet
    
        dt=f1['dt'][()]
        t=f1['t'][()][0]
        nstep=f1['nstep'][()][0]

        divb=np.memmap(mem_directory+'divb.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        divb[:]=f1['DivB'][()]; #del divb

        bu=f1['bcon'][()]
        bd=f1['bcov'][()]
        bu=csfn_hdf2harm(bu)
        bd=csfn_hdf2harm(bd)
        uu=f1['ucon'][()]
        ud=f1['ucov'][()]
        uu=csfn_hdf2harm(uu)
        ud=csfn_hdf2harm(ud)
        ug=np.memmap(mem_directory+'ug.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        ug[:]=f1['Energy'][()]
        pg =np.memmap(mem_directory+'pg.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        pg[:] = (gam-1.0)*ug[:]; #del ug

        ktot=np.memmap(mem_directory+'ktot.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        ktot[:]=pg/rho**gam; #del ktot,rho

        bsq=np.memmap(mem_directory+'bsq.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        beta=np.memmap(mem_directory+'beta.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        bsq[:]=mdot(bu,bd); 
        with np.errstate(divide='ignore'):
            beta[:]=2*pg[:]/bsq[:]
        #del bsq,pg,beta
        rhor=1+(1-a**2)**0.5 #the horizon

        # Nu-cooling quantities:
        Lambda_sim = np.memmap(mem_directory+'Lambda_sim.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        Lambda_sim[:] = f1['Lambda_sim'][()];
        P_dump = np.memmap(mem_directory+'P_dump.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        P_dump[:] = f1['P'][()];
        T_dump = np.memmap(mem_directory+'T_dump.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        T_dump[:] = f1['T'][()];
        Hd = np.memmap(mem_directory+'Hd.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        Hd[:] = f1['Hd'][()];
        Qnu = np.memmap(mem_directory+'Qnu.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz))
        Qnu[:] = f1['Qnu'][()];
        Tau = np.memmap(mem_directory+'Tau.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz)) #neutrino optical depth
        Tau[:] = f1['Tau'][()];
        Yee = np.memmap(mem_directory+'Yee.mymemmap', dtype=float_type, mode='w+', shape=(nx,ny,nz)) #electron fraction
        Yee[:] = f1['Yee'][()];

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
        M_BLH = 2.0 * MSUN

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
        
    else:
        print()
        print("ERROR: The file of the data <%s> does not exist!"%(dump_directory + "/" + dump + ".h5"))
        print()
    gc.collect()
    return(t,nstep,rho,ug,vu,B,pg,divb,uu,ud,bu,bd,bsq,ktot,rhor,beta,Lambda_sim,P_dump,T_dump,Hd,Qnu,Tau,Yee)
    #return(t,uu,ud,bu,bd)
    return()

def csfn_psicalc(phi_j=None,mem_directory=None):
    """
    Computes the field vector potential among a specific phi=phi_j plane
    """
    #global B,gdet,_dx2

    if mem_directory==None:
        mem_directory="./mems/"

    if phi_j==None:
        phi_j = 0
    
    B=np.memmap(mem_directory+'B.mymemmap', dtype=float_type, mode='r', shape=(4,nx,ny,nz))
    gdet=np.memmap(mem_directory+'gdet.mymemmap', dtype=float_type, mode='r', shape=(nx,ny,nz))
    B1 = B[1,:,:,phi_j]
    gdet_j=gdet[:,:,phi_j]
    del B,gdet;
    daphi = -(gdet_j*B1)*_dx2 #we removed the average of Br along the phi direction.   
    aphi=daphi[:,::-1].cumsum(axis=1)[:,::-1] #the [:,::-1] interchanges the order of the elements starting from end (th=\pi) to begining (th=0). Cumsum adds all the values of daphi allong the specific r=const surface
    aphi-=0.5*daphi #correction for half-cell shift between face and center in theta
    return(aphi)

def csfn_spher2cart(r,th,ph,A):
    '''
    Performs the simple (newtonian) transformation of a vector A from spherical to Cartesian Coordinates.
    comment expressions is to use the BLAS function instead for speed. It can be improved by
    a. calculate the transformation matrix once at the metric reading and setting it as global
    b. if we constrain the values through rmax before
    '''
    if ((r.ndim!=3) or (th.ndim!=3) or (ph.ndim!=3)):
        print("Error, your coordinates are not a 3D scalar field")
        return
    elif A.ndim!=4:
        print("Error, you must provide a 4D vector field")
        return
    else:
        nx,ny,nz=r.shape
        Acart=np.empty((3,nx,ny,nz));Acart[:]=np.nan
        Lsph2cart=np.empty((3,3),dtype=type(A[1,1,1,1]));
        Aspher=np.empty((1,3),dtype=type(A[1,1,1,1]));                   
        #gemm_dot = linalg.get_blas_funcs("gemm", arrays=(Aspher,Lsph2cart))
        for i in range(0,nx):
            for j in range(0,ny):
                for k in range(0,nz):
                    #thi=th[i,j,k];
                    #phii=ph[i,j,k]
                    sinth=np.sin(th[i,j,k]);costh=np.cos(th[i,j,k]);
                    sinphi=np.sin(ph[i,j,k]);cosphi=np.cos(ph[i,j,k]);
                    Lsph2cart=np.array([[sinth*cosphi,sinth*sinphi,costh],
                              [costh*cosphi,costh*sinphi,-sinth],
                              [-sinphi,cosphi,0]],dtype=type(A[1,1,1,1]))#Transfomation matrix
                    Aspher=np.array([A[0,i,j,k],A[1,i,j,k],A[2,i,j,k]],dtype=type(A[1,1,1,1]))
                    #Ares=gemm_dot(alpha=1.0, a=Aspher, b=Lsph2cart, trans_a=Aspher.flags.c_contiguous, trans_b=Lsph2cart.flags.c_contiguous)
                    Ares=np.dot(Aspher,Lsph2cart)
                    Acart[:,i,j,k]=Ares[:]
        return(Acart)
    
def csfn_remove_inner_bh(A):
    '''
        removes the values of A that they are inside the horizon
    '''
    global lapse
    
    inds_ins_hor=np.argwhere(np.isnan(alpha))
    
    if inds_ins_hor.size!=0:
        if (A.ndim==3):
            for i in range(0,inds_ins_hor.shape[0]):
                nnx,nny,nnz=inds_ins_hor[i,:]
                A[nnx,nny,nnz]=np.nan
        elif (A.ndim==4):
            for i in range(0,inds_ins_hor.shape[0]):
                nnx,nny,nnz=inds_ins_hor[i,:]
                A[:,nnx,nny,nnz]=np.nan
    return(A)

def csfn_psicalc3D(Br,phi_j,gdetg,_dx2):
    """
    Computes the field phi vector potential on a particular slice r,th of ph=phi1, at present phi1 is given as an indice
    of the original coordinates. If this procedure is succesful we might use interpolation
    """
    
   
    B1 = np.squeeze(Br[:,1:-1,phi_j])
    gdet=np.squeeze(gdetg[:,1:-1,phi_j])
    aphi=np.empty((gdetg.shape[0],gdetg.shape[1]));aphi[:]=np.nan
    daphi = -(gdet*B1)*_dx2 #that's the mean value of the (-gdet*B1*_dx2) array, where gdet,B1 arrays and _dx2 scalar,  taken the average from the last dimension. ie it averages Br along the phi direction.   
    print(daphi.shape)
    aphi2=daphi[:,::-1].cumsum(axis=1)[:,::-1] #the [:,::-1] interchanges the order of the elements starting from end (th=\pi) to begining (th=0). Cumsum adds all the values of daphi allong the specific r=const surface
    aphi2-=0.5*daphi #correction for half-cell shift between face and center in theta
    aphi[:,1:-1]=aphi2
    aphi[:,0]=aphi[:,1]
    aphi[:,-1]=aphi[:,-2]
    return(aphi)
