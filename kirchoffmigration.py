#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.fftpack

# cython for more speed
import mod_cyMigrate

def PlotSpectrum(trace):
    '''
    calculate and plot frequency spectrum of a single trace
    '''
    # use unfiltered data for that
    f = open('SEIS-orig', 'r')
    # read data in ieee754 format
    inputs = np.fromfile(f, dtype=np.float32)
    # reshape
    data = inputs.reshape(noff,ntrc,nsmp)

    FFT = abs(scipy.fft(data[1,trace,:]))
    freqs = scipy.fftpack.fftfreq(data[1,1,:].size, dt)

    half = len(freqs)/2
    # plot
    plt.figure()
    plt.plot(freqs[0:half],FFT[0:half],'-')
    #plt.xlim([0,200])
    plt.xlabel('Frequency in [Hz]')

    plt.savefig('./figures/spectrum.eps', bbox_inches=0)
    print freqs

def PlotWiggle(offset):
    '''
    Do a wiggle plot of selected offset
    '''
    x = np.linspace(0,nsmp*dt,nsmp)
    plt.figure()
    for n in xrange(1, ntrc):
        plt.plot(data[offset,n,:]+n,x,'k-')

    plt.gca().invert_yaxis()
    plt.xlim([0,ntrc])
    plt.ylim([nsmp*dt,0])
    plt.ylabel('Time in [s]')
    plt.title('Trace')
    plt.gca().xaxis.tick_top()


def PlotImg(data2d):
    data2d = inputs.reshape(noff*ntrc,nsmp) # reshape into 2d array with all offsets

    plt.figure()
    imgplot = plt.imshow(data2d.T, extent=[0, ntrc*noff/100,nsmp*dt,0],cmap=colormap)
    plt.axes().set_aspect('auto')                       # set aspect ratio to auto
    plt.xlabel('Offset')
    plt.ylabel('Time in [s]')
    plt.savefig('./figures/unprocessed.eps', bbox_inches=0)
    plt.show()


def Migrate(data,nx,dx,nz,dz,dt,dcdp,v,offsets,nsmp,ntrc,noff):
    
    R = np.zeros((nx, nz))
    for ix in xrange(0, nx):        #loop over discreticized undergroundpoints in x
        x = dx*ix
        for iz in xrange(1, nz):    #loop over discreticized undergroundpoints in z
            z = dz*iz               #(depth)
            
            for itrc in xrange(0, ntrc-1):     #loop over all traces
                for ioff in xrange(0, noff):  #loop over all offsets
                    ksi = dcdp * itrc           # cdp point
                    h = offsets[ioff]/2         # half offset
                    rs = np.sqrt( (x - (ksi-h))**2 + z**2)     # distance point<->source
                    rr = np.sqrt( (x - (ksi+h))**2 + z**2)     # distance point<->reciever
    
                    wco = ( z/rs * np.sqrt(rs/rr) + z/rr * np.sqrt(rr/rs) ) /v
    
                    t = (rs + rr)/v             # resulting time
                    it = np.floor(t/dt)         # nearest neighbor for timesample
    
                    #print rs, rr, ix, iz, t, it
    
                    if (it <= nsmp-1):
                        R[ix,iz] = R[ix,iz] + data[ioff,itrc,it] * wco /np.sqrt(2*np.pi)
    
    return R


def benchmark():
    #R = Migrate(data,nx,dx,nz,dz,dt,dcdp,3000,offsets,nsmp,ntrc,noff)
    R = mod_cyMigrate.cyMigrate(data,nx,dx,nz,dz,dt,dcdp,3000,offsets,nsmp,ntrc,noff)





def v_analysis(vmin, vmax):
    '''
    Velocity analysis due to testing a range of velocities (from vmin to vmax)
    plotting them together in one plot.
    '''

    nvels = 100                         # number of velocities to be tested
    V = np.linspace(vmin,vmax,nvels)
    nx=1                                # only compute one trace
    M = np.zeros((nz, nvels))
    for n in xrange(0, len(V)-1):
        R = mod_cyMigrate.cyMigrate(data,nx,dx,nz,dz,dt,dcdp,V[n],offsets,nsmp,ntrc,noff)
        M[:,n] = R[0,:]

    plt.figure()
    imgplot = plt.imshow(M, extent=[vmin, vmax, nz*dz, dz])
    plt.ylabel('Depth in [m]')
    plt.xlabel('Velocity in [m/s]')
    plt.savefig('./figures/v_analysis.eps', bbox_inches=0)
    plt.show()

    return M

def full_migration():
    nx = ntrc
    V = 2950
    R = mod_cyMigrate.cyMigrate(data,nx,dx,nz,dz,dt,dcdp,V,offsets,nsmp,ntrc,noff)
    plt.figure()
    mgplot = plt.imshow(R.T, extent=[0,nx*dx ,nz*dz,0],cmap=colormap)
    plt.ylabel('Depth z in [m]')
    plt.xlabel('x in [m]')
    plt.savefig('./figures/full_migration.eps', bbox_inches=0)
    plt.show()

    return R


def damp(data):

    traces=10
    for off in xrange(0,noff):
        for n in xrange(0,traces):
            data[off,n,:] = 0
            data[off,ntrc-1-n,:] = 1

    PlotImg(data)

##################################################################

ntrc = 101
noff = 5
nsmp = 1001

dt = 0.002


offsets = np.array([0, 500, 1000, 1500, 2000])

nz=300
dx=20
dz=10
nx=1
dcdp=20



colormap = 'gist_yarg'

f = open('SEIS-filt', 'r')
# read data in ieee754 format
inputs = np.fromfile(f, dtype=np.float32)
# reshape
data = inputs.reshape(noff,ntrc,nsmp)






#PlotWiggle(1)
#PlotSpectrum(59)
#plt.show()
#v_analysis(2000, 5000)
# v_analysis ergibt v=2950 bei z=1970
#PlotImg(data)
#plt.show()
#benchmark()
#full_migration()
damp(data)




