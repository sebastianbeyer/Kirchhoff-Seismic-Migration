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
    plt.show()


def PlotImg(data2d,name):
    data2d = inputs.reshape(noff*ntrc,nsmp) # reshape into 2d array with all offsets

    filepath = "./figures/" + name + ".eps"

    plt.figure()
    imgplot = plt.imshow(data2d.T, extent=[0, ntrc*noff/100,nsmp*dt,0],cmap=colormap)
    plt.axes().set_aspect('auto')                       # set aspect ratio to auto
    plt.xlabel('Offset')
    plt.ylabel('Time in [s]')
    plt.savefig(filepath, bbox_inches=0)
    plt.show()


def Migrate(data,nx,dx,nz,dz,dt,dcdp,v,offsets,nsmp,ntrc,noff):
    
    R = np.zeros((nx, nz))
    for ix in xrange(0, nx):        #loop over discreticized undergroundpoints in x
        x = dx*ix
        for iz in xrange(1, nz):    #loop over discreticized undergroundpoints in z
            z = dz*iz               #(depth)
            
            for itrc in xrange(0, ntrc):     #loop over all traces
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

def benchmarkfunc(mode):
    nx=1
    if (mode == 'python'):
        R = Migrate(data,nx,dx,nz,dz,dt,dcdp,3000,offsets,nsmp,ntrc,noff)
    if (mode == 'cython'):
        R = mod_cyMigrate.cyMigrate(data,nx,dx,nz,dz,dt,dcdp,3000,offsets,nsmp,ntrc,noff)



def benchmark():
    import timeit
    nx = 1
    print "Benchmarking python vs cython:"
    print "Migrating 10 times with nx=",nx
    print "python code:"
    py_time = timeit.timeit("kirchhoffmigration.benchmarkfunc('python')",setup="import kirchhoffmigration", number=10)
    print py_time

    print "cython code:"
    cy_time = timeit.timeit("kirchhoffmigration.benchmarkfunc('cython')",setup="importkirchhoffmigration", number=10)
    print cy_time



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

def full_migration(data,name):
    nx = ntrc
    V = 2950
    filepath = "./figures/" + name + ".eps"
    R = mod_cyMigrate.cyMigrate(data,nx,dx,nz,dz,dt,dcdp,V,offsets,nsmp,ntrc,noff)

    # image
    plt.figure()
    imgplot = plt.imshow(R.T, extent=[0,nx*dx ,nz*dz,0],cmap=colormap)
    plt.ylabel('Depth z in [m]')
    plt.xlabel('x in [m]')
    plt.savefig(filepath, bbox_inches=0)
    plt.show()

    return R


def taper(data,traces=20):
    '''
    Tapers the data to reduce migration artifacts via sinus curve.
    Optional Parameter traces: Number of traces to be affected.
    '''

    amount = list()
    i = 0
    for x in np.linspace(0, np.pi/2,traces):
        amount.append(np.sin(x))
        i=i+1

    for off in xrange(0,noff):
        for n in xrange(0,traces):
            data[off,n,:] = data[off,n,:]*amount[n]
            data[off,ntrc-1-n,:] = data[off,ntrc-1-n,:]*amount[n]

    return data

def check_amplitudes(data):
    maxamp = list()
    for i in xrange(0,data.shape[0]):
        maxamp.append(max(abs(data[i,:])))

    #plot
    plt.figure()
    plt.plot(range(len(maxamp)),maxamp,'-')
    plt.xlabel('Trace')
    plt.ylabel('Maximum amplitude')

    plt.savefig('./figures/amplitudes.eps', bbox_inches=0)


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






#full_migration(taperdata,'tapered')

#PlotWiggle(1)
#PlotSpectrum(59)
#plt.show()
#v_analysis(2000, 5000)
# v_analysis ergibt v=2950 bei z=1970
#PlotImg(data,'unprocessed')
#plt.show()
#benchmark()
#full_migration(data)
#taperdata = taper(data)
#PlotImg(taperdata,'tapered_data')
#full_migration(taperdata,'tapered')



