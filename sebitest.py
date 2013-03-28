#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.fftpack


def PlotSpectrum(trace):
    '''
    calculate and plot frequency spectrum of a single trace
    '''
    FFT = abs(scipy.fft(data[1,trace,:]))
    freqs = scipy.fftpack.fftfreq(data[1,1,:].size, dt)

    # plot
    plt.figure()
    plt.plot(freqs,FFT,'-')
    plt.xlim([0,200])
    plt.xlabel('Frequency in Hz')

def PlotData(offset):
    '''
    Plot data
    '''
    x = np.linspace(0,nsamples*dt,nsamples)
    plt.figure()
    for n in xrange(1, ntraces):
        plt.plot(data[offset,n,:]+n,x,'k-')

    plt.gca().invert_yaxis()
    plt.xlim([0,ntraces])
    plt.ylim([nsamples*dt,0])
    plt.ylabel('Time in seconds')
    plt.title('Trace')
    plt.gca().xaxis.tick_top()


nsamples = 1001
ntraces = 101
noffsets = 5

dt = 0.002


f = open('SEIS-filt', 'r')

# read data in ieee754 format
inputs = np.fromfile(f, dtype=np.float32)

# reshape
data = inputs.reshape(noffsets,ntraces,nsamples)

######################################
#PlotSpectrum(59)

#PlotData(1)
#PlotData(4)
#plt.show()
nx = 100
dx = 20
nz = 101
dz = 100
ntrc = ntraces
noff =noffsets

v=9000

offsets = [0, 500, 1000, 1500, 2000]

dcdp = 20

R = np.zeros((20, nz))

for ix in xrange(0, 19):        #loop over discretiziced undergroundpoints in x
    x = dx*ix
    for iz in xrange(1, nz):    #loop over discretiziced undergroundpoints in z
        z = dz*iz               #(depth)
        
        for itrc in xrange(0, ntrc-1):     #loop over all traces
            for ioff in xrange(0, noff-1):  #loop over all offsets
                ksi = dcdp * itrc           # cdp point
                h = offsets[ioff]/2         # half offset
                rs = np.sqrt( (x - (ksi-h))**2 + z**2)     # distance point-source
                rr = np.sqrt( (x - (ksi+h))**2 + z**2)     # distance point-reviever

                wco = ( z/rs * np.sqrt(rs/rr) + z/rr * np.sqrt(rr/rs) ) /v

                t = (rs + rr)/v
                it = np.floor(t/dt)

                #print rs, rr, ix, iz, t, it

                if (it <= nsamples-1):
                    #print 'data {0:.0e}'.format(data[ioff,itrc,it]*wco/np.sqrt(2*3.14))
                    R[ix,iz] = R[ix,iz] + data[ioff,itrc,it] * wco /np.sqrt(2*np.pi)


#print R[0,:]

x = np.linspace(0,nsamples*dt,nz)
plt.figure()
for n in xrange(0, 19):
    plt.plot(R[n,:]*100+n,x,'k-')

plt.show()
