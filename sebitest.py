#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


nsamples = 1001
ntraces = 101
noffsets = 5

dt = 0.002


f = open('SEIS-orig', 'r')

# read data in ieee754 format
inputs = np.fromfile(f, dtype=np.float32)


# reshape
#data = inputs.reshape(nsamples,ntraces,noffsets)
data = inputs.reshape(noffsets,ntraces,nsamples)

x = np.linspace(0,nsamples*dt,nsamples)

print x




#plt.plot(data[1,1,:],x,'-')
#plt.plot(data[1,2,:]+1,x,'-')

for n in xrange(1, ntraces):
    plt.plot(data[1,n,:]+n,x,'k-')

plt.gca().invert_yaxis()
plt.show()



