#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt


nsamples = 1001
ntraces = 101
noffsets = 5



f = open('SEIS-filt', 'r')

# read data in ieee754 format
inputs = np.fromfile(f, dtype=np.float32)


# reshape
#data = inputs.reshape(nsamples,ntraces,noffsets)
data = inputs.reshape(noffsets,ntraces,nsamples)

x = np.linspace(1,nsamples,nsamples)


print data.shape
print x.size

#plt.plot(x,data[:,1,1],'-')

plt.plot(data[1,2,:],x,'-')
plt.show()



