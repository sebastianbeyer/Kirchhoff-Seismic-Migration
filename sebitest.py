#!/usr/bin/env python

import struct
import numpy as np

print "hello"

nsamples = 1001
ntraces = 101
noffsets = 5



f = open('SEIS-filt', 'r')
nread = nsamples * ntraces * noffsets


# read data in ieee754 format
inputs = struct.unpack('f'*nread, f.read(4*nread))

print inputs
