cimport numpy as np
import numpy as np


def cyMigrate(np.ndarray[float, ndim=3] data,int nx, int dx,int nz,int dz,float dt,int dcdp,int v,list offsets,int nsamples,int ntrc,int noff):
    cdef np.ndarray[float, ndim=2] R = np.zeros((nx, nz),dtype=np.float32)
    cdef int ix, iz, itrc, ioff
    cdef int x, z
    cdef int ksi
    cdef float h, rs, rr, wco, t
    cdef int it
    for ix in range(0, nx):        #loop over discretiziced undergroundpoints in x
        x = dx*ix
        for iz in range(1, nz):    #loop over discretiziced undergroundpoints in z
            z = dz*iz               #(depth)
            
            for itrc in range(0, ntrc-1):     #loop over all traces
                for ioff in range(0, noff-1):  #loop over all offsets
                    ksi = dcdp * itrc           # cdp point
                    h = offsets[ioff]/2         # half offset
                    rs = np.sqrt( (x - (ksi-h))**2 + z**2)     # distance point-source
                    rr = np.sqrt( (x - (ksi+h))**2 + z**2)     # distance point-reviever
    
                    wco = ( z/rs * np.sqrt(rs/rr) + z/rr * np.sqrt(rr/rs) ) /v
    
                    t = (rs + rr)/v             # resulting time
                    it = np.floor(t/dt)         # nearest neighbor for timesample
    
                    #print rs, rr, ix, iz, t, it
    
                    if (it <= nsamples-1):

                        R[ix,iz] = R[ix,iz] + data[ioff,itrc,it] * wco /np.sqrt(2*np.pi)
    
    return R



