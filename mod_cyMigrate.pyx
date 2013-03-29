cimport numpy as np
import numpy as np

def cyMigrate(data,int nx, int dx,int nz,int dz,double dt,int dcdp,v,offsets,int nsamples,int ntrc,int noff):
    R = np.zeros((nx, nz))
    for ix in xrange(0, nx):        #loop over discretiziced undergroundpoints in x
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
    
                    t = (rs + rr)/v             # resulting time
                    it = np.floor(t/dt)         # nearest neighbor for timesample
    
                    #print rs, rr, ix, iz, t, it
    
                    if (it <= nsamples-1):

                        R[ix,iz] = R[ix,iz] + data[ioff,itrc,it] * wco /np.sqrt(2*np.pi)
    
    return R



