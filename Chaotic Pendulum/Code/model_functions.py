import numpy as np 
import scipy.signal as sig
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt
from IPython.display import HTML

t0 = 0
tmax = 500
h = 0.125e-3
N = int((tmax - t0) / h)

def rk4(**params):

    R  = params.get("R", 47.0)
    Rv = params.get("Rv", 70.0)
    R0 = params.get("R0", 157.0)
    Vc = params.get("Vc", 0.25)
    k = params.get("k", 6)
    C  = params.get("C", 1e-3)
    x0  = params.get("x0", np.array([0,0,R/R0 * Vc]))

    def f(t,V):
        ftV = np.array([V[1],
                        V[2],
                -R/Rv * V[2]-V[1] - R/R0*Vc
        ])
        if V[0] < 0:
            ftV[2] -= k *V[0]        
        return ftV
    
    m = np.size(x0)                   
    h = (tmax - t0) / N            
    t = np.arange(t0, tmax + h/2, h)  
    x = np.zeros((N+1, m))            

    x[0,:] = x0                    
    for n in range(N):             
        k1 = f(t[n],       x[n,:]             ,**params)
        k2 = f(t[n] + h/2, x[n,:] + (h/2) * k1,**params)
        k3 = f(t[n] + h/2, x[n,:] + (h/2) * k2,**params)
        k4 = f(t[n] + h,   x[n,:] +     h * k3,**params)
        x[n+1,:] = x[n,:] + (h/6) * ( k1 + 2*k2 + 2*k3 + k4 )  
    t_real = t * R * C
    return t_real, x

