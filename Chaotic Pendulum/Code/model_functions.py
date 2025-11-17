import numpy as np 
from numba import jit

R = 47
C = 1e-6
scale_factor = C * R * 1e3
k = 6
t0 = 0
tmax = 10 / scale_factor
h = 0.125e-3 / scale_factor
N = int((tmax - t0) / h)
half = N // 2

@jit(nopython=True)
def rk4_core(Rv, R0, Vc, R, k, t0, tmax, h, N,x0,scale_factor,renorm_interval):
    alpha = R / Rv
    beta  = R / R0 * Vc

    t = np.linspace(t0, tmax, N+1)
    x = np.zeros((N+1, len(x0)))
    x[0] = x0

    Phi = np.eye(3)
    S = np.zeros(3)
    t_accum = 0.0

    num_renorms = N // renorm_interval

    lam_inst = np.zeros((num_renorms, 3))

    m = 0 

    for n in range(N):
        tn = t[n]
        xn = x[n]

        if xn[0] >= 0:
            J = np.array([[0.0,1.0,0.0],
                          [0.0,0.0,1.0],
                          [0.0,-1.0,-alpha]])
        else:
            J = np.array([[0.0,1.0,0.0],
                          [0.0,0.0,1.0],
                          [-k,-1.0,-alpha]])

        ftV = np.empty_like(xn)
        k1 = np.empty_like(xn)
        k2 = np.empty_like(xn)
        k3 = np.empty_like(xn)
        k4 = np.empty_like(xn)

        # k1
        ftV[0] = xn[1]
        ftV[1] = xn[2]
        ftV[2] = -alpha * xn[2] - xn[1] - beta
        if xn[0] < 0:
            ftV[2] -= k * xn[0]
        k1[:] = ftV

        # k2
        V2 = xn + (h/2) * k1
        ftV[0] = V2[1]
        ftV[1] = V2[2]
        ftV[2] = -alpha * V2[2] - V2[1] - beta
        if V2[0] < 0:
            ftV[2] -= k * V2[0]
        k2[:] = ftV

        # k3
        V3 = xn + (h/2) * k2
        ftV[0] = V3[1]
        ftV[1] = V3[2]
        ftV[2] = -alpha * V3[2] - V3[1] - beta
        if V3[0] < 0:
            ftV[2] -= k * V3[0]
        k3[:] = ftV

        # k4
        V4 = xn + h * k3
        ftV[0] = V4[1]
        ftV[1] = V4[2]
        ftV[2] = -alpha * V4[2] - V4[1] - beta
        if V4[0] < 0:
            ftV[2] -= k * V4[0]
        k4[:] = ftV

        x[n+1] = xn + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

        Phi = Phi + h * (J @ Phi)

        if (n+1) % renorm_interval == 0:
            Q, Rmat = np.linalg.qr(Phi)
            diagR = np.abs(np.diag(Rmat))
            logs  = np.log(diagR)

            S += logs

            dt_block = h * renorm_interval
            lam_inst[m, :] = logs / dt_block

            Phi = Q
            t_accum += dt_block
            m += 1

    lambdas = S / t_accum
    t_real = t * scale_factor
    return t_real, x, lambdas, lam_inst

def rk4(**params):
    Rv = params.get("Rv", 70.0)
    R0 = params.get("R0", 157.0)
    Vc = params.get("Vc", 0.25)
    x0 = params.get("x0", np.array([0.0,0.0,-R / R0 * Vc]))
    renorm_interval = params.get("renorm_interval",10)

    t, x, lambdas,lam_inst = rk4_core(Rv, R0, Vc, R, k, t0, tmax, h, N,x0, scale_factor,renorm_interval)

    return t,x,lambdas,lam_inst,Rv

@jit(nopython=True)
def rk4_pendulum(q,omega_D,g):
    f0 = np.array([0.0,0.2,0.0])
    t = np.linspace(t0,tmax,N+1)
    f = np.zeros((N+1,len(f0))) # f = [omega, theta, phi]
    f[0] = f0
    for n in range(N):
        tn = t[n]
        fn = f[n]

        func = np.zeros_like(fn)
        k1 = np.zeros_like(fn)
        k2 = np.zeros_like(fn)
        k3 = np.zeros_like(fn)
        k4 = np.zeros_like(fn)

        #k1 
        func[0] = -1/q * fn[0] - np.sin(fn[1]) + g * np.cos(fn[2])
        func[1] = fn[0]
        func[2] = omega_D
        k1[:] = func 

        #k2 
        F2 = fn + (h/2) * k1 
        func[0] = -1/q * F2[0] - np.sin(F2[1]) + g * np.cos(F2[2])
        func[1] = F2[0]
        func[2] = omega_D 
        k2[:] = func

        #k3 
        F3 = fn + (h/2) * k2 
        func[0] = -1/q * F3[0] - np.sin(F3[1]) + g * np.cos(F3[2])
        func[1] = F3[0]
        func[2] = omega_D 
        k3[:] = func

        F4 = fn + h * k3 
        func[0] = -1/q * F4[0] - np.sin(F4[1]) + g * np.cos(F4[2])
        func[1] = F4[0]
        func[2] = omega_D 
        k4[:] = func

        f[n+1] = fn + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
    return t, f