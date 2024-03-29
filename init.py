#from data import GPEdata,CNdata
from numpy import loadtxt,reshape,pi,sqrt,exp
import random
def initialize(g,nstp)-> None:
    qr_al = sqrt(sqrt(g.gamma * g.nu) / pi)
    random.seed(13)
    if nstp == 0:
        re,im = loadtxt("input_wave.txt",unpack=True)
        psi1d = re + 1j * im
        g.psi = reshape(psi1d,[g.nx,g.ny])
    else:
        for i in range(g.ny):
            for j in range(g.nx):
                g.r[i,j] = sqrt(g.x2[i] + g.y2[j])
                k = i+1 * j+1
                tmp = (g.gamma * g.x2[i] + g.nu * g.y2[j]) / 2.0
                g.psi[i,j] = qr_al * exp(-tmp) * (g.x[i] + 1j * g.y[j]) * exp(2.0 * pi * 1j)

def calc_trap(g) -> None:
     gamma2 = g.gamma * g.gamma
     nu2 = g.nu * g.nu
     for i in range(g.ny):
         for j in range(g.nx):
             g.v[i,j] = g.xop * (gamma2 * g.x2[i] + nu2 * g.y2[j])/2.0

def coef(g,cn) -> None:
    dx2 = g.dx * g.dx
    dy2 = g.dy * g.dy
    dxx = 1.0 / dx2
    dyy = 1.0 / dy2
    cdt = g.cij * g.dt
    ca0 = 1.0 + cdt * dxx
    cn.ca0r = 1.0 - cdt * dxx
    cb0 = 1.0 + cdt * dyy
    cn.cb0r = 1.0 - cdt * dyy
    
    c0 = cdt * 1j * g.xop * g.omega / 2.0
    
    ctmpx = c0 / (2.0 * g.dx)
    ct0 = cdt * dxx / 2.0
    cn.cam = -(ct0 - ctmpx * g.y)
    cn.cap = -(ct0 + ctmpx * g.y)

    for j in range(g.ny):
        cn.cala[g.nxx-1,j] = 0.0
        cn.cgaa[g.nxx-1,j] = -1.0 / ca0
        for i in range(g.nxx-1,1,-1):
            cn.cala[i-1,j] = cn.cam[j] * cn.cgaa[i,j]
            cn.cgaa[i-1,j] = -1.0 / (ca0 + cn.cap[j] * cn.cala[i-1,j])

    ctmpy = c0 / (2.0 * g.dy)
    ct0 = cdt * dyy / 2.0
    cn.cbm = -(ct0 + ctmpy * g.x)
    cn.cbp = -(ct0 - ctmpy * g.x)

    for i in range(g.nx):
        cn.calb[i,g.nyy-1] = 0.0
        cn.cgab[i,g.nyy-1] = -1.0 / cb0
        for j in range(g.nyy-1,1,-1):
            cn.calb[i,j-1] = cn.cbm[i] * cn.cgab[i,j]
            cn.cgab[i,j-1] = -1.0 / (cb0 + cn.cbp[i] * cn.calb[i,j-1])

