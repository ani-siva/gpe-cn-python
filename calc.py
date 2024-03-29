from numpy import exp,dot,sqrt,zeros_like
from util import diff,integrate
from scipy.integrate import simpson


def calc_nu(g)->None:
    p2 = (g.psi * g.psi.conjugate()).real
    tmp = g.v + g.g * p2
    g.psi *= exp(-g.cij * g.dt * tmp)
def lux(g,cn)->None:
    for j in range(0,g.ny):
        if g.orm == 1:
            cn.cbex[g.nxx-1,j] = 0.
        else:
            cn.cbex[g.nxx-1,j] = g.psi[g.nx-1,j] 
        for i in range(g.nxx-1,1,-1):
            cxx = -cn.cap[j] * g.psi[i+1,j] + cn.ca0r * g.psi[i,j] - cn.cam[j] * g.psi[i-1,j]
            cn.cbex[i-1,j] = cn.cgaa[i,j] * (cn.cap[j] * cn.cbex[i,j] - cxx)
        g.psi[0,j] = 0.0
        for i in range(0,g.nxx):
            g.psi[i+1,j] = cn.cala[i,j] * g.psi[i,j] + cn.cbex[i,j] 
        g.psi[g.nx-1,j] = 0.0
def luy(g,cn)->None:
    for i in range(0,g.nx):
        if g.orm == 1:
            cn.cbey[i,g.nyy-1] = 0.0
        else:
            cn.cbey[i,g.nyy-1] = g.psi[i,g.ny-1]
        for j in range(g.nyy-1,1,-1):
            cyy = -cn.cbp[i] * g.psi[i,j+1] + cn.cb0r * g.psi[i,j] - cn.cbm[i] * g.psi[i,j-1]
            cn.cbey[i,j-1] = cn.cgab[i,j] * (cn.cbp[i] * cn.cbey[i,j] - cyy)
        g.psi[i,0] = 0.0
        for j in range(0,g.nyy):
            g.psi[i,j+1] = cn.calb[i,j] * g.psi[i,j] + cn.cbey[i,j]
        g.psi[i,g.ny-1] = 0.0

def rad(g)->None:
    psi2 = (g.psi * g.psi.conjugate()).real
    g.rms = sqrt(integrate(g.r * g.r * psi2,g.dx,g.dy))

def renorm(g)->None:
    psi2 = (g.psi * g.psi.conjugate()).real
    g.znorm = sqrt(integrate(psi2,g.dx,g.dy))
    if g.orm == 1:
        g.psi /= g.znorm

def chem(g,cd):
    psi_c = g.psi.imag
    psi_r = g.psi.real
    psi2 = (g.psi * g.psi.conjugate()).real
    for i in range(0,g.nx):
        cd.dpyr[i,:] = diff(psi_r[i,:],g.dy)
        cd.dpyi[i,:] = diff(psi_c[i,:],g.dy)
    for j in range(0,g.ny):
        cd.dpxr[:,j] = diff(psi_r[:,j],g.dx)
        cd.dpxi[:,j] = diff(psi_c[:,j],g.dx)
    dp2 = cd.dpxr * cd.dpxr + cd.dpyr * cd.dpyr

    for j in range(0,g.ny):
        for i in range(0,g.nx):
            cd.dplz[i,j] = psi_r[i,j] * (g.x[i] * cd.dpyi[i,j]  - g.y[j] * cd.dpxi[i,j]) 
    p2r = psi_r * psi_r
    gp2 = g.g * psi2
    tmp2d = (g.v + gp2) * p2r + dp2 - g.xop * (g.omega * cd.dplz)
    emp2d = (g.v + gp2/2.0) * p2r + dp2 - g.xop * (g.omega * cd.dplz)
    znorm = integrate(psi_r * psi_r, g.dx,g.dy)
    g.mu = integrate(tmp2d,g.dx,g.dy) / znorm
    g.en = integrate(emp2d,g.dx,g.dy) / znorm



