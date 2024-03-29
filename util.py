from numpy import savetxt,reshape,stack,zeros_like
from scipy.integrate import simpson

def write_psi(file,psi):
    wf = psi.reshape(-1)
    u = stack((wf.real,wf.imag),axis=-1)
    savetxt(file,u,fmt='%1.8e %1.8e')

def integrate(p,dx,dy):
    tmp1 = simpson(p,dx = dx,axis=0)
    res = simpson(tmp1,dx = dy)
    return res

def diff(p,dx):
    n = len(p)
    dp = zeros_like(p)
    dp[0] = 0.
    dp[1] = 0.5 * (p[2]-p[0]) / dx
    for i in range(2,n-2):
        dp[i] = (p[i-2] - 8. * p[i-1] + 8. * p[i+1] - p[i+2])/(12. * dx)
    dp[-2] = 0.5 * (p[-1] - p[-3]) / dx
    dp[-1] = 0.
    return dp
