import time
from numpy import dot,allclose,zeros,stack 
from data import GPEdata,CNdata,Chemdata
from init import initialize,calc_trap,coef
from util import write_psi
from calc import calc_nu,lux,luy,rad,renorm,chem



def main():
    nstp = 1
    npas = 10
    snap = 10
    it = npas / snap
    nrun = 1
    g = GPEdata(n=257,d=0.05,dt=0.001)
    cn = CNdata(g)
    cd = Chemdata(g)
    
    g.xop = 2 #Select to solve equation with double the time  
    g.orm = 1 #Select 1 for imaginary time and 2 for real time

    initialize(g,nstp)
    if g.orm == 1:
        outfile = open('im2d-out.txt','w')
        outfile.write('# Imaginary time propagation\n')
        g.cij = 1 #Imaginary time iterations 
    elif g.orm == 2:
        g.cij = 1j #Real time iterations
        outfile = open('re2d-out.txt')
        outfile.write('# Real time propagation\n')
    else:
        print('Invalid option for real or imaginary time iteration')
    outfile.write('# Option = %d\n'%g.xop)
    outfile.write('# Number of atoms = %d\n'%g.natoms)
    outfile.write('# Nonlineariy strength g2d = %3.4f\n'%g.g2d)
    outfile.write('# Trap aspect ratios gamma = %1.4f, nu = %1.4f\n'%(g.gamma,g.nu))
    outfile.write('# Rotation frequency = %1.4f\n'%g.omega)
    outfile.write('#\n')
    outfile.write('# Space grids nx = %d, ny = %d\n'%(g.nx,g.ny))
    outfile.write('# Space step dx = %1.4f, dy = %1.4f\n'%(g.dx,g.dy)) 
    outfile.write('# Time step dt = %1.4f\n'%(g.dt))

    calc_trap(g)
    coef(g,cn)

    dashe = 74 * '-'
    space = 10 * ' '
    space2 = 7 * ' '
    spaceinit = 2 * ' '
    spacenstp = 5 * ' '
    spacenpas = 11 * ' '
    outfile.write('#'+space+dashe+'\n')
    outfile.write('#'+space+'Iter'+space2+'Time'+space2+'Norm'+space2+'Chem'+space2+'Ener'+space2+'Radi'+space2+'Amom\n')
    outfile.write('#'+space+dashe+'\n')

    

    t = 0.0
    step = 0
    #---- NSTP iterations ----#
    if nstp!=0:
        gstp = g.xop * g.g2d / float(nstp)
        g.g = 0.
        renorm(g)
        chem(g,cd)
        rad(g)
        str = '#Initial:'+spaceinit+'{0:4d}    {1:9.4f} {2:9.4f}    {3:9.4f}  {4:9.4f}{5:9.4f}  {6:9.4f}\n'
        outfile.write(str.format(step,t,g.znorm,g.mu,g.en,g.rms,g.omega))
        print(g.psi.shape)
        if g.orm == 1:
            write_psi('im2d-initial-wave.txt',g.psi)
        else:
            write_psi('re2d-initial-wave.txt',g.psi)
        for i in range(nstp):
            g.g += gstp
            calc_nu(g)
            lux(g,cn)
            luy(g,cn)
            if g.orm == 1:
                renorm(g)
        if g.orm == 2:
            renorm(g)
        chem(g,cd)
        rad(g)
        str = '#NSTP:'+spacenstp+'{0:4d}    {1:9.4f} {2:9.4f}    {3:9.4f}  {4:9.4f}{5:9.4f}  {6:9.4f}\n'
        outfile.write(str.format(step,t,g.znorm,g.mu,g.en,g.rms,g.omega))
    else:
        g.g = g.xop * g.g2d

    #---- NPAS iterations ----#
    step = 1
    for k in range(npas):
        calc_nu(g)
        lux(g,cn)
        luy(g,cn)
        if g.orm == 1:
            renorm(g)
        if k%it == 0:
            if g.orm == 1:
                filename = 'im2d-wave-%d.txt'%step
            else:
                filename = 're2d-wave-%d.txt'%step
            write_psi(filename,g.psi)
            chem(g,cd)
            rad(g)
            if g.orm == 2:
                renorm(g)
            str = spacenpas+'{0:4d}    {1:9.4f} {2:9.4f}    {3:9.4f}  {4:9.4f}{5:9.4f}  {6:9.4f}\n'
            outfile.write(str.format(step,t,g.znorm,g.mu,g.en,g.rms,g.omega))
            step = step + 1
    
        t = t + g.dt
    #---- NRUN iteration -----#
    if nrun!=0:
        for k in range(nrun):
            calc_nu(g)
            lux(g,cn)
            luy(g,cn)
            if g.orm == 1:
                renorm(g)
        chem(g,cd)
        rad(g)
        if g.orm == 2:
            renorm(g)

        str = '#NRUN:'+spacenstp+'{0:4d}    {1:9.4f} {2:9.4f}    {3:9.4f}  {4:9.4f}{5:9.4f}  {6:9.4f}\n'
        outfile.write(str.format(step,t,g.znorm,g.mu,g.en,g.rms,g.omega))




    outfile.close()
if __name__=="__main__":
    timefile = open('time.txt','w')
    start = time.process_time()
    main()
    timefile.write('Clock time: ', time.process_time()-start,'seconds')
    timefile.close()
