from numpy import pi,sqrt,arange,zeros
class GPEdata:
    def __init__(self, n=256, d=0.1, dt=0.001):
        self.nx = n ; self.ny = n #grid size
        self.dx = d ; self.dy = d #space step
        self.nxx = self.nx - 1 ; self.nyy = self.ny - 1
        self.nx2 = self.nx // 2 ; self.ny2 = self.ny // 2 
        #--- Trap Params ---#
        self.gamma = 1.0 ; self.nu = 1.0 #Trap aspect ratios γ and ν 
        self.d_z = 0.1 #Axial gaussian width 
        #--- Condensate Params ---#
        aho = 1e-6 ; bohr_a0 = 5.2917720859e-11/aho
        self.natoms = 10000 # Number of atoms N
        self.scatter_len = 3.769458264*bohr_a0 #scattering length 
        self.g2d = 100.0#4.0 * pi * self.scatter_len * self.natoms / sqrt(2 * pi * self.d_z) #Interaction strength g 
        self.omega =  0.8 #angular frequency Ω 
        self.mu = 0 #Chemical potential
        self.rms = 0 #RMS radius 
        self.en = 0 #Energy per atom
        self.g = 0 #Temporary g value
        self.lz = 0. #Angular momentum
        #--- Space arrays ---#
        xmax = (self.nx-1) // 2 * self.dx ; xmin = -(self.nx-1) // 2 * self.dx
        ymax = (self.ny-1) // 2 * self.dy ; ymin = -(self.ny-1) // 2 * self.dy
        self.x = arange(xmin, xmax+0.001,self.dx)
        self.y = arange(ymin, ymax+0.001,self.dx)
        self.x2 = self.x * self.x ; self.y2 = self.y * self.y
        #--- Other matrices ---#
        self.psi = zeros([self.nx,self.ny],dtype='complex128') #wavefunction ψ 
        self.v = zeros([self.nx,self.ny]) # potential
        self.r = zeros([self.nx,self.ny]) #sqrt(x^2 + y^2)
        #--- Other quantities ---#
        # DON'T CHANGE THESE HERE. CHANGE IN MAIN PROGRAM #
        self.dt = dt #Time step
        self.znorm = 1 #Normalization constant
        self.cij = zeros(1,dtype='complex128') #parameter for real or imaginary time
        self.xop = 2 #Option to solve different equations
        self.orm = 0 #Select 1 for  imaginary time and 2 for real time



class CNdata:
    def __init__(self,g: GPEdata):
        self.cbp = zeros(g.nx, dtype='complex128')
        self.cbm = zeros(g.nx, dtype='complex128')
        self.cap = zeros(g.ny, dtype='complex128')
        self.cam = zeros(g.ny, dtype='complex128')
        
        self.cala = zeros([g.nx,g.ny],dtype='complex128')
        self.calb = zeros([g.nx,g.ny],dtype='complex128')
        self.cgaa = zeros([g.nx,g.ny],dtype='complex128')
        self.cgab = zeros([g.nx,g.ny],dtype='complex128')

        self.cbex = zeros([g.nx,g.ny],dtype='complex128')
        self.cbey = zeros([g.nx,g.ny],dtype='complex128')

        self.ca0r = zeros(1,dtype='complex128')
        self.cb0r = zeros(1,dtype='complex128')


class Chemdata:
    def __init__(self,g: GPEdata):
        self.dpyr = zeros([g.nx,g.ny]) 
        self.dpxr = zeros([g.nx,g.ny]) 
        self.dpyi = zeros([g.nx,g.ny]) 
        self.dpxi = zeros([g.nx,g.ny]) 

        self.dplz = zeros([g.nx,g.ny]) 



        
