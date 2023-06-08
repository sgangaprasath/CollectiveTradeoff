import numpy as np
import scipy as sp
# Parameters in the model
class params():
    def __init__(self) -> None:
        super().__init__()
        self.nSteps = 1500 #number of steps of each run
        self.nAgnts = 10 #number of agents
        self.nForag = 1 #number of foragers
        # self.nEpchs = 250 #number of epochs
        self.a = 1.0 #radius of semi-circle
        self.nSc = int(0.1*self.nSteps) #number of steps of intrinsic policy
        
        # Pheromone parameters
        self.nMsh = 500 #number of points in 
        self.kM = 1 #rate of pheromone decay
        self.kP = 0.1 #rate of pheromone generation
        
        # Agent dynamics parameters
        self.pt = 0.08 #initial pheromone trail thickness
        self.sz = 0.01 #agent size (radius)
        self.dt = 1e-2 #time-step size
        self.l = 0.5e-2 #length travelled in 1 time-step
        self.vo = self.l/self.dt #effective speed of motion
        self.nu = 50 #orientation relaxation rate (s^-1)
        self.diff = 50e-1 #noise diffusion coefficient (s^-1)
        self.diffCst = 50 #noise diffusion coefficient (s^-1)
        
        # Reward/Learning parameters
        self.alpha = 0.7 #learning rate
        self.phiSt = np.pi/8 #reward goes as exp(-\phi/\phiSt)
        self.sigma = 0.1 #radius of region near target within which agent succeeds
        self.nPtn = 40 #number of division of \phi over which value function is define
        self.epsilon = 0.9 #constant for epsilon-greedy strategy
        
        # Initial conditions for workers
        phiInit = np.random.uniform(0, np.pi, self.nAgnts)
        self.rInit = self.a*np.array([np.cos(phiInit), np.sin(phiInit)])
        self.tgt = np.array([-self.a, 0.0]) #target location of agent
        self.strt = np.array([self.a, 0.0]) #initial location of agent
        self.thetInit = np.random.uniform(0, 2*np.pi, self.nAgnts)
        # Initial conditions for foragers
        self.rFgInit = self.strt
        self.thetFgInit = np.pi*0.5