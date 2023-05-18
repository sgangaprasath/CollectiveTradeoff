import numpy as np
import scipy as sp
# Parameters in the model
class params():
    def __init__(self):
        super().__init__()
        self.nSteps = 800 #number of steps of each run
        self.nEpchs = 10000 #number of epochs
        self.a = 1.0 #radius of semi-circle
        self.nSc = int(0.1*self.nSteps) #number of steps of intrinsic policy
        
        # Pheromone parameters
        self.nMsh = 500 #number of points in
        self.kM = 1 #rate of pheromone decay
        self.kP = 0.1 #rate of pheromone generation
        
        # Agent dynamics parameters
        self.pt = 0.04 #initial pheromone trail half-thickness
        self.sz = 0.02 #agent size (radius)
        self.dt = 1e-2 #time-step size
        self.l = 5e-3 #length travelled in 1 time-step
        self.vo = self.l/self.dt #effective speed of motion
        self.nu = 50 #orientation relaxation rate (s^-1)
        # self.diff = 50e-1 #noise diffusion coefficient (s^-1)
        # self.diff = 10/(self.nSc*self.l) #noise diffusion coefficient (s^-1)
        self.diffCst = 50 #noise diffusion coefficient (s^-1)
        self.diff = self.diffCst #dynamic diffusion coefficient (s^-1)
        self.pe = self.diff/self.nu #peclet number
        
        # Reward/Learning parameters
        self.alpha = 0.7 #learning rate
        self.phiSt = np.pi/8 #reward goes as exp(-\phi/\phiSt)
        self.sigma = 0.1 #radius of region near target within which agent succeeds
        self.nPtn = 40 #number of division of \phi over which value function is define
        self.epsilon = 0.9 #constant for epsilon-greedy strategy
        
        # Initial conditions
        self.rInit = np.array([self.a + np.random.uniform(-self.pt, self.pt), 0.0]) #initial location of agent
        self.tgt = np.array([-self.a, 0.0]) #target location of agent
        self.thetInit = 0.5*np.pi #initial orientation of agent
        # self.rng = random.PRNGKey(1) #initial random number
        self.tgtPhi = 0.0 #target orientation from current location