import numpy as np
  
class lbNodes:
  def __init__(self, Nx, Ny, Obstacle=None):
    self.Nx = Nx
    self.Ny = Ny
    self.e = np.array([[0,0], [1,0], [-1,0], [0,1], [0,-1], [1,1], [-1,1], [-1,-1],  [1,-1]])
    self.w = np.array([4/9,1/9,1/9,1/9,1/9,1/36,1/36,1/36,1/36])
    self.f = np.ones((9, Nx, Ny))
    self.f_eq = np.ones((9, Nx, Ny))    
    self.u = np.zeros((2, Nx, Ny))
    self.rho = np.ones((Nx, Ny))
    self.obstacle = Obstacle
    self.bound_f = None
    
  def init(self, tau, rho, u):
    self.tau = tau
    self.nu = (tau-0.5)/3
    self.omega = 1/tau
    self.rho[:] = rho
    self.u[0] = u[0]
    self.u[1] = u[1]
    self.equilibrium(self.f)
                   
  def macro_vars(self):
    self.rho[:,:] = np.sum(self.f, axis=0)
    self.u[0,:,:] = np.tensordot(self.e[:,0], self.f, axes=(0,0))/self.rho
    self.u[1,:,:] = np.tensordot(self.e[:,1], self.f, axes=(0,0))/self.rho    
    
  def equilibrium(self, f_eq):
    u2 = 1.5*(self.u[0]**2 + self.u[1]**2)
    for i in range(9):
      eu = 3*(self.e[i,0]*self.u[0] + self.e[i,1]*self.u[1])
      f_eq[i] = self.rho*self.w[i]*(1 + eu + 0.5*eu**2 - u2)
          
  def collision(self):
    self.f -= self.omega*(self.f - self.f_eq)
        
  def streaming(self):
    self.f[1,1:,:] = self.f[1,:-1,:]
    self.f[2,:-1,:] = self.f[2,1:,:]
    self.f[3,:,1:] = self.f[3,:,:-1]
    self.f[4,:,:-1] = self.f[4,:,1:]
    self.f[5,1:,1:] = self.f[5,:-1,:-1]
    self.f[6,:-1,1:] = self.f[6,1:,:-1]
    self.f[7,:-1,:-1] = self.f[7,1:,1:]
    self.f[8,1:,:-1] = self.f[8,:-1,1:]     
#    for i in range(1,9):
#      self.f[i] = np.roll(self.f[i], self.e[i,0], axis=0)
#      self.f[i] = np.roll(self.f[i], self.e[i,1], axis=1)
          
  def bounce_back(self):
    self.f[:,self.obstacle] = self.bound_f[[0,2,1,4,3,7,8,5,6],:] 

  def set_boco(self, boco_macro, boco_f):
    self.boco_macro = boco_macro
    self.boco_f = boco_f
    
  def iteration(self):
    self.macro_vars()
    if self.obstacle is not None: self.u[:, self.obstacle] = 0.0
    self.boco_macro(self.f, self.rho, self.u)
    self.equilibrium(self.f_eq)
    if self.obstacle is not None: self.bound_f = self.f[:, self.obstacle]
    self.collision()
    if self.obstacle is not None: self.bounce_back()
    self.boco_f(self.f, self.rho, self.u, self.f_eq)
    self.streaming()

