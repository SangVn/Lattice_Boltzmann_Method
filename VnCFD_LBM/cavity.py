import sys
sys.path.insert(1,'lib_v1')

from lbBoco import *
from lbSolver import *

def boco_macro_1(f, rho, u):
  left_velocity_macro([0,0], f, rho, u)
  right_velocity_macro([0,0], f, rho, u)
  bottom_velocity_macro([0,0], f, rho, u)
  top_velocity_macro([ux_top,0], f, rho, u)
  #corner_macro(f, rho, u)
  
# u = const, gradient_p = 0
def boco_macro_2(f, rho, u):
  u[:,0,:] = 0.0
  u[:,-1,:] = 0.0
  u[:,:,0] = 0.0
  u[0,:,-1] = ux_top
  u[1,:,-1] = 0.0
  rho[0,:] = rho[1,:]
  rho[-1,:] = rho[-2,:]
  rho[:,0] = rho[:,1]
  rho[:,-1] = rho[:,-2]

def plot_streamline():
  ux = np.load('result/ux.npy')
  uy = np.load('result/uy.npy')
  V = (ux**2+uy**2)**0.5
  
  Nx = ux.shape[0]
  Ny = uy.shape[1]
  x = np.linspace(0, L, Nx)
  y= np.linspace(0, H, Ny)

  fig, ax = plt.subplots(figsize=(6,6))
  shw = plt.imshow(V.T, cmap='viridis')      
  plt.streamplot(x, y, ux.T, uy.T, density=[2, 2])
  plt.xlim(0, Nx)
  plt.ylim(0, Ny)  
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)  
  plt.savefig('result/cavity_Re1000.png', bbox_inches='tight')
  plt.show()  
            
if __name__ == "__main__":
  Re, ux_top = 1000, 0.1
  Nx, Ny = 101, 101
  L, H = Nx-1, Ny-1
  tau = 3*ux_top*L/Re + 0.5
  print('Re: ', Re, ', Mach: ', ux_top*3**0.5, ', tau: ', tau) # cs = 1/sqrt(3)
  
  rho0, u0 = 1.0, [0.0, 0.0]
  #rho0 = np.load('result/rho.npy')
  #u0 = [np.load('result/ux.npy'), np.load('result/uy.npy')] 
  
  # You guys check the combination of two set [boco_macro_1, boco_macro_2] and [boco_f, boco_neq_f, boco_gradient_neq_f]
  # with Re = [100, 1000, 4000]
  solver(Nx, Ny, rho0, u0, tau, 50001, boco_macro_1, boco_neq_f, plot_field='velocity')  
  plot_streamline()
  
