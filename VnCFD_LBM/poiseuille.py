import sys
sys.path.insert(1,'lib_v1')

from lbBoco import *
from lbSolver import *

def boco_macro(f, rho, u):
  left_pressure_macro(rho_left, 0.0, f, rho, u)
  right_pressure_macro(rho_right, 0.0, f, rho, u)
  bottom_velocity_macro([0,0], f, rho, u)
  top_velocity_macro([0,0], f, rho, u)
  #corner_macro(f, rho, u)
    
def plot_velocity():
  ux = np.load('result/ux.npy')
  uy = np.load('result/uy.npy')

  y = np.arange(Ny)
  dp = drho/3 # p=rho*cs^2, cs=1/sqrt(3)
  ux_exact = -dp*y*(y-H)/(2*L*nu*rho0)
  uy_exact = np.zeros(Ny)
  
  fig, ax = plt.subplots(1,2)
  ax[0].plot(ux[Nx//4,:], y, ux[Nx//2,:], y, ux[3*Nx//4,:], y, ux_exact, y)
  ax[1].plot(uy[Nx//4,:], y, uy[Nx//2,:], y, uy[3*Nx//4,:], y, uy_exact, y)
  ax[0].legend(['Nx/4', 'Nx/2', '3Nx/4', 'exact'])
  ax[0].grid(True)
  ax[1].legend(['Nx/4', 'Nx/2', '3Nx/4', 'exact'])
  ax[1].grid(True)
  plt.savefig('result/poiseuille_velocity.png', bbox_inches='tight')    
  plt.show()
    
    
if __name__ == "__main__":
  Nx, Ny = 101, 41
  L, H = Nx-1, Ny-1
  rho0, drho = 1.0, 0.03
  u0 = [0.0, 0.0]
  rho_left = rho0 + 0.5*drho
  rho_right = rho0 - 0.5*drho
  Re = 20.0
  
  # u(x,y, t) = dp*y*(H-y)/(2*mu*L), dp=drho/3 
  # u_max = drho*H^2/(24*nu*rho*L)
  # Re = u_max*L/nu = drho*H^2/(24*nu^2*rho)
  nu = np.sqrt(drho*H*H/(24*rho0*Re))
  tau = 3*nu+0.5
  u_max = drho*H*H/(24*nu*rho0*L)
  print('Re: ', u_max*L/nu, 'M:  ', u_max*np.sqrt(3))
  print('nu: ', nu, 'u_max: ', u_max, ', tau: ', tau)
  
  #rho0 = np.load('result/rho.npy')
  #u0 = [np.load('result/ux.npy'), np.load('result/uy.npy')]
  
  solver(Nx, Ny, rho0, u0, tau, 10001, boco_macro, boco_neq_f, plot_field='velocity')
  plot_velocity()
