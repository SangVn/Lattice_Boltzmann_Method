import sys
sys.path.insert(1,'lib_v1')

from lbBoco import *
from lbSolver import *
import matplotlib.path as mplPath

def boco_macro(f, rho, u):
  left_velocity_macro(u0_left, f, rho, u)
  bottom_velocity_macro([0,0], f, rho, u)
  top_velocity_macro([0,0], f, rho, u) 
  
def art(Nx, Ny, L, H):
  art = np.loadtxt(fname)
  poly_path = mplPath.Path(art)
  points = np.zeros((Nx,Ny, 2))
  x = np.linspace(0, L, Nx)
  y = np.linspace(0, H, Ny)
  for i in range(Nx):
    for j in range(Ny):
      points[i,j,0] = x[i]
      points[i,j,1] = y[j]

  list_points = points.reshape(Nx*Ny,2)
  is_inside = poly_path.contains_points(list_points, radius=-0.0).reshape(Nx, Ny)
  return is_inside

def plot_art():
  art = np.loadtxt(fname)
  art = np.concatenate((art, [art[0]]), axis=0)  
  dy = H/(Ny-1)
  art /= dy
  plt.plot(art[:,0], art[:,1], 'm')
           
if __name__ == "__main__":
  L, H, D = 8.0, 2.0, 1.0 #(m)
  Re, U_max = 200, 1.0 
  Ny = 161
  y, dy = np.linspace(0, H, Ny, retstep=True)
  Nx = int(L/dy)+1
  Nd = int(D/dy)
  u_max = U_max/10  
  u_avg = 2.*u_max/3.  
  U = 4*u_max*y*(H-y)/H**2 
  nu = u_avg*Nd/Re
  tau = 3*nu+0.5
  
  print('Re: ', Re, ', Mach: ', u_max*np.sqrt(3)) # cs = 1/sqrt(3)
  print('nu: ', nu, 'u_max: ', u_max, ', tau: ', tau)

  u0_left = [U, 0]
  fname = 'result/face.txt'
  #fname = 'result/whale.txt'
  obs = art(Nx, Ny, L, H)
  
  rho0, u0 = 1.0, [0.0, 0.0]
  #rho0 = np.load('result/rho.npy')
  #u0 = [np.load('result/ux.npy'), np.load('result/uy.npy')]
  
  solver(Nx, Ny, rho0, u0, tau, 5001, boco_macro, boco_right_open_neq_f, obstacle=obs, plot_field='vorticity', plot_obs=plot_art)
