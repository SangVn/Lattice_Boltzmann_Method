import sys
sys.path.insert(1,'lib_v1')

from lbBoco import *
from lbSolver import *
import matplotlib.path as mplPath

def boco_macro(f, rho, u):
  left_velocity_macro(u0_left, f, rho, u)
  bottom_velocity_macro([0,0], f, rho, u)
  top_velocity_macro([0,0], f, rho, u) 
  
def naca23012(Nx, Ny, L, H):
  naca23012 = np.loadtxt('result/naca23012.txt')
  poly_path = mplPath.Path(naca23012)
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
         
def plot_streamline():
  naca23012 = np.loadtxt('result/naca23012.txt')
  ux = np.load('result/ux.npy')
  uy = np.load('result/uy.npy')
  L, H = 6.0, 2.0
  Nx = ux.shape[0]
  Ny = uy.shape[1]
  x = np.linspace(0,L,Nx)
  y = np.linspace(0,H,Ny)

  fig, ax = plt.subplots(figsize=(5,5))
  plt.plot(naca23012[:,0], naca23012[:,1], 'r')
  plt.streamplot(x[60:250], y, ux[60:250].T, uy[60:250].T, density=[2, 4])
  plt.ylim(0,H)
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)  
  plt.savefig('result/naca_streamline.png', bbox_inches='tight')
  plt.show()
           
if __name__ == "__main__":
  '''
  U(0,y) = 4Um*y(H-y)/H^2, V = 0
  Re = Ua*D/nu, Ua = 2U(0,H/2)/3 = 2Um/3 -> Re = 2Um*D/3nu
  latice Boltzmann nondimensional velocity: u_lb = U/U0, U0 - the reference speed
  Let U0 = 10 -> Um_lb = Um/U0
  '''
  L, H, D = 6.0, 2.0, 1.0 #(m)

  # a) Um = 0.3 m/s, Re = 20
  #Re, U_max = 20, 0.3

  # b) Um = 1.5 m/s, Re = 400
  Re, U_max = 400, 1.5
 
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
  naca = naca23012(Nx, Ny, L, H)
  
  rho0, u0 = 1.0, [0.0, 0.0]
  #rho0 = np.load('result/rho.npy')
  #u0 = [np.load('result/ux.npy'), np.load('result/uy.npy')]
  
  solver(Nx, Ny, rho0, u0, tau, 5001, boco_macro, boco_right_open_neq_f, obstacle=naca, plot_field='vorticity')
  plot_streamline()
