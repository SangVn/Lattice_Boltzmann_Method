import sys
sys.path.insert(1,'lib_v1')

from lbBoco import *
from lbSolver import *

def cylinder(Nx, Ny, cx, cy, r):
  Y, X = np.meshgrid(range(Ny), range(Nx))
  wall = (X-cx)**2 + (Y-cy)**2 <= r**2
  return wall
  
def boco_macro(f, rho, u):
  left_velocity_macro(u0_left, f, rho, u)
  bottom_velocity_macro([0,0], f, rho, u)
  top_velocity_macro([0,0], f, rho, u)  
  
def plot_streamline():
  ux = np.load('result/ux.npy')
  uy = np.load('result/uy.npy')
  
  L, H = 2.2, 0.41
  Nx = ux.shape[0]
  Ny = uy.shape[1]
  x = np.linspace(0, L, Nx)
  y, dy = np.linspace(0, H, Ny, retstep=True)
  xr = 1.0
  Nr = int(xr/dy)

  fig, ax = plt.subplots(figsize=(10,10*H/xr))
  art = plt.Circle((0.2, 0.2), radius=0.05, linewidth=2, color='r', fill=True)
  ax.add_artist(art)
 
  plt.streamplot(x[:Nr], y, ux[:Nr].T, uy[:Nr].T, density=[4, 2])
  plt.xlim(0, xr)
  plt.ylim(0, H)
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)  
  plt.savefig('result/cylinder_streamline_Re20.png', bbox_inches='tight')
  plt.show()  
      
if __name__ == "__main__":
  '''
  1996, M. Schafer, Benchmark Computations of Laminar Flow Around a Cylinder
  rho = 1.0 kg/m^3, nu = 10e-3 m^2/2
  U(0,y) = 4Um*y(H-y)/H^2, V = 0
  Re = Ua*D/nu, Ua = 2U(0,H/2)/3 = 2Um/3 -> Re = 2Um*D/3nu
  latice Boltzmann nondimensional velocity: u_lb = U/U0, U0 - the reference speed
  Let U0 = 10 -> Um_lb = Um/U0
  '''
  L, H, D = 2.2, 0.41, 0.1 #(m)

  # a) Um = 0.3 m/s, Re = 20
  #Re, U_max = 20, 0.3

  # b) Um = 1.5 m/s, Re = 100
  Re, U_max = 100, 1.5
 
  Ny = 165  
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
  Cyl = cylinder(Nx, Ny, int(0.2/dy), int(0.2/dy), Nd/2)
  
  rho0, u0 = 1.0, [0.0, 0.0]
  #rho0 = np.load('result/rho.npy')
  #u0 = [np.load('result/ux.npy'), np.load('result/uy.npy')]
  
  solver(Nx, Ny, rho0, u0, tau, 501, boco_macro, boco_right_open_neq_f, obstacle=Cyl, plot_field='vorticity')
  plot_streamline()
