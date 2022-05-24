import time
from tqdm import tqdm
from lbData import *
import matplotlib.pyplot as plt

def show_field(block, field='velocity', plot_obstacle=None):
  plt.cla()
  if plot_obstacle is not None: plot_obstacle()    
  if field == 'rho':
    shw = plt.imshow(block.rho.T, cmap='viridis')      
  elif field == 'velocity':
    V = (block.u[0]**2+block.u[1]**2)**0.5
    shw = plt.imshow(V.T, cmap='viridis')#, aspect='auto')     
  else:# field=='vorticity':
    vorticity = (np.roll(block.u[1],-1,axis=0)-np.roll(block.u[1],1,axis=0)) - \
                (np.roll(block.u[0],-1,axis=1)-np.roll(block.u[0],1,axis=1))
    if block.obstacle is not None: vorticity[block.obstacle] = np.nan
    shw = plt.imshow(vorticity.T, norm=plt.Normalize(-0.025, 0.025), cmap='jet')
  ax = plt.gca()
  ax.invert_yaxis()
  ax.get_xaxis().set_visible(False)
  ax.get_yaxis().set_visible(False)
  plt.pause(0.001)
  return shw

def solver(Nx, Ny, rho0, u0, tau, Nt, boco_macro, boco_f, obstacle=None, plot_field='velocity', plot_obs=None, fname='result.png'):
  t0 = time.time()
  show = plot_field is not None
  
  block = lbNodes(Nx, Ny, obstacle)
  block.init(tau, rho0, u0)
  block.set_boco(boco_macro, boco_f)
  
  if show: plt.figure(figsize=(8,8/(Nx/Ny)))
  for n in tqdm(range(Nt)):
    block.iteration()
    if(show and n%100 == 0): 
      shw = show_field(block, field=plot_field, plot_obstacle=plot_obs)
      # save images to create gif file by ImageMagick:
      # convert -delay 20 $(ls -v r_*.png) name.gif
      #plt.savefig('result/r_'+str(n)+'.png', bbox_inches='tight')
      
  print("Time to execute: ", time.time()-t0)
  np.save('result/ux.npy', block.u[0])
  np.save('result/uy.npy', block.u[1])
  np.save('result/rho.npy', block.rho)
   
  if show:
    plt.savefig('result/'+fname, bbox_inches='tight')    
    plt.show()
