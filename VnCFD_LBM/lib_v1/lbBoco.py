import numpy as np

cl=[2,6,7]
cv=[0,3,4]
cr=[1,5,8]
ct=[3,5,6]
ch=[0,1,2]
cb=[4,7,8]

def left_velocity_macro(ul, f, rho, u):
  u[0,0,:] = ul[0]
  u[1,0,:] = ul[1]
  rv = np.sum(f[cv,0,:], axis=0)
  rl = np.sum(f[cl,0,:], axis=0)
  rho[0,:] = (2*rl + rv)/(1 - ul[0])

def right_velocity_macro(ur, f, rho, u):
  u[0,-1,:] = ur[0]
  u[1,-1,:] = ur[1]
  rv = np.sum(f[cv,-1,:], axis=0)
  rr = np.sum(f[cr,-1,:], axis=0)
  rho[-1,:] = (2*rr + rv)/(1 + ur[0])

def bottom_velocity_macro(ub, f, rho, u):
  u[0,:,0] = ub[0]
  u[1,:,0] = ub[1]  
  rh = np.sum(f[ch,:,0], axis=0)
  rb = np.sum(f[cb,:,0], axis=0)
  rho[:,0] = (2*rb + rh)/(1 - ub[1])
    
def top_velocity_macro(ut, f, rho, u):
  u[0,:,-1] = ut[0]
  u[1,:,-1] = ut[1]
  rh = np.sum(f[ch,:,-1], axis=0)
  rt = np.sum(f[ct,:,-1], axis=0)
  rho[:,-1] = (2*rt + rh)/(1 + ut[1])
  
def left_pressure_macro(rhol, uyl, f, rho, u):
  rho[0,:] = rhol
  rv = np.sum(f[cv,0,:], axis=0)
  rl = np.sum(f[cl,0,:], axis=0)
  u[0,0,:] = 1.0 - (2*rl + rv)/rhol
  u[1,0,:] = uyl

def right_pressure_macro(rhor, uyr, f, rho, u):
  rho[-1,:] = rhor
  rv = np.sum(f[cv,-1,:], axis=0)
  rr = np.sum(f[cr,-1,:], axis=0)
  u[0,-1,:] = (2*rr + rv)/rhor - 1.0
  u[1,-1,:] = uyr
  
def left_f(f, rho, u):
  f[1,0,:] = f[2,0,:] + 2/3.*rho[0,:]*u[0,0,:]
  f[5,0,:] = f[7,0,:] - 1/2.*(f[3,0,:] - f[4,0,:]) + rho[0,:]*(1/6.*u[0,0,:] + 1/2.*u[1,0,:])
  f[8,0,:] = f[6,0,:] + 1/2.*(f[3,0,:] - f[4,0,:]) + rho[0,:]*(1/6.*u[0,0,:] - 1/2.*u[1,0,:])  
    
def right_f(f, rho, u):  
  f[2,-1,:] = f[1,-1,:] - 2/3.*rho[-1,:]*u[0,-1,:]
  f[7,-1,:] = f[5,-1,:] + 1/2.*(f[3,-1,:] - f[4,-1,:]) - rho[-1,:]*(1/6.*u[0,-1,:] + 1/2.*u[1,-1,:])
  f[6,-1,:] = f[8,-1,:] - 1/2.*(f[3,-1,:] - f[4,-1,:]) - rho[-1,:]*(1/6.*u[0,-1,:] - 1/2.*u[1,-1,:])  
     
def bottom_f(f, rho, u):
  f[3,:,0] = f[4,:,0] + 2/3.*rho[:,0]*u[1,:,0]
  f[5,:,0] = f[7,:,0] - 1/2.*(f[1,:,0] - f[2,:,0]) + rho[:,0]*(1/2.*u[0,:,0] + 1/6.*u[1,:,0])
  f[6,:,0] = f[8,:,0] + 1/2.*(f[1,:,0] - f[2,:,0]) - rho[:,0]*(1/2.*u[0,:,0] - 1/6.*u[1,:,0])  

def top_f(f, rho, u):  
  f[4,:,-1] = f[3,:,-1] - 2/3.*rho[:,-1]*u[1,:,-1]
  f[7,:,-1] = f[5,:,-1] + 1/2.*(f[1,:,-1] - f[2,:,-1]) - rho[:,-1]*(1/2.*u[0,:,-1] + 1/6.*u[1,:,-1])  
  f[8,:,-1] = f[6,:,-1] - 1/2.*(f[1,:,-1] - f[2,:,-1]) + rho[:,-1]*(1/2.*u[0,:,-1] - 1/6.*u[1,:,-1])  
        
def left_neq_f(f, f_eq):
  f[cr,0,:] = f_eq[cr,0,:] + f[cl,0,:] - f_eq[cl,0,:]
  
def right_neq_f(f, f_eq):
  f[cl,-1,:] = f_eq[cl,-1,:] + f[cr,-1,:] - f_eq[cr,-1,:]  
  
def bottom_neq_f(f, f_eq):
  f[ct,:,0] = f_eq[ct,:,0] + f[cb,:,0] - f_eq[cb,:,0]
  
def top_neq_f(f, f_eq):
  f[cb,:,-1] = f_eq[cb,:,-1] + f[ct,:,-1] - f_eq[ct,:,-1]
  
def left_gradient_neq_f(f, f_eq):
  f[:,0,:] = f_eq[:,0,:] + f[:,1,:] - f_eq[:,1,:]
  
def right_gradient_neq_f(f, f_eq):
  f[:,-1,:] = f_eq[:,-1,:] + f[:,-2,:] - f_eq[:,-2,:]  
  
def bottom_gradient_neq_f(f, f_eq):
  f[:,:,0] = f_eq[:,:,0] + f[:,:,1] - f_eq[:,:,1]
  
def top_gradient_neq_f(f, f_eq):
  f[:,:,-1] = f_eq[:,:,-1] + f[:,:,-2] - f_eq[:,:,-2]
  
def right_open_f(f):
  f[cl,-1,:] = f[cl,-2,:]

def bottom_open_f(f):
  f[ct,:,0] = f[ct,:,1]
  
def top_open_f(f):
  f[cb,:,-1] = f[cb,:,-2]
    
def boco_f(f, rho, u, f_eq):
  left_f(f, rho, u)
  right_f(f, rho, u)
  bottom_f(f, rho, u)
  top_f(f, rho, u)
  #corner_f(f, rho, u)  
  
def boco_neq_f(f, rho, u, f_eq):
  left_neq_f(f, f_eq)
  right_neq_f(f, f_eq)
  bottom_neq_f(f, f_eq)
  top_neq_f(f, f_eq)
  #corner_f(f, rho, u)
  
def boco_right_open_neq_f(f, rho, u, f_eq):
  left_neq_f(f, f_eq)
  right_open_f(f)
  bottom_neq_f(f, f_eq)
  top_neq_f(f, f_eq)
    
def boco_gradient_neq_f(f, rho, u, f_eq):
  left_gradient_neq_f(f, f_eq)
  right_gradient_neq_f(f, f_eq)
  bottom_gradient_neq_f(f, f_eq)
  top_gradient_neq_f(f, f_eq)
  #corner_f(f, rho, u)
  
def left_bottom_f(f, rho, u):
  f[1,0,0] = f[2,0,0] + 2/3.*rho[0,0]*u[0,0,0]
  f[3,0,0] = f[4,0,0] + 2/3.*rho[0,0]*u[1,0,0]
  f[5,0,0] = f[7,0,0] + 1/6.*rho[0,0]*(u[0,0,0] + u[1,0,0])
  f[7,0,0] = 0.0
  f[8,0,0] = 0.0
  f[0,0,0] = rho[0,0] - np.sum(f[1:,0,0])
  
def left_top_f(f, rho, u):
  f[1,0,-1] = f[2,0,-1] + 2/3.*rho[0,-1]*u[0,0,-1]
  f[4,0,-1] = f[3,0,-1] - 2/3.*rho[0,-1]*u[1,0,-1]
  f[8,0,-1] = f[6,0,-1] + 1/6.*rho[0,-1]*(u[0,0,-1] - u[1,0,-1])
  f[5,0,-1] = 0.0
  f[7,0,-1] = 0.0
  f[0,0,-1] = rho[0,-1] - np.sum(f[1:,0,-1])
    
def right_bottom_f(f, rho, u):
  f[2,-1,0] = f[1,-1,0] - 2/3.*rho[-1,0]*u[0,-1,0]
  f[3,-1,0] = f[4,-1,0] + 2/3.*rho[-1,0]*u[1,-1,0]
  f[6,-1,0] = f[8,-1,0] - 1/6.*rho[-1,0]*(u[0,-1,0] - u[1,-1,0])
  f[5,-1,0] = 0.0
  f[7,-1,0] = 0.0
  f[0,-1,0] = rho[-1,0] - np.sum(f[1:,-1,0])
  
def right_top_f(f, rho, u):
  f[2,-1,-1] = f[1,-1,-1] - 2/3.*rho[-1,-1]*u[0,-1,-1]
  f[4,-1,-1] = f[3,-1,-1] - 2/3.*rho[-1,-1]*u[1,-1,-1]
  f[7,-1,-1] = f[5,-1,-1] - 1/6.*rho[-1,-1]*(u[0,-1,-1] + u[1,-1,-1])
  f[6,-1,-1] = 0.0
  f[8,-1,-1] = 0.0
  f[0,-1,-1] = rho[-1,-1] - np.sum(f[1:,-1,-1])  
   
def corner_macro(f, rho, u):
  u[:,0,0] = u[:,1,0]
  rho[0,0] = rho[1,0]
  u[:,0,-1] = u[:,1,-1]
  rho[0,-1] = rho[1,-1]
  u[:,-1,0] = u[:,-2,0]
  rho[-1,0] = rho[-2,0]
  u[:,-1,-1] = u[:,-2,-1]
  rho[-1,-1] = rho[-2,-1]
        
def corner_f(f, rho, u):
  left_bottom_f(f, rho, u)
  left_top_f(f, rho, u)
  right_bottom_f(f, rho, u)
  right_top_f(f, rho, u)  
                        
