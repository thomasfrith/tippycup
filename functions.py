import matplotlib.pyplot as plt
import math
import numpy as np
import pandas as pd

def slide_bottle(r,u,h,MU):
    # define constants of the bottle
    
    L = 0.3 # height of bottle
    RHO_B = 0.01 # mass per unit length for bottle
    RHO_L = 0.1 # mass per unit area for the liquid
    M_B = RHO_B*2*L*(1+r) # mass of bottle
    M_L = RHO_L*h*r*(L**2) # mass of liquid
    M = M_B+M_L # total mass
    X = L*(M_L*h + M_B)/(M_L+M_B) # 2 times the height of the COM (to make equations work nicely)
    
    I_B = 2*(r**3)*(L**3)*RHO_B/12 + r*L*RHO_B*(X/2)**2 + r*L*RHO_B*(L-X/2)**2 + 2*(L**3)*RHO_B/12 + 2*L*RHO_B*((r*L)**2 + (L/2-X/2)**2) # moment of inertia of bottle about COM
    I_L = (1/12)*M_L*(L**2)*(r**2+h**2) + M_L*((X/2)-(h*L/2))**2 # moment of inertia of liquid about COM
    I_T = I_B + I_L # total moment of inertia about COM
    
    g = 10 # gravity
    phi = math.atan2(X,(L*r)) # starting angle of the COM
    alpha = math.atan(1/MU)
    
    # define constants of the simulation
    
    GRID_LENGTH = 1000 # length of simulation
    dt = 0.001 # time increment
    
    # define the lists for the differential equation
    
    th = [0 for i in range(GRID_LENGTH)] # angle of rotation
    d_th = [0 for i in range(GRID_LENGTH)]
    v = [0 for i in range(GRID_LENGTH)] # linear velocity
    t = [i*dt for i in range(GRID_LENGTH)] # time
    tipped = False  # has it tipped 
    
    # define initial conditions
    
    th_i = phi # initial angle
    th_dot_i = 0 # initial angular velocity
    
    th[0] = th_i
    d_th[1] = th_dot_i
    
    v[0] = u
    v[1] = u
    
    # backwards difference method
    
    for i in range(GRID_LENGTH-1):
        
        d_th[i+1] = d_th[i] + dt*((M*g-(0.5)*M*math.sqrt(X**2+(r*L)**2)*( - math.sin(th[i])*((d_th[i]))**2))*(0.5)*(math.sqrt((X**2+(r*L)**2)*(1+MU**2)))*(math.sin(th[i]-alpha)))\
        /(I_T + (0.5)*M*math.sqrt(X**2+(r*L)**2)*math.cos(th[i])*(0.5)*(math.sqrt((X**2+(r*L)**2)*(1+MU**2)))*(math.sin(th[i]-alpha)))
        
        th[i+1] = th[i]+dt*d_th[i+1]
        
        if v[i]<=0:
            v[i+1]=0
            alpha = math.pi/2
            MU = 0
        else:
            v[i+1] = u-MU*g*t[i]+MU*(0.5)*math.sqrt(X**2+(r*L)**2)*math.cos(th[i])*(d_th[i+1])
        
        if th[i+1]>=math.pi/2:
            th[i+1] = math.pi/2
            tipped = True
            break
        if th[i+1]<=phi:
            th[i+1] = phi
            break
        
    df = pd.DataFrame(list(zip(th,v)),columns=["theta","velocity"],index=t)
    
    return df,tipped

#df,tipped=slide_bottle(0.25,0.8,0.4,0.3)
#plt.plot(df.theta)
#plt.plot(df.velocity)

