# =============================================================================
# Tippy Cup simulation models a bottle sliding across a table. The conditions 
# of whether it tips or not due to the friction from the table is the objective
# of the simulation.
# 2D
# v1 for a block of dimensions L and rL
# v2 for a bottle filled to a height h of liquid
# =============================================================================

import matplotlib.pyplot as plt
import pandas as pd

from functions import slide_bottle

# variable to change

MU = 0.9

RANGE_u = 20 # number of initial velocity values to check
RANGE_h = 100 # number of heights of liquid in the bottle as a fraction of total height (MAX 1)
RANGE_r = 100 # number of heights of liquid in the bottle as a fraction of total height (MAX 1)

MAX_u = 5

du = MAX_u/RANGE_u
dh = 1/RANGE_h
dr = MU/RANGE_r


#
#df,tipped=slide_bottle(0.2,0.8,0.4,0.3)
#plt.plot(df.theta)
#plt.plot(df.velocity)

#u = []
#h = []
#stable = []
#
#for i in range(RANGE_u):
#    for j in range(RANGE_h):
#        u.append(i*du)
#        h.append(j*dh)
#        df,tipped = slide_bottle(r,i*du,j*du)
#        if tipped == True:
#            stable.append(False)
#        else:
#            stable.append(True)
#            
#df = pd.DataFrame(list(zip(u,h,stable)),columns=["u","h","stable"])
#df.to_json("parameter_space_" + str(RANGE_u) + "_" + str(RANGE_h) + ".json")
#
#
#fig = plt.figure()
#ax1 = fig.add_subplot(111)
#
#ax1.scatter(df[df.stable==False].h,df[df.stable==False].u, label='Unstable',color = "red")
#ax1.scatter(df[df.stable==True].h,df[df.stable==True].u, label='Stable', color = "blue")
#plt.legend(loc='upper right')
#ax1.set_xlabel("Height of liquid as a fraction of the total height of the bottle")
#ax1.set_ylabel("Initial velocity in m/s")
#ax1.set_title("width/height="+str(r) +" coefficient of friction=0.3")
#plt.show()

u = [[None for i in range(RANGE_r)] for i in range(RANGE_h)]
h = [i*dh for i in range(RANGE_h)]
r = [i*dr for i in range(RANGE_r)]


for i in range(RANGE_h):
    for j in range(RANGE_r):
        for k in range(RANGE_u):
            if j*dr<MU:
                df,tipped = slide_bottle(j*dr,k*du,i*dh,MU)
                if tipped == True:
                    u[i][j] = k*du
                    break
            
df = pd.DataFrame(list(zip(u,h,r)),columns=["u","h","r"])
df.to_json("parameter_space_" + str(RANGE_r) + "_" + str(RANGE_h) + ".json")

fig,ax1=plt.subplots(1,1)
cp = ax1.contourf(r, h, u)

cbar = plt.colorbar(cp)
cbar.set_label("Velocity needed to tip bottle", labelpad = 10, rotation=270)

ax1.set_xlabel("ratio of height to width of bottle")
ax1.set_xlim(right=MU)
ax1.set_ylabel("ratio of height of liquid to height of bottle")
ax1.set_title("Coefficient of friction="+str(MU))
plt.show()



fig.savefig("parameter_space_" + str(MU) + "_" + str(RANGE_r) + "_" + str(RANGE_h) + ".pdf",bbox_inches='tight')


