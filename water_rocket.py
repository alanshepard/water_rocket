from scipy.integrate import ode
from math import sqrt, pi, cos, sin
import numpy as np
import matplotlib.pyplot as plt

#Using SI units all around.

#### CONSTANTS ####
rhoWater=998
rhoAir=1.2754
gamma=1.4
nozzleArea=22e-3**2*3.1416/4
totalVolume=2e-3
Cd=.75
m0 = 50e-3 
S=.008 #Surface area exposed to the flow
g=9.81
########################

#### INITIAL VALUES ####
theta = 90 * pi/180 
t = [ 0 ]
P = [ 500e3 ]
V = [ 1e-3 ] #initial air volume
vx = [1e-6* cos(theta)]
vz = [1e-6* sin(theta)]
x = [0]
z = [0]
#########################

def signal(x): return 1 if x>=0 else -1

def f(t, y, with_water):
    """System of differential equations, returns dy/dt
    with_water --- True or false, determine if there is water in the bottle
    """

    P = y[0] #manometric pressure inside the rocket
    V = y[1] #air volume
    vx = y[2] #velocity
    vz = y[3]
    x = y[4] 
    z  = y[5]#height
    v=sqrt(vx**2+vz**2)

    dydt=np.zeros(6)

 

    if with_water: 
        Q = signal(P)*sqrt(2*abs(P)/rhoWater)*nozzleArea #Q: volumetric flow rate
        waterVolume=totalVolume-V
        F = 2*P*nozzleArea #thrust
    else:
        #neglect thrust
        F = 0#thrust
        Q = 0
        P = 0
        waterVolume = 0

    dydt[0] = -(P * gamma * Q) / V
    dydt[1] = Q

    D = Cd*S*v**2 #drag
    m = m0+waterVolume*rhoWater
    ax = (F-D)*vx/(m*v)
    az = (F-D)*vz/(m*v) - g #acceleration
    
    dydt[2] = ax
    dydt[3] = az
    dydt[4] = vx #dx/dt = vx
    dydt[5] = vz

    return dydt


rocket = ode(f)
#rocket.set_integrator("dopri5")


solution = [P[0],V[0], vx[0], vz[0], x[0], z[0]]
rocket.set_initial_value(solution,0)

step = .01 #integration step

#The problem must be solved in two steps due to the discontinuity when it runs out of water
#Solve with water in the bottle
rocket.set_f_params(True)
while totalVolume-V[-1]>0 and z[-1]>=0:
    if(not rocket.successful()): 
        print ("ERROR: ode solver failed. Aborting.")
        break
    solution=rocket.integrate(rocket.t+step)
    t.append(rocket.t)
    P.append(solution[0])
    V.append(solution[1])
    vx.append(solution[2])
    vz.append(solution[3])
    x.append(solution[4])
    z.append(solution[5])

#remove last values, they are not valid
if(len(t)>1):
    t.pop()
    P.pop()
    V.pop()
    vx.pop()
    vz.pop()
    x.pop()
    z.pop()

solution = [P[-1],V[-1], vx[-1], vz[-1],  x[-1], z[-1]]

#solve with the bottle empty
rocket.set_f_params(False)
rocket.set_initial_value(solution,rocket.t-step)
while z[-1]>0:
    if(not rocket.successful()): 
        print ("ERROR: ode solver failed. Aborting.")
        break
    solution=rocket.integrate(rocket.t+step)
    t.append(rocket.t)
    P.append(solution[0])
    V.append(solution[1])
    vx.append(solution[2])
    vz.append(solution[3])
    x.append(solution[4])
    z.append(solution[5])

v=np.sqrt(np.power(vx,2)+np.power(vz,2))

print('Velocidade Máxima: ' )
print(np.max(v))
print('Altura Máxima: ')
print(np.max(z))
print('Alcance: ')
print(np.max(x))

plt.grid()
plt.plot(t, v)
plt.show()
