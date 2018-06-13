"""

 This is a code that combines Pennes.py with CoolPennes.py, i.e both scaling and time computations are performed.
 TODO: implement the code for combining plans and implement non-linear parameters for the perfusion
 (Also: Improve time stepping, not working correctly at the moment)
    
"""

import h5py

from dolfin import *
import numpy as np

# Load .mat parameter
# The optional input argument 'degree' is FEniCS internal interpolation type.
# The loaded data will additionally be fit with trilinear interpolation.
def load_data(filename, degree=0):
    
    # Load the .mat file
    f = h5py.File(filename, "r")
    data = np.array(list(f.items())[0][1], dtype=float)
    f.close()
    
    # Load the intepolation c++ code
    f = open('TheGreatInterpolator.cpp', "r")
    code = f.read()
    f.close()
    
    # Amend axis ordering to match layout in memory
    size = tuple(reversed(np.shape(data)))
    
    # Add c++ code to FEniCS
    P = Expression(code, degree=degree)
    
    # Add parameters about the data
    P.stridex  = size[0]
    P.stridexy = size[0]*size[1]
    P.sizex = size[0]
    P.sizey = size[1]
    P.sizez = size[2]
    P.sidelen = 1.0/1000
    
    # As the last step, add the data
    P.set_data(data)
    return P

# Load mesh
print("Reading and unpacking mesh...")
mesh = Mesh('../Input_to_FEniCS/mesh.xml')

# Define material properties
# ----------------------------------------------------------
# T_b:      blood temperature [K relative body temp]
# P:        power loss density [W/m^3]
# k_tis:    thermal conductivity [W/(m K)]
# w_c_b:    volumetric perfusion times blood heat capacity [W/(m^3 K)]
# alpha:    boundary heat transfer constant [W/(m^2 K)]
# T_out_ht  alpha times ambient temperature [W/(m^2)]

print('Importing material properties...')
# Load P matrices, either just one or several depending on how many HT plans one wants to combine. TODO make it possible to combine plans
P        = load_data("../Input_to_FEniCS/P.mat")
#P2        = load_data("../Input_to_FEniCS/P2.mat")
#P3        = load_data("../Input_to_FEniCS/P3.mat")

T_b = Constant(0.0) # Blood temperature relative body temp
k_tis    = load_data("../Input_to_FEniCS/thermal_cond.mat")

# Load the w_c_b, depending on whether one wants to use linear perfusion data or non-linear perfusion data. TODO create the non-linear perfusion data and implement it in the Matlab code
w_c_b    = load_data("../Input_to_FEniCS/perfusion_heatcapacity.mat") # This is the "standard" perfusion matrix with linear values
#w_c_b   = load_data("../Input_to_FEniCS/perfusion_heatcapacity_nonlinear.mat") # TODO This should be chosen if a non-linear scaling of the perfusion is wanted, not created yet though
alpha    = load_data("../Input_to_FEniCS/bnd_heat_transfer.mat", 0)
T_out_ht = load_data("../Input_to_FEniCS/bnd_temp_times_ht.mat", 0)


# Read current amplitudes and the amplitude limit, generated in MATLAB
with open("../Input_to_FEniCS/amplitudes.txt") as file:
    amplitudes = []
    for line in file:
        amplitudes.append(line.rstrip().split(","))

with open("../Input_to_FEniCS/ampLimit.txt") as file:
    ampLimit = []
    for line in file:
        ampLimit.append(line.rstrip().split(","))

print("Done loading.")

# Set parameters
#-----------------------
Tmax= 5 # 0 = 37C, 8 if head and neck, 5 if brain
Tmin= 4.5 # 0 = 37C
Time=60*10
dt=120
numSteps=Time/dt

#Change type of data
al=ampLimit[0][0]
ampLimit=float(al)
maxAmpl=max(amplitudes)
maxAmp=maxAmpl[0][0]
maxAmp=int(maxAmpl[0][0])
maxAmp=float(maxAmp)

# Define function space and test/trial functions needed for the variational formulation
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

# Perform the scaling iteratively ------------------------------------------------------

scaleTot=1;
nbrIter=0;
T=0
done=False

while (((np.max(T)<Tmin or np.max(T)>Tmax) and nbrIter<=maxIter) or maxAmp>ampLimit):
    
    #If amplitude is too high, maxAmp is set to amlitude limit
    if (maxAmp>ampLimit):# and np.max(T)<Tmax):
        print(np.max(T))
        scaleAmp=(ampLimit/maxAmp)**2
        maxAmp=ampLimit
        scaleTot=scaleTot*(scaleAmp)
        P=P*scaleAmp
        
        V = FunctionSpace(mesh, "CG", 1)
        u = TrialFunction(V)
        v = TestFunction(V)
        # Variation formulation of Pennes heat equation
        a = v*u*alpha*ds + k_tis*inner(grad(u), grad(v))*dx + w_c_b*v*u*dx
        L = T_out_ht*v*ds + P*v*dx # + w_c_b*T_b*v*dx not needed due to T_b = 0
        
        u = Function(V)
        solve(a == L, u, solver_parameters={'linear_solver':'gmres'}) #gmres is fast
        T =u.vector().array()
        print("Tmax:")
        print(np.max(T))
        print("Scale:")
        print(scaleTot)
        if (np.max(T)<Tmax):
            done = True # exit loop

    elif (maxAmp<=ampLimit):
        V = FunctionSpace(mesh, "CG", 1)
        u = TrialFunction(V)
        v = TestFunction(V)
        # Variation formulation of Pennes heat equation
        a = v*u*alpha*ds + k_tis*inner(grad(u), grad(v))*dx + w_c_b*v*u*dx
        L = T_out_ht*v*ds + P*v*dx # + w_c_b*T_b*v*dx not needed due to T_b = 0
        u = Function(V)
        solve(a == L, u, solver_parameters={'linear_solver':'gmres'}) #gmres is fast
        
        #Use T to find the scale for P and maxAmp (increase T)
        if(np.max(T)<=Tmin):
            if(np.max(T)<0.4*Tmin):
                P=P*(1.6)
                scaleTot=scaleTot*1.6
                maxAmp=sqrt(1.6)*maxAmp
            elif (np.max(T)>=0.4*Tmin and np.max(T)<0.8*Tmin):
                P=P*(1.3)
                scaleTot=scaleTot*1.3
                maxAmp=sqrt(1.3)*maxAmp
            elif (np.max(T)>=0.8*Tmin):
                P=P*1.05
                scaleTot=1.05*scaleTot
                maxAmp=sqrt(1.05)*maxAmp
        #Use T to find the scale for P and maxAmp (decrease T)
        if(np.max(T)>=Tmax):
            if(np.max(T)>1.4*Tmax):
                P=P*0.5
                scaleTot=scaleTot*0.5
                maxAmp=sqrt(0.5)*maxAmp
            elif(np.max(T)>1.2*Tmax and np.max(T)<=1.4*Tmax):
                P=P*0.7
                scaleTot=scaleTot*0.7
                maxAmp=sqrt(0.7)*maxAmp
            elif (np.max(T)<=1.2*Tmax):
                P=P*0.85
                scaleTot=scaleTot*(0.85)
                maxAmp=sqrt(0.85)*maxAmp

T =u.vector().array()

nbrIter=nbrIter+1
    print("Tmax:")
    print(np.max(T))
    print("Scale:")
    print(scaleTot)
    print("MaxAmp:")
    print(maxAmp)
    
    if(done):
        break

# Plot and save
if ((np.max(T)>Tmin and np.max(T)<Tmax and maxAmp<=ampLimit) or maxAmp==ampLimit):
    
    # Plot solution and mesh
    plot(u)
    plot(mesh)
    
    # Save data in a format readable by matlab
    Coords = mesh.coordinates()
    Cells  = mesh.cells()
    
    f = h5py.File('../FEniCS_results/temperature.h5','w')
    
    f.create_dataset(name='Temp', data=T)
    f.create_dataset(name='P',    data=Coords)
    f.create_dataset(name='T',    data=Cells)
    # Need a dof(degree of freedom)-map to permutate Temp
    f.create_dataset(name='Map',  data=dof_to_vertex_map(V))
    
    f.close()
    
    #Scale amplitudes and save in a new file
    amplitudeVec=[]
    fileAmp=open('../FEniCS_results/scaledAmplitudes.txt','w')
    for x in amplitudes:
        a=float(x[0])
        a=(round(a*100))*sqrt(scaleTot)/100
        amplitudeVec.append(a)
        fileAmp.write(str(a) + " ")

    fileAmp.close()

# Save the scale factor in a file
fileScale=open('../FEniCS_results/scale_factor.txt','w')
fileScale.write(str(scaleTot))
fileScale.close()

#Print parameters
print("Tmax:")
    print(np.max(T))
    print("Scale:")
    print(scaleTot)
    print("Nbr of iterations:")
    print(nbrIter)
    print("MaxAmp:")
    print(maxAmp)
    
    if (np.max(T)>Tmax and ampLimit==maxAmp):
        print(" High temperature. Try to increase the interval [Tmin,Tmax] or try a higher maxIter.")

else:
    print("Not enough iterations for the scaling")

print("Scaling finished")
#-----------------------------------------------------------------------------------------

# Perform the time calculations as in CoolPennes------------------------------------------

scale=Constant(scaleTotal)

# Define function space and test/trial functions needed for the variational formulation
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

#Initial condition
u_IC= Expression("0", t=0, degree=0) # degree=1?
u_n=interpolate(u_IC,V)

P=P*scale # Scale P according to previous calculations
F=alpha*u*v*ds + v*u*dx + dt*k_tis*dot(grad(u), grad(v))*dx - (u_n + dt*(P-w_c_b*u))*v*dx + T_out_ht*v*ds
#alpha*u*v*ds + v*u*dx + dt*k_tis*dot(grad(u), grad(v))*dx - (u_n + dt*(P-w_c_b*u))*v*dx + T_out_ht*v*ds
#alpha*u*v*dx + dt*k_tis*dot(grad(u), grad(v))*dx - (u_n + dt*(P-w_c_b*u))*v*dx + T_out_ht*v*ds
a=lhs(F)
L=rhs(F)

u=Function(V)

# Now take steps in time and estimate the temperature for each time step, until the full scaling is made.
t=0
for n in range(int(numSteps)):
    # Update time
    t += dt
    u_IC.t=t
    
    # Solve the system
    solve(a == L, u, solver_parameters={'linear_solver':'gmres'})   #might need to change from gmres to other solver?
    T =u.vector().array()
    
    # Print the highest temperature
    print("Tmax for time step number " + str(t/dt) + ":")
    print(np.max(T))
    
    u_n.assign(u)
    
    # If okay temperature then save data for each time step in format readable by MATLAB
    """
        if (np.max(T)<Tmax and np.max(T)>Tmin):
        Coords = mesh.coordinates()
        Cells  = mesh.cells()
        
        # Index for this time step should be included in the name for the temperature file
        index=t/dt
        f = h5py.File('../FEniCS_results/temperature_'+ str(index) + '.h5','w')
        f.create_dataset(name='Temp', data=T)
        f.create_dataset(name='P',    data=Coords)
        f.create_dataset(name='T',    data=Cells)
        # Need a dof(degree of freedom)-map to permutate Temp
        f.create_dataset(name='Map',  data=dof_to_vertex_map(V))
        f.close()
        print("saved T for step: ")
        print(index)
        """

print("Time iteration finished")
#-------------------------------------------------------------------------------------------















