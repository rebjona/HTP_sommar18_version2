#--------------
#Pseudokod:

#Kör pennes.py för varje behandlingsplan och få tillhörande amplituder. (step 1)

#SAR(x) är nu känt då amplituder är bestämda. Stega med dt i dT=(T_0-T)*c*dt +SAR(x)*dt. Uppdatera T i varje iteration.

# Lägg till data för värmepermittivitet och/eller kapacitet för varje rumskoordinat x. Prio 2 efter att ha löst homogenfall?

#Bryt om T>5 i någon punkt (även in tumören?) Spara T-matris.

#Byt plan, utgå ifrån gammal T-matris. Stega med dt i dT=(T_0-T)*c*dt +SAR(x)*dt.

# Eventuellt pausa för att låta allting kylas ned? Om T inte får vara högre än 5 i tumören kan detta ge bättre T50 samtidigt som eventuella hot spots motverkas.

#----------------

import h5py

from dolfin import *
import numpy as np

# This function calculates temperature from Pld. A loop is used to scale amplitudes, in order to make sure that the goal temperature is reached inside the tumor.
#  Data needed: P-matrix, original amplitude settings in a text file, amplitude limit in a text file, bnd_heat_transfer.mat, bnd_temp_times_ht.mat, bnd_temp.mat, mesh.xml, perfusion_heatcapacity.mat, thermal_cond.mat and tumor_mesh.obj.
#  Output(saved): scaled amplitude settings and temperature.h5 file for further calculations and plots in MATLAB. Temperaturefiles as .pvd and .vtu can be generated as well (for ParaView).



# Load in .mat parameter
# The optional input argument 'degree' is FEniCS internal interpolation type.
# The loaded data will additionally be fit with trilinear interpolation.
def load_data(filename, degree=0):
    
    # Load the .mat file
    f = h5py.File(filename, "r")
    data = np.array(list(f.items())[0][1], dtype=float) # might need to remove list() when using Python 2?
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
# -------------------------
# T_b:      blood temperature [K relative body temp]
# P:        power loss density [W/m^3]
# k_tis:    thermal conductivity [W/(m K)]
# w_c_b:    volumetric perfusion times blood heat capacity [W/(m^3 K)]
# alpha:    boundary heat transfer constant [W/(m^2 K)]
# T_out_ht  alpha times ambient temperature [W/(m^2)]

print('Importing material properties...')
T_b = Constant(0.0) # Blood temperature relative body temp
P        = load_data("../Input_to_FEniCS/P.mat")
k_tis    = load_data("../Input_to_FEniCS/thermal_cond.mat")
w_c_b    = load_data("../Input_to_FEniCS/perfusion_heatcapacity.mat")
alpha    = load_data("../Input_to_FEniCS/bnd_heat_transfer.mat", 0)
T_out_ht = load_data("../Input_to_FEniCS/bnd_temp_times_ht.mat", 0)

f = h5py.File('../FEniCS_results/temperature.h5','w')

    f.create_dataset(name='Temp', data=T)
    f.create_dataset(name='P',    data=Coords)
    f.create_dataset(name='T',    data=Cells)
    # Need a dof(degree of freedom)-map to permutate Temp
    f.create_dataset(name='Map',  data=dof_to_vertex_map(V))

    f.close()

print('Importing temperature matrixes...')
with h5py.File("../FEniCS_results/temperature.h5",'r') as hdf:
	print('List of datasets; \n')
	data=hdf.get('Temp')
	Temp=np.array(data)
	
	
	
	
#-----------------------
Tmax= 5 # 0 = 37C, 8 if head and neck, 5 if brain
Tmin= 4.5 # 0 = 37C
#scaleLocal=1
Time=10 # just to try something, change this later
dt=1 # just to try something, change this later
numSteps=Time/dt
#-----------------------

# Read current amplitudes and the amplitude limit, generated in MATLAB
with open("../Input_to_FEniCS/amplitudes.txt") as file:
    amplitudes = []
    for line in file:
        amplitudes.append(line.rstrip().split(","))

with open("../Input_to_FEniCS/ampLimit.txt") as file:
    ampLimit = []
    for line in file:
        ampLimit.append(line.rstrip().split(","))

# Read old T, which is the initial condition used TODO


# Behövs ens skalningen? Om inte: ta bort den kodbiten i pennes.py som sparar skalningen i txt-fil
# load the scale of amplitudes and P found when running the original version of Pennes (step 1, i.e pennes.py)
with open("../FEniCS_results/scale_factor.txt") as file:
    scaleTotal = file.read()
    print(scaleTotal)


print("Done loading.")

#Change type of data
al=ampLimit[0][0]
ampLimit=float(al)
maxAmpl=max(amplitudes)
maxAmp=maxAmpl[0][0]
maxAmp=int(maxAmpl[0][0])
maxAmp=float(maxAmp)

# Get a first estimate of the temperature, using just the original P, i.e without any scaling at all.
# TODO: Find equation to solve, then uncomment it all
V = FunctionSpace(mesh, "CG", 1)
u = TrialFunction(V)
v = TestFunction(V)

# Variational formulation but using steps of time instead
#a= insert LHP here
#L= insert RHP here
u=Function(V)

#Initial condition TODO
#u_IC= Expression('insert eq here etc', t=0)
u_n= interpolate(u_IC,V)

#Define variational formulation to solve
#F= insert eq here
#a,L= lhs(F), rhs(F)

# Now take steps in time and estimate the temperature for each time, until the full scaling is made.
t=0
for n in range(numSteps)
    # Update time
    t +=dt
    u_IC.t=t

    # Solve the system
    #solve(a == L, u, solver_parameters={'linear_solver':'gmres'}) #might need to change from gmres to other solver?

    # Print the highest temperature obtained in this time step
    #T =u.vector().array()
    #print("Tmax for time step number " + t/dt)
    #print(np.max(T))

    # If okay temperature then save data for each time step in format readable by MATLAB
    #if (np.max(T)<Tmax and np.max(T)>Tmin):
    #Coords = mesh.coordinates()
    #Cells  = mesh.cells()
    # Index for this time step, should be included in the name of the file
    #index = t/dt
    #f = h5py.File('../FEniCS_results/temperature_' + index +'.h5','w')
    #f.create_dataset(name='Temp', data=T)
    #f.create_dataset(name='P',    data=Coords)
    #f.create_dataset(name='T',    data=Cells)
    # Need a dof(degree of freedom)-map to permutate Temp
    #f.create_dataset(name='Map',  data=dof_to_vertex_map(V))
    #f.close()
    #print("saved T for step: ")
    #print(index)








