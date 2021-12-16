# 
import numpy as np
from matplotlib import pyplot as plt

from pic_functions import update_charge_density,solve_potential,solve_potentialfft,update_E_field,update_E_field2,update_E_field,update_particle_motion,generate_particles
from pic_classes import Particle,Cell,Node,DSTpoisson #import classes from other script

#------------------------------------------------------------------
# FLAGS
animate_particle_motion = False
animate_potential = False

show_gridlines = True
#------------------------------------------------------------------

#------------------------------------------------------------------
# GLOBAL CONSTANTS
q_e = -1.602e-19 #charge of an electron
m_e = 9.10938356e-31 #kg, mass of an electron
m_p = 1.6726219e-27 #kg, mass of a proton
epsilon_0 = 8.85418782e-12 # permitivity of free space [m^-3 kg^-1 s^4 A^2]
dt_half = 0.01 #leapfrog integration time interval (half of it) [s]
n_of_timesteps = 10000 #number of integration steps to make (full dt's)
#------------------------------------------------------------------

#------------------------------------------------------------------
# DOMAIN SETUP:
a = 1 #aspect ratio
nx = 81 #number of nodes in x-direction
ny = a*nx #number of nodes in y-direction
n_nodes = nx*ny #number of nodes total
n_cells = (nx-1)*(ny-1) #total number of cells


Lx = 1. #length of the x-side of the box [m]
Ly = a*Lx
x  = np.linspace(0, Lx, nx)
y  = np.linspace(0, Ly, ny)
h  = x[1] - x[0]
X, Y = np.meshgrid(x, y)
#------------------------------------------------------------------

#------------------------------------------------------------------
# INITIAL PARTICLE DISTRIBUTION
Nprotons = 1000000
Nelectrons = 0

mu_x, sigma_x = Lx/2, Lx/20 # mean and standard deviation of initial particle distribution (x-axis)
mu_y, sigma_y = Ly/2, Ly/20 # mean and standard deviation of initial particle distribution (y-axis)
#parts_x_batch1 = np.random.normal(Lx/4, sigma_x, N0) # generate initial particle d istribution of particle positions (x,y)
#parts_y_batch1 = np.random.normal(mu_y, sigma_y, N0)
#parts_x_batch2 = np.random.normal(Lx*3/4, sigma_x, N0) # generate initial particle distribution of particle positions (x,y)
#parts_y_batch2 = np.random.normal(mu_y, sigma_y, N0)

parts_x_protons = np.random.uniform(Lx/4,3*Lx/4,Nprotons)
parts_y_protons = np.random.uniform(Ly/4,3*Ly/4,Nprotons)
parts_x_electrons = np.random.uniform(Lx/4,3*Lx/4,Nelectrons)
parts_y_electrons = np.random.uniform(Ly/4,3*Ly/4,Nelectrons)
#------------------------------------------------------------------

#------------------------------------------------------------------
# GENERATING PARTICLES
particles = [] #a 1D list that will contain all particles (each being an object)
#for i in range(int(N0/2)): #generate particles, assign an initial position, name, charge, etc. to each, and append them all in a list
#    particles.append(Particle(ID=i, x=parts_x_batch1[i], y=parts_y_batch1[i], name="electron", charge=q_e, mass=m_e,velocity_x=400,velocity_y=0))
#for i in range(int(N0/2),N0):
#    particles.append(Particle(ID=i, x=parts_x_batch2[i], y=parts_y_batch2[i], name="electron", charge=q_e, mass=m_e,velocity_x=-400,velocity_y=0))
# particles.append(Particle(ID=1 + N0, x= 0.5, y= 0.5, name="electron", charge=q_e, mass=m_e,velocity_x=0,velocity_y=0))
for i in range(Nprotons):
    particles.append(Particle(ID=i, x=parts_x_protons[i], y=parts_y_protons[i], name="proton", charge=-q_e, mass=m_p,velocity_x=0,velocity_y=0))
for i in range(Nelectrons):
    particles.append(Particle(ID=i + Nprotons, x=parts_x_electrons[i], y=parts_y_electrons[i], name="electron", charge=q_e, mass=m_e,velocity_x=0,velocity_y=0))

#------------------------------------------------------------------

#------------------------------------------------------------------
# GENERATING NODES
node_ind_matrix = np.linspace(0,n_nodes-1,n_nodes,dtype=int).reshape(ny, -1) #a matrix of node IDs with m=ny, n=nx, starting from 0 in upper left corner
nodes = [] # a 1D list that will contain all nodes (each being an object)
nodes_x = np.linspace(0,Lx,nx) #find x and y coordinates of all nodes in the domain
nodes_y = np.linspace(0,Ly,ny)
node_count = 0
for i in range(ny): # go over the domain, from left upper corner, row by row, and generate nodes (objects), and store them in the list
    for j in range(nx):
        nodes.append(Node(ID=node_count, x=j*h, y=Ly-i*h, i=i, j=j, charge_dens=0,potential=0,inner=True,inner_ID=-1,E_x=0,E_y=0))
        node_count += 1
#------------------------------------------------------------------

#------------------------------------------------------------------
# GENERATING CELLS
cell_ind_matrix = np.linspace(0,n_cells-1,n_cells,dtype=int).reshape(ny-1, -1) #a matrix of cell ID's with m=ny-1, n=nx-1, starting from 0 in upper left corner
cells = [] # a 1D list that will contain all cells (each one being an object)
cell_count = 0
for i in range(ny-1): # go over the domain, from left upper row by row, and generate cells, and store them in the list
    for j in range(nx-1):
        d_ID = nx*i + j #compute what ID the cell's node "d" (upper left) has (among all other nodes)
        a_ID = nx*(i+1) + j #compute what ID the cell's node "a" (lower left) has 
        c_ID = d_ID + 1 #compute what ID the cell's node "c" (upper right) has
        b_ID = a_ID + 1 #compute what ID the cell's node "b" (lower right) has
        new_cell = Cell(ID=cell_count, N=0, x_min=j*h, x_max=(j+1)*h, y_min=Ly-(i+1)*h, y_max=Ly-i*h, node_IDs=[a_ID,b_ID,c_ID,d_ID]) #generate the cell
        cells.append(new_cell) #store the new cells
        cell_count += 1
#-----------------------------------------------------------------

#------------------------------------------------------------------
# SETTING BOUNDARY CONDITIONS
V_boundary = 0 #placeholder boundary potential
for j in range(0,nx):
    nodes[node_ind_matrix[0,j]].set_potential(V_boundary) #upper row
    nodes[node_ind_matrix[0,j]].set_inner(False) #upper row nodes are not considered as inner nodes anymore
    nodes[node_ind_matrix[ny-1,j]].set_potential(V_boundary) #lower row
    nodes[node_ind_matrix[ny-1,j]].set_inner(False) #lower row nodes are not considered as inner nodes anymore
for i in range(1,ny-1):
    nodes[node_ind_matrix[i,0]].set_potential(V_boundary) #left side
    nodes[node_ind_matrix[i,0]].set_inner(False) #left side nodes are not considered as inner nodes anymore
    nodes[node_ind_matrix[i,ny-1]].set_potential(V_boundary) #right side
    nodes[node_ind_matrix[i,ny-1]].set_inner(False) #right side nodes are not considered as inner nodes anymore
# ------------------------------------------------------------------
#------------------------------------------------------------------
# MAIN LOOP HERE (leapfrog), following the "kick-drift-kick" form (https://en.wikipedia.org/wiki/Leapfrog_integration)
step_counter = 0 #counter used for displaying current time in an animation

energy = []
for i in range(n_of_timesteps):
    # generate_particles(particles=particles,loc_x=Lx/4,  loc_y=Ly/2,sigma_x=Lx/20,sigma_y=Ly/20,name="electron",charge=q_e,mass=m_e,v_x=40, sigma_v_x=2,v_y=0,sigma_v_y=5,N=1)
    # generate_particles(particles=particles,loc_x=Lx/2,loc_y=Ly/2,sigma_x=0,sigma_y=0,name="proton",charge= - q_e,mass=m_p,v_x= 1,sigma_v_x=0,v_y=1,sigma_v_y=0,N=10)
    #calculate and update particle positions (--> update_pos flag set as True) at i+1/2
    part_pos = update_particle_motion(particles=particles,cells=cells,nodes=nodes,Lx=Lx,Ly=Ly,nx=nx,ny=ny,h=h,cell_ind_matrix=cell_ind_matrix,dt=dt_half,update_pos=True)
    if animate_particle_motion: #if plotting turned on (flag among globals)
        plt.clf()
        plt.figure(figsize=(6, 6), dpi=80)
        plt.suptitle('Particles')
        plt.xlim(0, Lx)
        plt.ylim(0, Lx)
        for pos in part_pos:
            if "electron" in pos[2]:
                particle_color = "red"
            elif "proton" in pos[2]:
                particle_color = "blue"
            else:
                particle_color = "black"
            plt.plot(pos[0],pos[1],'.',markersize = 3,color=particle_color)
        if show_gridlines: #if a flag among globals is raised, gridlines representing the mesh will be plotted in the animation
            plt.minorticks_on()
            # Customize the major grid
            plt.grid(which='major', linestyle='-', linewidth='0.5', color='grey')
            # Customize the minor grid
            plt.grid(which='minor', linestyle=':', linewidth='0.5', color='grey')
            legend_elements = [
                   plt.Line2D([0], [0], marker='o', color= 'w', label='Proton',
                          markerfacecolor='b', markersize=5),
                   plt.Line2D([0], [0], marker='o', color= 'w', label='Electron',
                          markerfacecolor='r', markersize=5)]
        plt.legend(handles = legend_elements, loc = 'upper right')
        plt.savefig("./Motion_plots/Fig-{:05d}.png".format(i))
    
    #calculate charge density, potential and E-field at time i:
    rho = update_charge_density(particles,cells,nodes,Lx,Ly,nx,ny,h,cell_ind_matrix) #calculate charge density at nodes (at i)
    V = solve_potentialfft(nodes, h, epsilon_0, nx, ny) #solve for the potential (at time i)
    update_E_field2(nodes,h,node_ind_matrix) #updating electric field at nodes (at time i)
    if animate_potential: #if animating the potential is turned on:
        plt.clf #clears the previously drawn figure
        potential = V.reshape(ny, -1)
        density = rho.reshape(ny, -1)
        fig, axes = plt.subplots(nrows = 1, ncols = 2, figsize=(14, 5))
        plt.suptitle("Time = {:10.5f} s".format(2*i*dt_half))
        c1 = axes[1].contourf(X, Y, potential, cmap = 'RdBu_r')
        axes[1].set_title('Potential')
        axes[1].set_xlabel('x')
        axes[1].set_ylabel('y')
        plt.colorbar(c1, ax=axes[1]);
        c2 = axes[0].contourf(X, Y, density, cmap = 'Spectral')
        axes[0].set_title('Charge density')
        axes[0].set_xlabel('x')
        axes[0].set_ylabel('y')
        plt.colorbar(c2, ax=axes[0]);
        plt.savefig("./Potential_Density_plots/Fig-{:05d}.png".format(i))

    #UPDATING EACH PARTICLE'S VELOCITY AT i+1 (not position though --> update_pos boolean is False!)
    energies = update_particle_motion(particles=particles,cells=cells,nodes=nodes,Lx=Lx,Ly=Ly,nx=nx,ny=ny,h=h,cell_ind_matrix=cell_ind_matrix,dt=dt_half,update_pos=False)
    energy.append(np.mean(energies)/(-q_e))
    step_counter += 1
    print("Solved step",step_counter,"/",n_of_timesteps)
    
#------------------------------------------------------------------
# 
plt.clf()
plt.figure(figsize=(8, 6), dpi=80)
plt.title('Plotting the mean of the overall energy over time')
plt.minorticks_on()
# Customize the major grid
plt.grid(which='major', linestyle='-', linewidth='0.5', color='black')
# Customize the minor grid
plt.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
plt.plot(np.linspace(0, n_of_timesteps*dt_half*2, n_of_timesteps), energy, color = 'red')
plt.xlabel('$t$')
plt.ylabel('$\overline{E}$')
plt.savefig("./Motion_plots/Energy.png")