
import numpy as np
from matplotlib import pyplot as plt
from celluloid import Camera #celluloid package for making an animation
from mpl_toolkits import mplot3d
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import ticker

from pic_functions import update_charge_density,solve_potential,update_E_field,update_particle_motion,generate_particles
from pic_classes import Particle,Cell,Node #import classes from other script

#------------------------------------------------------------------
# FLAGS
animate_particle_motion = True
animate_potential = True
show_gridlines = False #show node gridlines in animation..?
#------------------------------------------------------------------

#------------------------------------------------------------------
# GLOBAL CONSTANTS
q_e = -1.602e-19 #charge of an electron
m_e = 9.10938356e-31 #kg, mass of an electron
m_p = 1.6726219e-27 #kg, mass of a proton
epsilon_0 = 8.85418782e-12 # permitivity of free space [m^-3 kg^-1 s^4 A^2]
dt_half = 0.000005 #leapfrog integration time interval (half of it) [s]
n_of_timesteps = 75 #number of integration steps to make (full dt's)
#------------------------------------------------------------------

#------------------------------------------------------------------
# DOMAIN SETUP:
nx = 81 #number of nodes in x-direction
ny = 81 #number of nodes in y-direction
n_nodes = nx*ny #number of nodes total
n_cells = (nx-1)*(ny-1) #total number of cells

Lx = 1 #length of the x-side of the box [m]
h = Lx/(nx-1) #side length of a square cell
Ly = h*(ny-1) #length of the y-side of the box (domain side ratio given by nx, ny!)
#------------------------------------------------------------------

#------------------------------------------------------------------
# GENERATING INITIAL PARTICLE DISTRIBUTION
N0 = 1 #initial number of 

mu_x, sigma_x = Lx/2, Lx/20 # mean and standard deviation of initial particle distribution (x-axis)
mu_y, sigma_y = Ly/2, Ly/20 # mean and standard deviation of initial particle distribution (y-axis)

parts_x_protons = np.random.uniform(0,Lx,N0)
parts_y_protons = np.random.uniform(0,Ly,N0)

particles = [] #a 1D list that will contain all particles (each being an object)

for i in range(N0): #generate protons
    particles.append(Particle(ID=i, x=parts_x_protons[i], y=parts_y_protons[i], name="electron", charge=q_e, mass=m_p,velocity_x=0,velocity_y=0,in_domain=True))
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
#------------------------------------------------------------------

#------------------------------------------------------------------
# SETTING BOUNDARY CONDITIONS
V_boundary = -25e-6 #placeholder boundary potential
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
for i in range (25,55):
    for j in range(1,30):
        nodes[node_ind_matrix[i,j]].set_potential(V_boundary) #left side
        nodes[node_ind_matrix[i,j]].set_inner(False) #left side nodes are not considered as inner nodes anymore
#------------------------------------------------------------------

#------------------------------------------------------------------
# MAIN LOOP HERE (leapfrog), following the "kick-drift-kick" form (https://en.wikipedia.org/wiki/Leapfrog_integration)
#------------------------------------------------------------------
if animate_potential:
    potential_data = []
if animate_particle_motion:
    particle_motion_data = []

step_counter = 0 #counter used for displaying current time in an animation
time = np.arange(0,n_of_timesteps*dt_half*2,dt_half) #create an array of times
for t in time:
    #generate new particles:
    if step_counter == 0:
        generate_particles(particles=particles,loc_x=3/4*0.98*Lx,loc_y=Ly/2,sigma_x=Lx/10,sigma_y=Ly/10,name="electron",charge=q_e,mass=m_e,v_x=0, sigma_v_x=0,v_y=0,sigma_v_y=0,N=2000)
    
        
    #calculate charge density, potential and E-field at time i:
    update_charge_density(particles,cells,nodes,Lx,Ly,nx,ny,h,cell_ind_matrix) #calculate charge density at nodes (at i)
    solve_potential(nodes,h,epsilon_0,node_ind_matrix) #solve for the potential (at time i)

    if animate_potential: #if animating the potential is turned on:
        X,Y,V = [],[],[]
        for node in nodes: #go over nodes and read potential values, store them in a list for later plotting
            V.append(node.read_potential())
            X.append(node.read_x())
            Y.append(node.read_y())
        potential_data.append([X,Y,V,t])

    update_E_field(nodes,h,node_ind_matrix) #updating electric field at nodes (at time i)

    #calculate and update particle positions (--> update_pos flag set as True) at i+1/2
    part_pos = update_particle_motion(particles=particles,cells=cells,nodes=nodes,Lx=Lx,Ly=Ly,nx=nx,ny=ny,h=h,cell_ind_matrix=cell_ind_matrix,dt=dt_half*2,update_pos=True)
    if animate_particle_motion: #if plotting turned on (flag among globals)
        particle_motion_data.append([part_pos,t+dt_half])


    #calculate charge density, potential and E-field at time i+1
    update_charge_density(particles,cells,nodes,Lx,Ly,nx,ny,h,cell_ind_matrix) #calculate charge density at nodes (at i+1)
    solve_potential(nodes,h,epsilon_0,node_ind_matrix) #solve for the potential (at i+1)
    update_E_field(nodes,h,node_ind_matrix) #updating electric field at nodes (at i+1)

    #UPDATING EACH PARTICLE'S VELOCITY AT i+1 (not position though --> update_pos boolean is False!)
    update_particle_motion(particles=particles,cells=cells,nodes=nodes,Lx=Lx,Ly=Ly,nx=nx,ny=ny,h=h,cell_ind_matrix=cell_ind_matrix,dt=dt_half*2,update_pos=False)
    step_counter += 1
    print("Solved step",step_counter,"/",n_of_timesteps*2)
#------------------------------------------------------------------

#------------------------------------------------------------------
# MAKING THE PARTICLE MOTION ANIMATION
#------------------------------------------------------------------
if animate_particle_motion: #flag (animation on/off)
    motion_fig = plt.figure() #initialize a plot
    motion_camera = Camera(motion_fig) # "point" a camera onto the plot

    frame_counter = 0
    for position_dataset in particle_motion_data: #plot particle positions at time i+1/2
        if show_gridlines: #if a flag among globals is raised, gridlines representing the mesh will be plotted in the animation
            for ver_ax in np.linspace(0,Lx,nx): #plot vertical axes representing grid of nodes
                plt.axvline(x=ver_ax,linestyle='--',linewidth=0.2,color='black')
            for hor_ax in np.linspace(0,Ly,ny): #plot horizontal axes representing grid of nodes
                plt.axhline(y=hor_ax,linestyle='--',linewidth=0.2,color='black')
        positions = position_dataset[0]
        t = position_dataset[1]
        for pos in positions:
            if "electron" in pos[2]:
                particle_color = "blue"
            elif "proton" in pos[2]:
                particle_color = "red"
            else:
                particle_color = "black"
            plt.plot(pos[0],pos[1],'.',markersize=1,color=particle_color)
        plt.text(0.05, 0.05, "t = {:10.5f} s".format(t))
        plt.xlabel("x [m]")
        plt.ylabel("y [m]")
        plt.xlim([0,Lx]) #set limits on the plot
        plt.ylim([0,Ly])
        motion_camera.snap() #take a shot of the particle positions

        frame_counter += 1
        print("animating particle motion: step",frame_counter,"/",round(n_of_timesteps*2))

    motion_animation = motion_camera.animate(blit=False, interval=50) #collect the snaps and make an animation
    motion_animation.save('animations/TEST_particle_motion.gif', writer = 'imagemagick') #save the animation

#------------------------------------------------------------------
# MAKING THE POTENTIAL ANIMATION
if animate_potential:
    potential_fig = plt.figure()
    potential_camera = Camera(potential_fig)
    ax = plt.axes(projection='3d')
    my_cmap = plt.get_cmap('hsv')

    V_max = max(np.array(potential_data)[:,2])

    frame_counter = 0 
    for dataset_at_t in potential_data:
        X = np.array(dataset_at_t[0])
        Y = np.array(dataset_at_t[1])
        V = np.array(dataset_at_t[2])
        t = np.array(dataset_at_t[3])

        sctt = ax.scatter3D(X, Y, V,c = V,cmap = my_cmap)

        ax.set_xlabel('x [m]')
        ax.set_ylabel('y [m]')
        ax.set_zlabel('V')
        ax.text2D(0.05, 0.95, "t = {:10.5f} s".format(t), transform=ax.transAxes)
        ax.set_title('potential')
        niceMathTextForm = ticker.ScalarFormatter(useMathText=True)
        ax.w_zaxis.set_major_formatter(niceMathTextForm)

        potential_camera.snap() #take a snapshot into the animation

        frame_counter += 1
        print("animating potential: step",frame_counter,"/",round(n_of_timesteps*2))

    potential_animation = potential_camera.animate(blit=False, interval=50) #produce the animation with all snaps made
    potential_animation.save('animations/TEST_potential_animation.gif', writer = 'imagemagick') #save the animation
#------------------------------------------------------------------

#------------------------------------------------------------------
# MAKING THE PHASE SPACE ANIMATION
if animate_particle_motion:
    phase_space_fig = plt.figure() #initialize a plot
    phase_space_camera = Camera(phase_space_fig) # "point" a camera onto the plot

    frame_counter = 0
    for position_dataset in particle_motion_data: #plot particle positions at time i+1/2
        positions = position_dataset[0]
        t = position_dataset[1]
        for pos in positions:
            if "electron" in pos[2]:
                particle_color = "blue"
                plt.plot(pos[0],pos[3],'.',markersize=1,color=particle_color)
        plt.text(0.05, 0.05, "t = {:10.5f} s".format(t))
        plt.xlabel("position x [m]")
        plt.ylabel("velocity along x [m/s]")
        plt.xlim([0,Lx]) #set limits on the plot
        plt.ylim([-800,800])

        phase_space_camera.snap() #take a shot of the particle positions

        frame_counter += 1
        print("animating phase space: step",frame_counter,"/",round(n_of_timesteps*2))

    phase_space_animation = phase_space_camera.animate(blit=False, interval=50) #collect the snaps and make an animation
    phase_space_animation.save('animations/TEST_phase_space_animation.gif', writer = 'imagemagick') #save the animation

    #
