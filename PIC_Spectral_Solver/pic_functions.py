#this script contains functions that are called in the particle in a cell (PIC) code
import numpy as np
from matplotlib import pyplot as plt
from pic_classes import Particle,Cell,Node,DSTpoisson #import classes from other script
#from celluloid import Camera #celluloid package for making an animation

#------------------------------------------------------------------
# FUNCTION THAT TAKES POSITION OF A PARTICLE AND MATCHES IT TO A CELL
def match_particle_to_cell(x,y,nx,ny,Lx,Ly):
    #this function takes in coordinates (x,y) of a particle, dimensions of the domain (nx,ny,Lx,Ly), 
    # and returns back the indices of a cell it is in (upper left cell has indices [0,0])
    # if the particle is on the edge between two cells, it falls into the lower-left cell

    h = Lx/(nx-1) #length of x side of a cell

    x_ind = -1 #set initial cell indices as negative (assume the particle is not found within the domain)
    y_ind = -1

    # x-coordinate:
    found_x = False
    for j in range(nx-1):
        if x >= j*h and x < (j+1)*h :
            found_x = True
            x_ind = j #x-index of the cell
    # y-coordinate:
    found_y = False
    for i in range(ny-1):
        if y >= (Ly - (i+1)*h) and y < (Ly - i*h) :
            found_y = True
            y_ind = i #y-index of the cell
    return(y_ind,x_ind)
#------------------------------------------------------------------


#------------------------------------------------------------------
# FUNCTION THAT CALCULATES THE WEIGHTS FOR SPLITTING NODE CONTRIBUTIONS (OF CHARGE DENSITY,ELECTRIC FIELD, ETC.)
# https://www.particleincell.com/2010/es-pic-method/
def assign_weights(x,y,ax,ay,H):
    hx = (x-ax)/H
    hy = (y-ay)/H

    wa = (1-hx)*(1-hy)
    wb = hx*(1-hy)
    wd = (1-hx)*hy
    wc = hx*hy

    return(wa,wb,wc,wd)
#------------------------------------------------------------------

#------------------------------------------------------------------
# FUNCTION THAT UPDATES THE ELECTRIC CHARGE DENSITY OVER THE WHOLE DOMAIN
def update_charge_density(particles,cells,nodes,Lx,Ly,nx,ny,h,cell_ind_matrix):
    for node in nodes: #first, charge density must be deleted at each node, before new values will be assigned
        node.zero_charge_dens()
    for particle in particles: #go over all particles and for each:
        [x,y] = particle.locate() #find its coordinates x,y
        if x>=0 and x<=Lx and y>=0 and y<=Ly: #if particle in the domain at all
            [i,j] = match_particle_to_cell(x,y,nx,ny,Lx,Ly) #find which cell it belongs to (find the cell's indices)

            if i >= 0 or j >= 0: #if particle found within some cell:

                cell_ID = cell_ind_matrix[i,j] #get the respective cell ID from these indices

                a_ID = cells[cell_ID].read_a_ID() #find the IDs of nodes this cell is touching
                b_ID = cells[cell_ID].read_b_ID()
                c_ID = cells[cell_ID].read_c_ID()
                d_ID = cells[cell_ID].read_d_ID()

                ax = nodes[a_ID].read_x() #find the coordinates of the cell's node "a" (lower left corner)
                ay = nodes[a_ID].read_y()
                [wa,wb,wc,wd] = assign_weights(x,y,ax,ay,h) #calculate charge-splitting weights
                if wa < 0:
                    print(wa)

                nodes[a_ID].add_charge_dens(wa*particle.charge/h**2) #add charge density to all 4 nodes with corresponding weights
                nodes[b_ID].add_charge_dens(wb*particle.charge/h**2)
                nodes[c_ID].add_charge_dens(wc*particle.charge/h**2)
                nodes[d_ID].add_charge_dens(wd*particle.charge/h**2)
    rho = []
    for node in nodes:
        rho.append(node.read_charge_dens())
    return np.array(rho)
#------------------------------------------------------------------


#------------------------------------------------------------------
# FUNCTION THAT SETS UP THE LINEAR PROBLEM AND SOLVES FOR THE ELECTRIC POTENTIAL (POISSON EQUATION SOLVER)
def solve_potential(nodes,h,epsilon_0,node_ind_matrix):
    N_inner_nodes = 0 #number of nodes whose potential needs to be solved for (non-boundary nodes), initially 0
    for node in nodes: #iterate over all nodes and count how many of them are truly inner (after setting boundary conditions before calling  this function)
        if node.is_inner() == True:
            node.set_inner_ID(N_inner_nodes)
            N_inner_nodes += 1

    C_vect = np.zeros((N_inner_nodes)) #initializing a constant vector that is going to reflect boundary conditions (initially just zeros)
    B_vect = np.full((N_inner_nodes),h**2/epsilon_0) #initializing a right-hand side vector (for h^2*charge_dens/epsilon_0), so far without charge dens.
    A = np.zeros((N_inner_nodes,N_inner_nodes)) #initializing a matrix to be filled with Poisson coefficients (so far empty)

    inner_node_IDs = [] # a list that will contain IDs of all INNER nodes in ordered fashion (for backsubstitution after solving for potential)
    for node in nodes:
        if node.is_inner() == True: #if the node is NOT a boundary condition node
            A_row = np.zeros(N_inner_nodes) #create an empty array that will become a row of a matrix eventually

            central_inner_ID = node.read_inner_ID() #this node's ID among all inner nodes
            central_ID = node.read_ID() #this node's ID among all nodes (including boundary)
            i = node.read_i() #read i and j indices of this node
            j = node.read_j()

            inner_node_IDs.append(central_ID) #store this node's main ID (for updating after solving the linear system)

            A_row[central_inner_ID] = 4 #fill the corresponding position of the row with central difference coefficient
            B_vect[central_inner_ID] *= node.read_charge_dens() #multiply the corresponding entry of the right-hand side B vector by the node's charge density

            left_ID = node_ind_matrix[i,j-1] #get the indices of the other stencil nodes
            right_ID = node_ind_matrix[i,j+1]
            upper_ID = node_ind_matrix[i-1,j]
            lower_ID = node_ind_matrix[i+1,j]

            #for each of the stencil nodes (except the middle one):
            if nodes[left_ID].is_inner():
                left_inner_ID = nodes[left_ID].read_inner_ID()
                A_row[left_inner_ID] = -1
            if nodes[right_ID].is_inner():
                right_inner_ID = nodes[right_ID].read_inner_ID()
                A_row[right_inner_ID] = -1
            if nodes[upper_ID].is_inner():
                upper_inner_ID = nodes[upper_ID].read_inner_ID()
                A_row[upper_inner_ID] = -1
            if nodes[lower_ID].is_inner():
                lower_inner_ID = nodes[lower_ID].read_inner_ID()
                A_row[lower_inner_ID] = -1

            A[central_inner_ID] = A_row #pass the created row of coefficients into the matrix

    # SOLVING THE LINEAR PROBLEM AV+C=B for V, UPDATING POTENTIAL AT NODES
    B_tilde = B_vect - C_vect #subtracting C from B
    V = np.linalg.solve(A,B_tilde) #using linalg to solve for V
    for i in range(len(V)): #update all relevant nodes with new values that we solved for
        V_val = V[i]
        nodes[inner_node_IDs[i]].set_potential(V_val)
    return V
#------------------------------------------------------------------
def solve_potentialfft(nodes,h,epsilon_0,nx, ny):
    elliptic = DSTpoisson(h, nx, ny)
    rho1 = []
    for node in nodes:
        rho1.append(node.read_charge_dens())
    rho = np.array(rho1).reshape(ny, -1)
    V = np.zeros_like(rho)
    V = elliptic.solver(V, -rho/epsilon_0)
    V = np.ravel(V)
    for i in range(len(V)): #update all relevant nodes with new values that we solved for
        nodes[i].set_potential(V[i])
    return V

##########################################################################

def Laplacian9(x, h):
    """Returns the isotropic 9-point centered differences Laplacian,
assuming that the mesh size h is the same in the x and y direction.
Note that the size is two units less than the shape of V in each
direction."""

    Lapl = (
        x[:-2, :-2] + x[:-2, 2: ] + x[2: , 2:] + x[2:, :-2] + 
        4*(x[:-2, 1:-1] + x[2:, 1:-1] + x[1:-1, :-2] + x[1:-1, 2:])
        - 20*x[1:-1, 1:-1]
    ) / (6*h*h)
    return Lapl

######################################################################

def Laplacian5(x, h):
    """Returns the 5-point centered differences Laplacian, assuming that
the mesh size h is the same in the x and y direction.
Note that the size is two units less than the shape of V in each direction."""
    Lapl = (
        x[:-2, 1:-1] + x[2:, 1:-1] + x[1:-1, :-2] + x[1:-1, 2:]
        - 4*x[1:-1, 1:-1]
    ) / (h*h)
    return Lapl


#------------------------------------------------------------------
# FUNCTION THAT UPDATES THE ELECTRIC FIELD AT NODES
def update_E_field(nodes,h,node_ind_matrix):
    for node in nodes:
        if node.is_inner() == True: #if the node is NOT a boundary condition node

            i = node.read_i() #read i and j indices of this node
            j = node.read_j()

            left_ID = node_ind_matrix[i,j-1] #get the indices of the other stencil nodes
            right_ID = node_ind_matrix[i,j+1]
            upper_ID = node_ind_matrix[i-1,j]
            lower_ID = node_ind_matrix[i+1,j]

            E_x = (nodes[left_ID].read_potential() - nodes[right_ID].read_potential())/(2*h) #Electric field in x-direction
            E_y = (nodes[lower_ID].read_potential() - nodes[upper_ID].read_potential())/(2*h) #Electric field in y-direction

            node.set_E_x(E_x) #update the E-field component attributes in the node instance
            node.set_E_y(E_y)
    return
#------------------------------------------------------------------
# FUNCTION THAT UPDATES THE ELECTRIC FIELD AT NODES #2
def update_E_field2(nodes,h,node_ind_matrix):
    V = []
    for node in nodes:
        V.append(node.read_potential())
    V = np.array(V).reshape(81, -1)
    E = np.gradient(V)
    Ex = -np.ravel(E[0])
    Ey = -np.ravel(E[1])
    for i in range(len(nodes)):
        nodes[i].set_E_x(Ex[i]) #update the E-field component attributes in the node instance
        nodes[i].set_E_y(Ey[i])
    return
#-----------------------------------------

def update_particle_motion(particles=None,cells=None,nodes=None,Lx=None,Ly=None,nx=None,ny=None,h=None,cell_ind_matrix=None,dt=None,update_pos=bool):
    energy = []
    if update_pos:
        particle_positions = []
    for particle in particles: #go over all particles and for each:
        [x,y] = particle.locate() #find its coordinates x,y
        if x>=0 and x<=Lx and y>=0 and y<=Ly: #if particle in the domain at all
            [i,j] = match_particle_to_cell(x,y,nx,ny,Lx,Ly) #find which cell it belongs to (find the cell's indices)

            if i >= 0 or j >= 0: #if particle found within some cell:

                cell_ID = cell_ind_matrix[i,j] #get the respective cell ID from these indices

                a_ID = cells[cell_ID].read_a_ID() #find the IDs of nodes this cell is made of (corners)
                b_ID = cells[cell_ID].read_b_ID()
                c_ID = cells[cell_ID].read_c_ID()
                d_ID = cells[cell_ID].read_d_ID()

                ax = nodes[a_ID].read_x() #find the coordinates of the cell's node "a" (lower left corner)
                ay = nodes[a_ID].read_y()
                [wa,wb,wc,wd] = assign_weights(x,y,ax,ay,h) #calculate splitting weights

                a_E_x = nodes[a_ID].read_E_x() #read the local electric field at each node of the cell, x and y separately
                a_E_y = nodes[a_ID].read_E_y()

                b_E_x = nodes[b_ID].read_E_x()
                b_E_y = nodes[b_ID].read_E_y()

                c_E_x = nodes[c_ID].read_E_x()
                c_E_y = nodes[c_ID].read_E_y()

                d_E_x = nodes[d_ID].read_E_x()
                d_E_y = nodes[d_ID].read_E_y()

                q = particle.read_charge() #read the particle's charge
                m = particle.read_mass() #read the particle's mass

                #calculate acceleration (at i):
                acc_x = q/m*(wa*a_E_x + wb*b_E_x + wc*c_E_x + wd*d_E_x) #combide E-fields at nodes into particle's acceleration, using weighted contributions
                acc_y = q/m*(wa*a_E_y + wb*b_E_y + wc*c_E_y + wd*d_E_y)

                

                v_x = particle.read_velocity_x() + acc_x*dt #calculate new velocity
                v_y = particle.read_velocity_y() + acc_y*dt 

                particle.set_velocity_x(v_x) #update particle velocity
                particle.set_velocity_y(v_y)

                a_V = nodes[a_ID].read_potential() #read the local potentials

                b_V = nodes[b_ID].read_potential()

                c_V = nodes[c_ID].read_potential()

                d_V = nodes[d_ID].read_potential()

                #calculate potential energy for each particle
                potential_energy = q*(wa*a_V + wb*b_V + wc*c_V + wd*d_V)

                kinetic_energy = m * (v_x**2 + v_y**2)/2
            
                energy.append(kinetic_energy + potential_energy)

                if update_pos:
                    if particle.read_x() + 2*v_x*dt < Lx and particle.read_x() + 2*v_x*dt > 0:
                        x_new = particle.read_x() + 2*v_x*dt #new position
                        v_xnew = particle.read_velocity_x()
                    if particle.read_x() + 2*v_x*dt > Lx or particle.read_x() + 2*v_x*dt < 0:
                        x_new = particle.read_x() - 2*v_x*dt #new position
                        v_xnew = - particle.read_velocity_x()
                    if particle.read_y() + 2*v_y*dt < Ly and particle.read_y() + 2*v_y*dt > 0:
                        y_new = particle.read_y() + 2*v_y*dt #new position
                        v_ynew = particle.read_velocity_y()
                    if particle.read_y() + 2*v_y*dt > Ly or particle.read_y() + 2*v_y*dt < 0:
                        y_new = particle.read_y() - 2*v_y*dt #new position
                        v_ynew = - particle.read_velocity_y()

                    particle.set_x(x_new) #update particle position
                    particle.set_y(y_new)
                    particle.set_velocity_x(v_xnew)
                    particle.set_velocity_y(v_ynew)

                    particle_positions.append([x_new,y_new,particle.read_name()])
    if update_pos:
        return particle_positions
    else:
        return energy

def generate_particles(particles=None,loc_x=None,loc_y=None,sigma_x=None,sigma_y=None,name=None,charge=None,mass=None,v_x=None,sigma_v_x=None,v_y=None,sigma_v_y=None,N=None):
    #generate new particle:
    last_ID = particles[-1].read_ID()
    parts_x = np.random.normal(loc_x, sigma_x, N)
    parts_y = np.random.normal(loc_y, sigma_y, N)
    vel_x = np.random.normal(v_x, sigma_v_x, N)
    vel_y = np.random.normal(v_y, sigma_v_y, N)
    for i in range(N):
        particles.append(Particle(ID=int(last_ID+1), x=parts_x[i], y=parts_y[i], name=name, charge=charge, mass=mass,velocity_x=vel_x[i],velocity_y=vel_y[i]))
        last_ID += 1
    return