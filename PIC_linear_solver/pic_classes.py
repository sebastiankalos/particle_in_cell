class Particle(object):
    
    def __init__(self,  **kwargs):

        self._initialize( **kwargs)
        return

    def _initialize(self, **kwargs):
        self.ID = int(kwargs['ID'])
        self.x = float(kwargs['x'])
        self.y = float(kwargs['y']) 
        self.name = str(kwargs['name'])
        self.charge = float(kwargs['charge'])
        self.mass = float(kwargs['mass'])
        self.velocity_x = float(kwargs['velocity_x'])
        self.velocity_y = float(kwargs['velocity_y'])
        self.in_domain = bool(kwargs['in_domain'])
        return

    def locate(self):
        return(self.x,self.y)
    def set_x(self, x=None):
        self.x = float(x)
        return
    def set_y(self, y=None):   
        self.y = float(y) 
        return
    def set_name(self, name=None):
        self.name = str(name)
        return
    def set_charge(self, charge=None):
        self.charge = float(charge)
        return
    def set_ID(self, ID=None):
        self.ID = int(ID)
        return
    def set_mass(self, mass=None):
        self.mass = float(mass)
        return
    def set_velocity_x(self, velocity_x=None):
        self.velocity_x = float(velocity_x)
        return
    def set_velocity_y(self, velocity_y=None):
        self.velocity_y = float(velocity_y)
        return

    def read_name(self):
        return self.name
    def read_ID(self):
        return self.ID
    def read_x(self):
        return self.x
    def read_y(self):
        return self.y
    def read_charge(self):
        return self.charge
    def read_mass(self):
        return self.mass
    def read_velocity_x(self):
        return self.velocity_x
    def read_velocity_y(self):
        return self.velocity_y
    def is_in_domain(self):
        return self.in_domain

class Cell(object):
    
    def __init__(self, ID=None, N=None, x_min=None, x_max=None, y_min=None, y_max=None, node_IDs=None):

        self._initialize(ID=ID, N=N, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max, node_IDs=node_IDs)
        return

    def _initialize(self, ID=None, N=None, x_min=None, x_max=None, y_min=None, y_max=None, node_IDs=None):
        self.ID = int(ID)
        self.N = int(N)
        self.x_min = float(x_min)
        self.x_max = float(x_max) 
        self.y_min = float(y_min)
        self.y_max = float(y_max)
        self.a_ID = int(node_IDs[0])
        self.b_ID = int(node_IDs[1])
        self.c_ID = int(node_IDs[2])
        self.d_ID = int(node_IDs[3])
        return
    
    def info(self):
        print('Cell #',self.ID,', x_min=',self.x_min,', x_max=',self.x_max,', y_min=',self.y_min,', y_max=',self.y_max)
        return
    
    def set_x_min(self, x_min=None):
        self.x_min = float(x_min)
        return
    def set_x_max(self, x_max=None):   
        self.x_max = float(x_max) 
        return
    def set_y_min(self, y_min=None):
        self.y_min = float(y_min)
        return
    def set_y_max(self, y_max=None):
        self.y_max = float(y_max)
        return

    def set_ID(self, ID=None):
        self.ID = int(ID)
        return

    def set_N(self, N=None):
        self.N = int(N)
        return
    
    def add_particle(self):
        self.N = int(self.N+1)
        return

    def remove_particle(self):
        self.N = int(self.N-1)
        return

    def read_a_ID(self):
        return self.a_ID
    def read_b_ID(self):
        return self.b_ID
    def read_c_ID(self):
        return self.c_ID
    def read_d_ID(self):
        return self.d_ID

    def read_ID(self):
        return self.ID

class Node(object):
    
    def __init__(self, **kwargs):

        self._initialize(**kwargs)
        return

    def _initialize(self, **kwargs):
        self.ID = int(kwargs['ID'])
        self.x = float(kwargs['x'])
        self.y = float(kwargs['y']) 
        self.i = int(kwargs['i'])
        self.j = int(kwargs['j'])
        self.charge_dens = float(kwargs['charge_dens'])
        self.potential=float(kwargs['potential'])
        self.inner=bool(kwargs['inner'])
        self.inner_ID=int(kwargs['inner_ID'])
        self.E_x=float(kwargs['E_x'])
        self.E_y=float(kwargs['E_y'])
        return
    
    def set_x(self, x=None):
        self.x = float(x)
        return
    def set_y(self, y=None):   
        self.y = float(y) 
        return
    def set_i(self, i=None):
        self.i = int(i)
        return
    def set_j(self, j=None):
        self.j = int(j)
        return
    def set_ID(self, ID=None):
        self.ID = int(ID)
        return
    def set_inner_ID(self, inner_ID=None):
        self.inner_ID = int(inner_ID)
        return
    def set_charge_dens(self, charge_dens=None):
        self.charge_dens = float(charge_dens)
        return
    def set_potential(self, potential=None):
        self.potential = float(potential)
        return
    def set_inner(self, inner=None):
        self.inner = bool(inner)
        return
    def set_E_x(self, E_x=None):
        self.E_x = float(E_x)
        return
    def set_E_y(self, E_y=None):
        self.E_y = float(E_y)
        return
    def read_x(self):
        return self.x
    def read_y(self):
        return self.y
    def read_i(self):
        return self.i
    def read_j(self):
        return self.j
    def read_ID(self):
        return self.ID
    def read_inner_ID(self):
        return(self.inner_ID)
    def add_charge_dens(self,addition):
        self.charge_dens += addition
        return
    def zero_charge_dens(self):
        self.charge_dens = 0
        return
    def read_charge_dens(self):
        return(self.charge_dens)
    def read_potential(self):
        return(self.potential)
    def is_inner(self):
        return(self.inner)
    def read_E_x(self):
        return self.E_x
    def read_E_y(self):
        return self.E_y



