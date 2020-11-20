import numpy as np

from openmdao.api import Group, IndepVarComp

from fenics_motor_model.motor_mesh.PMSM_mesh import PMSMMeshGroup

class MotorMeshGroup(Group):

    def initialize(self):
        self.options.declare('shape', types = tuple)
        self.options.declare('mode', types = str)
    
    def setup(self):
        shape = self.options['shape']
        mode = self.options['mode']

        if mode == 'PMSM':
            comp = PMSMMeshGroup(shape = shape)
            self.add_subsystem('pmsm_comp', comp, promotes = ['*'])