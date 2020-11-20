import numpy as np

from fenics import *
from mshr import *

from openmdao.api import Group, IndepVarComp

class PMSMMeshGroup(Group):

    def initialize(self):
        self.options.declare('shape', types = tuple)
        self.options.declare('a', values = 10.) # inner radius ex. declaration from magnetostatics tutorial

    def setup(self):
        shape = self.options['shape']
        a = self.options['a']
        
        # just to test connections
        comp = IndepVarComp()
        comp.add_output('asdf')
        self.add_subsystem('test_comp', comp, promotes = ['*'])

        # import geometry and mesh from the magnetostatisc FEniCS tutorial and adjust for this motor
        # the motor inputs needs to be adjusted to fit variables as objects