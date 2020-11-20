'''This file models a Permanent Magnet Synchronous Motor (PMSM) 
for EVTOL aircrafts'''

import numpy as np

from openmdao.api import Problem, Group, IndepVarComp

from fenics_motor_model.motor_mesh.motor_mesh_group import MotorMeshGroup

n = 1
shape = (n,n)
mode = 'PMSM' # will attempt to include other motor models later

prob = Problem()

global_variable = IndepVarComp()
global_variable.add_output('motor_mass', val = 0.)
global_variable.add_output('normalized_torque', val = 0.)
global_variable.add_output('angular_speed', val = 0.)
global_variable.add_output('number_of_poles_per_phase', val = 0.)
global_variable.add_output('hysteresis_coeff', val = 0.)
global_variable.add_output('copper_resistivity', val = 0.) # or other material used for motor
global_variable.add_output('eta_slot', val = 0.)
global_variable.add_output('eta_fill', val = 0.)
prob.model.add_subsystem('global_variables_comp', global_variable, promotes = ['*'])

group = MotorMeshGroup(
    shape = shape,
    mode = mode
)
prob.model.add_subsystem('motor_mesh_group',group, promotes = ['*'])

prob.setup()
prob.run_model()