import gmsh
from math import sin, cos, tan, pi, floor, ceil
import numpy as np
import os

'''
Mesh Generation of 2D Radial Flux PMSM; only half is modeled (180 degrees) due to periodicity
'''

p = 8 # pole pairs
s = 9 # number of stator slots for windings per 180 degrees
m = 3 # number of phases for stator winding current

gmsh.initialize() # gmsh must be initialized in Python before using functions
gmsh.option.setNumber("General.Terminal",1) # for printing messages to Terminal

gmsh.model.add("Motor") # Creating new model called motor
motor = gmsh.model.occ

lc = 1e-2 # target mesh size

# Key Geometrical Parameters of Motor ----------------------------------------------------

# Key Rotor Angles for Geometry
theta_p = 2*pi/2/p # Angular sweep of each Rotor Slice
theta_m = .78 * theta_p # Angular sweep of magnet
theta_b = .0595 * theta_p # Angular sweep of magnet side piece separation (gap between each slice)
theta_g = (theta_p - theta_m - theta_b)/2 # Angular sweep of each magnet side piece

# Key Stator Angles for Geometry
theta_t = pi/s # Angular sweep of each Stator Tooth
theta_sso = .5 * theta_t # Angular sweep of total area of windings
theta_ssi = .3 * theta_sso # Angular sweep of tooth tip separation

# Need to define each of these
Rr = .80
Rtm = .79
Rtb = .775
Rbb = .75
Rbm = .745
Rin = .60 # Inner Radius of Rotor

Rout = 1.15
Rsy = 1.03
Rssi = .83
Rs = .81

D = 2 # Domain Radius
Ras = (Rr+Rs)/2 # Radius splitting air-gap mesh of rotor and stator
RR = (Rin+Rr)/2 # Midpoint to cut air-gap mesh in Rotor
RS = (Rsy+Rout)/2 # # Midpoint to cut air-gap mesh in Stator

# Points of first slice ----------------------------------------------------------------------------------
'''Origin'''
motor.addPoint(0,0,0,lc,1)

'''Rotor Section'''
# Inner surface of Rotor Core
motor.addPoint(Rin,0,0,lc,2) # Inner Radius of Rotor along "+" x-axis
motor.addPoint(0,Rin,0,lc,3)  # Inner Radius of Rotor along "+" y-axis
motor.addPoint(-Rin,0,0,lc,4)  # Inner Radius of Rotor along "-" x-axis

# Rotor End Points (Need so we don't have duplicates)
motor.addPoint(-Rtm,0,0,lc,5)

# Outer surface of Rotor Core
motor.addPoint(Rtm,0,0,lc,6)
motor.addPoint(Rtm * cos(theta_b/2),Rtm * sin(theta_b/2),0,lc,7)
motor.addPoint(Rr * cos(theta_b/2 + theta_g), Rr * sin(theta_b/2 + theta_g),0,lc,8)
motor.addPoint(Rr * cos(theta_b/2 + theta_g + theta_m), Rr * sin(theta_b/2 + theta_g + theta_m),0,lc,9)
motor.addPoint(Rtm * cos(theta_p - theta_b/2),Rtm * sin(theta_p - theta_b/2),0,lc,10)
# motor.addPoint(Rtm * cos(theta_p),Rtm * sin(theta_p),0,lc,10) # REMOVE

# Magnets and air slots
motor.addPoint(Rbb * cos(theta_b/2),Rbb * sin(theta_b/2),0,lc,11)
motor.addPoint(Rtb * cos(theta_b/2),Rtb * sin(theta_b/2),0,lc,12)
motor.addPoint(Rbm * cos(theta_b/2 + theta_g),Rbm * sin(theta_b/2 + theta_g),0,lc,13)
motor.addPoint(Rbb * cos(theta_b/2 + theta_g),Rbb * sin(theta_b/2 + theta_g),0,lc,14)
motor.addPoint(Rtm * cos(theta_b/2 + theta_g),Rtm * sin(theta_b/2 + theta_g),0,lc,15)
motor.addPoint(Rbm * cos(theta_b/2 + theta_g + theta_m),Rbm * sin(theta_b/2 + theta_g + theta_m),0,lc,16)
motor.addPoint(Rbb * cos(theta_b/2 + theta_g + theta_m),Rbb * sin(theta_b/2 + theta_g + theta_m),0,lc,17)
motor.addPoint(Rtm * cos(theta_b/2 + theta_g + theta_m),Rtm * sin(theta_b/2 + theta_g + theta_m),0,lc,18)
motor.addPoint(Rbb * cos(theta_p - theta_b/2),Rbb * sin(theta_p - theta_b/2),0,lc,19)
motor.addPoint(Rtb * cos(theta_p - theta_b/2),Rtb * sin(theta_p - theta_b/2),0,lc,20)

'''Stator Section'''
# Outer surface of Stator Core
motor.addPoint(Rout,0,0,lc,21) # Inner Radius of Stator along "+" x-axis
motor.addPoint(0,Rout,0,lc,22)  # Inner Radius of Stator along "+" y-axis
motor.addPoint(-Rout,0,0,lc,23)  # Inner Radius of Stator along "-" x-axis

# Stator End Points (Need so we don't have duplicates)
motor.addPoint(-Rsy,0,0,lc,24)

# Inner surface of Stator Core & Winding Definition
motor.addPoint(Rsy,0,0,lc,25)
motor.addPoint(Rsy*cos(theta_sso/2),Rsy*sin(theta_sso/2),0,lc,26)
motor.addPoint(Rssi*cos(theta_sso/2),Rssi*sin(theta_sso/2),0,lc,27)
motor.addPoint(Rssi*cos(theta_ssi/2),Rssi*sin(theta_ssi/2),0,lc,28)
motor.addPoint(Rs*cos(theta_ssi/2),Rs*sin(theta_ssi/2),0,lc,29)
motor.addPoint(Rs*cos(theta_t - theta_ssi/2),Rs*sin(theta_t - theta_ssi/2),0,lc,30)
motor.addPoint(Rssi*cos(theta_t - theta_ssi/2),Rssi*sin(theta_t - theta_ssi/2),0,lc,31)
motor.addPoint(Rssi*cos(theta_t - theta_sso/2),Rssi*sin(theta_t - theta_sso/2),0,lc,32)
motor.addPoint(Rsy*cos(theta_t - theta_sso/2),Rsy*sin(theta_t - theta_sso/2),0,lc,33)

'''Domain Boundary'''
motor.addPoint(D,0,0,lc,34) # Inner Radius of Rotor along "+" x-axis
motor.addPoint(0,D,0,lc,35)  # Inner Radius of Rotor along "+" y-axis
motor.addPoint(-D,0,0,lc,36)  # Inner Radius of Rotor along "-" x-axis

# Stator Winding Corner In Mid-Section (defines line between the stator windings)
motor.addPoint(Rssi,0,0,lc,37) # Same point on 'x' axis exists in the last index of array stator_winding_edge


# ----------------------------------------------------------------------------------
'''

What we want to do is minimize the number of domains we have to individually set in FEniCS.
Instead of slicing the domains for the rotor and stator (as done previously) we want to
rotate the defining points of each slice for the domains that are connected, which are:
- Rotor Core
- Stator Core
- Air Gap outside of stator
- Air Gap inside rotor
- Air gap between rotor & stator
After or during rotation, we can connect lines for each of these boundaries, and then 
connect as one whole domain

The areas that will need rotation of domains are:
- Magnets
- Windings
- Air Slots next to Magnets
(There is no way around this as these components are separated and not connected)

'''
# ----------------------------------------------------------------------------------
'''ROTOR CORE'''
rotor_points = [0 for i in range(p)] # 14 points define each rotor slice, with 8 pole pairs
rotor_out_bound_lines = [0 for i in range(5*p+3)]

rotor_points[0] = [(0,6),(0,7),(0,8),(0,9),(0,10)]
# rotor_out_bound_lines[0] = [(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)]

'''Rotor Point Rotation'''
for n in range(p-1):
    rotor_points[n+1] = motor.copy([(0,6),(0,7),(0,8),(0,9),(0,10),(0,11),(0,12),(0,13),(0,14),(0,15),(0,16),(0,17),(0,18),(0,19),(0,20)])
    for i in range(len(rotor_points[0])):
        motor.rotate([(rotor_points[n+1][i])],0,0,0,0,0,1,pi*(n+1)/p)
        # can't print lists so we need to loop through each index within each list 

'''Rotor Outer Boundary Lines'''
for n in range(p):
    if n < (p-1):
        rotor_end_line_slice = motor.addCircleArc(rotor_points[n][4][1],1,rotor_points[n+1][0][1],5*n + 5)
    else:
        rotor_end_line_slice = motor.addCircleArc(rotor_points[n][4][1],1,5,5*n + 5)
        # Point 5 is the limit of the outer rotor surface on the "-" x-axis
    # motor.addCircleArc(rotor_points[n][0][1],1,rotor_points[n][1][1],5*n + 1),
    # motor.addLine(rotor_points[n][1][1],rotor_points[n][2][1],5*n + 2),
    # motor.addCircleArc(rotor_points[n][2][1],1,rotor_points[n][3][1],5*n + 3),
    # motor.addLine(rotor_points[n][3][1],rotor_points[n][4][1],5*n + 4)

    rotor_out_bound_lines[5*n + 0] = motor.addCircleArc(rotor_points[n][0][1],1,rotor_points[n][1][1],5*n + 1)
    rotor_out_bound_lines[5*n + 1] = motor.addLine(rotor_points[n][1][1],rotor_points[n][2][1],5*n + 2)
    rotor_out_bound_lines[5*n + 2] = motor.addCircleArc(rotor_points[n][2][1],1,rotor_points[n][3][1],5*n + 3)
    rotor_out_bound_lines[5*n + 3] = motor.addLine(rotor_points[n][3][1],rotor_points[n][4][1],5*n + 4)
    rotor_out_bound_lines[5*n + 4] = rotor_end_line_slice

    # rotor_out_bound_lines[n] = [
    #     motor.addCircleArc(rotor_points[n][0][1],1,rotor_points[n][1][1],5*n + 1),
    #     motor.addLine(rotor_points[n][1][1],rotor_points[n][2][1],5*n + 2),
    #     motor.addCircleArc(rotor_points[n][2][1],1,rotor_points[n][3][1],5*n + 3),
    #     motor.addLine(rotor_points[n][3][1],rotor_points[n][4][1],5*n + 4),asdf
    #     ]
    # print(rotor_points[n][4][1])
    # print(rotor_points[n+1][0][1])
    # print(rotor_out_bound_lines[n])
    # print(rotor_out_bound_lines[n][4])
    # if n < (p-1):
    #     rotor_out_bound_lines[n][4] = motor.addCircleArc(rotor_points[n][4][1],1,rotor_points[n+1][0][1],p*n + 5)
        # asdf = motor.addCircleArc(rotor_points[n][4][1],1,rotor_points[n+1][0][1],p*n + 5)
    # else:
    #     rotor_out_bound_lines[n][4] = motor.addCircleArc(rotor_points[n][3][1],1,5,p*n + 5) 

'''Rotor Inner Boundary Lines'''
rotor_out_bound_lines[-3] = motor.addLine(5,4)
rotor_out_bound_lines[-2] = motor.addCircleArc(4,1,2)
rotor_out_bound_lines[-1] = motor.addLine(2,6)

'''Rotor Curve Loop & Surface'''
motor.addCurveLoop(rotor_out_bound_lines,1)
motor.addPlaneSurface([1],1)

# ----------------------------------------------------------------------------------

'''STATOR CORE'''
stator_points = [0 for i in range(s)] # 14 points define each rotor slice, with 9 pole pairs
stator_in_bound_lines = [0 for i in range(9*s+3)]

stator_points[0] = [(0,25),(0,26),(0,27),(0,28),(0,29),(0,30),(0,31),(0,32),(0,33)]

'''Stator Point Rotation'''
for n in range(s-1):
    stator_points[n+1] = motor.copy([(0,25),(0,26),(0,27),(0,28),(0,29),(0,30),(0,31),(0,32),(0,33)])
    for i in range(len(stator_points[0])):
        motor.rotate([(stator_points[n+1][i])],0,0,0,0,0,1,pi*(n+1)/s)

'''Stator Inner Boundary Lines'''
for n in range(s):
    if n < (s-1):
        stator_end_line_slice = motor.addCircleArc(stator_points[n][8][1],1,stator_points[n+1][0][1])
    else:
        stator_end_line_slice = motor.addCircleArc(stator_points[n][8][1],1,24)
        # Point 24 is the last point defining stator core on "-" y-axis to avoid duplicates

    stator_in_bound_lines[9*n + 0] = motor.addCircleArc(stator_points[n][0][1],1,stator_points[n][1][1])
    stator_in_bound_lines[9*n + 1] = motor.addLine(stator_points[n][1][1],stator_points[n][2][1])
    stator_in_bound_lines[9*n + 2] = motor.addCircleArc(stator_points[n][2][1],1,stator_points[n][3][1])
    stator_in_bound_lines[9*n + 3] = motor.addLine(stator_points[n][3][1],stator_points[n][4][1])
    stator_in_bound_lines[9*n + 4] = motor.addCircleArc(stator_points[n][4][1],1,stator_points[n][5][1])
    stator_in_bound_lines[9*n + 5] = motor.addLine(stator_points[n][5][1],stator_points[n][6][1])
    stator_in_bound_lines[9*n + 6] = motor.addCircleArc(stator_points[n][6][1],1,stator_points[n][7][1])
    stator_in_bound_lines[9*n + 7] = motor.addLine(stator_points[n][7][1],stator_points[n][8][1])
    stator_in_bound_lines[9*n + 8] = stator_end_line_slice

'''Stator Outer Boundary Lines'''
stator_in_bound_lines[-3] = motor.addLine(24,23)
stator_in_bound_lines[-2] = motor.addCircleArc(23,1,21)
stator_in_bound_lines[-1] = motor.addLine(21,25)

'''Stator Curve Loop & Surface'''
motor.addCurveLoop(stator_in_bound_lines,2)
motor.addPlaneSurface([2],2)

# ----------------------------------------------------------------------------------
'''AIR GAP'''
# Outer Air Gap
outer_air_gap = [
    motor.addLine(21,34),
    motor.addCircleArc(34,1,35),
    motor.addCircleArc(35,1,36),
    motor.addLine(36,23),
    stator_in_bound_lines[-2]
]
motor.addCurveLoop(outer_air_gap,3)
motor.addPlaneSurface([3],3)

# Inner Air Gap / Shaft
inner_air_gap = [
    motor.addLine(1,2),
    motor.addCircleArc(2,1,3),
    motor.addCircleArc(3,1,4),
    motor.addLine(4,1),
]
motor.addCurveLoop(inner_air_gap,4)
motor.addPlaneSurface([4],4)

# Mid Air Gap (first boundary)
stator_winding_edge = [0 for i in range(s+1)]
stator_winding_edge[0] = ([(0,37)])
stator_winding_edge_lines = [0 for i in range(2*(s))]

for n in range(s):
    stator_winding_edge[n+1] = motor.copy([(0,37)])
    motor.rotate(stator_winding_edge[n+1],0,0,0,0,0,1,pi*(n+1)/s)
    stator_winding_edge_lines[2*n] = motor.addCircleArc(stator_winding_edge[n][0][1],1,stator_points[n][3][1])
    stator_winding_edge_lines[2*n+1] = motor.addCircleArc(stator_points[n][6][1],1,stator_winding_edge[n+1][0][1])

outer_mid_air_gap_lines = [0 for i in range(5*s+2)]

# Computing BACKWARDS from 180 degrees to zero degrees bc outer edge of rotor is defined from 0 to 180
for n in range(s,0,-1): 
    outer_mid_air_gap_lines[5*(s-n)+1] = -stator_winding_edge_lines[2*(n-1)+1]
    outer_mid_air_gap_lines[5*(s-n)+2] = -stator_in_bound_lines[9*(n-1)+5]
    outer_mid_air_gap_lines[5*(s-n)+3] = -stator_in_bound_lines[9*(n-1)+4]
    outer_mid_air_gap_lines[5*(s-n)+4] = -stator_in_bound_lines[9*(n-1)+3]
    outer_mid_air_gap_lines[5*(s-n)+5] = -stator_winding_edge_lines[2*(n-1)]
    # We leave the first and last index of outer_mid_air_gap_lines open to connect the lines for the edge of the mid air gap

outer_mid_air_gap_lines[0] = motor.addLine(5,stator_winding_edge[-1][0][1])
outer_mid_air_gap_lines[-1] = motor.addLine(37,6)

outer_mid_air_gap_boundary = [0 for i in range(len(rotor_out_bound_lines[0:-3])+len(outer_mid_air_gap_lines))]

for n in range(len(rotor_out_bound_lines[0:-3])):
    outer_mid_air_gap_boundary[n] = rotor_out_bound_lines[n]
for n in range(len(outer_mid_air_gap_lines)):
    outer_mid_air_gap_boundary[n+len(rotor_out_bound_lines[0:-3])] = outer_mid_air_gap_lines[n]
# for n in range(len(outer_mid_air_gap_lines)):
#     outer_mid_air_gap_boundary[n] = outer_mid_air_gap_lines[n]
# for n in range(len(rotor_out_bound_lines[0:-4])):
#     outer_mid_air_gap_boundary[n+len(outer_mid_air_gap_lines)] = rotor_out_bound_lines[n]

'''Mid Air Gap Loop & Surface'''
motor.addCurveLoop(outer_mid_air_gap_boundary,5)
motor.addPlaneSurface([5],5)

# ----------------------------------------------------------------------------------
'''ROTOR MAGNET/AIR SLOTS'''
rotor_slot_points = [0 for i in range(p)]
rotor_slot_points[0] = [(0,11),(0,12),(0,13),(0,14),(0,15),(0,16),(0,17),(0,18),(0,19),(0,20)]

for n in range(p-1):
    rotor_slot_points[n+1] = motor.copy([(0,11),(0,12),(0,13),(0,14),(0,15),(0,16),(0,17),(0,18),(0,19),(0,20)])
    for i in range(len(rotor_slot_points[0])):
        motor.rotate([(rotor_slot_points[n+1][i])],0,0,0,0,0,1,pi*(n+1)/p)

magnet_lines = [0 for i in range(4*p)]
right_air_slot_lines = [0 for i in range(4*p)]
left_air_slot_lines = [0 for i in range(4*p)]

for i in range(p):
    magnet_lines[4*i+0] = motor.addLine(rotor_slot_points[i][2][1],rotor_slot_points[i][4][1]) # Point 13 to 15
    magnet_lines[4*i+1] = motor.addCircleArc(rotor_slot_points[i][4][1],1,rotor_slot_points[i][7][1]) # Point 15 to 18
    magnet_lines[4*i+2] = motor.addLine(rotor_slot_points[i][7][1],rotor_slot_points[i][5][1]) # Point 18 to 16
    magnet_lines[4*i+3] = motor.addCircleArc(rotor_slot_points[i][5][1],1,rotor_slot_points[i][2][1]) # Point 16 to 13
    ml_loop = magnet_lines[(4*i+0):(4*i+4)]
    motor.addCurveLoop(ml_loop,5+i+1) # Last index for curve loop is 5 + p
    motor.addPlaneSurface([5+i+1],5+i+1) # Last index for surface is 5 + p

for i in range(p): # uses points 11, 12, 15, 14
    right_air_slot_lines[4*i+0] = motor.addLine(rotor_slot_points[i][0][1],rotor_slot_points[i][1][1])
    right_air_slot_lines[4*i+1] = motor.addLine(rotor_slot_points[i][1][1],rotor_slot_points[i][4][1])
    right_air_slot_lines[4*i+2] = motor.addLine(rotor_slot_points[i][4][1],rotor_slot_points[i][3][1])
    right_air_slot_lines[4*i+3] = motor.addCircleArc(rotor_slot_points[i][3][1],1,rotor_slot_points[i][0][1])
    rasl_loop = right_air_slot_lines[(4*i+0):(4*i+4)]
    motor.addCurveLoop(rasl_loop,5+i+1+p) # Last index for curve loop is 5 + 2*p
    motor.addPlaneSurface([5+i+1+p],5+i+1+p) # Last index for surface is 5 + 2*p

for i in range(p): # uses points 18, 20, 19, 17
    left_air_slot_lines[4*i+0] = motor.addLine(rotor_slot_points[i][6][1],rotor_slot_points[i][7][1])
    left_air_slot_lines[4*i+1] = motor.addLine(rotor_slot_points[i][7][1],rotor_slot_points[i][9][1])
    left_air_slot_lines[4*i+2] = motor.addLine(rotor_slot_points[i][9][1],rotor_slot_points[i][8][1])
    left_air_slot_lines[4*i+3] = motor.addCircleArc(rotor_slot_points[i][8][1],1,rotor_slot_points[i][6][1])
    lasl_loop = left_air_slot_lines[(4*i+0):(4*i+4)]
    bbb = motor.addCurveLoop(lasl_loop,5+i+1+2*p) # Last index for curve loop is 5 + 3*p
    motor.addPlaneSurface([5+i+1+2*p],5+i+1+2*p) # Last index for surface is 5 + 3*p
    # print(bbb)


rotor_slot_bool = []
for i in range(3*p):
    rotor_slot_bool.append([])
    rotor_slot_bool[i] = (2,5+i+1)
    # rotor_slot_bool[p+i] = (2,5+i+1+p)
    # rotor_slot_bool[2*p+i] = (2,5+i+1+2*p)

# print(rotor_slot_bool)



# rotor_slot_lines = [0 for i in range(12*p)]
# First 10 spots represent boundary of rotor slots, last two represent the lines separating air gaps and magnet (right first, then left)
# for i in range(p):
#     rotor_slot_lines[12*i+0] = motor.addLine(rotor_slot_points[i][0][1],rotor_slot_points[i][1][1])
#     rotor_slot_lines[12*i+1] = motor.addLine(rotor_slot_points[i][1][1],rotor_slot_points[i][4][1])
#     rotor_slot_lines[12*i+2] = motor.addCircleArc(rotor_slot_points[i][4][1],1,rotor_slot_points[i][7][1])
#     rotor_slot_lines[12*i+3] = motor.addLine(rotor_slot_points[i][7][1],rotor_slot_points[i][9][1])
#     rotor_slot_lines[12*i+4] = motor.addLine(rotor_slot_points[i][9][1],rotor_slot_points[i][8][1])
#     rotor_slot_lines[12*i+5] = motor.addCircleArc(rotor_slot_points[i][8][1],1,rotor_slot_points[i][6][1])
#     rotor_slot_lines[12*i+6] = motor.addLine(rotor_slot_points[i][6][1],rotor_slot_points[i][5][1])
#     rotor_slot_lines[12*i+7] = motor.addCircleArc(rotor_slot_points[i][5][1],1,rotor_slot_points[i][2][1])
#     rotor_slot_lines[12*i+8] = motor.addLine(rotor_slot_points[i][2][1],rotor_slot_points[i][3][1])
#     rotor_slot_lines[12*i+9] = motor.addCircleArc(rotor_slot_points[i][3][1],1,rotor_slot_points[i][0][1])
#     rotor_slot_lines[12*i+10] = motor.addLine(rotor_slot_points[i][3][1],rotor_slot_points[i][4][1]) # Straight line separating right air slot (going up)
#     rotor_slot_lines[12*i+11] = motor.addLine(rotor_slot_points[i][6][1],rotor_slot_points[i][7][1]) # Straight line separating left air slot (going up)
#     print(rotor_slot_lines[(0+12*i):(10+12*i)])
#     asdf = rotor_slot_lines[(0+12*i):(10+12*i)]
#     motor.addCurveLoop(asdf,5+i+1) # Last index for curve loop is 5 + p
#     motor.addPlaneSurface([5+i+1],5+i+1) # Last index for surface is 5 + p

# rotor_slot_bool = []
# for i in range(p):
#     rotor_slot_bool.append([])
#     rotor_slot_bool[i] = (2,5+i+1)

# print(rotor_slot_bool)


# print(aaa)

# Stator Windings
right_stator_winding_lines = [0 for i in range(4*s)] # Right side of stator tooth
left_stator_winding_lines = [0 for i in range(4*s)] # Left side of stator tooth

for i in range(s):
    # print(i)
    right_stator_winding_lines[4*i+0] = motor.addLine(stator_winding_edge[i][0][1],stator_points[i][0][1])
    right_stator_winding_lines[4*i+1] = motor.addCircleArc(stator_points[i][0][1],1,stator_points[i][1][1])
    right_stator_winding_lines[4*i+2] = motor.addLine(stator_points[i][1][1],stator_points[i][2][1])
    right_stator_winding_lines[4*i+3] = motor.addCircleArc(stator_points[i][2][1],1,stator_winding_edge[i][0][1])

    left_stator_winding_lines[4*i+0] = motor.addCircleArc(stator_winding_edge[i+1][0][1],1,stator_points[i][7][1])
    left_stator_winding_lines[4*i+1] = motor.addLine(stator_points[i][7][1],stator_points[i][8][1])
    if i < (s-1):
        left_stator_winding_lines[4*i+2] = motor.addCircleArc(stator_points[i][8][1],1,stator_points[i+1][0][1])
        left_stator_winding_lines[4*i+3] = motor.addLine(stator_points[i+1][0][1],stator_winding_edge[i+1][0][1])
    else:
        left_stator_winding_lines[4*i+2] = motor.addCircleArc(stator_points[i][8][1],1,24)
        left_stator_winding_lines[4*i+3] = motor.addLine(24,stator_winding_edge[i+1][0][1])
    
    rswl_loop = right_stator_winding_lines[(4*i+0):(4*i+4)]
    lswl_loop = left_stator_winding_lines[(4*i+0):(4*i+4)]

    motor.addCurveLoop(rswl_loop,5+3*p+2+2*i)
    motor.addCurveLoop(lswl_loop,5+3*p+2+2*i+1) # Last index for curve loop is 5 + 3p + 2s + 1

    motor.addPlaneSurface([5+3*p+2+2*i],5+3*p+2+2*i) 
    motor.addPlaneSurface([5+3*p+2+2*i+1],5+3*p+2+2*i+1) # Last index for surface is 5 + 3p + 2s + 1

aaa = motor.cut([(2,1)],rotor_slot_bool,5 + 3*p + 2*s + 2,removeObject=True,removeTool=False)
    

# lasl_loop = left_air_slot_lines[(4*i+0):(4*i+4)]
#     motor.addCurveLoop(lasl_loop,5+i+1+2*p) # Last index for curve loop is 5 + 3*p
#     motor.addPlaneSurface([5+i+1+2*p],5+i+1+2*p) # Last index for surface is 5 + 3*p

gmsh.model.occ.removeAllDuplicates()

motor.synchronize()

gmsh.model.addPhysicalGroup(2,[5 + 3*p + 2*s + 2],1)
gmsh.model.setPhysicalName(2,1,"Rotor Core")

gmsh.model.addPhysicalGroup(1,[rotor_out_bound_lines],1)
gmsh.model.setPhysicalName(1,1,"Rotor Core Boundary")

gmsh.model.addPhysicalGroup(2,[2],2)
gmsh.model.setPhysicalName(2,2,"Stator Core")

gmsh.model.addPhysicalGroup(1,[stator_in_bound_lines],2)
gmsh.model.setPhysicalName(1,2,"Stator Core Boundary")

gmsh.model.addPhysicalGroup(2,[3],3)
gmsh.model.setPhysicalName(2,3,"Outer Motor Domain")

gmsh.model.addPhysicalGroup(1,[outer_air_gap],3)
gmsh.model.setPhysicalName(1,3,"Outer Motor Domain Boundary")

gmsh.model.addPhysicalGroup(2,[4],4)
gmsh.model.setPhysicalName(2,4,"Inner Motor Domain")

gmsh.model.addPhysicalGroup(1,[inner_air_gap],4)
gmsh.model.setPhysicalName(1,4,"Inner Motor Domain Boundary")

gmsh.model.addPhysicalGroup(2,[5],5)
gmsh.model.setPhysicalName(2,5,"Mid Air Gap")

gmsh.model.addPhysicalGroup(1,[outer_mid_air_gap_boundary],5)
gmsh.model.setPhysicalName(1,5,"Mid Air Gap Boundary")

# Need groups for magnets/air slots via for loop
# Air Slots (first 4 line-pairs of code)
for i in range(p):
    gmsh.model.addPhysicalGroup(2,[5+i+1+p],5+(i+1))
    gmsh.model.setPhysicalName(2,5+(i+1),"Right Air Slot Domain " + str(i+1))

    gmsh.model.addPhysicalGroup(1,right_air_slot_lines[(4*i+0):(4*i+4)],5+(i+1))
    gmsh.model.setPhysicalName(1,5+(i+1),"Right Air Slot Boundary " + str(i+1))

    gmsh.model.addPhysicalGroup(2,[5+i+1+2*p],5+(i+1)+p)
    gmsh.model.setPhysicalName(2,5+(i+1)+p,"Left Air Slot Domain " + str(i+1))

    gmsh.model.addPhysicalGroup(1,left_air_slot_lines[(4*i+0):(4*i+4)],5+(i+1)+p)
    gmsh.model.setPhysicalName(1,5+(i+1)+p,"Left Air Slot Boundary " + str(i+1))

    gmsh.model.addPhysicalGroup(2,[5+i+1],5+(i+1)+2*p)
    gmsh.model.setPhysicalName(2,5+(i+1)+2*p,"Magnet Domain " + str(i+1))

    gmsh.model.addPhysicalGroup(1,magnet_lines[(4*i+0):(4*i+4)],5+(i+1)+2*p)
    gmsh.model.setPhysicalName(1,5+(i+1)+2*p,"Magnet Boundary " + str(i+1))

# need groups for stator windings (unfortunately done manually)

for i in range(s):
    gmsh.model.addPhysicalGroup(2,[5+3*p+2+2*i],5+(i+1)+3*p)
    gmsh.model.setPhysicalName(2,5+(i+1)+3*p,"Right Stator Winding Domain " + str(i+1))

    gmsh.model.addPhysicalGroup(1,right_stator_winding_lines[(4*i+0):(4*i+4)],5+(i+1)+3*p)
    gmsh.model.setPhysicalName(1,5+(i+1)+3*p,"Right Stator Winding Boundary " + str(i+1))

    gmsh.model.addPhysicalGroup(2,[5+3*p+2+2*i+1],5+(i+1)+3*p+s)
    gmsh.model.setPhysicalName(2,5+(i+1)+3*p+s,"Left Stator Winding Domain " + str(i+1))

    gmsh.model.addPhysicalGroup(1,left_stator_winding_lines[(4*i+0):(4*i+4)],5+(i+1)+3*p+s)
    gmsh.model.setPhysicalName(1,5+(i+1)+3*p+s,"Left Stator Winding Boundary " + str(i+1))

gmsh.model.mesh.generate(2)

gmsh.write("motor_mesh.msh")

gmsh.finalize()

# -------------------------------------------------------------------------------------------
# os.system('python3 msh2xdmf.py -d 2 motor_mesh.msh')
# THIS LINE WORKS only when called using the general terminal window.
# The msh2xdmf call to terminal window in VS Code doesn't work here; errors 
# say that meshio and numpy are not available packages