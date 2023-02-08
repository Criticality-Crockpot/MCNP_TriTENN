# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 22:03:53 2023

@author: Nicholas Branam
@contact: nbranam@vols.utk.edu
@version: 1.0
"""
from math import pi
import numpy as np
import datetime
import matplotlib.pyplot as plt

filename = 'STA_mcnp_input.txt'
#%% Dimensions of STA

deck_title = 'STA MCNP Input Deck'
center_distance = 14.5 # distance from center of cylinders to center of assembly
# assuming container is cylindrical and has an h/d ratio of 1
liters_in_container = 2
volume_of_container = liters_in_container*1000 # cm^3
gpl = 100
pp = 2

# absorber: choose 1 for Boron-10, 2 for Aluminum, 3 for Air
absorber = 3

# Concentrations
h = 0.60689089
n = 0.02394325
o = 0.36916586
p239 = 1.20317020e-02
p240 = 5.30869186e-04
p241 = 3.62446221e-05
p242 = 1.99065079e-06

con_h =  (4*volume_of_container/pi)**(1./3) # height of container
con_r = con_h/2 # radius of container


rail_w = 22  # rail width
rail_l = 50  # rail length
rail_h = 10  # rail height
rail_edge = 10

plat_s = 18  # platform side length (square)
plat_h = 3  # platform height

ab_thick = 1 # absorber thickness
ab_w = plat_s # absorber width (same as platform width)
ab_h = con_h # absorber height (same as container height)

lip_h = 2 # fastener height
lip_r = con_r + 1

# k code
num_part = 1000 # number of particles per generation
num_skip = 10 #  number of generations skipped
num_tot_gen = 50 # number of total generations

#%% Header 

header = []
e = datetime.datetime.now()
x = str(e.strftime("%m-%d-%Y %H:%M:%S"))
date_and_time = 'c date and time: ' + x
header.append(deck_title)
header.append('c ╭━━━━╮╱╱╱╱╭━━━━┳━━━┳━╮╱╭┳━╮╱╭╮')
header.append('c ┃╭╮╭╮┃╱╱╱╱┃╭╮╭╮┃╭━━┫┃╰╮┃┃┃╰╮┃┃')
header.append('c ╰╯┃┃┣┻┳╮╱╱╰╯┃┃╰┫╰━━┫╭╮╰╯┃╭╮╰╯┃')
header.append('c ╱╱┃┃┃╭╋┫╭━━╮┃┃╱┃╭━━┫┃╰╮┃┃┃╰╮┃┃')
header.append('c ╱╱┃┃┃┃┃┃╰━━╯┃┃╱┃╰━━┫┃╱┃┃┃┃╱┃┃┃')
header.append('c ╱╱╰╯╰╯╰╯╱╱╱╱╰╯╱╰━━━┻╯╱╰━┻╯╱╰━╯')
header.append('c      Senior Design Project')
header.append('c Author: Nicholas Branam')
header.append(date_and_time)
header.append('c ')


#%% Surface Card

surface_values = []
surface_values.append('c Surface')
surface_values.append('1 px -'+str((np.round(rail_l+rail_edge,2)))+'    $ west edge of rail')
surface_values.append('2 px -' + str(np.round(center_distance+plat_s/2,2)) + '    $ west edge of platform')
surface_values.append('3 px -' + str(np.round(ab_thick+center_distance-plat_s/2,2)) + '    $ west edge of abs (west edge of platform minus abs thickness)')
surface_values.append('4 px -' + str(np.round(center_distance-plat_s/2,2)) + '    $ east edge of platform and of abs')
surface_values.append('5 px -'+ str(np.round(rail_edge,2)) +'   $ east edge of rail')
surface_values.append('6 py -'+ str(np.round(rail_w/2,2)) +'    $ south edge of rail')
surface_values.append('7 py -'+ str(np.round(plat_s/2,2)) +'     $ south edge of platform and abs')
surface_values.append('8 py '+ str(np.round(plat_s/2,2)) +'      $ north edge of platform and abs')
surface_values.append('9 py ' + str(np.round(rail_w/2,2)) + '     $ north edge of rail')
surface_values.append('10 pz 0     $ bottom of rail')
surface_values.append('11 pz ' + str(np.round(rail_h,2)) + '    $ top of rail, bottom of platform')
surface_values.append('12 pz ' + str(np.round(plat_h+rail_h,2)) + '    $ top of platform, bottom of lip and container')
surface_values.append('13 pz ' + str(np.round(plat_h+lip_h+rail_h,2)) + '    $ top of lip')
surface_values.append('14 pz ' + str(np.round(plat_h+con_h+rail_h,2)) + ' $ top of container and abs')
surface_values.append('15 c/z -' + str(np.round(center_distance,2)) + ' 0 '+str(np.round(con_r,2))+' $ center of container in xy and container radius')
surface_values.append('16 c/z -' + str(np.round(center_distance,2)) + ' 0 '+str(np.round(lip_r,2))+' $ center of lip in xy and lip radius')
surface_values.append('17 sph 0 0 0 100')
surface_values.append('')


#%% 
# listed in atomic fractions
# convert to WF, calculate the amount of Pu in the system to credible

# Convert to weight fractions here
# Adjust to credible concentrations
OG_AF = np.sum(np.array([h,n,o,p239,p240,p241,p242]))

denom = h*1 + o*16 + n*14 + p239*239 + p240*240 + p241*241 + p242*242
num = np.array([h*1, o*16, n*14, p239*239, p240*240, p241*241, p242*242])
weight_fractions = num/denom
total_plutonium_WF = np.sum([weight_fractions[3:7]])
density_plutonium = 19840 #g/L
grams_plutonium = total_plutonium_WF*density_plutonium*volume_of_container # grams
total_nitric_acid = np.sum([weight_fractions[0:3]])
density_nitric_acid = 1510 #g/L
grams_nitric_acid = total_nitric_acid*density_nitric_acid*volume_of_container # grams


Pu_atomic_fractions = np.array([p239, p240, p241, p242])
Pu_total_atomic_fractions = np.sum(Pu_atomic_fractions)
Pu_ratios = Pu_atomic_fractions / Pu_total_atomic_fractions

HNO3_atomic_fractions = np.array([h, n, o])
HNO3_total_atomic_fractions = np.sum(HNO3_atomic_fractions)
HNO3_ratios = HNO3_atomic_fractions /  HNO3_total_atomic_fractions
soln_den = 1

#%% USE MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
# assume 1L of solution
grams_plutonium_calc = gpl*liters_in_container # g in 1 L of solution
weight_fraction_plutonium_calc = grams_plutonium_calc / density_plutonium
weight_fraction_Pu_ratios = Pu_ratios * weight_fraction_plutonium_calc

# 1L of HNO3 is 1L*1510g/L
grams_HNO3_calc = 1510*liters_in_container # g of total solution
weight_fraction_HNO3_calc = grams_HNO3_calc / density_nitric_acid
weight_fraction_HNO3_ratios = HNO3_ratios * weight_fraction_HNO3_calc

MCNP_input_WF = np.append(weight_fraction_HNO3_ratios,weight_fraction_Pu_ratios)


#%% Data Card
data_values = []
data_values.append('c Data')
data_values.append('*tr1 0 0 0 120 210 90 30 120 90 90 90 0')
data_values.append('*tr2 0 0 0 240 330 90 150 240 90 90 90 0')
data_values.append('c Plutonium Nitrate Solution')   
data_values.append('m100 01001 -' + str(np.round(MCNP_input_WF[0],6)))
data_values.append('     07014 -' + str(np.round(MCNP_input_WF[1],6)))
data_values.append('     08016 -' + str(np.round(MCNP_input_WF[2],6)))
data_values.append('     94239 -' + str(np.round(MCNP_input_WF[3],6)))
data_values.append('     94240 -' + str(np.round(MCNP_input_WF[4],6)))
data_values.append('     94241 -' + str(np.round(MCNP_input_WF[5],6)))
data_values.append('     94242 -' + str(np.round(MCNP_input_WF[6],6)))
data_values.append('c Absorbing Material')
if absorber == 1:
    data_values.append('c Boron-10')
    data_values.append('m200 05010 -1.000')
    ab_den = 2.37
elif absorber == 2:
    data_values.append('c Aluminum')
    data_values.append('m200 13027 -1.000')
    ab_den = 2.7
elif absorber == 3:
    data_values.append('c Air (Dry, Near Sea Level)')
    data_values.append('m200 06000 -0.000124')
    data_values.append('     07014 -0.775268')
    data_values.append('     08016 -0.231781')
    data_values.append('     18000 -0.012827')
    ab_den = 0.001205
data_values.append('c Platform and Rail Material (steel)')
data_values.append('m300 06000 -0.0050')
data_values.append('     26000 -0.9950')
if absorber != 2:
    data_values.append('c Lip Material')
    data_values.append('m400 13027 -1.000')
data_values.append('mode n')
data_values.append('kcode ' + str(int(num_part))+ ' 1.0 '+str(int(num_skip))+' '+str(int(num_tot_gen)))
z = (np.round(plat_h+con_h+rail_h,2)+np.round(plat_h+rail_h,2))/2
k_pos1 = ' -'+str(np.round(center_distance,pp)) + ' 0 ' + str(np.round(z,pp))
k_pos2 = ' '+str(np.round(center_distance*np.cos(pi/3),pp))+' '+str(np.round(center_distance*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
k_pos3 = ' '+str(np.round(center_distance*np.cos(pi/3),pp))+' -'+str(np.round(center_distance*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
data_values.append('ksrc'+k_pos1+k_pos2+k_pos3)

extra = 5.5
Drawing_colored_circle1 = plt.Circle(( -(np.round(center_distance,pp)) , 0 ), con_r+extra,color='k' )
Drawing_colored_circle2 = plt.Circle(( np.round(center_distance*np.cos(pi/3),pp) , np.round(center_distance*np.cos(pi/6),pp) ), con_r+extra,color='k' )
Drawing_colored_circle3 = plt.Circle(( np.round(center_distance*np.cos(pi/3),pp) , -np.round(center_distance*np.cos(pi/6),pp) ), con_r+extra,color='k' )
Drawing_colored_circle4 = plt.Circle(( -(np.round(center_distance,pp)) , 0 ), con_r )
Drawing_colored_circle5 = plt.Circle(( np.round(center_distance*np.cos(pi/3),pp) , np.round(center_distance*np.cos(pi/6),pp) ), con_r )
Drawing_colored_circle6 = plt.Circle(( np.round(center_distance*np.cos(pi/3),pp) , -np.round(center_distance*np.cos(pi/6),pp) ), con_r )
fig, ax = plt.subplots()
dimen = 50
ax.set_xlim((-dimen, dimen))
ax.set_ylim((-dimen, dimen))
ax.add_patch(Drawing_colored_circle1)
ax.add_patch(Drawing_colored_circle2)
ax.add_patch(Drawing_colored_circle3)
ax.add_patch(Drawing_colored_circle4)
ax.add_patch(Drawing_colored_circle5)
ax.add_patch(Drawing_colored_circle6)
ax.set_aspect('equal', adjustable='box')
plt.show()
#%% Cell Card

cell_values = []
cell_values.append('c Cells')
cell_values.append('5 100 -'+str(np.round(soln_den,4)) +' -15 12 -14 imp:n=1         $ container')
cell_values.append('6 LIKE 5 BUT TRCL=1')
cell_values.append('7 LIKE 5 BUT TRCL=2')
if absorber == 2:
    cell_values.append('10 200 -2.7 15 -16 12 -13 imp:n=1     $ lip')
else:
    cell_values.append('10 400 -2.7 15 -16 12 -13 imp:n=1     $ lip')

cell_values.append('11 LIKE 10 BUT TRCL=1')
cell_values.append('12 LIKE 10 BUT TRCL=2')
cell_values.append('15 300 -7.82 2 -4 7 -8 11 -12 imp:n=1  $ platform')
cell_values.append('16 LIKE 15 BUT TRCL=1')
cell_values.append('17 LIKE 15 BUT TRCL=2')
cell_values.append('20 300 -7.82 1 -5 6 -9 10 -11 imp:n=1  $ rail')
cell_values.append('21 LIKE 20 BUT TRCL=1')
cell_values.append('22 LIKE 20 BUT TRCL=2')
cell_values.append('25 200 -' + str(np.round(ab_den,6)) + ' 3 -4 7 -8 12 -14 imp:n=1  $ absorber')
cell_values.append('26 LIKE 25 BUT TRCL=1')
cell_values.append('27 LIKE 25 BUT TRCL=2')
cell_values.append('100 0 -17 #5 #6 #7 #10 #11 #12')
cell_values.append('     #15 #16 #17 #20 #21 #22')
cell_values.append('     #25 #26 #27 imp:n=1 $ space around')
cell_values.append('999 0 17 imp:n=0')
cell_values.append('')


lines = [header,cell_values,surface_values,data_values]
with open(filename,'w',encoding="utf-8") as file:
    for inte,line in enumerate(lines):
        if type(line) == list:
            for i,j in enumerate(line):
                if type(j) == list:
                    for k,l in enumerate(j):
                        file.write(str(l))
                        file.write('\n')
                else:
                    file.write(j)
                    file.write('\n')
        else:
            file.write(line)
            file.write('\n')


#%%





