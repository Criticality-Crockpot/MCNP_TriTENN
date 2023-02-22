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
center_distance = 15 # distance from center of cylinders to center of assembly
# assuming container is cylindrical and has an h/d ratio of 1
liters_in_container = 2
h_to_d = 1
volume_of_container = liters_in_container*1000 # cm^3
gpl = 250 # Pu concentration
pp = 2
contact = 'yes'
walls = 'yes'

#if contact == 'yes':

# absorber: choose 1 for Boron-10, 2 for Aluminum, 3 for Air, 4 for Natural Boron,
# 5 for Cadmium, 6 for Gadolinium, 7 for Oak Ridge Concrete,
# 8 for Borated Polyethylene, 9 for Polyethylene
absorber = 6

# Concentrations
h = 0.60689089
n = 0.02394325
o = 0.36916586
p239 = 1.20317020e-02
p240 = 5.30869186e-04
p241 = 3.62446221e-05
p242 = 1.99065079e-06
soln_den = 1.24136


con_h = (h_to_d)*((4*volume_of_container)/(pi*h_to_d))**(1/3)
con_r = con_h/(h_to_d*2)
# con_h =  (4*volume_of_container/pi)**(1./3) # height of container
# con_r = con_h/2 # radius of container


rail_w = con_r*2+2  # rail width
rail_l = 50  # rail length
rail_h = 10  # rail height
rail_edge = 5

lip_h = 2 # fastener height
lip_r = con_r + 1

ex_th = 2
ab_thick = 1 # absorber thickness
plat_w = rail_w-1  # platform side length (square)
plat_l = lip_r*2 + ex_th + ab_thick
plat_h = 3  # platform height
ab_w = plat_w # absorber width (same as platform width)
ab_h = con_h # absorber height (same as container height)





# interaction check
threshold = 0.21
y_interaction_rail = -rail_edge*np.cos((pi/180)*30)+(rail_w/2)*np.sin((pi/180)*30)
if y_interaction_rail > -threshold:
    print('Error: Bounds of the rail overlap and MCNP will get very upset. :(')
    
y_interaction_plat = center_distance - plat_l/2
if (y_interaction_plat-threshold) < rail_edge:
    print('Error: Bounds of the platform overlap and MCNP will get very upset. :(')

# k code
num_part = 10000 # number of particles per generation
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
surface_values.append('2 px -' + str(np.round(center_distance+plat_l/2,2)) + '    $ west edge of platform')
surface_values.append('3 px -' + str(np.round(ab_thick+center_distance-plat_l/2,2)) + '    $ west edge of abs (west edge of platform minus abs thickness)')
surface_values.append('4 px -' + str(np.round(center_distance-plat_l/2,2)) + '    $ east edge of platform and of abs')
surface_values.append('5 px -'+ str(np.round(rail_edge,2)) +'   $ east edge of rail')
surface_values.append('6 py -'+ str(np.round(rail_w/2,2)) +'    $ south edge of rail')
surface_values.append('7 py -'+ str(np.round(plat_w/2,2)) +'     $ south edge of platform and abs')
surface_values.append('8 py '+ str(np.round(plat_w/2,2)) +'      $ north edge of platform and abs')
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

h = 0.60689089
n = 0.02394325
o = 0.36916586
p239 = 1.20317020e-02
p240 = 5.30869186e-04
p241 = 3.62446221e-05
p242 = 1.99065079e-06

#%% Concentrations Block
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
    filename = 'STA_mcnp_input_B10.txt'
elif absorber == 2:
    data_values.append('c Aluminum')
    data_values.append('m200 13027 -1.000')
    ab_den = 2.7
    filename = 'STA_mcnp_input_Al.txt'
elif absorber == 3:
    data_values.append('c Air (Dry, Near Sea Level)')
    data_values.append('m200 06000 -0.000124')
    data_values.append('     07014 -0.775268')
    data_values.append('     08016 -0.231781')
    data_values.append('     18000 -0.012827')
    ab_den = 0.001205
    filename = 'STA_mcnp_input_Air.txt'
elif absorber == 4:
    data_values.append('c Natural Boron')
    data_values.append('m200 05010 -0.20000')
    data_values.append('     05011 -0.80000')
    ab_den = 2.37
    filename = 'STA_mcnp_input_NatB.txt'
elif absorber == 5:
    data_values.append('c Cadmium')
    data_values.append('m200 48000 -1.000')
    ab_den = 8.65
    filename = 'STA_mcnp_input_Cd.txt'
elif absorber == 6:
    data_values.append('c Gadolinium')
    data_values.append('m200 64000 -1.000')
    ab_den = 7.9004
    filename = 'STA_mcnp_input_Gd.txt'
elif absorber == 7:
    data_values.append('c Oak Ridge Concrete')
    data_values.append('m200 01001 -0.006187')
    data_values.append('     06000 -0.175193')
    data_values.append('     08016 -0.410184')
    data_values.append('     11023 -0.000271')
    data_values.append('     11200 -0.032649')
    data_values.append('     13027 -0.010830')
    data_values.append('     14000 -0.034479')
    data_values.append('     19000 -0.001138')
    data_values.append('     20000 -0.321287')
    data_values.append('     26000 -0.007784')
    ab_den = 2.30
    filename = 'STA_mcnp_input_Concrete.txt'
elif absorber == 8:
    data_values.append('c Borated Polyethylene')
    data_values.append('m200 01001 -0.125355')
    data_values.append('     05010 -0.020000')
    data_values.append('     05011 -0.080000')
    data_values.append('     06000 -0.774645')
    ab_den = 1.00
    filename = 'STA_mcnp_input_BorPoly.txt'
elif absorber == 9:
    data_values.append('c Polyethylene')
    data_values.append('m200 06000 -0.240183')
    data_values.append('     01001 -0.143716')
    ab_den = 0.93
    filename = 'STA_mcnp_input_Poly.txt'
print('Your filename is:', filename)
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
# ax.add_patch(Drawing_colored_circle1)
# ax.add_patch(Drawing_colored_circle2)
# ax.add_patch(Drawing_colored_circle3)
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

# new_header = []
# new_cell_values = []
# new_surface_values = []
# new_data_values = []
# if contact == 'yes':
#     new_header = header
#     new_header.append('c CONTACTING CONTAINERS')    
#     new_cell_values.append('c Cells')
#     new_cell_values.append('5 100 -'+str(np.round(soln_den,4)) +' -15 12 -14 imp:n=1         $ container')
#     new_cell_values.append('6 LIKE 5 BUT TRCL=1')
#     new_cell_values.append('7 LIKE 5 BUT TRCL=2')
#     new_cell_values.
    
    
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
email = 'nbranam@vols.utk.edu'
header_qsub = []
header_qsub.append('#!/bin/bash')
header_qsub.append('#PBS -q fill')
header_qsub.append('#PBS -V')
header_qsub.append('#PBS -l nodes=1:ppn=8')
header_qsub.append('#PBS -m abe')
header_qsub.append('#PBS -N' + filename[0:-4])
header_qsub.append('#PBS -M ' + email)
header_qsub.append('')
header_qsub.append('cd $PBS_O_WORKDIR')
header_qsub.append('')
header_qsub.append('RTP=runtape_$(date "+%s%N")')
header_qsub.append('')
header_qsub.append('module load MCNP6/2.0')
header_qsub.append('')
header_qsub.append('mcnp6 TASKS 8 i=' + filename + ' o=' + filename[0:-4] + '_out.txt  runtpe=/tmp/$RTP')
header_qsub.append('rm /tmp/$RTP')

qsub_filename = filename[0:-4] + '.sh'
lines_for_qsub = [header_qsub]
with open(qsub_filename,'w') as file:
    for inte,line in enumerate(lines_for_qsub):
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
            
#%% Converts to Unix Format     

windows_line_ending = b'\r\n'
linux_line_ending = b'\n'
with open(qsub_filename, 'rb') as f:
    content = f.read()
    content = content.replace(windows_line_ending, linux_line_ending)

with open(qsub_filename, 'wb') as f:
    f.write(content)


