# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 15:10:41 2023

@author: Nicholas Branam
@contact: nbranam@vols.utk.edu
@ Version 2.0: This contact function allows for the addition
               of hands and walls around the container for
               criticality safety evaluation investigations.

"""
#%%

from math import pi
import numpy as np
import datetime
import matplotlib.pyplot as plt

#%%

# filename = 'STA_mcnp_input.txt'

def STA_jockey_contact(f,l,g,hd,flood,vr,d,w,h):
    print('         Welcome to the ')
    print(' ╭━━━━╮╱╱╱╱╭━━━━┳━━━┳━╮╱╭┳━╮╱╭╮')
    print(' ┃╭╮╭╮┃╱╱╱╱┃╭╮╭╮┃╭━━┫┃╰╮┃┃┃╰╮┃┃')
    print(' ╰╯┃┃┣┻┳╮╱╱╰╯┃┃╰┫╰━━┫╭╮╰╯┃╭╮╰╯┃')
    print(' ╱╱┃┃┃╭╋┫╭━━╮┃┃╱┃╭━━┫┃╰╮┃┃┃╰╮┃┃')
    print(' ╱╱┃┃┃┃┃┃╰━━╯┃┃╱┃╰━━┫┃╱┃┃┃┃╱┃┃┃')
    print(' ╱╱╰╯╰╯╰╯╱╱╱╱╰╯╱╰━━━┻╯╱╰━┻╯╱╰━╯')
    print('        File Generator')
    print()
    print('Generating File . . .')
    print()
    print('Filename: ',f)
    print('File description: ',d)
    
    # f - filename
    # l - liters of container
    # g - concentration in grams per liter
    # hd - h to d ratio
    # flood - ability to add water to room
    # vr - various runs, turns of .sh writer
    # d - description to write in file
    # w - walls (default thickness 0.5 cm)
    # h - hands (default thickness 2.5 cm)
    deck_title = 'STA MCNP Input Deck ' + str(d)
    # absorber: choose 1 for Boron-10, 2 for Aluminum, 3 for Air, 4 for Natural Boron,
    # 5 for Cadmium, 6 for Gadolinium, 7 for Oak Ridge Concrete,
    # 8 for Borated Polyethylene, 9 for Polyethylene
    # absorber = a
    volume_of_container = l*1000 # 1000 cm^3 per liter
    # Concentrations
    hy = 0.60689089
    n = 0.02394325
    o = 0.36916586
    p239 = 1.20317020e-02
    p240 = 5.30869186e-04
    p241 = 3.62446221e-05
    p242 = 1.99065079e-06
    soln_den = 1.24136
    
    pp = 2
    
    wall_thickness = 0.5 # cm
    hands_thickness = 2.5 # cm
    
    threshold = 0.1
    con_r = (((4*volume_of_container/pi)/(hd))**(1/3))/2
    con_h = con_r*2*hd

    
    # con_h = (hd)*((4*volume_of_container)/(pi*hd))**(1/3)
    # con_r = con_h/(hd*2)
    # con_h =  (4*volume_of_container/pi)**(1./3) # height of container
    # con_r = con_h/2 # radius of container
    print()
    if (w == 'yes') and (h == 'yes'):
        center_to_center = 2*con_r + 2*wall_thickness + 2*hands_thickness + threshold
        print('You have added walls and hands to the model.')
    elif (w == 'yes') or (h == 'yes'):
        if (w == 'yes'):
            center_to_center = 2*con_r + 2*wall_thickness + threshold
            print('You have added walls to the model.')
        else:
            center_to_center = 2*con_r + 2*hands_thickness + threshold
            print('You have added hands to the model.')
    else:
        center_to_center = 2*con_r + threshold
        print('This model has no walls or hands added.')
    
    print()
    cd = np.sqrt((center_to_center**2)/(2*(1-np.cos(120*pi/180))))
    
    
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
    
    plotting = 'no'
     
    
    
    # interaction check
    # threshold = 0.21
    # y_interaction_rail = -rail_edge*np.cos((pi/180)*30)+(rail_w/2)*np.sin((pi/180)*30)
    # if y_interaction_rail > -threshold:
    #     print('Error: Bounds of the rail overlap and MCNP will get very upset. :(')
        
    # y_interaction_plat = cd - plat_l/2
    # if (y_interaction_plat-threshold) < rail_edge:
    #     print('Error: Bounds of the platform overlap and MCNP will get very upset. :(')
    
    # k code
    num_part = 10000 # number of particles per generation
    num_skip = 10 #  number of generations skipped
    num_tot_gen = 50 # number of total generations
    
    print('K-code run selections: ')
    print('Number of particles per generation: ',num_part)
    print('Number of generations skipped: ',num_skip)
    print('Number of total generations: ',num_tot_gen)
    print()
    
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
    
    if (w == 'yes') and (h == 'yes'):
        surface_values = []
        surface_values.append('c Surface')
        surface_values.append('1 px -50    $ west edge of rail')
        surface_values.append('5 px 50   $ east edge of rail')
        surface_values.append('6 py -50    $ south edge of rail')
        surface_values.append('9 py 50     $ north edge of rail')
        surface_values.append('10 pz 0     $ bottom of rail')
        surface_values.append('11 pz ' + str(np.round(rail_h,2)) + '    $ top of rail, bottom of lower hands')
        surface_values.append('12 pz ' + str(np.round(rail_h+hands_thickness,2)) + ' $ top of lower hands, bottom of lower wall')
        surface_values.append('13 pz ' + str(np.round(rail_h+hands_thickness+wall_thickness,2)) + ' $ top of lower wall, bottom of solution')
        surface_values.append('14 pz ' + str(np.round(rail_h+hands_thickness+wall_thickness+con_h,2)) + ' $ top of solution, bottom of upper wall')
        surface_values.append('15 pz ' + str(np.round(rail_h+hands_thickness+2*wall_thickness+con_h,2)) + ' $ top of upper wall, bottom of upper hands')
        surface_values.append('16 pz ' + str(np.round(rail_h+2*hands_thickness+2*wall_thickness+con_h,2)) + ' $ top of upper hands')
        surface_values.append('17 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r,2))+' $ center of container in xy and container radius')
        surface_values.append('18 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r+wall_thickness,2))+' $ center of container wall in xy and wall radius')
        surface_values.append('19 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r+wall_thickness+hands_thickness,2))+' $ center of hands in xy and hands radius')
        surface_values.append('20 sph 0 0 0 150')
        surface_values.append('')
        
        
    elif w == 'yes' or h == 'yes':
        if w == 'yes':
            surface_values = []
            surface_values.append('c Surface')
            surface_values.append('1 px -50    $ west edge of rail')
            surface_values.append('5 px 50   $ east edge of rail')
            surface_values.append('6 py -50    $ south edge of rail')
            surface_values.append('9 py 50     $ north edge of rail')
            surface_values.append('10 pz 0     $ bottom of rail')
            surface_values.append('11 pz ' + str(np.round(rail_h,2)) + '    $ top of rail, bottom of lower wall')
            surface_values.append('12 pz ' + str(np.round(rail_h+wall_thickness,2)) + ' $ top of lower wall, bottom of solution')
            surface_values.append('13 pz ' + str(np.round(rail_h+wall_thickness+con_h,2)) + ' $ top of solution, bottom of upper wall')
            surface_values.append('14 pz ' + str(np.round(rail_h+2*wall_thickness+con_h,2)) + ' $ top of upper wall')
            surface_values.append('15 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r,2))+' $ center of container in xy and container radius')
            surface_values.append('16 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r+wall_thickness,2))+' $ center of container wall in xy and wall radius')
            surface_values.append('17 sph 0 0 0 150')
            surface_values.append('')
        else:
            surface_values = []
            surface_values.append('c Surface')
            surface_values.append('1 px -50    $ west edge of rail')
            surface_values.append('5 px 50   $ east edge of rail')
            surface_values.append('6 py -50    $ south edge of rail')
            surface_values.append('9 py 50     $ north edge of rail')
            surface_values.append('10 pz 0     $ bottom of rail')
            surface_values.append('11 pz ' + str(np.round(rail_h,2)) + '    $ top of rail, bottom of lower hands')
            surface_values.append('12 pz ' + str(np.round(rail_h+hands_thickness,2)) + ' $ top of lower hands, bottom of solution')
            surface_values.append('13 pz ' + str(np.round(rail_h+hands_thickness+con_h,2)) + ' $ top of solution, bottom of upper hands')
            surface_values.append('14 pz ' + str(np.round(rail_h+2*hands_thickness+con_h,2)) + ' $ top of upper hands')
            surface_values.append('15 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r,2))+' $ center of container in xy and container radius')
            surface_values.append('16 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r+hands_thickness,2))+' $ center of container hands in xy and hands radius')
            surface_values.append('17 sph 0 0 0 150')
            surface_values.append('')
        
    else:
        surface_values = []
        surface_values.append('c Surface')
        surface_values.append('1 px -50    $ west edge of rail')
        surface_values.append('5 px 50   $ east edge of rail')
        surface_values.append('6 py -50    $ south edge of rail')
        surface_values.append('9 py 50     $ north edge of rail')
        surface_values.append('10 pz 0     $ bottom of rail')
        surface_values.append('11 pz ' + str(np.round(rail_h,2)) + '    $ top of rail, bottom of platform')
        surface_values.append('12 pz 13    $ top of platform, bottom of lip and container')
        surface_values.append('14 pz ' + str(np.round(plat_h+con_h+rail_h,2)) + ' $ top of container and abs')
        surface_values.append('15 c/z -' + str(np.round(cd,2)) + ' 0 '+str(np.round(con_r,2))+' $ center of container in xy and container radius')
        surface_values.append('17 sph 0 0 0 100')
        surface_values.append('')
    
    
    #%% 
    
    hy = 0.60689089
    n = 0.02394325
    o = 0.36916586
    p239 = 1.20317020e-02
    p240 = 5.30869186e-04
    p241 = 3.62446221e-05
    p242 = 1.99065079e-06
    
    #%% USE MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    Pu_atomic_fractions = np.array([p239, p240, p241, p242])
    Pu_total_atomic_fractions = np.sum(Pu_atomic_fractions)
    Pu_ratios = Pu_atomic_fractions / Pu_total_atomic_fractions
    HNO3_atomic_fractions = np.array([hy, n, o])
    HNO3_total_atomic_fractions = np.sum(HNO3_atomic_fractions)
    HNO3_ratios = HNO3_atomic_fractions /  HNO3_total_atomic_fractions
    
    volume_solution = l # L
    concentration_plutonium = g # g/L
    density_plutonium = 19840 # g/L
    density_nitric_acid = 1510 #g/L
    
    volume_plutonium_solution = volume_solution*concentration_plutonium / density_plutonium # L
    volume_HNO3 = volume_solution - volume_plutonium_solution
    grams_HNO3_calc = volume_HNO3*density_nitric_acid # g of total solution
    grams_plutonium_calc = concentration_plutonium*volume_solution # g in X L of solution
    total_mass_solution = grams_HNO3_calc + grams_plutonium_calc
    weight_fraction_Pu_ratios = Pu_ratios * grams_plutonium_calc / total_mass_solution
    weight_fraction_HNO3_ratios = HNO3_ratios * grams_HNO3_calc / total_mass_solution
    MCNP_input_WF = np.append(weight_fraction_HNO3_ratios,weight_fraction_Pu_ratios)
    #print(np.sum(MCNP_input_WF))
    
    #%%
    
    # # listed in atomic fractions
    # # convert to WF, calculate the amount of Pu in the system to credible
    
    # # Convert to weight fractions here
    # # Adjust to credible concentrations
    # OG_AF = np.sum(np.array([h,n,o,p239,p240,p241,p242]))
    
    # denom = h*1 + o*16 + n*14 + p239*239 + p240*240 + p241*241 + p242*242
    # num = np.array([h*1, o*16, n*14, p239*239, p240*240, p241*241, p242*242])
    # weight_fractions = num/denom
    # total_plutonium_WF = np.sum([weight_fractions[3:7]])
    # density_plutonium = 19840 #g/L
    # grams_plutonium = total_plutonium_WF*density_plutonium*volume_of_container # grams
    # total_nitric_acid = np.sum([weight_fractions[0:3]])
    # density_nitric_acid = 1510 #g/L
    # grams_nitric_acid = total_nitric_acid*density_nitric_acid*volume_of_container # grams
    
    
    # Pu_atomic_fractions = np.array([p239, p240, p241, p242])
    # Pu_total_atomic_fractions = np.sum(Pu_atomic_fractions)
    # Pu_ratios = Pu_atomic_fractions / Pu_total_atomic_fractions
    
    # HNO3_atomic_fractions = np.array([h, n, o])
    # HNO3_total_atomic_fractions = np.sum(HNO3_atomic_fractions)
    # HNO3_ratios = HNO3_atomic_fractions /  HNO3_total_atomic_fractions
    
    
    # #%% USE MEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
    # # assume 1L of solution
    # grams_plutonium_calc = gpl*liters_in_container # g in 1 L of solution
    # weight_fraction_plutonium_calc = grams_plutonium_calc / density_plutonium
    # weight_fraction_Pu_ratios = Pu_ratios * weight_fraction_plutonium_calc
    
    # # 1L of HNO3 is 1L*1510g/L
    # grams_HNO3_calc = 1510*liters_in_container # g of total solution
    
    # weight_fraction_HNO3_calc = grams_HNO3_calc / density_nitric_acid
    # weight_fraction_HNO3_ratios = HNO3_ratios * weight_fraction_HNO3_calc
    
    # MCNP_input_WF = np.append(weight_fraction_HNO3_ratios,weight_fraction_Pu_ratios)
    
    
    #h = 'yes'
    
    #%% Data Card
    
    if (w == 'yes') and (h == 'yes'):
        # print('here')
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
        data_values.append('c Platform and Rail Material (steel)')
        data_values.append('m200 06000 -0.0050')
        data_values.append('     26000 -0.9950')
        mat_add = 300
        data_values.append('c Water')
        data_values.append('m' + str(int(mat_add)) + ' 01001 -0.111894')
        hands_mat_num = mat_add
        data_values.append('     08016 -0.888106')
        mat_add += 100
        data_values.append('c Wall of Container (Stainless Steel 304)')
        data_values.append('m'+str(mat_add)+' 06000 -0.000400')
        ss_mat_num = mat_add
        mat_add += 100
        data_values.append('     14000 -0.005000')
        data_values.append('     15031 -0.000230')
        data_values.append('     16000 -0.000150')
        data_values.append('     24000 -0.190000')
        data_values.append('     25055 -0.010000')
        data_values.append('     26000 -0.701730')
        data_values.append('     28000 -0.092500')
        data_values.append('mode n')
        data_values.append('kcode ' + str(int(num_part))+ ' 1.0 '+str(int(num_skip))+' '+str(int(num_tot_gen)))
        z = (np.round(con_h+2*wall_thickness+2*hands_thickness+rail_h,2)+np.round(rail_h,2))/2
        k_pos1 = ' -'+str(np.round(cd,pp)) + ' 0 ' + str(np.round(z,pp))
        k_pos2 = ' '+str(np.round((cd)*np.cos(pi/3),pp))+' '+str(np.round((cd)*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
        k_pos3 = ' '+str(np.round((cd)*np.cos(pi/3),pp))+' -'+str(np.round((cd)*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
        data_values.append('ksrc'+k_pos1+k_pos2+k_pos3)
    
    elif (w == 'yes') or (h == 'yes'):
        if w == 'yes':
            # print('5555')
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
            data_values.append('c Platform and Rail Material (steel)')
            data_values.append('m200 06000 -0.0050')
            data_values.append('     26000 -0.9950')
            mat_add = 300
            if flood == 'yes_steam' or flood == 'yes_flood':
                data_values.append('c Water')
                data_values.append('m' + str(int(mat_add)) + ' 01001 -0.111894')
                data_values.append('     08016 -0.888106')
                mat_add += 100
            data_values.append('c Wall of Container (Stainless Steel 304)')
            data_values.append('m'+str(mat_add)+' 06000 -0.000400')
            ss_mat_num = mat_add
            mat_add += 100
            data_values.append('     14000 -0.005000')
            data_values.append('     15031 -0.000230')
            data_values.append('     16000 -0.000150')
            data_values.append('     24000 -0.190000')
            data_values.append('     25055 -0.010000')
            data_values.append('     26000 -0.701730')
            data_values.append('     28000 -0.092500')
            data_values.append('mode n')
            data_values.append('kcode ' + str(int(num_part))+ ' 1.0 '+str(int(num_skip))+' '+str(int(num_tot_gen)))
            z = (np.round(con_h+2*wall_thickness+rail_h,2)+np.round(rail_h,2))/2
            k_pos1 = ' -'+str(np.round(cd,pp)) + ' 0 ' + str(np.round(z,pp))
            k_pos2 = ' '+str(np.round((cd)*np.cos(pi/3),pp))+' '+str(np.round((cd)*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
            k_pos3 = ' '+str(np.round((cd)*np.cos(pi/3),pp))+' -'+str(np.round((cd)*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
            data_values.append('ksrc'+k_pos1+k_pos2+k_pos3)
        else:
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
            data_values.append('c Platform and Rail Material (steel)')
            data_values.append('m200 06000 -0.0050')
            data_values.append('     26000 -0.9950')
            mat_add = 300
            hands_mat_num = mat_add
            data_values.append('c Water')
            data_values.append('m' + str(int(mat_add)) + ' 01001 -0.111894')
            data_values.append('     08016 -0.888106')
            mat_add += 100
            data_values.append('mode n')
            data_values.append('kcode ' + str(int(num_part))+ ' 1.0 '+str(int(num_skip))+' '+str(int(num_tot_gen)))
            z = (np.round(con_h+2*hands_thickness+rail_h,2)+np.round(rail_h,2))/2
            k_pos1 = ' -'+str(np.round(cd,pp)) + ' 0 ' + str(np.round(z,pp))
            k_pos2 = ' '+str(np.round((cd)*np.cos(pi/3),pp))+' '+str(np.round((cd)*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
            k_pos3 = ' '+str(np.round((cd)*np.cos(pi/3),pp))+' -'+str(np.round((cd)*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
            data_values.append('ksrc'+k_pos1+k_pos2+k_pos3)
    else:
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
        data_values.append('c Platform and Rail Material (steel)')
        data_values.append('m200 06000 -0.0050')
        data_values.append('     26000 -0.9950')
        if flood == 'yes_steam' or flood == 'yes_flood':
            data_values.append('c Water')
            data_values.append('m300 01001 -0.111894')
            data_values.append('     08016 -0.888106') 
        data_values.append('mode n')
        data_values.append('kcode ' + str(int(num_part))+ ' 1.0 '+str(int(num_skip))+' '+str(int(num_tot_gen)))
        z = (np.round(plat_h+con_h+rail_h,2)+np.round(plat_h+rail_h,2))/2
        k_pos1 = ' -'+str(np.round(cd,pp)) + ' 0 ' + str(np.round(z,pp))
        k_pos2 = ' '+str(np.round(cd*np.cos(pi/3),pp))+' '+str(np.round(cd*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
        k_pos3 = ' '+str(np.round(cd*np.cos(pi/3),pp))+' -'+str(np.round(cd*np.cos(pi/6),pp)) +' '+ str(np.round(z,pp))
        data_values.append('ksrc'+k_pos1+k_pos2+k_pos3)
    
    extra = 5.5
    
    if plotting == 'yes':
        Drawing_colored_circle1 = plt.Circle(( -(np.round(cd,pp)) , 0 ), con_r+extra,color='k' )
        Drawing_colored_circle2 = plt.Circle(( np.round(cd*np.cos(pi/3),pp) , np.round(cd*np.cos(pi/6),pp) ), con_r+extra,color='k' )
        Drawing_colored_circle3 = plt.Circle(( np.round(cd*np.cos(pi/3),pp) , -np.round(cd*np.cos(pi/6),pp) ), con_r+extra,color='k' )
        Drawing_colored_circle4 = plt.Circle(( -(np.round(cd,pp)) , 0 ), con_r )
        Drawing_colored_circle5 = plt.Circle(( np.round(cd*np.cos(pi/3),pp) , np.round(cd*np.cos(pi/6),pp) ), con_r )
        Drawing_colored_circle6 = plt.Circle(( np.round(cd*np.cos(pi/3),pp) , -np.round(cd*np.cos(pi/6),pp) ), con_r )
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
    
    
    if (w == 'yes') and (h == 'yes'):
        # print('here')
        cell_values = []
        cell_values.append('c Cells')
        cell_values.append('5 100 -'+str(np.round(soln_den,4)) +' -17 13 -14 imp:n=1         $ container')
        cell_values.append('6 LIKE 5 BUT TRCL=1')
        cell_values.append('7 LIKE 5 BUT TRCL=2')
        cell_values.append('10 '+str(ss_mat_num)+' -8.0 -18 17 13 -14 imp:n=1 $ walls')
        cell_values.append('11 LIKE 10 BUT TRCL=1')
        cell_values.append('12 LIKE 10 BUT TRCL=2')
        cell_values.append('15 '+str(ss_mat_num)+' -8.0 -18 12 -13 imp:n=1  $ wall bottom')
        cell_values.append('16 LIKE 15 BUT TRCL=1')
        cell_values.append('17 LIKE 15 BUT TRCL=2')
        cell_values.append('20 '+str(ss_mat_num)+' -8.0 -18 14 -15 imp:n=1  $ wall top')
        cell_values.append('21 LIKE 20 BUT TRCL=1')
        cell_values.append('22 LIKE 20 BUT TRCL=2')
        cell_values.append('30 '+ str(hands_mat_num) + ' -1.0 -19 18 12 -15 imp:n=1   $ hands')
        cell_values.append('31 LIKE 30 BUT TRCL=1')
        cell_values.append('32 LIKE 30 BUT TRCL=2')
        cell_values.append('35 '+ str(hands_mat_num) + ' -1.0 -19 11 -12 imp:n=1   $ hands bottom')
        cell_values.append('36 LIKE 35 BUT TRCL=1')
        cell_values.append('37 LIKE 35 BUT TRCL=2')
        cell_values.append('40 '+ str(hands_mat_num) + ' -1.0 -19 15 -16 imp:n=1    $ hands top')
        cell_values.append('41 LIKE 40 BUT TRCL=1')
        cell_values.append('42 LIKE 40 BUT TRCL=2')
        cell_values.append('50 200 -7.82 1 -5 6 -9 10 -11 imp:n=1  $ rail')
        if flood == 'yes_liquid':
            cell_values.append('100 300 -1.0 -20 #5 #6 #7 #10 #11 #12')
        elif flood== 'yes_steam':
            cell_values.append('100 300 -0.000756 -20 #5 #6 #7 #10 #11 #12')
        else:
            cell_values.append('100 0 -20 #5 #6 #7 #10 #11 #12')
        cell_values.append('     #15 #16 #17 #20 #21 #22')
        cell_values.append('     #30 #31 #32 #35 #36 #37')
        cell_values.append('     #40 #41 #42 #50	 imp:n=1')
        cell_values.append('999 0 20 imp:n=0')
        cell_values.append('')
        
    elif (w == 'yes') or (h == 'yes'):
        if w == 'yes':
            cell_values = []
            cell_values.append('c Cells')
            cell_values.append('5 100 -'+str(np.round(soln_den,4)) +' -15 12 -13 imp:n=1         $ container')
            cell_values.append('6 LIKE 5 BUT TRCL=1')
            cell_values.append('7 LIKE 5 BUT TRCL=2')
            cell_values.append('10 '+str(ss_mat_num)+' -8.0 -16 15 12 -13 imp:n=1 $ walls')
            cell_values.append('11 LIKE 10 BUT TRCL=1')
            cell_values.append('12 LIKE 10 BUT TRCL=2')
            cell_values.append('15 '+str(ss_mat_num)+' -8.0 -16 11 -12 imp:n=1  $ wall bottom')
            cell_values.append('16 LIKE 15 BUT TRCL=1')
            cell_values.append('17 LIKE 15 BUT TRCL=2')
            cell_values.append('20 '+str(ss_mat_num)+' -8.0 -16 13 -14 imp:n=1  $ wall top')
            cell_values.append('21 LIKE 20 BUT TRCL=1')
            cell_values.append('22 LIKE 20 BUT TRCL=2')
            cell_values.append('50 200 -7.82 1 -5 6 -9 10 -11 imp:n=1  $ rail')
            if flood == 'yes_liquid':
                cell_values.append('100 300 -1.0 -20 #5 #6 #7 #10 #11 #12')
            elif flood== 'yes_steam':
                cell_values.append('100 300 -0.000756 -20 #5 #6 #7 #10 #11 #12')
            else:
                cell_values.append('100 0 -17 #5 #6 #7 #10 #11 #12')
            cell_values.append('     #15 #16 #17 #20 #21 #22')
            cell_values.append('     #50	 imp:n=1')
            cell_values.append('999 0 17 imp:n=0')
            cell_values.append('')
        else:
            cell_values = []
            cell_values.append('c Cells')
            cell_values.append('5 100 -'+str(np.round(soln_den,4)) +' -15 12 -13 imp:n=1         $ container')
            cell_values.append('6 LIKE 5 BUT TRCL=1')
            cell_values.append('7 LIKE 5 BUT TRCL=2')
            cell_values.append('10 '+str(hands_mat_num)+' -1.0 -16 15 12 -13 imp:n=1 $ hands')
            cell_values.append('11 LIKE 10 BUT TRCL=1')
            cell_values.append('12 LIKE 10 BUT TRCL=2')
            cell_values.append('15 '+str(hands_mat_num)+' -1.0 -16 11 -12 imp:n=1  $ hands bottom')
            cell_values.append('16 LIKE 15 BUT TRCL=1')
            cell_values.append('17 LIKE 15 BUT TRCL=2')
            cell_values.append('20 '+str(hands_mat_num)+' -1.0 -16 13 -14 imp:n=1  $ hands top')
            cell_values.append('21 LIKE 20 BUT TRCL=1')
            cell_values.append('22 LIKE 20 BUT TRCL=2')
            cell_values.append('50 200 -7.82 1 -5 6 -9 10 -11 imp:n=1  $ rail')
            if flood == 'yes_liquid':
                cell_values.append('100 300 -1.0 -20 #5 #6 #7 #10 #11 #12')
            elif flood== 'yes_steam':
                cell_values.append('100 300 -0.000756 -20 #5 #6 #7 #10 #11 #12')
            else:
                cell_values.append('100 0 -17 #5 #6 #7 #10 #11 #12')
            cell_values.append('     #15 #16 #17 #20 #21 #22')
            cell_values.append('     #50	 imp:n=1')
            cell_values.append('999 0 17 imp:n=0')
            cell_values.append('')
    else:
        cell_values = []
        cell_values.append('c Cells')
        cell_values.append('5 100 -'+str(np.round(soln_den,4)) +' -15 12 -14 imp:n=1         $ container')
        cell_values.append('6 LIKE 5 BUT TRCL=1')
        cell_values.append('7 LIKE 5 BUT TRCL=2')    
        cell_values.append('20 200 -7.82 1 -5 6 -9 10 -12 imp:n=1  $ rail')
        if flood == 'yes_liquid':
            cell_values.append('100 300 -1.0 -17 #5 #6 #7 #10 #11 #12 imp:n=1')
        elif flood== 'yes_steam':
            cell_values.append('100 300 -0.000756 -17 #5 #6 #7 #10 #11 #12 imp:n=1')
        else:
            cell_values.append('100 0 -17 #5 #6 #7 #10 #11 #12 imp:n=1')
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
        
        
    with open(f,'w',encoding="utf-8") as file:
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
    print('File has been generated. Thank you.')
    
    #%%
    if vr == 'no':
        email = 'nbranam@vols.utk.edu'
        header_qsub = []
        header_qsub.append('#!/bin/bash')
        header_qsub.append('#PBS -q fill')
        header_qsub.append('#PBS -V')
        header_qsub.append('#PBS -l nodes=1:ppn=8')
        header_qsub.append('#PBS -m abe')
        header_qsub.append('#PBS -N' + f[0:-4])
        header_qsub.append('#PBS -M ' + email)
        header_qsub.append('')
        header_qsub.append('cd $PBS_O_WORKDIR')
        header_qsub.append('')
        header_qsub.append('RTP=runtape_$(date "+%s%N")')
        header_qsub.append('')
        header_qsub.append('module load MCNP6/2.0')
        header_qsub.append('')
        header_qsub.append('mcnp6 TASKS 8 i=' + f + ' o=' + f[0:-4] + '_out.txt  runtpe=/tmp/$RTP')
        header_qsub.append('rm /tmp/$RTP')
        
        qsub_filename = f[0:-4] + '.sh'
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