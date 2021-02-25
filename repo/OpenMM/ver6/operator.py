import os 
import numpy as np 
import yaml
import time as timer
import math


'''
There should be two version of calculations 
1. Slab structure version 
2. Ring structure version 
'''

###############################################
##               Command Center              ##
## Choose your model of simulation structure ##
##  1. Slab                                  ##
##  2. Square (radial distribution)          ##
###############################################
## Structure 
Structure = 'Slab' 
#Structure = 'Square' # or 'Square'
## Configuration based on structure specification 

############################################
##         SLAB STRUCTURE CASE            ##
############################################

if Structure == 'Slab':
    shape_slab = 'slab'
    unitlength = 50
    density_list1 = [0.004]
    density_list2 = [0.004,0.008,0.012]
    cellPerdead = 10
    DiffState = 'error'  # error (for error function), steady, or linear 

############################################
##        Square STRUCTURE CASE           ##
############################################
elif Structure == 'Square':
    shape_sq = 'square'
    NoParticleZone_sq = 50 # if Slab, then None
    R_of_inner_layer_sq = 150 # This can be used to assign the radial section.
    R_of_outer_layer_sq = 300
    density_list1 = [0.0001] # maximum 0.001 with current dimension in square/round geometry
    density_list2 = [0.0001]#,50,100]
    cellPerdead = 0
    #DiffState = 'error'
    DiffState = 'steady'  # error (for error function), steady, or linear 


## Other common physical properties of simulation 
### This may not change unless necessary
time = 1000000
DispScale = 35
cellconc = 0.005
oriConc = 1.0
DiffRate = 0.01
kd = oriConc*0.8

## A Number of Simulation Related
repeat = 1
restingratio1 = 1
restingratio2 = [0.1]#,0.5,0.9]
#density_list1 = [0.0005] # maximum 0.001 with current dimension in square/round geometry
#density_list2 = [0.00025]#,50,100]
no_of_series_sim = len(density_list1)*len(density_list2)


## Timer ON
start_whole = timer.time()
print('|-----SIMULATION INITIATED-----|')
count = 1

## Simulation Setup
# 1. Slab structure case
if Structure == 'Slab':
    shape = shape_slab
    highlim = unitlength
    Area = 5*highlim**2
    
    
    for i in np.arange(repeat):
        if i == 0:
            start_single = timer.time()
    
        for ratio2 in restingratio2:
            for nocell1 in density_list1:
                for nocell2 in density_list2:
                    
                    cell_count1 = round(nocell1*Area)
                    cell_count2 = round(nocell2*Area)
                    print('Number of Cells in the 1st Compartment = ',cell_count1)
                    print('Number of Cells in the 2nd Compartment = ',cell_count2)
        
                    inputname = 'test'+str(count)
                    outputname = 'test'+str(count)
        
                    count += 1
                    os.system('python3 runner.py -t '+str(time)+
                              ' -highlimBoxlen '+str(highlim)+
                              ' -numOfCells1 '+str(cell_count1)+
                              ' -numOfCells2 '+str(cell_count2)+
                              ' -restingRatio1 '+str(restingratio1)+
                              ' -restingRatio2 '+str(ratio2)+
                              ' -DispScale '+str(DispScale)+
                              ' -name '+outputname+
                              ' -input '+inputname+
                              ' -cellConc '+str(cellconc)+
                              ' -oriConc '+str(oriConc)+
                              ' -shape '+shape+
                              ' -DiffState '+DiffState+
                              ' -Diff '+str(DiffRate)+
                              ' -kd '+str(kd)+
                              ' -cellPerdead '+str(cellPerdead))

        if i == 0:
            currentg_single = timer.time()
            print('|-----Single Run Job Statistics-----|')
            print('|------ %s secs ------|' % (currentg_single-start_single))
            print('|------ %s mins ------|' % ((currentg_single-start_single)/60))
            print('|------ %s hrs ------|' % ((currentg_single-start_single)/3600))
            print('|--------------------------------|')
            print('|-----Estimated Job Duration-----|')
            estimated_sec = (currentg_single-start_single+1*repeat)*repeat
            estimated_min = estimated_sec/60
            estimated_hrs = estimated_min/60
            print('|------ %s secs ------|' % (estimated_sec))
            print('|------ %s mins ------|' % (estimated_min))
            print('|------ %s hrs ------|' % (estimated_hrs))
        
        else : 
            current_single = timer.time()
            current_sec = current_single - start_single
            current_min = current_sec/60
            current_hrs = current_min/60
            print('|-----Simulation has been excuted since... -----|')
            print('|------'+str(current_sec)+' secs ago and possibly '+str(estimated_sec-current_sec)+' secs remaining ------|')
            print('|------'+str(current_min)+' mins ago and possibly '+str(estimated_min-current_min)+' mins remaining ------|')
            print('|------'+str(current_hrs)+' hrs ago and possibly '+str(estimated_hrs-current_hrs)+' hrs remaining ------|')

            
# 2. Square structure case 
elif Structure == 'Square':
    shape = shape_sq
    NoParticleZone = NoParticleZone_sq # if Slab, then None
    R_of_inner_layer = R_of_inner_layer_sq # This can be used to assign the radial section.
    R_of_outer_layer = R_of_outer_layer_sq
    Area_inner = math.pi*(R_of_inner_layer**2 - NoParticleZone**2)
    Area_outer = math.pi*(R_of_outer_layer**2 - R_of_inner_layer**2)
    
    for i in np.arange(repeat):
        if i == 0:
            start_single = timer.time()
    
        for ratio2 in restingratio2:
            for nocell1 in density_list1:
                for nocell2 in density_list2:
        
                    cell_count1 = round(Area_inner*nocell1)
                    cell_count2 = round(Area_outer*nocell2)
                    print('Number of Cells in the 1st Compartment = ',cell_count1)
                    print('Number of Cells in the 2nd Compartment = ',cell_count2)
                
                    inputname = 'test'+str(count)
                    outputname = 'test'+str(count)
        
                    count += 1
                    os.system('python3 runner.py -t '+str(time)+
                              ' -highlimBoxlen '+str(R_of_outer_layer)+
                              ' -numOfCells1 '+str(cell_count1)+
                              ' -numOfCells2 '+str(cell_count2)+
                              ' -restingRatio1 '+str(restingratio1)+
                              ' -restingRatio2 '+str(ratio2)+
                              ' -CentoR '+str(R_of_inner_layer)+
                              ' -DispScale '+str(DispScale)+
                              ' -name '+outputname+
                              ' -input '+inputname+
                              ' -cellConc '+str(cellconc)+
                              ' -oriConc '+str(oriConc)+
                              ' -shape '+shape+
                              ' -DiffState '+DiffState+
                              ' -Diff '+str(DiffRate)+
                              ' -kd '+str(kd)+
                              ' -NoParticleZone '+str(NoParticleZone)+
                              ' -cellPerdead '+str(cellPerdead))

        if i == 0:
            currentg_single = timer.time()
            print('|-----Single Run Job Statistics-----|')
            print('|------ %s secs ------|' % (currentg_single-start_single))
            print('|------ %s mins ------|' % ((currentg_single-start_single)/60))
            print('|------ %s hrs ------|' % ((currentg_single-start_single)/3600))
            print('|--------------------------------|')
            print('|-----Estimated Job Duration-----|')
            estimated_sec = (currentg_single-start_single+1*repeat)*repeat
            estimated_min = estimated_sec/60
            estimated_hrs = estimated_min/60
            print('|------ %s secs ------|' % (estimated_sec))
            print('|------ %s mins ------|' % (estimated_min))
            print('|------ %s hrs ------|' % (estimated_hrs))
        
        else : 
            current_single = timer.time()
            current_sec = current_single - start_single
            current_min = current_sec/60
            current_hrs = current_min/60
            print('|-----Simulation has been excuted since... -----|')
            print('|------'+str(current_sec)+' secs ago and possibly '+str(estimated_sec-current_sec)+' secs remaining ------|')
            print('|------'+str(current_min)+' mins ago and possibly '+str(estimated_min-current_min)+' mins remaining ------|')
            print('|------'+str(current_hrs)+' hrs ago and possibly '+str(estimated_hrs-current_hrs)+' hrs remaining ------|')
            
print('|-----SIMULATION COMPLETED-----|')       
print('|------ %s secs ------|' % (timer.time()-start_whole))
print('|------ %s mins ------|' % ((timer.time()-start_whole)/60))
print('|------ %s hrs ------|' % ((timer.time()-start_whole)/3600))
