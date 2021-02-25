import os 
import numpy as np 
import yaml
import time as timer

time = 10000
DispScale = 7500 #2500
highlim = 150
the_len_of_1st_box = 0 #highlim*3
repeat = 1
cellconc = 0.0
oriConc = 5.0
restingratio1 = 0.5
restingratio2 = [0.1]#,0.5,0.9]
shape = 'slab'       # slab - work with SourceOfOrigin = None or square - work with SourceOfOrigin = Center 
DiffState = 'error'  # error (for error function), steady, or linear 
DiffRate = 5


density_list1 = [50]
density_list2 = [50]#,50,100]

no_of_series_sim = len(density_list1)*len(density_list2)

start_whole = timer.time()

print('|-----SIMULATION INITIATED-----|')


count = 1
for i in np.arange(repeat):
    if i == 0:
        start_single = timer.time()
    
    for ratio2 in restingratio2:
        for nocell1 in density_list1:
            for nocell2 in density_list2:
        
                inputname = 'test'+str(count)
                outputname = 'test'+str(count)
        
                count += 1
                os.system('python3 runner.py -t '+str(time)+
                          ' -highlimBoxlen '+str(highlim)+
                          ' -numOfCells1 '+str(nocell1)+
                          ' -numOfCells2 '+str(nocell2)+
                          ' -restingRatio1 '+str(restingratio1)+
                          ' -restingRatio2 '+str(ratio2)+
                          ' -CentoR '+str(the_len_of_1st_box)+
                          ' -DispScale '+str(DispScale)+
                          ' -name '+outputname+
                          ' -input '+inputname+
                          ' -cellConc '+str(cellconc)+
                          ' -oriConc '+str(oriConc)+
                          ' -shape '+shape+
                          ' -DiffState '+DiffState+
                          ' -Diff '+str(DiffRate))

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
