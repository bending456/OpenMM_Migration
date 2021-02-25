import numpy as np 
import yaml


def stats(Repeat, no_of_series_sim, prefix):
    with open(r'test'+str(prefix)+'_raw.yaml') as file:
        raw_list = yaml.load(file, Loader=yaml.FullLoader)
        
    no_of_r = int(raw_list['n_r'])
    no_of_a = int(raw_list['n_a'])
    
    rest = np.zeros([3,no_of_r*Repeat])
    act = np.zeros([3,no_of_a*Repeat])
    
    num1 = 0
    num2 = 0
    
    for j in np.arange(Repeat):
        with open(r'test'+str(prefix+no_of_series_sim*j)+'_raw.yaml') as file:
            raw_list = yaml.load(file, Loader=yaml.FullLoader)

        no_of_r = int(raw_list['n_r'])
        no_of_a = int(raw_list['n_a'])


        # resting
        for n in np.arange(no_of_r):
            rx = float(raw_list['x_r'][n])
            ry = float(raw_list['y_r'][n])
            rr = float(raw_list['r_r'][n])

            rest[0][num1] = rx
            rest[1][num1] = ry
            rest[2][num1] = rr
            num1 += 1
       

        # activated
        for m in np.arange(no_of_a):
            ax = float(raw_list['x_a'][m])
            ay = float(raw_list['y_a'][m])
            ar = float(raw_list['r_a'][m])

            act[0][num2] = ax
            act[1][num2] = ay
            act[2][num2] = ar
            num2 += 1

    r_x_avg = np.mean(rest[0])
    r_y_avg = np.mean(rest[1])
    r_r_avg = np.mean(rest[2])
    r_x_std = np.std(rest[0])
    r_y_std = np.std(rest[1])
    r_r_std = np.std(rest[2])
    r_x_sem = r_x_std/(num1**0.5)
    r_y_sem = r_y_std/(num1**0.5)
    r_r_sem = r_r_std/(num1**0.5)
    r_x_rms = np.sqrt(np.sum(rest[0]**2)/num1)
    r_y_rms = np.sqrt(np.sum(rest[1]**2)/num1)
    r_r_rms = np.sqrt(np.sum(rest[2]**2)/num1)

    a_x_avg = np.mean(act[0])
    a_y_avg = np.mean(act[1])
    a_r_avg = np.mean(act[2])
    a_x_std = np.std(act[0])
    a_y_std = np.std(act[1])
    a_r_std = np.std(act[2])
    a_x_sem = a_x_std/(num2**0.5)
    a_y_sem = a_y_std/(num2**0.5)
    a_r_sem = a_r_std/(num2**0.5)
    a_x_rms = np.sqrt(np.sum(act[0]**2)/num2)
    a_y_rms = np.sqrt(np.sum(act[1]**2)/num2)
    a_r_rms = np.sqrt(np.sum(act[2]**2)/num2)
    
    total = np.hstack((rest,act))
    totalnum = num1 + num2
    t_x_avg = np.mean(total[0])
    t_y_avg = np.mean(total[1])
    t_r_avg = np.mean(total[2])
    t_x_std = np.std(total[0])
    t_y_std = np.std(total[1])
    t_r_std = np.std(total[2])
    t_x_sem = t_x_std/(totalnum**0.5)
    t_y_sem = t_y_std/(totalnum**0.5)
    t_r_sem = t_r_std/(totalnum**0.5)
    t_x_rms = np.sqrt(np.sum(total[0]**2)/totalnum)
    t_y_rms = np.sqrt(np.sum(total[1]**2)/totalnum)
    t_r_rms = np.sqrt(np.sum(total[2]**2)/totalnum)
    
    

    calculated_statistics = {'no_of_r_total':num1,
                             'no_of_r_per':no_of_r,
                             'no_of_a_total':num2,
                             'no_of_a_per':no_of_a,
                             'r_x_a':r_x_avg,
                             'r_y_a':r_y_avg,
                             'r_r_a':r_r_avg,
                             'r_x_s':r_x_std,
                             'r_y_s':r_y_std,
                             'r_r_s':r_r_std,
                             'r_x_e':r_x_sem,
                             'r_y_e':r_y_sem,
                             'r_r_e':r_r_sem,
                             'r_x_r':r_x_rms,
                             'r_y_r':r_y_rms,
                             'r_r_r':r_r_rms,
                             'a_x_a':a_x_avg,
                             'a_y_a':a_y_avg,
                             'a_r_a':a_r_avg,
                             'a_x_s':a_x_std,
                             'a_y_s':a_y_std,
                             'a_r_s':a_r_std,
                             'a_x_e':a_x_sem,
                             'a_y_e':a_y_sem,
                             'a_r_e':a_r_sem,
                             'a_x_r':a_x_rms,
                             'a_y_r':a_y_rms,
                             'a_r_r':a_r_rms,
                             't_x_a':t_x_avg,
                             't_x_s':t_x_std,
                             't_x_e':t_x_sem,
                             't_x_r':t_x_rms,
                             't_y_a':t_y_avg,
                             't_y_s':t_y_std,
                             't_y_e':t_y_sem,
                             't_y_r':t_y_rms,
                             't_r_a':t_r_avg,
                             't_r_s':t_r_std,
                             't_r_e':t_r_sem,
                             't_r_r':t_r_rms,}

    return calculated_statistics


def distribution(Repeat, no_of_series_sim, prefix):
    
    with open(r'test'+str(prefix)+'_raw.yaml') as file:
        raw_list = yaml.load(file, Loader=yaml.FullLoader)
        
    no_of_r = int(raw_list['n_r'])
    no_of_a = int(raw_list['n_a'])
    
    rest = np.zeros([3,no_of_r*Repeat])
    act = np.zeros([3,no_of_a*Repeat])
    
    num1 = 0
    num2 = 0
    
    for j in np.arange(Repeat):
        with open(r'test'+str(prefix+no_of_series_sim*j)+'_raw.yaml') as file:
            raw_list = yaml.load(file, Loader=yaml.FullLoader)

        no_of_r = int(raw_list['n_r'])
        no_of_a = int(raw_list['n_a'])


        # resting
        for n in np.arange(no_of_r):
            rx = float(raw_list['x_r'][n])
            ry = float(raw_list['y_r'][n])
            rr = float(raw_list['r_r'][n])

            rest[0][num1] = rx
            rest[1][num1] = ry
            rest[2][num1] = rr
            num1 += 1
       

        # activated
        for m in np.arange(no_of_a):
            ax = float(raw_list['x_a'][m])
            ay = float(raw_list['y_a'][m])
            ar = float(raw_list['r_a'][m])

            act[0][num2] = ax
            act[1][num2] = ay
            act[2][num2] = ar
            num2 += 1


    return rest, act