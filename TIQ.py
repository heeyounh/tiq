# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 22:50:41 2020

@author: Heeyoun
"""
import numpy as np
import sys


def printProgress (iteration, total, prefix = '', suffix = '', decimals = 1, barLength = 100):
    formatStr = "{0:." + str(decimals) + "f}"
    percent = formatStr.format(100 * (iteration / float(total)))
    filledLength = int(round(barLength * iteration / float(total)))
    bar = '#' * filledLength + '-' * (barLength - filledLength)
    sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percent, '%', suffix)),
    if iteration == total:
        sys.stdout.write('\n')
    sys.stdout.flush()
    
rt_tol = 5.0
ppm = 20.0

def ms1_parsing(file_):
    
    #return x + x, a
    # for file_ in file_list:
    #     if ms1_file in file_:
    scan_dic = {}
    ms1_file_ = open(file_ + '.ms1')
    print(file_ + ' is opened.')
    ms1_file_lines = ms1_file_.readlines()
    
    

    scan = ''
    mz = []
    intensity = []
    charge = []
    RetTime = ''
    noise = []
    
    for line in ms1_file_lines:
        if line[0] == 'H':
            pass
        elif line[0] == 'S':
            if mz == []:
                pass
            else:
                #print noise
                scan_dic [RetTime] = [scan, mz, intensity, charge, np.average(noise)]
            mz = []
            intensity = []
            charge = []
            noise = []
            
            lines = line[:-1].split('\t')
            scan = lines[1]
            
        elif line[0] == 'I':
            lines = line[:-1].split('\t')
            if lines[1] == 'RetTime':
                RetTime = float(lines[2])
        else:
            lines = line[:-1].split(' ')
            
            
            if int(lines[2]) >= 2:
                mz.append(float(lines[0]))
                intensity.append(float(lines[1]))
                charge.append(int(lines[2]))
            else:
                noise.append(float(lines[1]))

    #dump_file = open(file_[:-4] + 'ms1.dmp', 'wb')
    #pickle.dump(scan_dic, dump_file)                        
#    if scan_dic == {}:
#        pass
#    else:
#                        
#        del scan_dic['']

    return scan_dic


def quantification(data):   ##[pep_Z_, ms1_dic[scanfile_list[i]], RT_dic]
    pep_Z_list = data[0]
    ms1_dic = data[1]
    RT_dic = data[2]
    
    TIQdata = []
    k = 0
    for pep_Z_ in pep_Z_list:
        k= k + 1
        printProgress(k, len(pep_Z_list), 'Progress:', 'complete', 1, 50)    
    

        TIQ = 0.0
        pep_Z = pep_Z_.split(':')
        peptide = pep_Z[1]
        Z =  float(pep_Z[2])
        mz = float(pep_Z[3])
        MH = mz * Z - 1.007825 * (Z - 1)
        M_data = []
        M_1_data = []
        M_2_data = []
        rt_data = []    
        rt_list = list(ms1_dic.keys())
        RT_list = RT_dic[pep_Z_]
        for rt_ in rt_list:
            if rt_ > min(RT_list) - rt_tol and rt_ < max(RT_list) + rt_tol:
                frac_data_list = ms1_dic [rt_]
                frac_mz_list = frac_data_list[1]
                frac_int_list = frac_data_list[2]
                frac_Z_list = frac_data_list[3]        
    
                for j in range(0, len(frac_mz_list)):
                    frac_MH = frac_mz_list[j] * frac_Z_list[j] - 1.007825 * (frac_Z_list[j] - 1)
                    
                    if abs(frac_MH - MH)/MH * 1000000 < ppm  and Z == frac_Z_list[j]:
                        if j == 0:
                            try:
                                frac_MH_1 = frac_mz_list[j+1] * frac_Z_list[j+1] - 1.007825 * (frac_Z_list[j+1] - 1)
                                frac_MH_2 = frac_mz_list[j+2] * frac_Z_list[j+2] - 1.007825 * (frac_Z_list[j+2] - 1)
                                
                                if abs(frac_MH + 1.0 - frac_MH_1)/frac_MH_1 * 1000000 < ppm and abs(frac_MH + 2.0 - frac_MH_2)/frac_MH_2 * 1000000 < ppm:
                                    M_data.append(frac_int_list[j])
                                    M_1_data.append(frac_int_list[j+1])
                                    M_2_data.append(frac_int_list[j+2])
                                    rt_data.append(rt_)
                                    #print exp_, rep_, pep_z_list[i], frac_int_list[j], rt_, [frac_mz_list[j], frac_mz_list[j+1], frac_mz_list[j+2]]
                            except:
                                pass
                                #print('no three isotope', exp_, rep_, pep_z_list[i], frac_int_list[j], rt_, frac_MH, file = log_file)
    
    
    
    
                            
                        else:
                            frac_MH_before = frac_mz_list[j-1] * frac_Z_list[j-1] - 1.007825 * (frac_Z_list[j-1] - 1)
                            
                            if abs(frac_MH - 1.0 - frac_MH_before)/frac_MH_before * 1000000 < ppm and frac_int_list[j] / 2 < frac_int_list[j-1]:
                                #print('not mono isotope', exp_, rep_, pep_z_list[i], frac_int_list[j], rt_, frac_MH, frac_MH_before, file = log_file)
                                pass
                            else:
                                try:
                                    frac_MH_1 = frac_mz_list[j+1] * frac_Z_list[j+1] - 1.007825 * (frac_Z_list[j+1] - 1)
                                    frac_MH_2 = frac_mz_list[j+2] * frac_Z_list[j+2] - 1.007825 * (frac_Z_list[j+2] - 1)
                                    
                                    if abs(frac_MH + 1.0 - frac_MH_1)/frac_MH_1 * 1000000 < ppm and abs(frac_MH + 2.0 - frac_MH_2)/frac_MH_2 * 1000000 < ppm:
                                    
                                        M_data.append(frac_int_list[j])
                                        M_1_data.append(frac_int_list[j+1])
                                        M_2_data.append(frac_int_list[j+2])
                                        rt_data.append(rt_)
                                        #print exp_, rep_, pep_z_list[i], frac_int_list[j], rt_, [frac_mz_list[j], frac_mz_list[j+1], frac_mz_list[j+2]]
                                except:
                                    pass
                                    #print('no three isotope', exp_, rep_, pep_z_list[i], frac_int_list[j], rt_, frac_MH, file = log_file)
        if M_data == []:
            # TIQ_data.append(np.nan)
            # SumQ.append(np.nan)
            # TSumQ.append(np.nan)
            # TIQ_3data.append(np.nan)
            # TIQ2_data.append(np.nan)
            # TIQ4_data.append(np.nan)
            # TIQ5_data.append(np.nan)
            TIQdata.append([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
            pass
        else:
            index = M_data.index(max(M_data))
            TIQ = max(M_data) + M_1_data[index] + M_2_data[index]
            TIQ_RT = rt_data [index]
            #TIQ_data.append(TIQ)
            #TIQ2_data.append(TIQ * delta_RT)
            TIQ4_data = TIQ * pow(np.log10(TIQ),2) * 3.14 * 2.0 / 3.0
            delta_RT = max(rt_data) - min(rt_data)
            TIQ5_data = TIQ * pow(delta_RT, 2) * 3.14 * 2.0 /3.0
            TIQ2_data = TIQ * delta_RT
            SumQ = np.sum(M_data)
            TSumQ = np.sum(M_data) + np.sum(M_1_data) + np.sum(M_2_data)
            try:
                TIQ_2data = M_data[index - 1] + M_1_data[index - 1] + M_2_data[index - 1]
            except:
                TIQ_2data = 0.0
            try:
                TIQ_2data2 = M_data[index + 1] + M_1_data[index + 1] + M_2_data[index + 1]
            except:
                TIQ_2data2 = 0.0
            
            TIQ_3data = TIQ + TIQ_2data + TIQ_2data2
            TIQdata.append([TIQ, SumQ,TSumQ, TIQ_3data, TIQ2_data, TIQ4_data, TIQ5_data, TIQ_RT])
    return TIQdata



