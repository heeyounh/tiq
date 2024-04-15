# -*- coding: utf-8 -*-
"""
Created on Sat Sep 12 22:09:00 2020

@author: Heeyoun
"""



import multiprocessing
import time
import TIQ
import pandas as pd
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

def protein_quant(TIQ_data_df, name):
    data_df_index = TIQ_data_df.index
    data = []
    for index in data_df_index:
        data_ = []
        series = TIQ_data_df.loc[index].tolist()
        index_ = index.split(':')
        data_ = index_ + series
        data.append(data_)
                    
    TIQ_data2_df = pd.DataFrame(data, columns = ['protein', 'peptide', 'charge', 'm/z'] + samples)
    
    protein_list = TIQ_data2_df['protein'].tolist()
    
    protein_list = list(set(protein_list))
    
    sample_df = pd.DataFrame()
    max_sample_df = pd.DataFrame()
    avg_sample_df = pd.DataFrame()
    for protein_ in protein_list:
        protein_df = TIQ_data2_df[TIQ_data2_df['protein'] == protein_]
        if len(protein_df) > 1:
            sum_df = protein_df[samples].sum(axis = 0)
            max_df = protein_df[samples].max(axis = 0)
            avg_df = protein_df[samples].mean(axis = 0)
            sum_df.name = protein_
            max_df.name = protein_
            avg_df.name = protein_
            sample_df = pd.concat([sample_df, sum_df], axis = 1)
            max_sample_df = pd.concat([max_sample_df, max_df], axis = 1)
            avg_sample_df = pd.concat([avg_sample_df, avg_df], axis = 1)
        
    protein_TIQ_df = np.log2(sample_df.T)
    protein_TIQ_df.to_excel(name + '_protein_MP.xlsx')
    
    protein_max_TIQ_df = np.log2(max_sample_df.T)
    protein_max_TIQ_df.to_excel('Max_' + name +'_protein_MP.xlsx')
    
    protein_avg_TIQ_df = np.log2(avg_sample_df.T)
    protein_avg_TIQ_df.to_excel('Avg_' + name +'_protein_MP.xlsx')
    
# def count(name):
#     name2 = ''
#     for i in range(1, 10):
#         #a = a + i
#         #print(name, "  :  ", i)
#         name2 = name2 + name
#     return name2

# num_list = ['p1', 'p2', 'p3', 'p4']

# a = 0

process = 12

if __name__ == "__main__":
    start_time = time.time()
    scanfile_df = pd.read_excel('scanfile_exp.xlsx',  index_col = 0)
    
    print('Multiprocessing = ' + str(process))
    
    scanfile_list = scanfile_df.index.tolist()

    
    #testfile = open('testfile.txt', 'w')
    pool = multiprocessing.Pool(processes = process)
    
    output = pool.map(TIQ.ms1_parsing, scanfile_list)
    print("--- %s seconds --- ms1 parsing" % (time.time() - start_time))
    #print (len(output))

    ms1_dic = {}
    for i in range(0, len(scanfile_list)):
        ms1_dic [scanfile_list[i]] = output[i]
        
    print('Dataframe file is opened.')
    start_time = time.time()       
    data_df = pd.read_excel('Peptide_Dataframe.xlsx')


        
    peplist = data_df ['pep'].tolist()
    
    peplist = list(set(peplist))
    
    peplist.sort()
    
    pep_z_list = []
    for peptide in peplist:
        pep_df = data_df[data_df ['pep'] == peptide] 
        charge_list = pep_df ['charge'].tolist()
        charge_list = list(set(charge_list))
        
        
        for charge in charge_list:
            data_ = []
            charge_df = pep_df[pep_df['charge'] == charge]
            cal_MZ_list = charge_df['cal_MZ'].tolist()
            protein_list = charge_df['protein'].tolist()
            
            
            pep_z_list.append(str(protein_list[0]) + ':' + str(peptide) + ':' + str(charge) + ':' + str(cal_MZ_list[0]))
            
            
    print("--- %s seconds --- peptide dataframe load" % (time.time() - start_time))
    start_time = time.time()
    TIQ_DATA = []
    SUMQ_DATA = []
    TSUMQ_DATA = []
    TIQ_3_DATA = []
    TIQ2_DATA = []
    TIQ4_DATA = []
    TIQ5_DATA = []
    RT_DATA = []
    k = 0
    pep_z_RT_dic = {}
    print (str(len(pep_z_list)) + ' peptides\t' + str(len(scanfile_list)) + ' files')
    for pep_Z_ in pep_z_list:
        

        
        
        #pep_Z_ = data[0]
        #ms1_dic = data[1]
        #TIQ = 0.0
        pep_Z = pep_Z_.split(':')
        peptide = pep_Z[1]
        Z =  float(pep_Z[2])
        mz = float(pep_Z[3])
        MH = mz * Z - 1.007825 * (Z - 1)
        
        process_data = []
        pep_frac_data_df = data_df[data_df['pep'] == peptide]
        z_pep_frac_data_df = pep_frac_data_df[pep_frac_data_df['charge'] == Z]
        RT_list = z_pep_frac_data_df['rt']


        if len(RT_list) == 0 :
            pass
        else:
            if max(RT_list) - min(RT_list) > 20.0:
                RT_list = [min(RT_list), min(RT_list) + 5.0]
                
        pep_z_RT_dic [pep_Z_] = RT_list
                
        
            
    for i in range(0, len(scanfile_list)):
        # k = k + 1
        # printProgress(k, len(scanfile_list), 'Progress:', 'submission', 1, 50)
        #tqdm(i)
        
        #frac_data_df = data_df[data_df['scanfile'] == scanfile_list[i]]            
        process_data.append([pep_z_list, ms1_dic[scanfile_list[i]], pep_z_RT_dic])

        
        
    pool = multiprocessing.Pool(processes = process)

    output2 = pool.map(TIQ.quantification, process_data)
    #TIQ, SumQ,TSumQ, TIQ_3data, TIQ2_data, TIQ4_data, TIQ5_data, TIQ_RT

    # TIQ_DATA1 = []
    # SUMQ_DATA1 = []
    # TSUMQ_DATA1 = []
    # TIQ_3_DATA1 = []
    # TIQ2_DATA1 = []
    # TIQ4_DATA1 = []
    # TIQ5_DATA1 = []
    # RT_DATA1 = []

    for i in range(0, len(scanfile_list)):
        TIQ_data = []
        SumQ_data = []
        TSumQ_data = []
        TIQ_3_data = []
        TIQ2_data = []
        TIQ4_data = []
        TIQ5_data = []
        TIQ_RT_data = []
        for j in range(0, len(pep_z_list)):
            TIQ_data.append(output2[i][j][0])
            SumQ_data.append(output2[i][j][1])
            TSumQ_data.append(output2[i][j][2])
            TIQ_3_data.append(output2[i][j][3])
            TIQ2_data.append(output2[i][j][4])
            TIQ4_data.append(output2[i][j][5])
            TIQ5_data.append(output2[i][j][6])
            TIQ_RT_data.append(output2[i][j][7])
        TIQ_DATA.append(TIQ_data)
        SUMQ_DATA.append(SumQ_data)
        TSUMQ_DATA.append(TSumQ_data)
        TIQ_3_DATA.append(TIQ_3_data)
        TIQ2_DATA.append(TIQ2_data)
        TIQ4_DATA.append(TIQ4_data)
        TIQ5_DATA.append(TIQ5_data)
        RT_DATA.append(TIQ_RT_data)
        
        
        
        

    samples = scanfile_list
    print("--- %s seconds --- TIQ calculation" % (time.time() - start_time))

    TIQ_array = np.array(TIQ_DATA)
    
    TIQ_data_df = pd.DataFrame(TIQ_array.T, columns = samples, index = pep_z_list)
    
    TIQ_data_df.to_excel('TIQ_data_MP.xlsx')
    
    
    Sumq_array = np.array(SUMQ_DATA)
    Sumq_data_df = pd.DataFrame(Sumq_array.T, columns = samples, index = pep_z_list)
    Sumq_data_df.to_excel('SUMQ_data_MP.xlsx')
    
    TSumq_array = np.array(TSUMQ_DATA)
    TSumq_data_df = pd.DataFrame(TSumq_array.T, columns = samples, index = pep_z_list)
    TSumq_data_df.to_excel('TSUMQ_data_MP.xlsx')
    
    TIQ_3_array = np.array(TIQ_3_DATA)
    TIQ_3_data_df = pd.DataFrame(TIQ_3_array.T, columns = samples, index = pep_z_list)
    
    TIQ_3_data_df.to_excel('TIQ_3data_MP.xlsx')
    
    TIQ2_array = np.array(TIQ2_DATA)
    TIQ2_data_df = pd.DataFrame(TIQ2_array.T, columns = samples, index = pep_z_list)
    TIQ2_data_df.to_excel('TIQ2_data_MP.xlsx')
    
    TIQ4_array = np.array(TIQ4_DATA)
    TIQ4_data_df = pd.DataFrame(TIQ4_array.T, columns = samples, index = pep_z_list)
    TIQ4_data_df.to_excel('TIQ4_data_MP.xlsx')
    
    TIQ5_array = np.array(TIQ5_DATA)
    TIQ5_data_df = pd.DataFrame(TIQ5_array.T, columns = samples, index = pep_z_list)
    TIQ5_data_df.to_excel('TIQ5_data_MP.xlsx')

    protein_quant(TIQ_data_df, 'TIQ')
    protein_quant(Sumq_data_df, 'SUMQ')
    protein_quant(TSumq_data_df, 'TSUMQ')
    protein_quant(TIQ_3_data_df, '3TIQ')
    protein_quant(TIQ2_data_df, 'TIQ2')
    protein_quant(TIQ4_data_df, 'TIQ4')
    protein_quant(TIQ5_data_df, 'TIQ5')


#testfile.close()    
    