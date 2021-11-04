import scipy
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib as mpl
import collections
from collections import defaultdict
import matplotlib
from matplotlib import pyplot as plt
import itertools
from scipy.stats import hypergeom
import pickle
import scipy as sci


matplotlib.rcParams['pdf.fonttype'] = 42

#Create an N-dimentional nested dictionary
def nested_dict(n, type):
    if n == 1:
        return defaultdict(type)
    else:
        return defaultdict(lambda: nested_dict(n-1, type))

#Function to define the turning point of a cumulative curve by using the minimum derivative method
def turn_point(sgRNA_name, df):
    sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
    sgRNA_cumsum = sgRNA_count.cumsum()

    #get the total cell number of this sgRNA
    cell_num = sum(sgRNA_count > 0)

    #calculate the turning point by using the max derivative
    if sgRNA_count.sum() == 0:
        turning_point = 0
    else:
        turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape
    
    return(sgRNA_count.iloc[turning_point])

#load sgRNA-barcode data
def load_data(input_file):
    input_fh  = open(input_file, 'r')
    #define parameters
    cell_bc_list   = []
    num_sgRNA_list = np.array([])
    sgRNAs         = []
    umis           = []
    #initiate a 2D dictionary
    data_dict = nested_dict(2, list)
    for line in input_fh:
        cell_bc    = line.strip().split('\t')[0] + '-1'
        num_sgRNA  = line.strip().split('\t')[2]
        sgRNA_list = line.strip().split('\t')[3].split(';')
        umi_list   = line.strip().split('\t')[5].split(';')
        for i in zip(sgRNA_list, umi_list):
            data_dict[cell_bc][i[0]] = i[1]
    #read the 2D dictionary into the pandas DataFrame
    df = pd.DataFrame(data_dict).fillna(0).astype(int)
    return df

#filter the sgRNA UMI count based on the cutoff values
def filter_umi (df_input, replace=True):
    if replace:
        df = df_input.copy()
    else:
        df = df_input
    #calculate the cutoff value for each sgRNA in the dataset
    cutoffs = []
    sgRNA_cutoff = [turn_point(i, df) for i in list(df.index)]
    cutoffs.append(sgRNA_cutoff)
    for i in range(0, len(sgRNA_cutoff)):
        df.iloc[i, :].loc[df.iloc[i, :] <= sgRNA_cutoff[i]] = 0
    return df, cutoffs


#### find_faulty_sg
def gen_sg_combos(target_dict, target):
    c_list = []
    c_str = []
    t_combos = []
    
    c_iters = list(range(1, len(target_dict[target]) + 1))
    for num in c_iters:
        iter_combos = itertools.combinations(list(range(1, len(target_dict[target]) + 1)), num)
        for com in iter_combos:
            curr_combos = []
            [curr_combos.append(target + ':' + str(i)) for i in com]
            c_list.append(list(com))
            c_str.append(str(com))
            t_combos.append(curr_combos)
    c_dict = dict(zip(c_str,c_list))
    
    return(c_list,c_str,t_combos)

def cluster_bias(target_cells, control_cells, t_matrix, direction = 'depletion'):
    clus_pvals = []
    com_cells = np.unique(target_cells + control_cells)
    t_total = len(np.unique(target_cells))
    back_total = len(np.unique(target_cells + control_cells))
    clus_list = t_matrix.obs.louvain.cat.categories.tolist()
    
    for cluster in clus_list:
        sg_cluster_count = sum(t_matrix[target_cells].obs.louvain == cluster)
        total_sg_cluster_count = sum(t_matrix[com_cells].obs.louvain == cluster)
        
        if direction == 'depletion':
            p_value = hypergeom.cdf(int(sg_cluster_count), int(back_total), int(t_total), int(total_sg_cluster_count))
        
        elif direction == 'enrichment':
            p_value = 1 - hypergeom.cdf(int(sg_cluster_count) - 1, int(back_total), int(t_total), int(total_sg_cluster_count))
        
        clus_pvals.append(p_value)
    return(clus_pvals, clus_list)

def CB_Filter_SG(pval_path, target, df_suffix, dict_suffix):
    pval_df_dep = pd.read_csv(pval_path + target + df_suffix + 'Depletion.csv', index_col = 0)
    pval_df_enr = pd.read_csv(pval_path + target + df_suffix + 'Enrichment.csv', index_col = 0)
    pval_df_dep_min = pval_df_dep.shape[0]/20
    pval_df_enr_min = pval_df_enr.shape[0]/20

    Dict_File = open(pval_path + target + dict_suffix + ".pkl", "rb")
    pval_dict = pickle.load(Dict_File)

    target_sg = len(pval_dict[list(pval_dict.keys())[-1]]) 
    cluster_list = pval_df_dep.columns

    ### Order pvals for each cluster. Then find turning point by finding the index value in which pval undergoes steepest change
    Biased_SG = []
    for cluster in cluster_list:
        Highest_Change = 0
        for i, value in enumerate(pval_df_dep.sort_values(cluster)[cluster]):
            if i <= pval_df_dep_min or i >= pval_df_dep.shape[0] - pval_df_dep_min :
                continue
            try:
                if abs(pval_df_dep.sort_values(cluster)[cluster][i+1] - pval_df_dep.sort_values(cluster)[cluster][i]) > Highest_Change:
                    Highest_Change = abs(pval_df_dep.sort_values(cluster)[cluster][i+1] - pval_df_dep.sort_values(cluster)[cluster][i])
                    Index_Pos = i + 1
            except:
                continue
        Turning_Combos = pval_df_dep.sort_values(cluster).iloc[Index_Pos:].index
        Selected_Combos = []
        [Selected_Combos.append(pval_dict[i]) for i in Turning_Combos]
        Frequency_Array = []
        for sg in list(range(1,target_sg+1)):
            sg_count = []
            for comb in Selected_Combos:
                if sg in comb:
                    sg_count.append(1)
                else:
                    sg_count.append(0)
            Frequency_Array.append((sum(sg_count)/len(Selected_Combos))*100)
        Initial_Biased_SG = []
        Frequency_Mean = np.mean(Frequency_Array)
        Frequency_STD = 2*(np.std(Frequency_Array))
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) > Frequency_Mean + Frequency_STD)]
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) < Frequency_Mean - Frequency_STD)] 

        Turning_Combos = pval_df_dep.sort_values(cluster).iloc[:Index_Pos].index
        Selected_Combos = []
        [Selected_Combos.append(pval_dict[i]) for i in Turning_Combos]
        Frequency_Array = []
        for sg in list(range(1,target_sg+1)):
            sg_count = []
            for comb in Selected_Combos:
                if sg in comb:
                    sg_count.append(1)
                else:
                    sg_count.append(0)
            Frequency_Array.append((sum(sg_count)/len(Selected_Combos))*100)
        Frequency_Mean = np.mean(Frequency_Array)
        Frequency_STD = 2*(np.std(Frequency_Array))
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) > Frequency_Mean + Frequency_STD)]
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) < Frequency_Mean - Frequency_STD)] 
        if len(Initial_Biased_SG) > 0:        
            [Biased_SG.append(np.unique(Initial_Biased_SG))]

        Highest_Change = 0
        for i, value in enumerate(pval_df_enr.sort_values(cluster)[cluster]):
            if i <= pval_df_enr_min or i >= pval_df_enr.shape[0] - pval_df_enr_min :
                continue
            try:
                if abs(pval_df_enr.sort_values(cluster)[cluster][i+1] - pval_df_enr.sort_values(cluster)[cluster][i]) > Highest_Change:
                    Highest_Change = abs(pval_df_enr.sort_values(cluster)[cluster][i+1] - pval_df_enr.sort_values(cluster)[cluster][i])
                    Index_Pos = i + 1
            except:
                continue
        Turning_Combos = pval_df_enr.sort_values(cluster).iloc[Index_Pos:].index
        Selected_Combos = []
        [Selected_Combos.append(pval_dict[i]) for i in Turning_Combos]
        Frequency_Array = []
        for sg in list(range(1,target_sg+1)):
            sg_count = []
            for comb in Selected_Combos:
                if sg in comb:
                    sg_count.append(1)
                else:
                    sg_count.append(0)
            Frequency_Array.append((sum(sg_count)/len(Selected_Combos))*100)
        Initial_Biased_SG = []
        Frequency_Mean = np.mean(Frequency_Array)
        Frequency_STD = 2*(np.std(Frequency_Array))
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) > Frequency_Mean + Frequency_STD)]
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) < Frequency_Mean - Frequency_STD)] 


        Turning_Combos = pval_df_enr.sort_values(cluster).iloc[:Index_Pos].index
        Selected_Combos = []
        [Selected_Combos.append(pval_dict[i]) for i in Turning_Combos]
        Frequency_Array = []
        for sg in list(range(1,target_sg+1)):
            sg_count = []
            for comb in Selected_Combos:
                if sg in comb:
                    sg_count.append(1)
                else:
                    sg_count.append(0)
            Frequency_Array.append((sum(sg_count)/len(Selected_Combos))*100)
        Frequency_Mean = np.mean(Frequency_Array)
        Frequency_STD = 2*(np.std(Frequency_Array))
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) > Frequency_Mean + Frequency_STD)]
        [Initial_Biased_SG.append(i[0]) for i in np.argwhere(np.array(Frequency_Array) < Frequency_Mean - Frequency_STD)] 
        if len(Initial_Biased_SG) > 0:        
            [Biased_SG.append(np.unique(Initial_Biased_SG))]



    faulty_sg = []
    if len(Biased_SG) > 0:
        faulty_sg = np.unique([i[0]+ 1 for i in Biased_SG]).tolist()
        
    return(faulty_sg)


def CB_Filter_SG_Stat(pval_path, target, df_suffix, dict_suffix, total_sg):
    pval_df_dep = pd.read_csv(pval_path + target + df_suffix + 'Depletion.csv', index_col = 0)
    pval_df_enr = pd.read_csv(pval_path + target + df_suffix + 'Enrichment.csv', index_col = 0)

    Dict_File = open(pval_path + target + dict_suffix + ".pkl", "rb")
    pval_dict = pickle.load(Dict_File)

    target_sg = len(pval_dict[list(pval_dict.keys())[-1]]) 
    cluster_list = pval_df_dep.columns

    ### Order pvals for each cluster. Then find turning point by finding the index value in which pval undergoes steepest change
    Biased_SG = []
    Enr_SG = []
    Dep_SG = []
    for cluster in cluster_list:
        sg_count = []
        for i in pval_dict.keys():
            if 1 in pval_dict[i]:
                sg_count.append(1)

        Index_Pos = sum(sg_count)
        Before_Turning_Index = pval_df_dep.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_dep.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos

        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))

            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)


        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]


        Before_Turning_Index = pval_df_enr.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_enr.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos

        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))

            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)


        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]


    faulty_enr_sg = []
    faulty_dep_sg = []
    faulty_sg = []

    if len(np.unique(Biased_SG)) < target_sg/2:
        
        if len(Enr_SG) > 0:
            faulty_enr_sg = np.unique(Enr_SG).tolist()

        if len(Dep_SG) > 0:
            faulty_dep_sg = np.unique(Dep_SG).tolist()

        if len(Biased_SG) > 0:
            faulty_sg = np.unique(Biased_SG).tolist()

    return(faulty_enr_sg, faulty_dep_sg, faulty_sg)

def CB_Filter_SG_Stat_ALL(pval_path, target, df_suffix, dict_suffix, total_sg):
    pval_df_dep = pd.read_csv(pval_path + target + df_suffix + 'Depletion_ALL.csv', index_col = 0)
    pval_df_enr = pd.read_csv(pval_path + target + df_suffix + 'Enrichment_ALL.csv', index_col = 0)

    Dict_File = open(pval_path + target + dict_suffix + ".pkl", "rb")
    pval_dict = pickle.load(Dict_File)

    target_sg = len(pval_dict[list(pval_dict.keys())[-1]]) 
    cluster_list = pval_df_dep.columns

    ### Order pvals for each cluster. Then find turning point by finding the index value in which pval undergoes steepest change
    Biased_SG = []
    Enr_SG = []
    Dep_SG = []
    for cluster in cluster_list:
        sg_count = []
        for i in pval_dict.keys():
            if 1 in pval_dict[i]:
                sg_count.append(1)

        Index_Pos = sum(sg_count)
        Before_Turning_Index = pval_df_dep.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_dep.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos


        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))

            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)


        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]


        Before_Turning_Index = pval_df_enr.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_enr.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos

        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))

            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)


        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]
    


    faulty_enr_sg = []
    faulty_dep_sg = []
    faulty_sg = []

    if len(np.unique(Biased_SG)) < target_sg/2:
        
        if len(Enr_SG) > 0:
            faulty_enr_sg = np.unique(Enr_SG).tolist()

        if len(Dep_SG) > 0:
            faulty_dep_sg = np.unique(Dep_SG).tolist()

        if len(Biased_SG) > 0:
            faulty_sg = np.unique(Biased_SG).tolist()


    return(faulty_enr_sg, faulty_dep_sg, faulty_sg)

def CB_Filter_NC(pval_df_dep, pval_df_enr, pval_dict):


    target_sg = len(pval_dict[list(pval_dict.keys())[-1]]) 
    cluster_list = pval_df_dep.columns

    ### Order pvals for each cluster. Then find turning point by finding the index value in which pval undergoes steepest change
    Biased_SG = []
    Enr_SG = []
    Dep_SG = []
    for cluster in cluster_list:
        sg_count = []
        for i in pval_dict.keys():
            if 1 in pval_dict[i]:
                sg_count.append(1)

        Index_Pos = sum(sg_count)
        Before_Turning_Index = pval_df_dep.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_dep.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos

        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))

            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)


        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]


        Before_Turning_Index = pval_df_enr.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_enr.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos

        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))

            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)


        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]


    faulty_enr_sg = []
    faulty_dep_sg = []
    faulty_sg = []

    if len(np.unique(Biased_SG)) < target_sg/2:
        
        if len(Enr_SG) > 0:
            faulty_enr_sg = np.unique(Enr_SG).tolist()

        if len(Dep_SG) > 0:
            faulty_dep_sg = np.unique(Dep_SG).tolist()

        if len(Biased_SG) > 0:
            faulty_sg = np.unique(Biased_SG).tolist()

    return(faulty_enr_sg, faulty_dep_sg, faulty_sg)