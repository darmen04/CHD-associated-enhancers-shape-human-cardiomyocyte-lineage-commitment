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


def CB_Filter_SG_Stat(pval_path, target, df_suffix, dict_suffix, total_sg):
    ### Load sgRNA Combo pvals for each target
    pval_df_dep = pd.read_csv(pval_path + target + df_suffix + 'Depletion.csv', index_col = 0)
    pval_df_enr = pd.read_csv(pval_path + target + df_suffix + 'Enrichment.csv', index_col = 0)

    ### Dictionary: Converts index name to guide combo
    Dict_File = open(pval_path + target + dict_suffix + ".pkl", "rb")
    pval_dict = pickle.load(Dict_File)

    target_sg = len(pval_dict[list(pval_dict.keys())[-1]]) 
    cluster_list = pval_df_dep.columns

    ### Order pvals for each cluster. Then set a cutoff at middle of pval list.
    Biased_SG = []
    Enr_SG = []
    Dep_SG = []
    for cluster in cluster_list:
        sg_count = []
        for i in pval_dict.keys():
            if 1 in pval_dict[i]:
                sg_count.append(1)
        
        ### For depletion in cluster pvals
        Index_Pos = sum(sg_count)
        Before_Turning_Index = pval_df_dep.sort_values(cluster).iloc[:Index_Pos].index
        Before_Combos = []
        [Before_Combos.append(pval_dict[i]) for i in Before_Turning_Index]

        M_tot = pval_df_dep.shape[0]
        N_draw = Index_Pos
        n_tot = Index_Pos

        Enriched_Biased_SG = []
        Depleted_Biased_SG = []

        ### for each sg in target, calculates hypergeometric
        for sg in list(range(1,target_sg+1)):
            Inc_Vals = []
            for com in Before_Combos:
                if sg in com:
                    Inc_Vals.append(1)

            X_Vals = (sum(Inc_Vals))
            ### Fewer sg in significant section
            if hypergeom.cdf(X_Vals, M_tot, n_tot, N_draw) < 0.05:
                Depleted_Biased_SG.append(sg)

            ### More sg in significant section
            if 1 - hypergeom.cdf(X_Vals - 1, M_tot, n_tot, N_draw) < 0.05:
                Enriched_Biased_SG.append(sg)

        ### Only consider a guide biased if no other sg follows trend,
        if len(Depleted_Biased_SG) == 1:        
            [Dep_SG.append(i) for i in Depleted_Biased_SG]
            [Biased_SG.append(i) for i in Depleted_Biased_SG]

        if len(Enriched_Biased_SG) == 1:        
            [Enr_SG.append(i) for i in Enriched_Biased_SG]
            [Biased_SG.append(i) for i in Enriched_Biased_SG]

        ### Repeat for enrichment in cluster pvals
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

    ### If more than half of sgRNA would be removed, remove none.
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