#%%

import numpy as np
import pandas as pd
from SERGIO.sergio import sergio

# Simulate Clean Data _ differentiation Simulation

#%%

df = pd.read_csv('data_sets/Gene_100_cell_2/bMat_cID100.tab', sep='\t', header=None, index_col=None)
bMat = df.values

sim = sergio(number_genes=100, number_bins = 2,
             number_sc = 100, noise_params = 0.00001,
             decays=0.8, sampling_state = 1,
             noise_params_splice = 0.00001, noise_type='sp', dynamics=True, bifurcation_matrix= bMat)
sim.build_graph(input_file_taregts ='data_sets/Gene_100_cell_2/Interaction_100.txt',
                input_file_regs='data_sets/Gene_100_cell_2/Regs_100.txt', shared_coop_state=2)
sim.simulate_dynamics()
exprU, exprS = sim.getExpressions_dynamics()
count_matrix_U1, count_matrix_S1 = sim.convert_to_UMIcounts_dynamics(exprU, exprS)
count_matrix_U1 = np.concatenate(count_matrix_U1, axis = 1)
count_matrix_S1 = np.concatenate(count_matrix_S1, axis = 1)

#%% md

# Add Technical Noise _ differentiation Simulations

#%%

"""
Add outlier genes
"""
exprU_O, exprS_O = sim.outlier_effect_dynamics(exprU, exprS, outlier_prob = 0.01, mean = 0.8, scale = 1)

"""

Add Library Size Effect
"""
libFactor, exprU_O_L, exprS_O_L = sim.lib_size_effect_dynamics(exprU_O, exprS_O, mean = 4.6, scale = 0.4)

"""
Add Dropouts
"""
binary_indU, binary_indS = sim.dropout_indicator_dynamics(exprU_O_L, exprS_O_L, shape = 6.5, percentile = 82)
exprU_O_L_D = np.multiply(binary_indU, exprU_O_L)
exprS_O_L_D = np.multiply(binary_indS, exprS_O_L)

"""
Convert to UMI count
"""
count_matrix_U, count_matrix_S = sim.convert_to_UMIcounts_dynamics(exprU_O_L_D, exprS_O_L_D)

"""
Make 2d spliced and unspliced expression matrices
"""
count_matrix_U = np.concatenate(count_matrix_U, axis = 1)
count_matrix_S = np.concatenate(count_matrix_S, axis = 1)
import pandas as pd
count_matrix_U = pd.DataFrame(count_matrix_U1)
writer = pd.ExcelWriter('./data_sets/Gene_100_cell_2/G100_C2_S100U.xlsx')		# 写入Excel文件
count_matrix_U.to_excel(writer, float_format='%.5f')		# ‘page_1’是写入excel的sheet名
writer.save()
writer.close()

count_matrix_S = pd.DataFrame(count_matrix_S1)
writer = pd.ExcelWriter('./data_sets/Gene_100_cell_2/G100_C2_S100S.xlsx')		# 写入Excel文件
count_matrix_S.to_excel(writer, float_format='%.5f')		# ‘page_1’是写入excel的sheet名
writer.save()
writer.close()

