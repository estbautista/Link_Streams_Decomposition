import networkx as nx
import math
import sys
sys.path.append('../../Library')
import decomposition_utils
import data_processing
import numpy as np
import matplotlib.pyplot as plt


########################
# REFERENCE GRAPH #
########################
num_com = 2
nod_com = 16
ref_graph = nx.planted_partition_graph(num_com, nod_com, p_in=1, p_out=0.05)
ref_graph = nx.to_edgelist( ref_graph )
reg_graph = [(e1, e2) for e1,e2,e3 in ref_graph] + [(e2, e1) for e1,e2,e3 in ref_graph] 

########################
# GRAPH DECOMPOSE DATA #
########################
mapping = decomposition_utils.mapping_edgespace(ref_graph, 'SVD')

########################
# GRAPH 1 			#
########################
graph_1 = nx.planted_partition_graph(num_com, nod_com, p_in=0.5, p_out=0.01, seed=1)
graph_1 = nx.to_edgelist( graph_1 )
graph_1 = [(e1, e2) for e1,e2,e3 in graph_1] + [(e2, e1) for e1,e2,e3 in graph_1] 
print('num_edges G1 = ', len(graph_1))


########################
# GRAPH 2 			#
########################
graph_2 = nx.planted_partition_graph(num_com, nod_com, p_in=0.5, p_out=0.01, seed=2)
graph_2 = nx.to_edgelist( graph_2 )
graph_2 = [(e1, e2) for e1,e2,e3 in graph_2] + [(e2, e1) for e1,e2,e3 in graph_2] 
print('num_edges G2 = ', len(graph_2))

########################
# DECOMPOSITION		#
########################
my_level = 8
Dec1 = decomposition_utils.graph_decomposition( graph_1, mapping, 'SVD', level=my_level )
Dec2 = decomposition_utils.graph_decomposition( graph_2, mapping, 'SVD', level=my_level )

coef1 = np.concatenate(Dec1).ravel()
coef2 = np.concatenate(Dec2).ravel()

coef_dif = np.power(coef1 - coef2, 2)

def plot_spectrum(spec, axis):	
	# levels of resolution
	bins = [2**(10 - my_level + i) for i in range(0, my_level)]
	# plot figure
	axis.plot(np.arange(1, spec.size), spec[1:], color='black')
	axis.vlines(x=bins, ymin=0, ymax=max(spec[1:]), colors='purple', ls='--', lw=1, label='levels')
	axis.set_xscale('log')
	axis.spines['right'].set_visible(False)
	axis.spines['top'].set_visible(False)
	return axis

plt.rcParams.update({'font.size': 10})
f, axs = plt.subplots(1, figsize=(8, 2.6))
axs = plot_spectrum(coef_dif, axs)
axs.set_xlabel('coefficient', fontsize=12)
axs.set_ylabel('squared difference', fontsize=12)
plt.title('Edit distance distribution', fontsize=12)
plt.savefig('edit_distance_fig.png', dpi=300, format='png', bbox_inches='tight')
plt.show()
