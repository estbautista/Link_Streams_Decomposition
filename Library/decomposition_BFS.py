import pywt
import numpy as np
import scipy.sparse as sps
from collections import deque
from python_utils import *

##########################################
####### FUNCTIONS BASED ON BFS ###########
##########################################

##########################################
# Partitioning of edge-space
##########################################

# Function : Mapping of edge-space to the integer line
def mapping(edgelist):
	'edgelist can be a set of or list with tuples'

	# Make a copy of the data that will be transformed
	tmp_edgelist = [{e for e in edgelist}]

	# Recursive partitioning of edge-space based on BFS
	while len(tmp_edgelist[0]) > 1:
		partition = []

		# E_i^{(-j)} -> E_{2i}^{(-j+1)} U E_{2i+1}^{(-j+1)}
		for curr_set in tmp_edgelist:
			if len(curr_set) > 1:
				partition += partition_edges( curr_set )
			else:
				partition += [curr_set]
		tmp_edgelist = partition

	# Mapping
	mapping = dict()
	for index, elem in enumerate(tmp_edgelist):
		mapping[ get_set_elem(elem) ] = index
	
	return mapping


# Function : Implementation of BFS to split edge-space 
def partition_edges( edgelist ):
	" The function should always receive a set that can be divided in two "

	# Get the target size
	m_size = next_power_of_2( len(edgelist) ) // 2

	# Transform edgelist into adjacency list
	adj_list = adjlist_from_edgelist( edgelist )
	
	# Get the list of source nodes to explore 
	to_explore = set(adj_list.keys())

	# set to store selected edges
	queue = deque()
	split = set()

	# add edges until we cover the quota
	while len(split) < m_size:		
		if len(queue) == 0:
			source_node = get_set_elem( to_explore )
			to_explore.remove( source_node )
			queue.append( source_node )

		curr_node = queue.popleft()
		for neighbor in adj_list[curr_node]:
			
			if len(split) < m_size:
				split.add( (curr_node, neighbor) )
				edgelist.remove( (curr_node, neighbor ) )
			else: break

			if neighbor in to_explore:
				queue.append( neighbor )
				to_explore.remove( neighbor )	
	
	return [split, edgelist]

##########################################
# Graph Decomposition
##########################################

# Function : graph decomposition by mapping graph to time-series
def graph_decomposition( edgelist, mapping, level='auto' ):	
	x = np.zeros(next_power_of_2(len(mapping)))
	for e in edgelist:
		x[mapping[(e[0], e[1])]] += 1
	
	max_level = pywt.dwt_max_level(len(x), 'haar')
	
	if level == 'auto': chosen_level = max_level
	else: chosen_level = min(max_level, level)
		
	coeffs = pywt.wavedec(x, 'haar', level=chosen_level)
	return coeffs

# Function : inverse transformation
def inv_graph_decomposition( coeffs, inv_mapping ):
	x_hat = pywt.waverec(coeffs, 'haar')
	inv_graph = set()
	for ix in range(len(inv_mapping)): 
		if abs(x_hat[ix]) > 1e-6:
			s, d = inv_mapping[ix]
			inv_graph.add( (s, d, x_hat[ix]) )	
	return inv_graph

##########################################
# Graph sequence decomposition 
##########################################

# Function : graph sequence decomposition
def graph_sequence_decomposition( graph_seq, mapping, level='auto', info='all'):
	
	# Get size of edgeset
	edgeset_size = next_power_of_2(len(mapping))

	# Choice of level
	max_level = pywt.dwt_max_level(edgeset_size, 'haar')	
	if level == 'auto': chosen_level = max_level
	else: chosen_level = min(max_level, level)

	# Parameters for decomposition in chunks
	num_snaps = len(graph_seq)
	window = 100
	x = np.zeros( [min(window, num_snaps),  edgeset_size] )

	# Store decomp chunks in sparse format
	sig_dict = dict()
	sig_dict['scaling'] = { chosen_level: [] }
	sig_dict['wavelet'] = { cl : [] for cl in range(chosen_level, 0, -1)}

	for ix, t in enumerate(graph_seq):

		# Fill time series
		for e in graph_seq[t]:
			x[ix%window, mapping[(e[0], e[1])]] += 1

		# Perform decomp
		if (ix%window == 0 and ix > 0) or ix == (num_snaps-1):	
			coeffs = pywt.wavedec(x, 'haar', level=chosen_level, axis=1)
			sig_dict['scaling'][chosen_level].append( coeffs[0] )
			if info == 'all':
				for cn, vec in enumerate(coeffs[1:]):
					sig_dict['wavelet'][chosen_level-cn].append( sps.csr_matrix(vec) )

			# Clean matrix	
			x = np.zeros( [min(window, num_snaps-ix), edgeset_size] )

	# Merge scaling coefficients
	sig_dict['scaling'][chosen_level] = np.vstack( sig_dict['scaling'][chosen_level] )

	# Merge wavelet coefficients
	if info == 'all':
		for cl in sig_dict['wavelet']:
			sig_dict['wavelet'][cl] = sps.vstack( sig_dict['wavelet'][cl] )
	
	return sig_dict
