import math
import pywt
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spl
from python_utils import *

##########################################
####### FUNCTIONS BASED ON SVD ###########
##########################################

##########################################
# Partitioning of edge-space
##########################################

# Function : Sort the vertices according to the SVD
def mapping( edgelist ):	
	# Transform edgelist to scipy matrix format
	data = [1]*len(edgelist)
	row, col = [], []
	for e in edgelist:
		row.append( e[0] )
		col.append( e[1] )
	max_size = max(max(row), max(col)) + 1

	# Build matrix
	A = sps.coo_matrix((data,(row, col)), shape=(max_size, max_size), dtype=float).tocsr()

	# We get the second largest then we sort it 
	sorted_indices = [[i for i in range(max_size)]]
	curr_max_length = max([len(indices) for indices in sorted_indices])
	while curr_max_length > 2:
		tmp_indices = []
		for indices in sorted_indices:
			if len(indices) > 2:
				tmp_indices += SVD_partitioning(A, indices)
			else:
				tmp_indices += [indices]
		sorted_indices = tmp_indices
		curr_max_length = max([len(indices) for indices in sorted_indices])

	clean_list = [item for sublist in sorted_indices for item in sublist]
	return { elem : index+1 for index, elem in enumerate(clean_list) }


# Function : Split specific sub-matrices using the SVD
def SVD_partitioning(A, indices):
	A_sub = A[indices, :]
	u, _, _ = spl.svds(A_sub, k=2, which='LM', return_singular_vectors="u")

	sort_indices = np.argsort(-u[:,0].flatten())
	set1 = [indices[i] for i in sort_indices[: sort_indices.size//2]]
	set2 = [indices[i] for i in sort_indices[sort_indices.size//2 :]]
	return [set1, set2]

##########################################
# Graph Decomposition
##########################################

# Function : graph decomposition by mapping graph to time-series
def graph_decomposition( edgelist, mapping, level='auto' ):	
	x = np.zeros(next_power_of_2(len(mapping))**2)
	for e in edgelist:
		x[ edge_time_map( mapping[e[0]], mapping[e[1]] ) - 1 ] += 1
	
	# Get the level if not specified
	max_level = pywt.dwt_max_level(len(x), 'haar')	
	if level == 'auto': chosen_level = max_level
	else: chosen_level = min(max_level, level)
		
	# Decomposition
	coeffs = pywt.wavedec(x, 'haar', level=chosen_level)
	return coeffs

# Function : map the relabelled edge (i, j) to a point in the positive integers
def edge_time_map( i, j ): 
	
	# Base case
	max_num = max(i, j)
	if max_num == 1: 
		return 1

	# Get the k parameter
	if is_power_of_two( max_num ): 
		k = previous_power_2( max_num - 1 )
	else:
		k = previous_power_2( max_num )

	# Apply the function 
	value = 0
	if i <= k and j > k:
		value += k**2 + edge_time_map( i, j - k)	
	elif i > k and j <= k:
		value += 2*(k**2) + edge_time_map( i-k, j )
	elif i > k and j > k:
		value += 3*(k**2) + edge_time_map( i-k, j-k )

	return value

# Function : inverse transformation of the graph decomposition
def inv_graph_decomposition( coeffs, inv_mapping ):
	
	# Recovered time series from mra
	x_rec = pywt.waverec(coeffs, 'haar')
	
	# Mapping the time series back to graph
	inv_graph = set()
	for index in range(len(x_rec)): 
		if abs(x_rec[index]) > 1e-6:
			x, y = inv_edge_time_map(index + 1)
			if (x in inv_mapping) and (y in inv_mapping):
				u = inv_mapping[x]
				v = inv_mapping[y]
				inv_graph.add( (u, v, x_rec[index]) )
	return inv_graph

def inv_edge_time_map( index ):
	# Base case
	if index == 1:
		return 1, 1
	
	# Get the k parameter
	exp = math.floor( math.log2(index) )
	if index == 2**exp:
		k = 2**((exp-1)//2)
	else:
		k = 2**(exp//2)	
	
	# Apply the inverse function 
	i = 0
	j = 0
	if index > k**2 and index <= 2*(k**2):
		v1, v2 = inv_edge_time_map( index - k**2)
		i = v1 
		j = k + v2
	elif index > 2*(k**2) and index <= 3*(k**2):
		v1, v2 = inv_edge_time_map( index - 2*(k**2) )
		i = k + v1
		j = v2 
	elif index > 3*(k**2) and index <= 4*(k**2):
		v1, v2 = inv_edge_time_map( index - 3*(k**2) )
		i = k + v1 
		j = k + v2 	
	return i, j

##########################################
# Graph sequence decomposition 
##########################################

# Function : graph sequence decomposition
def graph_sequence_decomposition( graph_seq, mapping, level='auto', info='all'):
	
	# Get size of edgeset
	edgeset_size = next_power_of_2(len(mapping))**2

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
			x[ix%window, edge_time_map( mapping[e[0]], mapping[e[1]] ) - 1] += 1

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
