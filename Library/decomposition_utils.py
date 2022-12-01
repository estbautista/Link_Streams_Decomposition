import numpy as np
import decomposition_BFS as g_dec_BFS
import decomposition_SVD as g_dec_SVD
from python_utils import *

# edge-space mapping 
def mapping_edgespace(graph, method):
	if method == 'BFS':
		mapping = g_dec_BFS.mapping( graph )
		return mapping
	elif method == 'SVD':
		mapping = g_dec_SVD.mapping( graph )
		return mapping
	else:
		print('Error. Please choose SVD or BFS as method')
		return 0

# General graph decomposition
def graph_decomposition( graph, mapping, method, level='auto' ): 
	if method == 'BFS':
		coef = g_dec_BFS.graph_decomposition( graph, mapping, level)
		return coef
	elif method == 'SVD':
		coef = g_dec_SVD.graph_decomposition( graph, mapping, level)
		return coef
	else:
		print('Error. Please choose SVD or BFS as method')
		return 0

# Decomposition of graph sequence
def graph_sequence_decomposition( graph_seq, mapping, method, level='auto', info='all'):	
	# Get the theoretical size of the edge-space
	if method == 'BFS':
		coef = g_dec_BFS.graph_sequence_decomposition( graph_seq, mapping, level, info )
		return coef
	elif method == 'SVD':
		coef = g_dec_SVD.graph_sequence_decomposition( graph_seq, mapping, level, info )
		return coef
	else:
		print('Error. Please choose SVD or BFS as method')
		return 0
