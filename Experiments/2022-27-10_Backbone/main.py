import sys
sys.path.append('../../Library')
import decomposition_utils
import data_processing
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import pywt
import math

def ts_to_graph_map( x_rec, inv_mapping ):
		
	# Mapping the time series back to graph
	inv_graph = set()
	for index in range(len(x_rec)): 
		if abs(x_rec[index]) > 1e-2:
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


def draw_graph(N, m):
	G = set()
	for e in range(m):
		# Chose edge
		orig = np.random.randint(N//2)
		dest = np.random.randint(N//2)
		# Avoid self loops
		while (dest == orig) or ((dest,orig) in G): dest = np.random.randint(N//2)
		# Chose community
		comm = int(np.round(np.random.rand()))
		# Add edge (undirected)
		G.add( (comm*(N//2) + orig, comm*(N//2) + dest))
		G.add( (comm*(N//2) + dest, comm*(N//2) + orig))
	return G

def X_to_L(X, num_coefs):
	L = np.zeros(X.shape)
	for i in range(X.shape[0]):
		x = [X[i, 0:num_coefs], X[i, num_coefs: 2*num_coefs]]
		step = 2*num_coefs
		last_val = 2*num_coefs
		while last_val < X.shape[1]:
			x.append( X[i, last_val : last_val + step] )
			last_val = last_val + step
			step = 2*step	
		L[i] = pywt.waverec(x, 'haar')
	L[abs(L) < 1e-5] = 0
	return L

#######################
### Create Raw Data ###
#######################
np.random.seed(10)
interval = 10
t_max = 320
N = 8
Link_stream = dict()

# Create a link stream with two comms and day-night activity
for t in range(t_max):	
	if (t//interval)%2 == 0:
		Link_stream[t] = draw_graph(N, 2)	
	else:
		Link_stream[t] = set()

# Fix some intra-edge communities
sampled_times = [i*interval + np.random.randint(interval) for i in range(t_max//interval) if i%2==0]

for curr_t in sampled_times:
	new_graph = draw_graph(N, 1)
	new_graph.add((3,5)) 
	new_graph.add((5,3))
	Link_stream[curr_t] = new_graph

# Fig the graphs from the figure
Link_stream[0] = {(0,2),(2,0), (5,7), (7,5)}
Link_stream[9] = {(0,3),(3,0),(6,7),(7,6)}
Link_stream[20] = {(0,1),(1,0),(1,3),(3,1)}
Link_stream[29] = {(3,5),(5,3),(4,5),(5,4)}

###################################
### DECOMPOSING THE LINK STREAM ###
###################################

# decomposing the edge-space
agg_graph = {e for t in Link_stream for e in Link_stream[t]}
mapping = decomposition_utils.mapping_edgespace(agg_graph, 'SVD')

# Taking the Time-Structural representation
X = [np.hstack(decomposition_utils.graph_decomposition(Link_stream[t], mapping, 'SVD', level=4)) for t in Link_stream]
X = np.vstack(X)

# Taking the Time-Structural representation
C = np.fft.rfft(X, norm='ortho', axis=0)
freq_axis = np.around(np.fft.rfftfreq(X.shape[0]), 3)

###############
#### FIGURE ###
###############
structures = ['$s_0^{(4)}$',
			  '$s_{1}^{(4)}$', 
			  '$s_{2}^{(4)}$', 
			  '$s_{3}^{(4)}$', 
			  '$w_{0}^{(4)}$', 
			  '$w_{1}^{(4)}$',
			  '$w_{2}^{(4)}$', 
			  '$w_{3}^{(4)}$', 
			  '$w_{0}^{(3)}$', 
			  '$w_{1}^{(3)}$', 
			  '$w_{2}^{(3)}$', 
			  '$w_{3}^{(3)}$',
			  '$w_{4}^{(3)}$', 
			  '$w_{5}^{(3)}$', 
			  '$w_{6}^{(3)}$', 
			  '$w_{7}^{(3)}$']

plt.figure(figsize=(6,2.5))
plt.imshow(np.power(np.abs(C),1), aspect='auto', cmap='binary')
plt.colorbar(orientation='horizontal', aspect=80)
plt.xlim([-1,15.5])
plt.ylim([29.5,-2.5])

ax = plt.gca()
ax.set_ylabel('Frequency', fontsize=13)
ax.set_xlabel('Structure', fontsize=13)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xticks(list(range(16)))
ax.set_xticklabels(structures, fontsize=11)
ax.set_yticks([i for i in range(30) if i%8==0])
ax.set_yticklabels(freq_axis[[i for i in range(30) if i%8==0]], fontsize=11)
plt.savefig('backbone_spectrum.png', dpi=300, bbox_inches='tight')

#################
#### BACKBONE ###
#################
C_hat = np.copy(C)
C_hat[:,4:] = 0
C_hat[17:,:] = 0

X_hat = np.fft.irfft(C_hat, norm='ortho', axis=0)
L_hat = X_to_L(X_hat, 4)

plt.figure(figsize=(6,2.5))
plt.imshow(abs(L_hat), aspect='auto', cmap='binary')
plt.colorbar(orientation='horizontal', aspect=80)
plt.ylim([40,-1])
plt.xticks(fontsize=11)
ax = plt.gca()
ax.set_ylabel('Time', fontsize=13)
ax.set_xlabel('Relations', fontsize=13)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_yticks([0,10,20,30,40])
ax.set_yticklabels([0,10,20,30,40])
ax.set_xticks([0,10,20,30,40,50,60])
ax.set_xticklabels(['$e_0$','$e_{10}$','$e_{20}$','$e_{30}$','$e_{40}$','$e_{50}$','$e_{60}$'], fontsize=11)
plt.savefig('backbone_ls.png', dpi=300, bbox_inches='tight')

########################
#### GRAPHS AT TIMES ###
########################
inv_mapping = {val:key for key,val in mapping.items()}
G_5 = ts_to_graph_map( np.abs(L_hat[5,:]), inv_mapping )

########################
#### GRAPHS AT TIMES ###
########################
max_val = np.max(abs(L_hat))
min_val = 0
norm = mpl.colors.Normalize(vmin=min_val, vmax=max_val)
cmap = cm.binary
m = cm.ScalarMappable(norm=norm, cmap=cmap)

for e in G_5:
	print('(',e[0], e[1],')', ' with code RGB = ', m.to_rgba(e[2]))
