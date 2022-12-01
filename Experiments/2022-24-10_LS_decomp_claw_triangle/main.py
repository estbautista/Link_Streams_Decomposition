import sys
sys.path.append('../../Library')
import decomposition_utils
import data_processing
import numpy as np
import matplotlib.pyplot as plt
import pywt


##############################
## GRAPHS IN ENCODED FORMAT ## 
##############################
claw = np.zeros(16)
claw[[2,3,4,5,6,7]] = 1

tria = np.zeros(16)
tria[[10,11,12,13,14,15]] = 1

#################
## LINK STREAM ## 
#################
t_max = 6
L = np.zeros([t_max, 16])
for t in range(t_max):
	if t%2==0:
		L[t]=claw
	else:
		L[t]=tria

################
## AXIS TICKS ##
################
relations = ['$e_1$', '$e_2$', '$e_3$', '$e_4$', '$e_5$', '$e_6$', '$e_7$', '$e_8$',
			 '$e_9$', '$e_{10}$', '$e_{11}$', '$e_{12}$', '$e_{13}$', '$e_{14}$', '$e_{15}$', '$e_{16}$']

structures = ['$s_0^{(3)}$', '$s_{1}^{(3)}$', '$w_{0}^{(3)}$', '$w_{1}^{(3)}$', '$w_{0}^{(2)}$', '$w_{1}^{(2)}$',
				'$w_{2}^{(2)}$', '$w_{3}^{(2)}$', '$w_{0}^{(1)}$', '$w_{1}^{(1)}$', '$w_{2}^{(1)}$', '$w_{3}^{(1)}$',
				'$w_{4}^{(1)}$', '$w_{5}^{(1)}$', '$w_{6}^{(1)}$', '$w_{7}^{(1)}$']

#####################
## TIME-RELAT REPR ## 
#####################
fig, ax = plt.subplots(figsize=(6,2))
ax.set_ylabel('Time', fontsize=13)
ax.set_xlabel('Relation', fontsize=13)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xticks(list(range(16)))
ax.set_xticklabels(relations)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.imshow(L, cmap='RdBu', vmin=-3, vmax=3, aspect='auto')
plt.savefig('link_stream.png', dpi=300, bbox_inches='tight')


######################
## TIME-STRUCT REPR ## 
######################
TS_coef = pywt.wavedec(L, 'haar', level=3, axis=1)
TS_coef = np.hstack(TS_coef)

fig, ax = plt.subplots(figsize=(6,2))
ax.set_ylabel('Time', fontsize=13)
ax.set_xlabel('Structure', fontsize=13)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xticks(list(range(16)))
ax.set_xticklabels(structures)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.imshow(TS_coef, cmap='RdBu', vmin=-3, vmax=3, aspect='auto')
plt.savefig('time_struct.png', dpi=300, bbox_inches='tight')

#####################
## FREQ-RELAT REPR ##
#####################
FR_coef = np.fft.rfft(L, norm='ortho', axis=0)
freqs = np.around(np.fft.rfftfreq(t_max), 2)
print(freqs)

fig, ax = plt.subplots(figsize=(6,2))
ax.set_ylabel('Frequency', fontsize=13)
ax.set_xlabel('Relation', fontsize=13)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xticks(list(range(16)))
ax.set_xticklabels(relations)
ax.set_yticks(list(range(0,t_max//2+1)))
ax.set_yticklabels(freqs)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.imshow(np.abs(FR_coef), cmap='RdBu', vmin=-3, vmax=3, aspect='auto')
plt.savefig('freq_relat.png', dpi=300, bbox_inches='tight')

######################
## FREQ-STRUCT REPR ##
######################
FS_coef = np.fft.rfft(TS_coef, norm='ortho', axis=0)

fig, ax = plt.subplots(figsize=(6,2))
ax.set_ylabel('Frequency', fontsize=13)
ax.set_xlabel('Structure', fontsize=13)
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')
ax.set_xticks(list(range(16)))
ax.set_xticklabels(structures)
ax.set_yticks(list(range(0,t_max//2+1)))
ax.set_yticklabels(freqs)
ax.tick_params(axis='x', labelsize=11)
ax.tick_params(axis='y', labelsize=11)
plt.imshow(np.abs(FS_coef), cmap='RdBu', vmin=-3, vmax=3, aspect='auto')
plt.savefig('freq_struct.png', dpi=300, bbox_inches='tight')
#plt.show()
