import sys
sys.path.append('../../Library')
import decomposition_utils
import data_processing
import numpy as np
import matplotlib.pyplot as plt

#################
## FREQ FILTER ## 
#################
h = np.array([1, 1, 0, 0, 0, 0])
lh = np.fft.rfft(h, norm='ortho')
freqs = np.around(np.fft.rfftfreq(6), 2)
freq_str = [str(f) + ' ($\chi_' + str(i) + '$)' for i,f in enumerate(freqs)]
print(freq_str)

fig, ax = plt.subplots(figsize=(5.5,1.2))
plt.stem(freqs, np.abs(lh))
ax.set_xticks(freqs)
ax.set_xticklabels(freq_str)
ax.set_xlabel('Frequency', fontsize=10)
ax.set_ylabel('Magnitude', fontsize=10)
ax.set_title('Aggregation (filter response)', fontsize=11)
ax.tick_params(axis='x', labelsize=9)
ax.tick_params(axis='y', labelsize=9)
ax.set_ylim([-0.1,1.1])
plt.savefig('freq_filt_response.png', dpi=300, bbox_inches='tight')

###################
## STRUCT FILTER ## 
###################
structures = ['$\sigma_0^{(3)}$', '$\sigma_{1}^{(3)}$', r'$\nu_{0}^{(3)}$', r'$\nu_{1}^{(3)}$', r'$\nu_{0}^{(2)}$', r'$\nu_{1}^{(2)}$',
				r'$\nu_{2}^{(2)}$', r'$\nu_{3}^{(2)}$', r'$\nu_{0}^{(1)}$', r'$\nu_{1}^{(1)}$', r'$\nu_{2}^{(1)}$', r'$\nu_{3}^{(1)}$',
				r'$\nu_{4}^{(1)}$', r'$\nu_{5}^{(1)}$', r'$\nu_{6}^{(1)}$', r'$\nu_{7}^{(1)}$']


lq = np.zeros(16)
lq[[0,1]] = 1

fig, ax = plt.subplots(figsize=(5.5,1.2))
plt.stem(lq)
ax.set_xticks(range(16))
ax.set_xticklabels(structures)
ax.set_xlabel('Structure', fontsize=10)
ax.set_ylabel('Magnitude', fontsize=10)
ax.set_title('Embedding (filter response)', fontsize=11)
ax.tick_params(axis='x', labelsize=9)
ax.tick_params(axis='y', labelsize=9)
ax.set_ylim([-0.1,1.1])
plt.savefig('struct_filt_response.png', dpi=300, bbox_inches='tight')
plt.show()
