import math

# Get next power of 2
def next_power_of_2(x):
    return 1 if x == 0 else 2**math.ceil(math.log2(x))

# Get element from set
def get_set_elem( s ):
	for e in s: break
	return e

# Adjacency list from edgelist
def adjlist_from_edgelist( edgelist ):
	adj_list = dict()
	for e in edgelist:
		if e[0] not in adj_list.keys(): adj_list[e[0]] = []
		adj_list[e[0]].append( e[1] )
	return adj_list

# Max number of levels
def max_number_bipartitions(len_set):
	return 0 if len_set == 0 else math.ceil(math.log2(len_set)) 

# Number of partitioned sets at a level
def num_partitions_at_level(len_set, level):
	full_set_size = next_power_of_2(len_set)
	partition_set_size = 2**level
	return full_set_size // partition_set_size

# Test if an integer is a power of two
def is_power_of_two(n):
    return (n != 0) and (n & (n-1) == 0)

# Get the largest power 2 before a given number
def previous_power_2(n):
    res = 0;
    for i in range(n, 0, -1):
        if ((i & (i - 1)) == 0):
            res = i;
            break;
    return res;

# Get the MCD to find sampling rate
def find_gcd(x, y):
    while(y):
        x, y = y, x % y
    return x
