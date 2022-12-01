import csv

def read_link_stream(edgelist_path, weighted=True):
	
	edgelist = dict()
	with open(edgelist_path, 'r') as f1:
		
		csv_reader = csv.reader(f1, delimiter=',')

		for line_edgelist in csv_reader:
				
			source = int(line_edgelist[0])
			destin = int(line_edgelist[1])
			time = int(line_edgelist[2])
	
			if weighted : 
				if time not in edgelist: edgelist[time] = []
				edgelist[time].append( (source, destin) )
			else :
				if time not in edgelist: edgelist[time] = set()
				edgelist[time].add( (source, destin) )

	return edgelist

def aggregate_link_stream( link_stream ):
	agg_graph = set()	
	for t in link_stream:
		for e in link_stream[t]:
			agg_graph.add( e )

	return agg_graph
