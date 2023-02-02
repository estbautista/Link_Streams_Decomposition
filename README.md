# Link_Streams_Decomposition
Repository for the preprint : A Joint Frequency-Structure Decomposition for Link Streams

The repository is structured in two folders:
  - Experiments: <br/>
    &ensp; Contains the code to reproduce the experiments and figures shown in the article. <br/>
    &ensp; Each experiment contain a main.py file that needs to be executed in order to reproduce the results. 
    
      
  - Library: <br/>
    &ensp; Implements the proposed graph decomposition via Filter-Banks.
    &ensp; The implementation is dependant on the Haar Wavelet method of the PyWavelets library. 
    &ensp; The decomposition_utils.py file provides an interface to the graph decomposition with a construction of basis elements via SVD or BFS-based procedures. It also provides a methods (graph_sequence_decomposition) that is optimized to transform a graph sequence. 
    &ensp; The files decomposition_BFS.py and decomposition_SVD.py contain all the implementation details of the graph decomposition (partition of the edge-space, mapping of the graph edge-space to a time series, forward and inverse transformations).
    
    
