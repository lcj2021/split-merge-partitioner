Edge Partitioning Algorithms for Large Graphs
=============================================

Acknowledgements
Parts of the implementation are based on the reference implementation provided by Qin Liu(NE): https://github.com/ansrlab/edgepart and Mayerrn(HEP): https://github.com/mayerrn/hybrid_edge_partitioner

This refers to the following classes: conversions.cpp/hpp, dense_bitset.hpp, edgepart.hpp, graph.cpp/hpp, min_heap.hpp, util.cpp/hpp

We added vertex & edge balance statistics to them.
1. NE: ne_partitioner.cpp/hpp dbh_partitioner.cpp/hpp 
2. HEP: hep_graph.cpp/hpp, hep_min_heap.hpp, hep_partitioner.cpp/hpp, hdrf_partitioner.cpp/hpp
