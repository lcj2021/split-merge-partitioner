FSM: A Fine-grained Splitting and Merging Framework for Dual-balanced Graph Partition
=============================================



## Compile

```sh
mkdir release 
cd release
cmake ..
make -j 4
```



## Run

Run with FSM-H (HEP-100) k=2, p=32

```shell
./main -p 32 -k 2 -method smp_hep -hdf 100 -filename ../../dataset/...
```

Run with FSM-N (NE) k=2, p=32

```shell
./main -p 32 -k 2 -method smp_ne -filename ../../dataset/...
```



If you want to write the partition, please add `-write true` to the command.

If you want to merge with *Fast Merge*, please add `-fastmerge true` to the command.



## Acknowledgements

Parts of the implementation are based on the reference implementation provided by Qin Liu(NE): https://github.com/ansrlab/edgepart and Mayerrn(HEP): https://github.com/mayerrn/hybrid_edge_partitioner

This refers to the following classes: conversions.cpp/hpp, dense_bitset.hpp, edgepart.hpp, graph.cpp/hpp, min_heap.hpp, util.cpp/hpp

We added vertex & edge balance statistics to them.
1. NE: ne_partitioner.cpp/hpp dbh_partitioner.cpp/hpp 
2. HEP: hep_graph.cpp/hpp, hep_min_heap.hpp, hep_partitioner.cpp/hpp, hdrf_partitioner.cpp/hpp

