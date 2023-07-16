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

