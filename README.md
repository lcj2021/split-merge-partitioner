FSM: A Fine-grained Splitting and Merging Framework for Dual-balanced Graph Partition
=============================================



## Compile

cmake, glog, gflags, boost are required:

`sudo apt-get install libgoogle-glog-dev libgflags-dev libboost-all-dev`

```sh
mkdir build 
cd build
cmake ..
make -j 4
```

## Run

### Prioritize Edge Balance

Run with FSM-H (HEP-100) k=2, p=32

```shell
./main -p 32 -k 2 -method fsm_hep -hdf 100 -filename ../../dataset/...
```

Run with FSM-N (NE) k=2, p=32

```shell
./main -p 32 -k 2 -method fsm_ne -filename ../../dataset/...
```

### Prioritize Vertex Balance

```shell
# -write true must be set
./main -p 8 -method fennel -filename ../../dataset/hollywood-2011.txt -write true

./main -p 8 -method e2a -filename ../../dataset/hollywood-2011.txt

./main -p 8 -method v2e_fennel -filename ../../dataset/hollywood-2011.txt
```



If you want to write the partition, please add `-write true` to the command.

If you want to merge with *Fast Merge*, please add `-fastmerge true` to the command.

## Dataset

### SNAP and Networkrepository

[SNAP](http://snap.stanford.edu/data/index.html)
[Networkrepository](https://networkrepository.com/networks.php)

### Webgraph
We prepare a sript for Webgraph environment setup. Run `setup_webgraph.sh` in `dataset` folder. 
Then you can download `{graphname}.graph` and `{graphname}.properties` from [webgraph](https://law.di.unimi.it/datasets.php). 

Take indochina graph as an example:

```shell
cd dataset/webgraph

wget http://data.law.di.unimi.it/webdata/indochina-2004/indochina-2004.graph
wget http://data.law.di.unimi.it/webdata/indochina-2004/indochina-2004.properties
```

Use `java -cp "webgraph_dep/*" it.unimi.dsi.webgraph.ArcListASCIIGraph indochina-2004 indochina-2004.txt` to get final graph. 

## FAQ

Q1: Why it fails on `./main -p 16 -method fsm_hep -hdf 100 -k 4 -filename ../dataset/snap/email-Enron.txt` with "hep_partitioner.hpp:76] Check failed: edge2bucket[edge_id] == -1 (17 vs. -1) " ? 
> In the *Split* Phase, HEP sample vertices in a sequential way. 
> Instead, NE sample randomly so that it might partition smoothly. 

Q2: Why it takes so long for `-method ne` to partition a large graph under `-p 512`? 
> When NE is close to finish, it needs to pull vertices into boundary in a random way much more frequently. However, it might take a long time to "lottery" a valid vertex. 
> HEP can avoid this phenomenon in a certain degree. 
