# Tag_Propagation_Partition

## Build
```shell
mkdir release
cd release && cmake ../
make
```

## Run
```shell
./main -p 8 -filename ../dataset/mini.txt
./main -p 8 -filename ../dataset/com-amazon.ungraph.txt
./main -p 16 -filename ../dataset/com-lj.ungraph.txt
```



# Statistic
- The datasets are from [snap](http://snap.stanford.edu/data/index.html)


[LJ](http://snap.stanford.edu/data/com-LiveJournal.html)

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    |
| ---- | ---------- | ----------- | ------- |
| 8    | 3997962    | 7827080     | 35.0114 |
| 16   | 3997962    | 9115827     | 46.0436 |
| 32   | 3997962    | 9053842     | 47.9411 |
| 128  | 3997962    | 9877578     | 68.9839 |
| 256  | 3997962    | 9732271     | 82.4134 |
| 512  | 3997962    | 10944786    | 123.017 |



[Amazon](http://snap.stanford.edu/data/com-Amazon.html)

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    |
| ---- | ---------- | ----------- | ------- |
| 8    | 334863     | 557979      | 1.01654 |
| 16   | 334863     | 651342      | 1.10861 |
| 32   | 334863     | 709462      | 1.43612 |
| 128  | 334863     | 762469      | 1.71044 |
| 256  | 334863     | 781182      | 2.2125  |
| 512  | 334863     | 794393      | 2.93649 |



[Orkut](http://snap.stanford.edu/data/com-Orkut.html)

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    |
| ---- | ---------- | ----------- | ------- |
| 8    | 3072441    | 6507625     | 131.083 |
| 16   | 3072441    | 7768115     | 178.266 |
| 32   | 3072441    |             |         |
| 128  | 3072441    | 8741467     | 255.466 |
| 256  | 3072441    | 8726848     | 286.928 |
| 512  | 3072441    | 9174221     | 383.242 |


# Results Comparison on Ubuntu 18.04 (Single 128Gb Memory Machine)

## NE 

- The results of Neighborhood Expansion (NE) algorithm.
[Graph Edge Partitioning via Neighborhood Heuristic](http://www.kdd.org/kdd2017/papers/view/graph-edge-partitioning-via-neighborhood-heuristic)


[LJ](http://snap.stanford.edu/data/com-LiveJournal.html)

| p    | VERTEX_CNT | Time    | RF      |
| ---- | ---------- | ------- | ------- |
| 8    | 3997962    | 24.3755 | 1.31405 |
| 16   | 3997962    | 25.3646 | 1.43483 |
| 32   | 3997962    | 31.9216 | 1.56744 |
| 128  | 3997962    | 76.6040 | 1.89865 |
| 256  | 3997962    | 146.308 | 2.08690 |
| 512  | 3997962    | 280.972 | 2.31126 |


[Orkut](http://snap.stanford.edu/data/com-Orkut.html)

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    | RF      |
| ---- | ---------- | ----------- | ------- | ------- |
| 8    | 3072441    | 6507625     | 75.3679 | 1.74282 |
| 16   | 3072441    | 7768115     | 84.4289 | 2.05309
| 32   | 3072441    |             | 114.714 | 2.51501 |
| 128  | 3072441    | 8741467     | 263.342 | 3.75545 |
| 256  | 3072441    | 8726848     | 457.818 | 4.56794 | 
| 512  | 3072441    | 9174221     | 871.835 | 5.60136 |

## TP

TBD