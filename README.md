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

LJ

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    |
| ---- | ---------- | ----------- | ------- |
| 8    |            |             |         |
| 16   | 3997962    | 9115827     | 46.0436 |
| 32   |            |             |         |
| 128  | 3997962    | 9877578     | 68.9839 |
| 256  | 3997962    | 9732271     | 82.4134 |

Amazon

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    |
| ---- | ---------- | ----------- | ------- |
| 8    | 334863     | 557979      | 1.01654 |
| 16   | 334863     | 651342      | 1.10861 |
| 32   | 334863     | 709462      | 1.43612 |
| 128  | 334863     | 762469      | 1.71044 |

Orkut

| p    | VERTEX_CNT | ALL_TAG_CNT | Time    |
| ---- | ---------- | ----------- | ------- |
| 8    |            |             |         |
| 16   | 3072441    | 7768115     | 178.266 |
| 32   | 3072441    |             |         |
| 128  | 3072441    | 8741467     | 255.466 |
| 256  | 3072441    | 8726848     | 286.928 |

