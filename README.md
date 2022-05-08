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