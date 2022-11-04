cd release
make
./main -p   8 -filename ../../dataset/com-lj.ungraph.txt
./main -p  16 -filename ../../dataset/com-lj.ungraph.txt
./main -p  32 -filename ../../dataset/com-lj.ungraph.txt
./main -p  64 -filename ../../dataset/com-lj.ungraph.txt
./main -p 128 -filename ../../dataset/com-lj.ungraph.txt
./main -p 256 -filename ../../dataset/com-lj.ungraph.txt
./main -p 512 -filename ../../dataset/com-lj.ungraph.txt