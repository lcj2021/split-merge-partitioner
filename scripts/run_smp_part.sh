#!/bin/bash
file_list=(
    "arabic-2005.edge" 
    "com-orkut.ungraph.txt"
    "com-lj.ungraph.txt"
    "flickr.edges"
    "hollywood-2011.edge"
    "indochina-2004.edge"
    "soc-sinaweibo.mtx"
    "twitter-2010.edge"
    "uk-2005.edge"
    "webbase-2001.edge"
    "web-it-2004-all.mtx"
    "web-wikipedia_link_en13-all.edges"
            )

file_list_len=${#file_list[@]}
k=2
echo "file_list_len =" $file_list_len
echo "k =" $k
cd release
for ((i = 0; i < file_list_len; ++i))
do
    filename=${file_list[i]}
    echo $filename"_k"$k".log"
    ./main -p 8 -k $k  -method smp -filename "../../dataset/"$filename > $filename"_k"$k".log"
    pid=$!
    echo $pid
    wait
done
