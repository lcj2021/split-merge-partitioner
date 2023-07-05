declare -A map=(
    ["amazon"]="com-amazon.ungraph.txt"
    ["arabic"]="arabic-2005.edge"
    ["hollywood"]="hollywood-2011.edge"
    ["indo"]="indochina-2004.edge"
    ["it"]="web-it-2004-all.mtx"
    ["lj"]="com-lj.ungraph.txt"
    ["orkut"]="com-orkut.ungraph.txt"
    ["sina"]="soc-sinaweibo.mtx"
    ["sk"]="web-sk-2005-all.mtx"
    ["twitter"]="twitter-2010.edge"
    ["uk"]="uk-2005.edge"
    ["webbase"]='webbase-2001.edge'
    ["wiki"]="web-wikipedia_link_en13-all.edges"
)

# arg1: shortname of data file

short=$1
filename=${map[$short]}

# arg2: p

partition_info_path=../partition_log
p=$2

cd ../release

# Partition SMP
for ((k = 2; k <= 4; ++k))
do
    logname=$filename"_smp_k_"$k".log"
    echo $logname
    nohup ./main -p $p -k $k -method smp -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
    pid=$!
    echo $pid
    wait
done

# Partition SMP-HEP
for ((k = 2; k <= 3; ++k))
do
    logname=$filename"_smp_hep_k_"$k".log"
    echo $logname
    nohup ./main -p $p -k $k -method smp_hep -hdf 100 -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
    pid=$!
    echo $pid
    wait
done

# Partition HEP
hdf_list=(1 10 100)
for ((i = 0; i < ${#hdf_list[@]}; ++i))
do
    hdf=${hdf_list[i]}
    logname=$filename"_hep_hdf_"$hdf".log"
    echo $logname
    nohup ./main -p $p -hdf $hdf -method hep -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
    pid=$!
    echo $pid
    wait
done

# Partition NE
logname=$filename"_ne.log"
echo $logname
nohup ./main -p $p -method ne -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
pid=$!
echo $pid
wait

# Partition EBV
logname=$filename"_ebv.log"
echo $logname
nohup ./main -p $p -method ebv -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
pid=$!
echo $pid
wait

# Partition HDRF
logname=$filename"_hdrf.log"
echo $logname
nohup ./main -p $p -method hdrf -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
pid=$!
echo $pid
wait

# Partition DBH
logname=$filename"_dbh.log"
echo $logname
nohup ./main -p $p -method dbh -filename "../../dataset/"$filename -write true > $partition_info_path/p_$p/$short/$logname &
pid=$!
echo $pid
wait

ls -lh ../../dataset/$filename.edgepart.*
wc -l ../../dataset/$filename.edgepart.*
rm -f ../../dataset/$filename.edgepart.*