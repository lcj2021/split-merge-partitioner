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

# arg2: exp round
# arg3: p
# arg4: $4 = 0 means no need to partition

partition_info_path=/home/lcj_graph_partitioning/partition_info
exp_res_path=/home/lcj_graph_partitioning/exp_res_$2
p=$3

# Tasks
cd /home/lcj_graph_partitioning/cloud/PowerGraph_lcj/release/toolkits/graph_analytics
data_path=/home/lcj_graph_partitioning/dataset
res_path=/home/lcj_graph_partitioning/exp_res_$2

# Tasks on SMP
for ((k = 2; k <= 4; ++k))
do
    part_name=$filename".edgepart.smp_k_"$k"."$p
    echo $part_name

    nohup mpiexec -n $p -f ~/machines ./pagerank --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name --iterations=100 > $res_path/p_$p/$short/pagerank/smp_k_$k &
    pid=$!
    echo $pid
    wait

    nohup mpiexec -n $p -f ~/machines ./connected_component --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/cc/smp_k_$k &
    pid=$!
    echo $pid
    wait

    # nohup mpiexec -n $p -f ~/machines ./approximate_diameter --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/diameter/smp_k_$k &
    # pid=$!
    # echo $pid
    # wait
done

# Tasks on HEP
hdf_list=(1 10 100)
for ((i = 0; i < ${#hdf_list[@]}; ++i))
do
    hdf=${hdf_list[i]}
    part_name=$filename".edgepart.hep_hdf_"$hdf"."$p
    echo $part_name
    
    nohup mpiexec -n $p -f ~/machines ./pagerank --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name --iterations=100 > $res_path/p_$p/$short/pagerank/hep_hdf_$hdf &
    pid=$!
    echo $pid
    wait

    nohup mpiexec -n $p -f ~/machines ./connected_component --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/cc/hep_hdf_$hdf &
    pid=$!
    echo $pid
    wait

    # nohup mpiexec -n $p -f ~/machines ./approximate_diameter --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/diameter/hep_hdf_$hdf &
    # pid=$!
    # echo $pid
    # wait
done

# Tasks on NE
part_name=$filename".edgepart.ne".$p
echo $part_name

nohup mpiexec -n $p -f ~/machines ./pagerank --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name --iterations=100 > $res_path/p_$p/$short/pagerank/ne &
pid=$!
echo $pid
wait

nohup mpiexec -n $p -f ~/machines ./connected_component --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/cc/ne &
pid=$!
echo $pid
wait

# nohup mpiexec -n $p -f ~/machines ./approximate_diameter --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/diameter/ne &
# pid=$!
# echo $pid
# wait

# Tasks on EBV
part_name=$filename".edgepart.ebv".$p
echo $part_name

nohup mpiexec -n $p -f ~/machines ./pagerank --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name --iterations=100 > $res_path/p_$p/$short/pagerank/ebv &
pid=$!
echo $pid
wait

nohup mpiexec -n $p -f ~/machines ./connected_component --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/cc/ebv &
pid=$!
echo $pid
wait

# nohup mpiexec -n $p -f ~/machines ./approximate_diameter --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/diameter/ebv &
# pid=$!
# echo $pid
# wait


# Tasks on HDRF
part_name=$filename".edgepart.hdrf".$p
echo $part_name

nohup mpiexec -n $p -f ~/machines ./pagerank --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name --iterations=100 > $res_path/p_$p/$short/pagerank/hdrf &
pid=$!
echo $pid
wait

nohup mpiexec -n $p -f ~/machines ./connected_component --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/cc/hdrf &
pid=$!
echo $pid
wait

# nohup mpiexec -n $p -f ~/machines ./approximate_diameter --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/diameter/hdrf &
# pid=$!
# echo $pid
# wait

# Tasks on DBH
part_name=$filename".edgepart.dbh".$p
echo $part_name

nohup mpiexec -n $p -f ~/machines ./pagerank --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name --iterations=100 > $res_path/p_$p/$short/pagerank/dbh &
pid=$!
echo $pid
wait

nohup mpiexec -n $p -f ~/machines ./connected_component --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/cc/dbh &
pid=$!
echo $pid
wait

# nohup mpiexec -n $p -f ~/machines ./approximate_diameter --format=snap_dist --graph_opts ingress=my --graph=$data_path/$part_name > $res_path/p_$p/$short/diameter/dbh &
# pid=$!
# echo $pid
# wait