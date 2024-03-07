
#!/bin/bash

hosts=(
    slave0 
    slave1
    slave2
    slave3
    slave4
    slave5
    slave6
    slave7
    # slave8
    # slave9
)

for host in "${hosts[@]}"
do
    echo "Shutting down $host"
    ssh $host "systemctl stop firewalld"
done


