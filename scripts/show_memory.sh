
#!/bin/bash

hosts=(slave0 
        slave1
        slave2
        slave3
        slave4
        slave5
        slave6
        slave7)

for host in "${hosts[@]}"
do
    echo "Memory of $host"
    ssh $host "free -g"
done


