#!/bin/bash

path=workload8

gen="./datagen"

for i in {1..16}
do 
    rm dataset.dat
    rm /mnt/pmem/*tree.pool
    ${gen} -l ${i} 
    ./main 2 ${i} >> static.txt
done


# # ${gen} -o 32000 -n preRead.txt

# for i in {1..4} 
# do 
#     ${gen} -z -o ${i} -r 0 -i 1.0 -n ${path}3/insert${i}
#     ${gen} -z -o ${i} -r 0 -u 1.0 -n ${path}3/update${i}
#     ${gen} -z -o ${i} -r 0 -d 1.0 -n ${path}3/delete${i}
#     ${gen}  -o ${i} -r 0 -i 1.0 -n ${path}2/insert${i}
#     ${gen}  -o ${i} -r 0 -u 1.0 -n ${path}2/update${i}
#     ${gen}  -o ${i} -r 0 -d 1.0 -n ${path}2/delete${i}
# done




# gen="./datagen"
# path=./workload5

# for i in {5..9}
# do 
#     ${gen} -p ${i} -n ${path}/hotspot${i}
# done
# # ${gen} -n ${path}/random
# # ${gen} -z -s 0.95 -n ${path}/95
# # ${gen} -z -s 0.99 -n ${path}/99



#gen="./datagen -z 0.85 -o 32"




# ${gen} -n workload4/RU.txt -r 0.5 -u 0.5
# ${gen} -n workload4/Ru.txt -r 0.95 -u 0.05

# echo -n "datagen ReadOnly"
# ${gen} -n workload4/R.txt 


# ${gen} -n workload4/Ri.txt -r 0.95 -i 0.05
# echo "datagen ReadInsert"
# ${gen} -n workload4/RI.txt -r 0.5 -i 0.5


# echo "datagen Mixed"
# ${gen} -n workload4/Mixed.txt -r 0.6 -i 0.2 -u 0.2


# echo "datagen hotspot8";
#${gen} -n workload4/hotspot8.txt -p 8

# echo "datagen random";
# ./datagen -o 64 -n workload2/random.txt 

# cp -a workload2 ../lbtree




