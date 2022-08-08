#!/bin/bash

path="./workload"
files=$(ls ${path})
nvm=/mnt/pmem1/lgc/*

nvmfile="/mnt/pmem1/lbtree.pool"
cmdinit="../../lbtree/lbtree thread 1 mempool 128 nvmpool ${nvmfile} 1024"
task="taskset -c 0-17"
workloadSize=33554432
scalr=4
millon=`expr 1024 \* 1024`
preLoadSize=`expr $scalr \* $millon`
preLoadFile=dataset.dat
fillrate=0.7
#workloadFile=workload2/Mixed.txt
output="./output"
preRead=preRead.txt


rm $nvm
${task} ${cmdinit} bulkload 1 ../../lbtree/keygen-8B/1.txt $fillrate$ insert $preLoadSize$ dataset.dat  workload workload.txt 1

for i in 1 2 3 4
do
    ${task} ./main -c $i 
    # echo 3 > /proc/sys/vm/drop_caches
done

# for filename in $files
# do
#     rm $nvm
#     # echo 3 > /proc/sys/vm/drop_caches
#     ${task} ${cmdinit} bulkload 1 ../lbtree/keygen-8B/1.txt $fillrate$ insert $preLoadSize$ dataset.dat  workload $path/$filename 1 
    
#     # for i in {2, 3, 11}
#     # do
#     # ${task} ./main $i 32 $path/$filename >> ${output}/$filename
#     # done
#     # for i in {11..13}
#     # do
#     # ${task} ./main $i 32 $path/$filename >> ${output}/$filename
#     # done
#     # ${task} ./main 3 $path/$filename >> ${output}/$filename
#     # ${task} ./main 3 $path/$filename ${preRead} >> ${output}/$filename
# done


# path=./workload
# for i in {2..3}
# do
#     path=./workload${i}
#     files=$(ls ${path})
#     for filename in $files 
#     do
#         rm $nvm
#         # echo 3 > /proc/sys/vm/drop_caches
#         ${task} ${cmdinit} bulkload 1 ../lbtree/keygen-8B/1.txt $fillrate$ insert $preLoadSize$ dataset.dat  workload $path/$filename >> ${output}$i/$filename
#         for j in {1..2}
#         do
#         ${task} ./main $j $path/$filename >> ${output}$i/$filename
#         done
#         ${task} ./main 3 $path/$filename ${preRead} >> ${output}$i/$filename
#     done
# done



