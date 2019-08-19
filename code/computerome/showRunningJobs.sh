#!/bin/bash
IDs=($(qstat -r | grep -oP '.*?(?=\.risoe)' | tr '\n' ' '))
numberOfJobs=${#IDs[@]}

for (( i=0; i<${numberOfJobs}; i++ ));
do 
  echo ${IDs[$i]}
  checkjob ${IDs[$i] }
done 
