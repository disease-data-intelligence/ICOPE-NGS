#! /usr/bin/env bash
(qstat | tail -n +3 | sed 's/ \+/\t/g' | cut -f1,5) > temp 
while read line; do 
  status=$(echo $line | awk '{print $2}')
  if [ "$status" == "Q" ]
  then 
    jobid=$(echo $line | awk '{print $1}')
    showstart $jobid; 
    fi 
  done < temp
rm temp 
