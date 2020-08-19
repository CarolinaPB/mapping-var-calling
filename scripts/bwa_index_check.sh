#!/bin/bash
FILE=$1

if [ -f "$FILE.amb" ]; then
    echo "$FILE.amb exists."
else 
    echo "$FILE.amb does not exist."
    module load bwa
    bwa index $FILE
fi