#!/bin/bash
FILE=$1

if [ -f "$FILE.amb" ]; then
    echo "$FILE.amb exists." > $2
else 
    echo "$FILE.amb does not exist." > $2
    module load bwa
    bwa index $FILE
fi

# added so snakemake doesn't complain that the file is not there when it is
sleep 10