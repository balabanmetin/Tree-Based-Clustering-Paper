#! /bin/bash

export PASTA_DEBUG=TRUE
export PASTA_LOGGING_LEVEL=DEBUG
export PASTA_LOGGING_FORMAT=RICH

config=$1
input=$2
temp=`mktemp -d`

run_pasta.py $config -i $input --temporaries $temp 
