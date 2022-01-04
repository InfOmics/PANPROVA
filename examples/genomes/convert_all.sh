#!/bin/bash

for f in `ls *.gb`
do
b=`echo $f | sed s/\.gb/\.peg/g`
echo $f $b
python3 ../../gb2peg.py $f $b
done