#!/bin/bash

cd test
echo "--parse--test--"

ulimit -s unlimited
ulimit -f unlimited


fn="mytest"
echo "realy prcoess $fn.txt"

# Must add tag <s>
sed 's/^/<s> /g;s/$/ <\/s>/g' $fn.txt > $fn.tmp

# This is only for parseStdin (must clean-text first)
sed -f clean-text.sed $fn.txt | sed 's/^/<s> /g;s/$/ <\/s>/g' > $fn.cut 

echo "realy prcoess $fn.1.PAR"
../parseIt -l399 -d5 "../DATA/" $fn.tmp > $fn.1.PAR  2> $fn.1.log

echo "realy prcoess $fn.2.PAR"
../parseStdin -l399 -d5 "../DATA/" < $fn.cut > $fn.2.PAR 2> $fn.2.log


echo "--DONE--"
cat $fn.txt

echo "[file]---"
cat $fn.1.PAR

echo "[pipe]---"
cat $fn.2.PAR

echo 
echo
cd ..
