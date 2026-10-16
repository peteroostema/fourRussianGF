#!/bin/bash

#for (( i = 6; i < 300; i = i + 33));
for i in 4 6 8 9 10 16 128 256
do
for (( j = 1; j < 300; j = j + 33));
do
for (( k = 1; k < 300; k = k + 33));
do
echo "./a.out $i $j $k"
./a.out $i $j $k
done
done
done