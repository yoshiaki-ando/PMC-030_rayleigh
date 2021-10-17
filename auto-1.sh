#!/bin/bash

find ./ -name "C-1.dat" -exec rm {} \;
for ((i=60;i<82;i++));do
  ./main $i -1
done
