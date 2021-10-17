#!/bin/bash

find ./ -name "C1.dat" -exec rm {} \;
for ((i=38;i<82;i++));do
  ./main $i 1
done
