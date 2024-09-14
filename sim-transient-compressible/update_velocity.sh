#!/bin/bash

num_procs="$1"

for i in `seq 0 $num_procs`;
do
        #cp 0.orig/* processor$i/0/
        cd processor$i/0/
        #perl -pi -e 's/massFlowRate    constant 9.971e-09;/massFlowRate    constant 0.9971e-4;/g' U
        perl -pi -e 's/massFlowRate    constant 0.0001;/massFlowRate    constant 0.01;/g' U  #For 0.4 m_A
        cd ../../
done

