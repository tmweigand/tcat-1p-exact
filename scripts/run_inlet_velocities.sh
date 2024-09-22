#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

canCompile || exit 0    # Dynamic code

rm log.tcat-1p-steady*

restore0Dir

rm log.decomposePar

rm -rf processor*

runApplication decomposePar

rm -rf tcat/

mkdir tcat/

mkdir tcat/local

# Get number of procs
base_dir="."
num_procs=-1
# Loop through all directories that match the pattern 'processor*'
for dir in "${base_dir}/processor"*; do
    # Extract the number from the directory name
    num=$(basename "$dir" | grep -o '[0-9]\+')
    
    # Check if the number is greater than the current max_num
    if [[ $num -gt $num_procs ]]; then
        num_procs=$num
    fi
done
num_procs=$((num_procs + 1))


echo $num_procs

velocity=(1e-07 1e-06 5e-06 1e-05 5e-05 0.0001 0.0025 0.0005 0.00075 0.001 0.0025 0.005 0.01 0.025 0.05 0.075 0.1)
length=${#velocity[@]}
i=1
while [ $i -ne $length ]
do
    # Run code 
    # mpirun -np $num_procs tcat-1p-steady -parallel > log.tcat-1p-steady

    # Update velocity
    ./update_velocity.py -n $num_procs -in ${velocity[i-1]} -out ${velocity[i]} -f 0
    
    mv log.tcat-1p-steady log.tcat-1p-steady${velocity[i-1]}

    echo $i
    
    i=$((i + 1))

done