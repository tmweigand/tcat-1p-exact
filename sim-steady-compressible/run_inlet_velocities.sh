#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

canCompile || exit 0    # Dynamic code

rm log.tcat-1p-steady*

rm -rf processor*

rm -rf 0/

rm log.decomposePar

restore0Dir

runApplication decomposePar

rm -rf tcat/

mkdir tcat/

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

echo $num_procs

velocity=(0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1 5 10 50 100 0.001)
length=${#velocity[@]}
i=1
while [ $i -ne $length ]
do
    # Run code 
    runParallel tcat-1p-steady

    # Update velocity
    largest_folder=$(./get_largest_folder.sh)
    ./update_velocity.sh $num_procs ${velocity[i-1]} ${velocity[i]} $largest_folder
    
    mv log.tcat-1p-steady log.tcat-1p-steady${velocity[i-1]}
    
    i=$((i + 1))

done
