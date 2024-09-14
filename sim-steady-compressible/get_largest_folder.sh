#!/bin/bash

# Initialize variables
largest_num=-1
largest_folder=""

# Loop through the directories
for dir in processor*/[0-9]*; do
    # Extract the second part (Y in processorX/Y)
    folder_name=$(basename "$dir")
    
    # Check if it's a number and larger than the current largest
    if [[ "$folder_name" =~ ^[0-9]+$ ]]; then
        if (( folder_name > largest_num )); then
            largest_num=$folder_name
            largest_folder=$dir
        fi
    fi
done

# Output the largest folder
if [[ -n "$largest_folder" ]]; then
    echo "$largest_num"
else
    echo "No valid folders found."
fi