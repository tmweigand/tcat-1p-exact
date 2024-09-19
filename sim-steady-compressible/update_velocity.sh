#!/bin/bash

num_procs="$1"
in_vel="$2"
out_vel="$3"
folder="$4"

echo $folder

old_line="massFlowRate    constant ${in_vel};"
new_line="massFlowRate    constant ${out_vel};" 

for i in `seq 0 $num_procs`;
do
        file="processor$i/$folder/U"

	if [[ "$OSTYPE" == "darwin"* ]]; then
    		# macOS
    		sed -i '' "/$old_line/c\\
    		$new_line
    		" "$file"
	else
    		# Linux
    		sed -i "/$old_line/c\\
    		$new_line
    		" "$file"
	fi

done


file="constant/transportProperties"


## As number
old_line="v_in v_in \[0 1 -1 0 0 0 0\] ${in_vel};"
new_line="v_in v_in [0 1 -1 0 0 0 0] ${out_vel};"
if [[ "$OSTYPE" == "darwin"* ]]; then
	sed -i '' "s|$old_line|$new_line|" "$file"
else
        sed -i "s|$old_line|$new_line|" "$file"
fi
## Output File - as string
old_line="file_out tcat_velocity_${in_vel};"
new_line="file_out tcat_velocity_${out_vel};"
if [[ "$OSTYPE" == "darwin"* ]]; then
	sed -i '' "s|$old_line|$new_line|" "$file"
else
	sed -i "s|$old_line|$new_line|" "$file"
fi
