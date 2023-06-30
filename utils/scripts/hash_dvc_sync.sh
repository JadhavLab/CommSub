#!/bin/bash
# This script links the matfiles from the hash directory to the code directory
# and adds them to dvc if they are not already there

# Specify MATLAB command to be run
MATLAB_CMD="hfolder=hashdefine();disp(hfolder);"

# Run MATLAB with specified command
folder=/Volumes/Ark/commsubspace/hash/
codefolder=/Volumes/MATLAB-Drive/Shared/hash/

echo "Linking files from $folder to /Volumes/MATLAB-Drive/Shared/hash/"
echo "$folder/*.mat"
echo "ln -sf $folder/*.mat $codefolder"
read -p "Press enter to continue"
ln -sf $folder/*.mat $codefolder

for file in /Volumes/MATLAB-Drive/Shared/hash/*.mat; do
	if [ ! -f "${file}.dvc" ]; then
	    echo "Adding $file to dvc"
	    dvc add ${file}
    	fi
done
