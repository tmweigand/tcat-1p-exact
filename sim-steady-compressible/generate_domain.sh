#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

canCompile || exit 0    # Dynamic code

restore0Dir

mkdir tcat

touch tcat-1p-steady.foam

runApplication blockMesh

runApplication decomposePar

runParallel snappyHexMesh -overwrite

runApplication reconstructParMesh -constant

#rm -rf processor*

runApplication renumberMesh -overwrite

# runApplication transformPoints -scale 1.e-5
#------------------------------------------------------------------------------
