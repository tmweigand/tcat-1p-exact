#!/bin/sh
cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------

canCompile || exit 0    # Dynamic code

restore0Dir

touch tcat-1p-steady.foam

rm tcat/*

rm -rf processor*

rm log.decomposePar

runApplication decomposePar

rm log.tcat-1p-steady

runParallel tcat-1p-steady  

#------------------------------------------------------------------------------
