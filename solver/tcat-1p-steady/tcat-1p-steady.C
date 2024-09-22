/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    buoyantSimpleFoam

Group
    grpHeatTransferSolvers

Description
    Steady-state solver for buoyant, turbulent flow of compressible fluids,
    including radiation, for ventilation and heat-transfer.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
// #include "rhoThermo.H"
// #include "turbulentFluidThermoModel.H"
// #include "radiationModel.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"
#include "TCAT_compressible.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Steady-state solver for buoyant compressible fluid flow"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    int count = 0;

    // MacroscaleCompressible tcat_global(false,file_out,mesh,count,mu);
    // MacroscaleCompressible tcat_proc(true,file_out,mesh,count,mu);

    // double mom_residual = 1.312;
    // double _residual = 1.e-10;

    // while (mom_residual > _residual)
    while (simple.loop())
    {

        // Take over some of simple's job
        // p_rgh.storePrevIter();
        // U.storePrevIter();
        // rho.storePrevIter();


        Info<< "Time = " << runTime.timeName() << nl << endl;
        // Info << "Tim! " << mom_residual << '\t' << _residual << endl;

        // Pressure-velocity SIMPLE corrector
        {
            #include "UEqn.H"
            #include "pEqn.H"
        }

        runTime.write();

        runTime.printExecutionTime(Info);

        // mom_residual = tcat_global.residual(U,p,p_rgh,rho,phi,g,gh,ghf);

        // if (mom_residual < _residual)
        // {
        //     tcat_global.update(count,runTime.value(),U,p,p_rgh,rho,phi,chem_potential,grav_potential,g,gh,ghf);
        //     tcat_proc.update(count,runTime.value(),U,p,p_rgh,rho,phi,chem_potential,grav_potential,g,gh,ghf);
        //     count = count + 1;

        // }
    }

    return 0;
}


// ************************************************************************* //
