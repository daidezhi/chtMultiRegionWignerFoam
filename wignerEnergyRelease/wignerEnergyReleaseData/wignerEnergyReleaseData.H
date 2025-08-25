
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2025 Dezhi Dai, Argonne National Laboratory (ANL)
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

Class
    Foam::functionObjects::wignerEnergyReleaseData

Description
    Function object for post-processing Wigner energy release and solid-fluid
    heat transfer in conjugate heat transfer cases.

    The object evaluates:
      - minimum, maximum and average solid temperature
      - Wigner energy release power [W]
      - accumulated released Wigner energy [J]
      - solid-to-fluid heat transfer power [W]
      - accumulated transferred energy from solid to fluid [J]

    Output is written to the postProcessing directory at run-time according
    to the selected write controls. An example is given below:

    wignerEnergyReleaseData1
    {
        type            wignerEnergyReleaseData;
        libs            (fieldFunctionObjects wignerEnergyRelease);
        region          fuelZone;
        patches         (fuelToHelium);

        writeFields     false;

        // executeControl timeStep = 1 to ensure accurately transferred Wigner
        // energy calculation.
        executeControl  timeStep;
        executeInterval 1;
        writeControl    adjustableRunTime;
        writeInterval   1;
    }

SourceFiles
    wignerEnergyReleaseData.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_wignerEnergyReleaseData_H
#define functionObjects_wignerEnergyReleaseData_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "HashSet.H"
#include "volFields.H"
#include "surfaceFields.H"

#include "wignerEnergyRelease.H"

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                Class wignerEnergyReleaseData Declaration
\*---------------------------------------------------------------------------*/

class wignerEnergyReleaseData
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // Protected data

        //- Total volume
        scalar vol_;

        //- Min T
        scalar TMin_;

        //- Max T
        scalar TMax_;

        //- Average T
        scalar TAvg_;

        //- Wigner energy release power [W]
        scalar WignerEnergyReleasePower_;

        //- Released Wigner energy [J]
        scalar releasedWignerEnergy_;

        //- Solid-to-fluid heat transfer power [W]
        scalar energyTransferPower_;

        //- Transferred energy from solid to fluid [J]
        scalar transferredEnergy_;


        // Data for calculating transferred energy from solid to fluid

            //- Patch name(s) of solid-to-fluid interface(s) (wall(s))
            wordRes patchNames_;

            //- Optional list of wall patches to process
            labelHashSet patchSet_;

            //- Total patch area
            scalar area_;

            //- Wall conduction heat transfer rate
            scalar QDotc_;


    // Protected Member Functions

        //- File header information
        virtual void writeFileHeader(Ostream& os) const;

        //- Calculate the Wigner energy data
        void calcWignerEnergyData(const volScalarField& rho);

        //- Calculate the heat transfer rate (power)
        void calcHeatPower
        (
            const volScalarField& alpha,
            const volScalarField& he
        );


public:

    //- Runtime type information
    TypeName("wignerEnergyReleaseData");


    // Constructors

        //- Construct from Time and dictionary
        wignerEnergyReleaseData
        (
            const word& name,
            const Time& runTime,
            const dictionary&
        );

        //- No copy construct
        wignerEnergyReleaseData(const wignerEnergyReleaseData&) = delete;

        //- No copy assignment
        void operator=(const wignerEnergyReleaseData&) = delete;


    //- Destructor
    virtual ~wignerEnergyReleaseData() = default;


    // Member Functions

        //- Read the wallHeatPower data
        virtual bool read(const dictionary& dict);

        //- Calculate the wall heat-flux
        virtual bool execute();

        //- Write the wall heat-flux
        virtual bool write();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //