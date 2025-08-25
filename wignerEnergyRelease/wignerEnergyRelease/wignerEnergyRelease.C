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

\*---------------------------------------------------------------------------*/

#include "wignerEnergyRelease.H"

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::Wigner::wignerEnergyRelease::typeName
    = "wignerEnergyRelease";


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::Wigner::wignerEnergyRelease::wignerEnergyRelease
(
    const Foam::volScalarField& T,
    const Foam::fvMesh& mesh
)
:
    // External data references
    mesh_(mesh),
    T_(T),
    TIn_(T.internalField()),

    // Wigner energy release rate calculator
    releaseRateCalculator_(mesh_),

    // Wigner energy fields
    S_(nullptr),
    SDot_(nullptr),

    // Release Wigner energy or not
    isReleasingEnergy_(false),

    // Step update guard
    updated_(false),
    lastTimeIndex_(-1)
{
    Info<< "\n========================== Wigner Energy Release ==========================="
        << endl;

    if (releaseRateCalculator_.hasDataset())
    {
        isReleasingEnergy_ = true;

        // Allocate S_
        S_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "S",
                    mesh_.time().timeName()/mesh_.name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimEnergy/dimMass, scalar(0.0)),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        // Allocate SDot_
        SDot_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "SDot",
                    mesh_.time().timeName()/mesh_.name(),
                    mesh_,
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimEnergy/dimMass/dimTime, scalar(0.0)),
                zeroGradientFvPatchScalarField::typeName
            )
        );

        // Print information
        Info<< "Region name    : " << mesh_.name() << endl;

        Info<< "Total volume   : " << gSum(mesh_.V()) << " [m^3]" << endl;

        Info<< "Stored energy  : "
            << releaseRateCalculator_.dataset().last().first()
            << " [J/kg]" << nl << endl;

        Info<< "Current time   : " << mesh_.time().timeName() << " [s]" << endl;

        Info<< "Released energy: min/max/avg = "
            << min(S_()).value() << "/"
            << max(S_()).value() << "/"
            << S_->weightedAverage(mesh_.V()).value()
            << " [J/kg]" << endl;

        Info<< "Release rate   : min/max/avg = "
            << min(SDot_()).value() << "/"
            << max(SDot_()).value() << "/"
            << SDot_->weightedAverage(mesh_.V()).value()
            << " [J/kg/s]" << endl;
    }
    else
    {
        Info<< "No Wigner PLWC discretized dataset found for region "
            << mesh_.name()
            << endl;

        // Allocate SDot_ as a zero field
        SDot_.reset
        (
            new volScalarField
            (
                IOobject
                (
                    "SDot",
                    mesh_.time().timeName()/mesh_.name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh_,
                dimensionedScalar(dimEnergy/dimMass/dimTime, scalar(0.0)),
                zeroGradientFvPatchScalarField::typeName
            )
        );
    }

    Info<< "============================================================================"
        << nl << endl;
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

void Foam::Wigner::wignerEnergyRelease::update()
{
    if (isReleasingEnergy_)
    {
        const label currentTimeIndex(mesh_.time().timeIndex());

        if (currentTimeIndex != lastTimeIndex_)
        {
            updated_ = false;
            lastTimeIndex_ = currentTimeIndex;
        }

        if (updated_) return;

        Info<< "\nUpdating Wigner energy release rates and released energy ..."
            << endl;

        scalarField& SIn(S_->ref());
        scalarField& SDotIn(SDot_->ref());

        forAll(TIn_, cellI)
        {
            SDotIn[cellI]
          = 
            releaseRateCalculator_.releaseRate(TIn_[cellI], SIn[cellI]);
        }

        SIn += (SDotIn * mesh_.time().deltaTValue());

        S_->correctBoundaryConditions();
        SDot_->correctBoundaryConditions();

        Info<< "Min/max/avg Wigner energy released: "
            << min(S_()).value() << "/"
            << max(S_()).value() << "/"
            << S_->weightedAverage(mesh_.V()).value()
            << " J/kg" << endl;

        Info<< "Min/max/avg Wigner energy release rate: "
            << min(SDot_()).value() << "/"
            << max(SDot_()).value() << "/"
            << SDot_->weightedAverage(mesh_.V()).value()
            << " J/kg/s" << nl << endl;

        updated_ = true;
    }
    else
    {
        return;
    }
}

// ************************************************************************* //