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

#include "wignerEnergyReleaseData.H"

#include "solidThermo.H"
#include "wallPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(wignerEnergyReleaseData, 0);
    addToRunTimeSelectionTable(functionObject, wignerEnergyReleaseData, dictionary);
}
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::functionObjects::wignerEnergyReleaseData::writeFileHeader
(
    Ostream& os
) const
{
    writeCommented(os, "Region: ");
    os << obr_.name() << nl;

    writeCommented(os, "Total solid volume: ");
    os << vol_ << " m^3" << nl;

    writeCommented(os, "Total solid-to-fluid interface area: ");
    os << area_ << " m^2" << nl;

    writeCommented(os, "Time");
    writeTabbed(os, "TMin [K]");
    writeTabbed(os, "TMax [K]");
    writeTabbed(os, "TAvg [K]");
    writeTabbed(os, "QDot_Wigner [W]");
    writeTabbed(os, "QDot_inter [W]");
    writeTabbed(os, "Q_Wigner [J]");
    writeTabbed(os, "Q_inter [J]");

    os << endl;
}


void Foam::functionObjects::wignerEnergyReleaseData::calcWignerEnergyData
(
    const volScalarField& rho
)
{
    const volScalarField& T(mesh_.lookupObject<volScalarField>("T"));
    const volScalarField& SDot(mesh_.lookupObject<volScalarField>("SDot"));
    const volScalarField& S(mesh_.lookupObject<volScalarField>("S"));

    const scalarField& TIn(T.internalField());
    const scalarField& SDotIn(SDot.internalField());
    const scalarField& SIn(S.internalField());
    
    const scalarField& rhoIn(rho.internalField());
    const scalarField& volIn(mesh_.V().field());
    const scalarField  massIn(rhoIn * volIn);


    TMin_ = gMin(TIn);
    TMax_ = gMax(TIn);
    TAvg_ = gSum(TIn * volIn) / vol_;

    WignerEnergyReleasePower_ = gSum(SDotIn * massIn);
    releasedWignerEnergy_ = gSum(SIn * massIn);

    energyTransferPower_ = QDotc_;
    transferredEnergy_ += (QDotc_ * mesh_.time().deltaTValue());
}


void Foam::functionObjects::wignerEnergyReleaseData::calcHeatPower
(
    const volScalarField& alpha,
    const volScalarField& he
)
{
    QDotc_ = 0.0;

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();
    const volScalarField::Boundary& heBf = he.boundaryField();
    const volScalarField::Boundary& alphaBf = alpha.boundaryField();

    for (const label patchi : patchSet_)
    {
        QDotc_ += gSum(magSf[patchi] * (alphaBf[patchi] * heBf[patchi].snGrad()));
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::wignerEnergyReleaseData::wignerEnergyReleaseData
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(obr_, name, typeName, dict),
    vol_(gSum(mesh_.V().field())),
    TMin_(0.0),
    TMax_(0.0),
    TAvg_(0.0),
    WignerEnergyReleasePower_(0.0),
    releasedWignerEnergy_(0.0),
    energyTransferPower_(0.0),
    transferredEnergy_(0.0),

    // Data for calculating transferred Wigner energy
    patchNames_(dict.getOrDefault<wordRes>("patches", wordRes())),
    patchSet_(),
    area_(0.0),
    QDotc_(0.0)
{
    read(dict);
    writeFileHeader(file());
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::wignerEnergyReleaseData::read
(
    const dictionary& dict
)
{
    fvMeshFunctionObject::read(dict);
    writeFile::read(dict);

    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();

    patchSet_ = pbm.patchSet(patchNames_);

    const surfaceScalarField::Boundary& magSf = mesh_.magSf().boundaryField();

    if (not patchSet_.empty())
    {
        for (const label patchi : patchSet_)
        {
            area_ += gSum(magSf[patchi]);
        }
    }

    return true;
}


bool Foam::functionObjects::wignerEnergyReleaseData::execute()
{
    if (foundObject<solidThermo>(solidThermo::dictName))
    {
        const solidThermo& thermo =
            lookupObject<solidThermo>(solidThermo::dictName);

        calcHeatPower(thermo.alpha(), thermo.he());

        calcWignerEnergyData(thermo.rho());
    }
    else
    {
        FatalErrorInFunction
            << "Unable to find solid thermal model in the database"
            << exit(FatalError);
    }

    return true;
}


bool Foam::functionObjects::wignerEnergyReleaseData::write()
{
    if (Pstream::master())
    {
        writeCurrentTime(file());

        file()
            << tab << setw(charWidth()) << TMin_
            << tab << setw(charWidth()) << TMax_
            << tab << setw(charWidth()) << TAvg_
            << tab << setw(charWidth()) << WignerEnergyReleasePower_
            << tab << setw(charWidth()) << energyTransferPower_
            << tab << setw(charWidth()) << releasedWignerEnergy_
            << tab << setw(charWidth()) << transferredEnergy_
            << endl;
    }

    return true;
}


// ************************************************************************* //