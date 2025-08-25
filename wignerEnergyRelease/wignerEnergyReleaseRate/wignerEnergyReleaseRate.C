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

#include "wignerEnergyReleaseRate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::Wigner::wignerEnergyReleaseRate::typeName =
"wignerEnergyReleaseRate";


// * * * * * * * * * * * * * * * * Constructor * * * * * * * * * * * * * * * //

Foam::Wigner::wignerEnergyReleaseRate::wignerEnergyReleaseRate
(
    const fvMesh& mesh
)
:
    IOdictionary
    (
        IOobject
        (
            "plwcDiscretizedDataset",
            mesh.time().rootPath()/mesh.time().globalCaseName()/"constant"/mesh.name(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        )
    ),
    dataset_(),
    hasDataset_(false)
{
    if (headerOk())  // file was present and parsed
    {
        if (this->found("plwcDataset"))
        {
            dataset_ = interpolationTable<Wigner::plwcDataset>(this->subDict("plwcDataset"));

            hasDataset_ = !dataset_.empty();
        }
    }
}


// * * * * * * * * * * * * Public Member Functions * * * * * * * * * * * * * //

Foam::scalar Foam::Wigner::wignerEnergyReleaseRate::releaseRate
(
    const Foam::scalar T,
    const Foam::scalar S
) const
{
    const Foam::Wigner::plwcDataset intePlwcDataset(dataset_.operator()(S));

    //Info << component(intePlwcDataset, componentIndex("Ek")) << endl;
    //Info << intePlwcDataset << endl;

    const scalar& EkS(intePlwcDataset[0]);
    const scalar& T1S(intePlwcDataset[1]);
    const scalar& T2S(intePlwcDataset[2]);
    const scalar& SDot1S(intePlwcDataset[3]);
    const scalar& SDot2S(intePlwcDataset[4]);

    const scalar SDot1(SDot1S * Foam::exp(-EkS * (1.0/T - 1.0/T1S)));
    const scalar SDot2(SDot2S * Foam::exp(-EkS * (1.0/T - 1.0/T2S)));

    return (0.5 * (SDot1 + SDot2));
}


// ************************************************************************* //