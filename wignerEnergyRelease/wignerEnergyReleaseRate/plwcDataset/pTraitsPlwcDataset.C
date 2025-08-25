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

#include "pTraitsPlwcDataset.H"
#include "zero.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::pTraits<Foam::Wigner::plwcDataset>::typeName
    = "plwcDataset";


const char* const Foam::pTraits<Foam::Wigner::plwcDataset>::componentNames[] =
{
    "Ek","T1","T2","SDot1","SDot2"
};


const Foam::Wigner::plwcDataset Foam::pTraits<Foam::Wigner::plwcDataset>::zero
    = Foam::Wigner::plwcDataset(Zero);


const Foam::Wigner::plwcDataset Foam::pTraits<Foam::Wigner::plwcDataset>::one
    = Foam::Wigner::plwcDataset(1);


// ************************************************************************* //