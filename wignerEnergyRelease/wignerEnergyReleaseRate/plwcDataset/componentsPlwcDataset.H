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

Description
    Add component functions of plwcDataset.

\*---------------------------------------------------------------------------*/

#ifndef componentsPlwcDataset_H
#define componentsPlwcDataset_H

#include "plwcDataset.H"
#include "pTraitsPlwcDataset.H"
#include "direction.H"
#include "word.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
// Return i-th component (read-only)
inline scalar component(const Wigner::plwcDataset& v, direction i)
{
    return v[i];
}

// Return i-th component (assignable reference)
inline scalar& setComponent(Wigner::plwcDataset& v, direction i)
{
    return v[i];
}

// Return index by name
inline direction componentIndex(const word& name)
{
    constexpr direction N = pTraits<Wigner::plwcDataset>::nComponents;

    wordList valid(N);

    for (direction i=0; i<N; ++i)
    {
        valid[i] = pTraits<Wigner::plwcDataset>::componentNames[i];
        if (name == valid[i]) return i;
    }

    FatalErrorInFunction
        << "Unknown component name: " << name << nl
        << "Valid: " << valid << nl
        << exit(FatalError);

    return 0;   // never reached
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
