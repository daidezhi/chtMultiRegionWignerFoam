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
    Foam::Wigner::pTraits<plwcDataset>

Description
    Specialised version of plwcDataset for Foam::pTraits<>.

SourceFiles
    pTraitsPlwcDataset.C

\*---------------------------------------------------------------------------*/

#ifndef pTraitsPlwcDataset_H
#define pTraitsPlwcDataset_H

#include "plwcDataset.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class pTraits<plwcDataset> Declaration
\*---------------------------------------------------------------------------*/

template<>

class pTraits<Wigner::plwcDataset>
{
private:

    //- Private data
    Wigner::plwcDataset p_;


public:
    
    //- Type define
    typedef scalar      cmptType;

    //- Static data members
    enum { rank = 0 };
    static const direction nComponents = 5;

    //- Optional, but useful for diagnostics and CSV headers
    static const char* const typeName;
    static const char* const componentNames[];

    static const Wigner::plwcDataset zero;
    static const Wigner::plwcDataset one;

    //- Standard pTraits ctors/conversions
    pTraits()
    {
    }

    pTraits(const Wigner::plwcDataset& val) : p_(val)
    {
    }

    operator const Wigner::plwcDataset&() const
    {
        return p_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //