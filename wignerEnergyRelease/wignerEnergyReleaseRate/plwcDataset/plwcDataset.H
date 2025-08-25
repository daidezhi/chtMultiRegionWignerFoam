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

Typedef
    Foam::Wigner::plwcDataset

Description
    A plwcDataset is {E_k, T_1, T_2, S_dot_1, S_dot_2}.

\*---------------------------------------------------------------------------*/

#ifndef plwcDataset_H
#define plwcDataset_H

#include "FixedList.H"
#include "direction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Wigner
{
    typedef FixedList<scalar,5> plwcDataset; // E_k, T_1, T_2, S_dot_1, S_dot_2

    // a - b
    inline plwcDataset operator-(const plwcDataset& a, const plwcDataset& b)
    {
        plwcDataset r(Zero);
        for (direction i=0; i<5; ++i)
        {
            r[i] = a[i] - b[i];
        }

        return r;
    }

    // a + b
    inline plwcDataset operator+(const plwcDataset& a, const plwcDataset& b)
    {
        plwcDataset r(Zero);
        for (direction i=0; i<5; ++i)
        {
            r[i] = a[i] + b[i];
        }
        return r;
    }

    // a * s  (right scalar)
    inline plwcDataset operator*(const plwcDataset& a, const scalar s)
    {
        plwcDataset r(Zero);
        for (direction i=0; i<5; ++i)
        {
            r[i] = a[i] * s;
        }
        return r;
    }

    // s * a  (left scalar)
    inline plwcDataset operator*(const scalar s, const plwcDataset& a)
    {
        plwcDataset r(Zero);
        for (direction i=0; i<5; ++i)
        {
            r[i] = s * a[i];
        }
        return r;
    }
    
    // a / s
    inline plwcDataset operator/(const plwcDataset& a, const scalar s)
    {
        plwcDataset r(Zero);
        const scalar inv = 1.0/s;
        for (direction i=0; i<5; ++i)
        {
            r[i] = a[i] * inv;
        }
        return r;
    }
}
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //