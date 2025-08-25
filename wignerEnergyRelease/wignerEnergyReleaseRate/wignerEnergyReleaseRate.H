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
    Foam::Wigner::wignerEnergyReleaseRate

Description
    Evaluates Wigner energy release rate as a function of released energy
    and temperature using the Piecewise Linear Wigner energy release rate
    Calculation (PLWC) method. The class also manages reading the discretized
    PLWC dataset.

SourceFiles
    wignerEnergyReleaseRate.C

\*---------------------------------------------------------------------------*/

#ifndef wignerEnergyReleaseRate_H
#define wignerEnergyReleaseRate_H

#include "IOdictionary.H"
#include "interpolationTable.H"
#include "fvMesh.H"

#include "plwcDataset.H"
#include "pTraitsPlwcDataset.H"
#include "componentsPlwcDataset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

using Wigner::operator-;
using Wigner::operator+;
using Wigner::operator*;
using Wigner::operator/;

namespace Wigner
{

/*---------------------------------------------------------------------------*\
                Class wignerEnergyReleaseRate Declaration
\*---------------------------------------------------------------------------*/

class wignerEnergyReleaseRate : public IOdictionary
{
private:

    // Private data

        //- PLWC discretized dataset {S, E_k, T_1, T_2, S_dot_1, S_dot_2}
        //    S is the reference column
        //    {E_k, T_1, T_2, S_dot_1, S_dot_2} = plwcDataset = FixedList<scalar,5>
        interpolationTable<plwcDataset> dataset_;

        //- Loaded properly or not
        bool hasDataset_;


public:

    //- Static data member
    static const char* const typeName;


    //- Constructor
    wignerEnergyReleaseRate(const fvMesh& mesh);

    //- Destructor
    virtual ~wignerEnergyReleaseRate() = default;


    // Member functions

        //- Calculate release rate with given temperature T and released energy S
        scalar releaseRate(const scalar T, const scalar S) const;


        // Access functions

            //- Return dataset_
            const interpolationTable<plwcDataset>& dataset() const noexcept
            {
                return dataset_;
            }

            //- Return dataset_
            interpolationTable<plwcDataset>& dataset()
            {
                return dataset_;
            }

            //- Return hasDataset_
            const bool& hasDataset() const noexcept
            {
                return hasDataset_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Wigner
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
