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
    Foam::Wigner::wignerEnergyRelease

Description
    Provides the evaluation and storage of Wigner energy release in solid
    regions. The class manages the calculation of stored energy (S) and
    energy release rate (SDot) fields.

    The class controls allocation of S and SDot fields, ensures they are
    only updated once per time step, and supports optional creation when
    no dataset is present.  It also supplies run-time access to S and SDot
    as registered volScalarFields for coupling to energy equations or for
    post-processing.

SourceFiles
    wignerEnergyRelease.C

\*---------------------------------------------------------------------------*/

#ifndef wignerEnergyRelease_H
#define wignerEnergyRelease_H

#include "volFields.H"

#include "wignerEnergyReleaseRate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace Wigner
{

/*---------------------------------------------------------------------------*\
                    Class wignerEnergyRelease Declaration
\*---------------------------------------------------------------------------*/

class wignerEnergyRelease
{
private:

    // Private data

        // External data references

            //- Reference to the mesh
            const fvMesh& mesh_;

            //- Reference to the temperature field
            const volScalarField& T_;

            //- Often used reference to T_ internal field
            const scalarField& TIn_;


        // Wigner energy fields

            //- Released Wigner energy field
            autoPtr<volScalarField> S_;

            //- Wigner energy release rate field
            autoPtr<volScalarField> SDot_;


        //- Wigner energy release rate calculator
        wignerEnergyReleaseRate releaseRateCalculator_;

        //- Release Wigner energy or not
        bool isReleasingEnergy_;

        // Step update guard
        bool  updated_;
        label lastTimeIndex_;


public:

    //- Static data member
    static const char* const typeName;


    //- Construct with given temperature field and mesh
    wignerEnergyRelease(const volScalarField& T, const fvMesh& mesh);

    //- Destructor
    virtual ~wignerEnergyRelease() = default;


    // Member functions

        //- Update once per time step
        void update();


        // Access functions

            //- Return S_
            const volScalarField& S() const
            {
                return S_();
            }

            //- Return SDot_
            const volScalarField& SDot() const
            {
                return SDot_();
            }

            //- Return S_
            volScalarField& S()
            {
                return S_();
            }

            //- Return SDot_
            volScalarField& SDot()
            {
                return SDot_();
            }

            //- Return isReleasingEnergy_
            inline bool isReleasingEnergy() const
            {
                return isReleasingEnergy_;
            }

            //- Return updated_
            inline bool updated() const
            {
                return updated_;
            }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Wigner
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //