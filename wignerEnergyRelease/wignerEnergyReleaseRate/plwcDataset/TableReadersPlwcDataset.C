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

#include "tableReader.H"
#include "tableReaders.H"
#include "FixedList.H"
#include "addToRunTimeSelectionTable.H"
#include "csvTableReader.H"

#include "plwcDataset.H"
#include "pTraitsPlwcDataset.H"
#include "componentsPlwcDataset.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    typedef Wigner::plwcDataset WignerPlwcDataset;
    
    // Instantiate the base selection table for plwcDataset
    defineNamedTemplateTypeNameAndDebug(tableReader<WignerPlwcDataset>, 0);
    defineTemplateRunTimeSelectionTable(tableReader<WignerPlwcDataset>, dictionary);

    // Register concrete readers for plwcDataset
    makeTableReaderType(csvTableReader, WignerPlwcDataset);
}

// ************************************************************************* //