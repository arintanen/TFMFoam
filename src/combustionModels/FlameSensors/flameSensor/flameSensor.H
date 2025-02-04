/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
     \\/     M anipulation  |
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

Namespace
    Foam::flameSensorModels

Description
    Namespace for laminar flame speed models


Class
    Foam::flameSensor

Description
    Abstract class for laminar flame speed

SourceFiles
    flameSensor.C
    flameSensorNew.C

\*---------------------------------------------------------------------------*/

#ifndef flameSensor_H
#define flameSensor_H

#include "dictionary.H"
#include "volFields.H"
#include "fvMesh.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class flameSensor Declaration
\*---------------------------------------------------------------------------*/

class flameSensor
{
protected:

    // Protected data

        const fvMesh& mesh_;
        volScalarField S_;

public:

    //- Runtime type information
    TypeName("flameSensor");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            flameSensor,
            dictionary,
            (
                const dictionary& dict,
                const fvMesh& mesh
            ),
            (dict, mesh)
        );


    // Constructors

        //- Construct from dictionary and psiuReactionThermo
        flameSensor
        (
            const dictionary&,
            const fvMesh& mesh
        );

        //- Disallow default bitwise copy construction
        flameSensor(const flameSensor&) = delete;


    // Selector

        static autoPtr<flameSensor> New
        (
            const dictionary& dict,
            const fvMesh& mesh
            
        );

    // Member functions
    
    virtual void correct() = 0;
    
    const volScalarField& S() const{
    	return S_;
    }
     
    //- Destructor
    virtual ~flameSensor();

        //- Disallow default bitwise assignment
        void operator=(const flameSensor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
