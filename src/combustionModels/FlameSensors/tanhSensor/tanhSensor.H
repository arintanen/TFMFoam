/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2013-2020 OpenFOAM Foundation
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

Class
    Foam::laminarFlameSpeedModels::tanhSensor

Description
    Laminar flame speed obtained from Ravi and Petersen's correlation.

    The correlation for the laminar flame speed \f$Su\f$ is of the following
    form:
    \f[
        Su = \left( \sum \alpha_i \phi^i \right)
        \left( \frac{T}{T_{ref}} \right)^{\left( \sum \beta_j \phi^j \right)}
    \f]

    Where \f$\phi\f$ is the equivalence ratio, and \f$\alpha\f$ and \f$\beta\f$
    are polynomial coefficients given for a number of pressure and equivalence
    ratio points.

SourceFiles
    tanhSensor.C

\*---------------------------------------------------------------------------*/

#ifndef tanhSensor_H
#define tanhSensor_H

#include "flameSensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace flameSensorModels
{

/*---------------------------------------------------------------------------*\
                        Class tanhSensor Declaration
\*---------------------------------------------------------------------------*/

class tanhSensor
:
    public flameSensor
{
    // Private Data
    dictionary coeffsDict_;
    scalar beta_;
    word fieldName_;
    

public:

    //- Runtime type information
    TypeName("tanh");


    // Constructors

        //- Construct from dictionary and psiuReactionThermo
        tanhSensor
        (
            const dictionary&,
            const fvMesh&
        );

        //- Disallow default bitwise copy construction
        tanhSensor(const tanhSensor&) = delete;


    //- Destructor
    virtual ~tanhSensor();
    void correct() override;

    // Member Operators
    
        //- Disallow default bitwise assignment
        void operator=(const tanhSensor&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End flameSensorModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
