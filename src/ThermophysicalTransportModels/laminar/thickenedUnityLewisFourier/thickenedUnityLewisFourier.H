/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2022 OpenFOAM Foundation
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
    Foam::laminarThermophysicalTransportModels::thickenedUnityLewisFourier

Description
    thickenedUnityLewisFourier's energy gradient heat flux model for laminar flow.
    Specie fluxes are computed assuming a unity turbulent Lewis number.

SourceFiles
    thickenedUnityLewisFourier.C

\*---------------------------------------------------------------------------*/

#ifndef thickenedUnityLewisFourier_H
#define thickenedUnityLewisFourier_H

#include "laminarThermophysicalTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                           Class thickenedUnityLewisFourier Declaration
\*---------------------------------------------------------------------------*/

template<class BasicThermophysicalTransportModel>
class thickenedUnityLewisFourier
:
    public unityLewisFourier
    <
        BasicThermophysicalTransportModel
    >
    
{
const fvMesh& mesh_;
public:

    typedef typename BasicThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename BasicThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename BasicThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("thickenedUnityLewisFourier");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        thickenedUnityLewisFourier
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Construct from a type name, a momentum transport model and a thermo
        //  model
        thickenedUnityLewisFourier
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );


    //- Destructor
    virtual ~thickenedUnityLewisFourier()
    {}


        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const
        {
        const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
            return volScalarField::New
            (
                "DEff",
                this->thermo().kappa()/this->thermo().Cp()*EF
            );
        }


        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thickenedUnityLewisFourier.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
