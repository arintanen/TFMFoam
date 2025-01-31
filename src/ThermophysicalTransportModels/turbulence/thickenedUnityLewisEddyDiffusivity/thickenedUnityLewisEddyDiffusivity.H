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
    Foam::turbulenceThermophysicalTransportModels::thickenedUnityLewisEddyDiffusivity

Description
    Eddy-diffusivity based energy gradient heat flux model for RAS or LES
    of turbulent flow.  Specie fluxes are computed assuming a unity turbulent
    Lewis number.

Usage
    \verbatim
    LES
    {
        model           thickenedUnityLewisEddyDiffusivity;
        Prt             0.85;
    }
    \endverbatim

SourceFiles
    thickenedUnityLewisEddyDiffusivity.C

\*---------------------------------------------------------------------------*/

#ifndef thickenedUnityLewisEddyDiffusivity_H
#define thickenedUnityLewisEddyDiffusivity_H
#include "unityLewisEddyDiffusivity.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace turbulenceThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                       Class thickenedUnityLewisEddyDiffusivity Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceThermophysicalTransportModel>
class thickenedUnityLewisEddyDiffusivity
:
    public unityLewisEddyDiffusivity<TurbulenceThermophysicalTransportModel>
{

const fvMesh& mesh_;

public:

    typedef typename TurbulenceThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        TurbulenceThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename TurbulenceThermophysicalTransportModel::thermoModel
        thermoModel;


    //- Runtime type information
    TypeName("thickenedUnityLewisEddyDiffusivity");


    // Constructors

        //- Construct from a momentum transport model and a thermo model
        thickenedUnityLewisEddyDiffusivity
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        //- Construct from a type name, a momentum transport model and a thermo
        //  model, and whether to default the turbulent Prandtl number to one
        //  if it is not specified
        thickenedUnityLewisEddyDiffusivity
        (
            const word& type,
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo,
            const bool allowDefaultPrt = false
        );


    //- Destructor
    virtual ~thickenedUnityLewisEddyDiffusivity()
    {}


    // Member Functions

        //- Effective thermal turbulent diffusivity for temperature
        //  of mixture [W/m/K]
        virtual tmp<volScalarField> kappaEff() const
        {
        const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
        const volScalarField& S = mesh_.lookupObject<volScalarField>("S");
            return this->thermo().kappaEff((1-S)/EF*this->alphat())*EF;
        }


        //- Effective thermal turbulent diffusivity of mixture [kg/m/s]
        virtual tmp<volScalarField> alphaEff() const
        {
        const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
        const volScalarField& S = mesh_.lookupObject<volScalarField>("S");
            return this->thermo().alphaEff((1-S)*this->alphat()/EF)*EF;
        }


        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const
        {
        const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
        const volScalarField& S = mesh_.lookupObject<volScalarField>("S");
            return volScalarField::New
            (
                "DEff",
                this->thermo().kappa()/this->thermo().Cp()*EF + (1-S)*this->alphat()
            );
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace turbulenceThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "thickenedUnityLewisEddyDiffusivity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
