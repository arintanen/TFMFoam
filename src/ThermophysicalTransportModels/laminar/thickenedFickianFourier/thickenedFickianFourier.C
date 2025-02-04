/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020-2021 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "thickenedFickianFourier.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



template<class laminarThermophysicalTransportModel>
thickenedFickianFourier<laminarThermophysicalTransportModel>::
thickenedFickianFourier
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    Fickian<unityLewisFourier<laminarThermophysicalTransportModel>>
    (
        typeName,
        momentumTransport,
        thermo
    ),
    mesh_(momentumTransport.mesh())
    
{
    read();
    this->correct();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
bool
thickenedFickianFourier<laminarThermophysicalTransportModel>::read()
{
    return Fickian
    <
        unityLewisFourier<laminarThermophysicalTransportModel>
    >::read();
}

template<class laminarThermophysicalTransportModel>
tmp<volScalarField> thickenedFickianFourier<laminarThermophysicalTransportModel>::DEff(const volScalarField& Yi) const {
const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
 return Fickian<unityLewisFourier<laminarThermophysicalTransportModel>>::DEff(Yi)*EF;
}

template<class laminarThermophysicalTransportModel>
tmp<volScalarField> thickenedFickianFourier<laminarThermophysicalTransportModel>::kappaEff() const {
const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
 return Fickian<unityLewisFourier<laminarThermophysicalTransportModel>>::kappaEff()*EF;
}

template<class laminarThermophysicalTransportModel>
tmp<volScalarField> thickenedFickianFourier<laminarThermophysicalTransportModel>::alphaEff() const {
const volScalarField& EF = mesh_.lookupObject<volScalarField>("EF");
 return Fickian<unityLewisFourier<laminarThermophysicalTransportModel>>::alphaEff()*EF;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
