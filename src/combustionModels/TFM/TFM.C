/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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

#include "TFM.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace combustionModels
{
    defineTypeNameAndDebug(TFM, 0);
    addToRunTimeSelectionTable(combustionModel, TFM, dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::combustionModels::TFM::TFM
(
    const word& modelType,
    const fluidReactionThermo& thermo,
    const compressibleMomentumTransportModel& turb,
    const word& combustionProperties
)
:
    laminar(modelType, thermo, turb, combustionProperties),
    Fmax_(this->coeffs().template lookup<scalar>("Fmax")),
    F_
     (
         IOobject
         (
             "F",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    EF_
     (
         IOobject
         (
             "EF",
             this->mesh().time().timeName(),
             this->mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->mesh(),
         dimensionedScalar(dimless, 1)
     ),
    flameSensor_(flameSensor::New(coeffs(),turb.mesh())),
    efficiencyFunction_(efficiencyFunction::New(turb.mesh(),coeffs(),turb,Fmax_))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::combustionModels::TFM::~TFM()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::combustionModels::TFM::correct()
{
    laminar::correct();
    flameSensor_->correct();
    efficiencyFunction_->correct();
    F_ = (Fmax_-1)*flameSensor_->S() +1;
    EF_ = efficiencyFunction_->E()*F_;

}


Foam::tmp<Foam::fvScalarMatrix>
Foam::combustionModels::TFM::R(volScalarField& Y) const
{
    return efficiencyFunction_->E()/F_*laminar::R(Y);
}


Foam::tmp<Foam::volScalarField>
Foam::combustionModels::TFM::Qdot() const
{
    return volScalarField::New
    (
        this->thermo().phasePropertyName(typeName + ":Qdot"),
        efficiencyFunction_->E()/F_*laminar::Qdot()
    );
}


bool Foam::combustionModels::TFM::read()
{
    if (laminar::read())
    {
        this->coeffs().lookup("Fmax") >> Fmax_;
        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
