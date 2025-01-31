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

\*---------------------------------------------------------------------------*/

#include "powerLaw.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcCurl.H"
#include "fvcLaplacian.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace efficiencyFunctionModels
{
    defineTypeNameAndDebug(powerLaw, 0);

    addToRunTimeSelectionTable
    (
        efficiencyFunction,
        powerLaw,
        dictionary
    );
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::efficiencyFunctionModels::powerLaw::powerLaw
(
    const dictionary& dict,
    const fvMesh& mesh,
    const compressibleMomentumTransportModel& turb,
    scalar Fmax
)
:
    efficiencyFunction(dict, mesh, turb, Fmax),
    coeffsDict_(dict.subDict("powerLawCoeffs")),
    SL_(coeffsDict_.lookup<scalar>("SL")),
    delta_L_(coeffsDict_.lookup<scalar>("delta_L")),
    n_filters_(coeffsDict_.lookup<scalar>("n_filters")),
    sFilter_(mesh_)
{
}

void Foam::efficiencyFunctionModels::powerLaw::correct() {
    volVectorField U_hat(turb_.U());
    for(int i = 0; i < n_filters_; i++) {
      U_hat = sFilter_(U_hat);
    }
    volScalarField u_d = mag(fvc::laplacian(fvc::curl(U_hat)));
    const scalar b = 1.4;
    const scalar C_k = 1.5; 
    const scalar alpha = 0.5;
    const scalar pi34 = pow(Foam::constant::mathematical::pi,4.0/3.0);
    const scalarField& V = E_.mesh().V();
    
    scalar delta_e = F_*delta_L_; 
    scalar f_d = sqrt(27*C_k*pi34/110.0*(pow(delta_e/delta_L_,4.0/3.0)-1.0));
    
    forAll(E_,celli) {
      scalar u_dot = 2*V[celli]*u_d[celli];
      scalar Re_d = 4*delta_e*u_dot/delta_L_/SL_; 
      if (Re_d > 0) {
        scalar tp = exp(-3.0/2.0*C_k*pi34/Re_d);
        scalar f_Re = sqrt(9.0/55.0*tp*Re_d);
        if(f_Re > 0) {
          scalar a = 0.6 + 0.2*exp(-0.1*u_dot/SL_) - 0.2*exp(-0.01*delta_e/delta_L_);
          scalar f_u = 4*sqrt(27.0*C_k/110.0)*18*C_k/55.0*pow(u_dot/SL_,2);
          scalar f1 = pow(1/f_u,a);
          scalar f2 = pow(1/f_d,a);
          scalar f3 = pow(f1 + f2 ,b/a);
          scalar f4 = pow(1/f_Re,b);
          scalar Gamma = pow(f3 + f4,-1.0/b);
          E_[celli] = pow(1 + min(delta_e/delta_L_,Gamma),alpha);
        }
      }
    }
 }


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::efficiencyFunctionModels::powerLaw::~powerLaw()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// ************************************************************************* //
