/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2026 A. N. Nuriev, Kazan Federal University
-------------------------------------------------------------------------------
License
    This file is an extension library for OpenFOAM; it is not part of OpenFOAM.

    This file is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    It is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with this file.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "linearJointSpring.H"
#include "rigidBodyModel.H"
#include "rigidBodyModelState.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{
    defineTypeNameAndDebug(linearJointSpring, 0);

    addToRunTimeSelectionTable
    (
        restraint,
        linearJointSpring,
        dictionary
    );
}
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::RBD::restraints::linearJointSpring::linearJointSpring
(
    const word& name,
    const dictionary& dict,
    const rigidBodyModel& model
)
:
    restraint(name, dict, model),
    refAngle_(0),
    stiffness_(0),
    damping_(0)
{
    read(dict);
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::RBD::restraints::linearJointSpring::restrain
(
    scalarField& tau,
    Field<spatialVector>& fx,
    const rigidBodyModelState& state
) const
{
    // Joint connecting the restrained body to its parent
    // (the joint list is indexed by body, cf. rigidBodyModel::forwardDynamics)
    const joint& jnt = model_.joints()[bodyIndex_];

    if (jnt.nDoF() != 1)
    {
        FatalErrorInFunction
            << "Restraint " << name() << ": the joint of body "
            << model_.name(bodyIndex_) << " has " << jnt.nDoF()
            << " degrees of freedom; linearJointSpring requires a"
            << " single-DoF joint (e.g. Rz)." << nl
            << exit(FatalError);
    }

    const label qi = jnt.qIndex();

    // Joint coordinate relative to the reference angle and its rate
    const scalar theta = state.q()[qi] - refAngle_;
    const scalar thetaDot = state.qDot()[qi];

    // Hinge moment (generalised force of the joint)
    const scalar moment = -(stiffness_*theta + damping_*thetaDot);

    tau[qi] += moment;

    if (model_.debug)
    {
        Info<< " " << name()
            << " body " << model_.name(bodyIndex_)
            << " theta " << theta
            << " thetaDot " << thetaDot
            << " moment " << moment
            << endl;
    }
}


bool Foam::RBD::restraints::linearJointSpring::read
(
    const dictionary& dict
)
{
    restraint::read(dict);

    coeffs_.readEntry("stiffness", stiffness_);
    coeffs_.readEntry("damping", damping_);
    refAngle_ = coeffs_.getOrDefault<scalar>("referenceAngle", 0);

    return true;
}


void Foam::RBD::restraints::linearJointSpring::write
(
    Ostream& os
) const
{
    restraint::write(os);

    os.writeEntry("referenceAngle", refAngle_);
    os.writeEntry("stiffness", stiffness_);
    os.writeEntry("damping", damping_);
}


// ************************************************************************* //
