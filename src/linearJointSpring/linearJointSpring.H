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

Class
    Foam::RBD::restraints::linearJointSpring

Group
    grpRigidBodyDynamicsRestraints

Description
    Linear spring-damper acting on the single degree of freedom of the joint
    that connects the given body to its parent, e.g. a revolute joint Rz.
    The generalised joint force is

        tau = -stiffness*(q - referenceAngle) - damping*qDot

    where q is the joint coordinate (rotation of the body relative to its
    parent) and qDot is its rate. Since tau is applied in joint space, the
    moment acts on the body and, with the opposite sign, on its parent, and
    the damping term depends on the relative angular velocity.

    The restraint is intended for approximating an Euler-Bernoulli beam by a
    chain of rigid segments, segment i being attached to segment i-1 by an
    Rz joint with

        stiffness = 2*EI/(dx_{i-1} + dx_i),  damping = eta*stiffness,

    where dx_{-1} = 0 for the first segment (clamped end, stiffness 2*EI/dx_0).

Usage
    \verbatim
    restraints
    {
        hinge_segment0_segment1
        {
            type            linearJointSpring;
            body            segment1;     // restrained body (joint to its parent)
            stiffness       100;          // [N m/rad]
            damping         0.1;          // [N m s/rad]
            referenceAngle  0;            // [rad], optional (default 0)
        }
    }
    \endverbatim

SourceFiles
    linearJointSpring.C

\*---------------------------------------------------------------------------*/

#ifndef RBD_restraints_linearJointSpring_H
#define RBD_restraints_linearJointSpring_H

#include "rigidBodyRestraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RBD
{
namespace restraints
{

/*---------------------------------------------------------------------------*\
                      Class linearJointSpring Declaration
\*---------------------------------------------------------------------------*/

class linearJointSpring
:
    public restraint
{
    // Private data

        //- Joint coordinate at which the moment vanishes [rad]
        scalar refAngle_;

        //- Spring stiffness [N m/rad]
        scalar stiffness_;

        //- Damping coefficient [N m s/rad]
        scalar damping_;


public:

    //- Runtime type information
    TypeName("linearJointSpring");


    // Constructors

        //- Construct from components
        linearJointSpring
        (
            const word& name,
            const dictionary& dict,
            const rigidBodyModel& model
        );

        //- Construct and return a clone
        virtual autoPtr<restraint> clone() const
        {
            return autoPtr<restraint>
            (
                new linearJointSpring(*this)
            );
        }


    //- Destructor
    virtual ~linearJointSpring() = default;


    // Member Functions

        //- Add the restraint moment to the joint force field tau
        virtual void restrain
        (
            scalarField& tau,
            Field<spatialVector>& fx,
            const rigidBodyModelState& state
        ) const;

        //- Update properties from given dictionary
        virtual bool read(const dictionary& dict);

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace restraints
} // End namespace RBD
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
