/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | Unsupported Contributions for OpenFOAM
 \\    /   O peration     |
  \\  /    A nd           | Copyright (C) <year> <author name(s)>
   \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is a derivative work of OpenFOAM.

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

#include "chemPointISAT.H"
#include "scalarField.H"
#include "binaryNode.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
template<class CompType, class ThermoType>
binaryNode<CompType, ThermoType>::binaryNode
(
)
:
    elementLeft_(NULL),
    elementRight_(NULL),
    left_(NULL), 
    right_(NULL),
    parent_(NULL)
{}


template<class CompType, class ThermoType>
binaryNode<CompType, ThermoType>::binaryNode
(
    chemPointISAT<CompType, ThermoType>* elementLeft,
    chemPointISAT<CompType, ThermoType>* elementRight,
    binaryNode<CompType, ThermoType>* parent
)
:
    elementLeft_(elementLeft),
    elementRight_(elementRight),
    left_(NULL),
    right_(NULL),
    parent_(parent),
    v_(elementLeft->spaceSize(),0.0)
{
    calcV(elementLeft, elementRight, v_);
    a_ = calcA(elementLeft, elementRight);
}

template<class CompType, class ThermoType>
binaryNode<CompType, ThermoType>::binaryNode
(
    binaryNode<CompType, ThermoType> *bn
)
:
    elementLeft_(bn->elementLeft()),
    elementRight_(bn->elementRight()),
    left_(bn->left()), 
    right_(bn->right()),
    parent_(bn->parent()),
    v_(bn->v()),
    a_(bn->a())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class CompType, class ThermoType>
void
binaryNode<CompType, ThermoType>::calcV
(
    chemPointISAT<CompType, ThermoType>*& elementLeft,
    chemPointISAT<CompType, ThermoType>*& elementRight,
    scalarField& v
)
{
    //LT is the transpose of the L matrix
    List<List<scalar> >& LT = elementLeft->LT();
    label dim = elementLeft->spaceSize();
    if (elementLeft->DAC()) dim = elementLeft->NsDAC()+2;
    scalarField phiDif = elementRight->phi() - elementLeft->phi();
    const scalarField& scaleFactor = elementLeft->scaleFactor();
    const scalar epsTol = elementLeft->epsTol();

    for (label i=0; i<elementLeft->spaceSize(); i++)
    {
        label si = i;
        bool outOfIndexI = true;
        if (elementLeft->DAC())
        {
            if (i<elementLeft->spaceSize()-2)
            {
                si = elementLeft->completeToSimplifiedIndex(i);
                outOfIndexI = (si==-1);
            }
            else
            {
                outOfIndexI = false;
                label dif = i-(elementLeft->spaceSize()-2);
                si = elementLeft->NsDAC()+dif;
            }
        }
        if (!(elementLeft->DAC()) || (elementLeft->DAC() && !(outOfIndexI)))
        {
            v[i]=0.0;
            for (label j=0; j<elementLeft->spaceSize(); j++)
            {
                label sj = j;
                bool outOfIndexJ = true;
                if (elementLeft->DAC())
                {
                    if (j<elementLeft->spaceSize()-2)
                    {
                        sj = elementLeft->completeToSimplifiedIndex(j);
                        outOfIndexJ = (sj==-1);
                    }
                    else
                    {
                        outOfIndexJ = false;
                        label dif = j-(elementLeft->spaceSize()-2);
                        sj = elementLeft->NsDAC()+dif;
                    }
                }
                if (!(elementLeft->DAC()) || (elementLeft->DAC() && !(outOfIndexJ)))
                {
                    //since L is a lower triangular matrix k=0->min(i,j)
                    for (label k=0; k<=min(si,sj); k++)
                    {
                        v[i] += LT[k][si]*LT[k][sj]*phiDif[j];
                    }
                }
            }
        }
        else
        {
            //when it is an inactive species the diagonal element of LT is
            //  1/(scaleFactor*epsTol)
            scalar div = scaleFactor[i]*scaleFactor[i]*epsTol*epsTol;
            v[i] = phiDif[i]/div;
        }
    }
}


template<class CompType, class ThermoType>
scalar binaryNode<CompType, ThermoType>::calcA
(
    chemPointISAT<CompType, ThermoType>* elementLeft,
    chemPointISAT<CompType, ThermoType>* elementRight
)
{
    scalar a = 0.0;
    scalarField phih = (elementLeft->phi()+elementRight->phi())/2;
    label spaceSize = elementLeft->spaceSize();
    const scalarField& V = v();
    for (label i=0; i<spaceSize; i++)
    {
        a += V[i]*phih[i]; 
    }
    return a;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //