/*---------------------------------------------------------------------------*\
 =========                 |
 \\      /  F ield         | Unsupported Contributions for OpenFOAM
  \\    /   O peration     |
   \\  /    A nd           | Copyright (C) 2014 F. Contino, S. Backaert,
    \\/     M anipulation  |                    N. Bourgeois, T. Lucchini
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
 
Class
    Foam::mechanismReduction

Description
    An abstract class for reducing chemical mechanism

SourceFiles
    mechanismReduction.C

\*---------------------------------------------------------------------------*/



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CompType, class ThermoType>
inline const Foam::List<bool>&
Foam::mechanismReduction<CompType,ThermoType>::activeSpecies() const
{
    return activeSpecies_;
}

template<class CompType, class ThermoType>
inline Foam::label
Foam::mechanismReduction<CompType,ThermoType>::NsSimp()
{
    return NsSimp_;
}


template<class CompType, class ThermoType>
inline Foam::label
Foam::mechanismReduction<CompType,ThermoType>::nSpecie()
{
    return nSpecie_;
}


template<class CompType, class ThermoType>
inline Foam::scalar
Foam::mechanismReduction<CompType,ThermoType>::tolerance() const
{
    return tolerance_;
}


template<class CompType, class ThermoType>
inline Foam::Switch
Foam::mechanismReduction<CompType,ThermoType>::active() const
{
    return active_;
}



// ************************************************************************* //
