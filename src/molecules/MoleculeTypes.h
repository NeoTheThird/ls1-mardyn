/*
 * MoleculeTypes.h
 *
 *  Created on: May 10, 2011
 *      Author: eckhardw
 */

#ifndef MOLECULETYPES_H_
#define MOLECULETYPES_H_

#include "molecules/BasicMolecule.h"
#include "molecules/CachingMolecule.h"

typedef BasicMolecule Molecule;
//typedef CachingMolecule Molecule;

typedef CachingMolecule HandlerMoleculeType;
//typedef BasicMolecule HandlerMoleculeType;

typedef utils::DynamicArray<HandlerMoleculeType, false, false> HandlerMoleculeTypeArray;
typedef utils::DynamicArray<Molecule, false, false> MoleculeArray;

/**
 * The following code can be used to determine if the types Molecule and HandlerMoleculeType
 * are the same:
 *
 * \verbatim
 * if ( IsSame<FirstClass,SecondClass>::Result::value ) {
 *   ...
 * }
 * \endverbatim
 */

struct FalseType { enum { value = false }; };
struct TrueType { enum { value = true }; };


template <typename T1, typename T2>
struct IsSame
{
typedef FalseType Result;
};


template <typename T>
struct IsSame<T,T>
{
typedef TrueType Result;
};


#endif /* MOLECULETYPES_H_ */
