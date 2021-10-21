// BasicFunctions.h
#ifndef BASIC_FUNCTIONS_H
#define BASIC_FUNCTIONS_H

#include "TVector3.h"

TVector3 CalcEField(const TVector3 sourcePosition, const TVector3 ePosition,
		    const TVector3 eVelocity, const TVector3 eAcceleration);

#endif
