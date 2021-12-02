/*
  BaseField.h
 
  Abstract base class for magnetic fields used electron trackin.
  Transcription of T. Goffrey's original python code 
*/

#ifndef BASE_FIELD_H
#define BASE_FIELD_H

#include "TVector3.h"

namespace rad
{
  class BaseField {
    
  public:
    virtual TVector3 evaluate_field_at_point(const TVector3 vec) = 0;

    double evaluate_field_magnitude(const TVector3 vec);
  };
  
}

#endif
