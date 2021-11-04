/*
  Signal.h
  Class for a signal detected on an antenna or antennas
*/

#ifndef SIGNAL_H
#define SIGNAL_H

#include "TGraph.h"

namespace rad
{
  /// Signal base class
  /// Class members are a voltage signal in the time domain

  class Signal
  {
  private:
    TGraph* grSignalTime;

  public:
    Signal();
    ~Signal();
  };
}

#endif
