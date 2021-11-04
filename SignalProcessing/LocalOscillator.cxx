// LocalOscillator.cxx

#include "SignalProcessing/LocalOscillator.h"

rad::LocalOscillator::LocalOscillator(double angFreq) {
  angularFreq = angFreq;
}

rad::LocalOscillator::LocalOscillator() {
  angularFreq = 0;
}

rad::LocalOscillator::~LocalOscillator() { }
