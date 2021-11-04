// Signal.cxx

#include "SignalProcessing/Signal.h"

#include "TGraph.h"

rad::Signal::Signal() {
  grSignalTime = new TGraph();
}

rad::Signal::~Signal() {
  delete grSignalTime;
}


