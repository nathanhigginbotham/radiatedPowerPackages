// FieldClasses.cxx

#include "FieldClasses.h"

#include "TGraph.h"
#include "TVector3.h"

rad::FieldPoint::FieldPoint() {
  Ex = new TGraph();
  Ey = new TGraph();
  Ez = new TGraph();
  Bx = new TGraph();
  By = new TGraph();
  Bz = new TGraph();
  antennaPoint = TVector3(0.0, 0.0, 0.0);
}

rad::FieldPoint::~FieldPoint() {
  delete Ex;
  delete Ey;
  delete Ez;
  delete Bx;
  delete By;
  delete Bz;
}
