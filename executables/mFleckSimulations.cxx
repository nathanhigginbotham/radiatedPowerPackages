/*
    mFleckSimulations.cxx
*/

#include "ElectronDynamics/ComsolFields.h"
#include "ElectronDynamics/QTNMFields.h"
#include "ElectronDynamics/TrajectoryGen.h"
#include "Antennas/IAntenna.h"
#include "Antennas/HalfWaveDipole.h"
#include "SignalProcessing/InducedVoltage.h"
#include "BasicFunctions/BasicFunctions.h"
#include "BasicFunctions/Constants.h"
#include "FieldClasses/FieldClasses.h"

#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TString.h"
#include "TVector3.h"

using namespace rad;

int main(int argc, char *argv[])
{
    TString outputFile{argv[1]};
    double desiredField{atof(argv[2])};
    double scaleFactor{desiredField / 0.94};

    TFile *fout = new TFile(outputFile, "RECREATE");

    const double coilRadius{0.05}; // metres
    const double trapDepth{5e-3};  // Tesla
    const double coilCurrent{2 * trapDepth * coilRadius / MU0};

    std::string fieldFileStem{"/home/sjones/work/qtnm/radiatedPowerPackages/files/"};
    std::string file110mm{"fieldmap_50A_110mmx800mm.csv"};

    ComsolField *field_110mm = new ComsolField(fieldFileStem + file110mm,
                                               scaleFactor);
    ComsolHarmonicField *harm_110mm =
        new ComsolHarmonicField(coilRadius, coilCurrent,
                                fieldFileStem + file110mm, scaleFactor);

    const double centralField{harm_110mm->evaluate_field_at_point(TVector3(0, 0.001, 0)).Mag()};

    // Electron kinematics
    const double electronKE{18600}; // eV
    const double electronSpeed{GetSpeedFromKE(electronKE, ME)};
    const double tau = 2 * R_E / (3 * TMath::C());
    TVector3 v0(electronSpeed, 0, 0);
    const double gyroradius{GetGyroradius(v0, harm_110mm->evaluate_field_at_point(TVector3(0, 0.001, 0)), ME)};
    TVector3 x0(0, -gyroradius, 0);

    const double centralFreq{CalcCyclotronFreq(electronKE, centralField)};
    const double centralPeriod{1.0 / centralFreq};
    const double centralLambda{TMath::C() / centralFreq};
    std::cout << "Electron frequency, wavelength = " << centralFreq / 1e9 << " GHz, " << centralLambda * 1e2 << " cm" << std::endl;

    // Setup the track writing
    const double simTime{1e-7};
    const double simStepSize{1e-12};
    TString trackFile{"~/work/qtnm/outputs/mFleckSimulations/track.root"};
    ElectronTrajectoryGen traj(trackFile, harm_110mm, x0,
                               v0, simStepSize, simTime, 0.0, tau);
    traj.GenerateTraj();
    std::cout << "Generated trajectory\n";

    const double antennaRadius{28.7e-3};
    const double loadResistance{73.0};

    std::cout << "Central field = " << centralField << "T\n";

    int nMeasuredFieldPoints{201};
    double zMin{-0.4};
    double zMax{0.4};
    TGraph *grFieldMag = new TGraph();
    setGraphAttr(grFieldMag);
    grFieldMag->SetTitle("110mm bore; z [m]; |B| [T]");

    TGraph *grHarmFieldMag = new TGraph();
    setGraphAttr(grHarmFieldMag);
    grHarmFieldMag->SetTitle("110mm bore; z [m]; |B| [T]");

    for (int iPnt{0}; iPnt < nMeasuredFieldPoints; iPnt++)
    {
        double thisZ{zMin + (zMax - zMin) * double(iPnt) / double(nMeasuredFieldPoints - 1)};
        TVector3 thisPoint(0, 0.002, thisZ);
        grFieldMag->SetPoint(iPnt, thisZ,
                             field_110mm->evaluate_field_magnitude(thisPoint));
        grHarmFieldMag->SetPoint(iPnt, thisZ,
                                 harm_110mm->evaluate_field_magnitude(thisPoint));
    }
    fout->cd();
    grFieldMag->Write("grFieldMag");
    grHarmFieldMag->Write("grHarmFieldMag");

    const double radiusWavelengths{antennaRadius / centralLambda};

    TGraph *grPowerAntennas = new TGraph();
    setGraphAttr(grPowerAntennas);
    grPowerAntennas->SetTitle(Form("110mm bore B = %d mT; N_{antennas}; Collected Power [fW]", int(desiredField * 1000)));
    grPowerAntennas->SetMarkerStyle(20);

    for (int iAnt{1}; iAnt < 5; iAnt++)
    {
        std::vector<IAntenna *> antennaArray;

        for (int ii{0}; ii < iAnt; ii++)
        {
            double antennaAngle{2 * TMath::Pi() * double(ii) / double(iAnt)};
            double timeShift{centralPeriod * (1.0 - double(ii) / double(iAnt))};
            std::cout << "Dipole " << ii << ": Angle = " << antennaAngle * 180 / TMath::Pi() << " degrees. Required time shift is " << timeShift * 1e12 << " ps" << std::endl;
            TVector3 antennaPoint(antennaRadius * cos(antennaAngle),
                                  antennaRadius * sin(antennaAngle), 0.0);
            TVector3 antennaDirZ(-1 * sin(antennaAngle), cos(antennaAngle), 0.0);
            TVector3 antennaDirX(cos(antennaAngle), sin(antennaAngle), 0.0);

            HalfWaveDipole *antenna = new HalfWaveDipole(antennaPoint, antennaDirX,
                                                         antennaDirZ, centralFreq,
                                                         timeShift);
            antennaArray.push_back(antenna);
        }

        InducedVoltage iv(trackFile, antennaArray, true);
        iv.GenerateVoltage(0, 5e-8);
        TGraph *grV{iv.GetVoltageGraph()};
        TGraph *grPower{iv.GetPowerPeriodogram(loadResistance)};

        fout->cd();
        grV->Write(Form("grV%d", iAnt));
        grPower->Write(Form("grPower%d", iAnt));

        // Integrate the total power
        double collectedPower{0};
        for (int n{0}; n < grPower->GetN(); n++)
        {
            if (grPower->GetPointX(n) < centralFreq - 2e9 ||
                grPower->GetPointX(n) > centralFreq + 2e9)
                continue;
                
            collectedPower += grPower->GetPointY(n) * 1e15;
        }

        grPowerAntennas->SetPoint(grPowerAntennas->GetN(), iAnt, collectedPower);

        delete grV;
        delete grPower;

        for (auto &ant : antennaArray)
        {
            delete ant;
        }
        antennaArray.clear();
    }
    fout->cd();
    grPowerAntennas->Write("grPowerAntennas");

    fout->Close();
    delete fout;
}