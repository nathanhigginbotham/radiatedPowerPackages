/// ComsolFields.cxx

#include "ElectronDynamics/ComsolFields.h"

#include "TGraph2D.h"
#include "TVector3.h"

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

rad::ComsolField::ComsolField(std::string fieldFile)
{
    fieldValues = new TGraph2D();

    // Need to read in field data first of all
    std::ifstream file110mm(fieldFile);
    if (file110mm.is_open())
    {
        std::string line;
        int nLines{0};
        int ix{1};
        int iy{1};

        // Read in line by line
        while (std::getline(file110mm, line))
        {
            nLines++;
            // Exclude the lines containing the headers and stuff
            if (nLines >= 10)
            {
                std::stringstream ss{line};
                std::string r;
                std::string z;
                std::string magB;
                std::getline(ss, r, ',');
                std::getline(ss, z, ',');
                std::getline(ss, magB, ',');
                double magBDouble{std::stod(magB)};
                double zDouble{std::stod(z)};
                double rDouble{std::stod(r)};
                // Convert from mm to m
                zDouble *= 1e-3;
                rDouble *= 1e-3;
                fieldValues->SetPoint(fieldValues->GetN(),
                                      zDouble, rDouble, magBDouble);
            }
        }
    }
    else
    {
        std::cout << "Unable to open field map file. Exiting." << std::endl;
        exit(1);
    }
}

TVector3 rad::ComsolField::evaluate_field_at_point(const TVector3 vec)
{
    double r{vec.Perp()};
    double z{vec.Z()};
    // Consider B field to only B in z direction
    return TVector3(0, 0, fieldValues->Interpolate(z, r));
}