/*
 * Example code: how to read BK solution and extract dipole amplitude
 *
 * Heikki Mäntysaari <heikki.mantysaari@jyu.fi>, 2019
 */

#include <pybind11/pybind11.h>
#include "amplitudelib.hpp"
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <ctime>
#include <unistd.h>

using namespace std;
namespace py = pybind11;

//void executeDipoleAmplitude(const string& datafile, double xbj) {
std::string DipoleAmplitude(string datafile, double xbj) {    
    // Read data
    AmplitudeLib N(datafile);
    
    // Initialize interpolation for fast evaluation of dipole
    N.InitializeInterpolation(xbj);
    // Create a stringstream to capture the output
    std::stringstream output;
    
    //cout << "# " << N.GetString() << endl;  

    //cout << datafile << endl;

    //const double Ns =  1.0-std::exp(-0.5);
    //cout << "# Q_s^2(x=0.01) = " << N.SaturationScale(0.01, Ns) << " GeV^2" <<endl;
    //cout << "# Q_s^2(x=0.0001) = " << N.SaturationScale(0.0001, Ns) << " GeV^2" << endl;

    // add the evolution parameter C^2 in the beginning:
    output << "# " << N.getC2() << endl;
   
    //cout << "# r [1/GeV]   N(r,x=" << xbj << ")" << endl;
    double minr = N.MinR()*1.01; double maxr=N.MaxR()*0.99;
    for (double r=minr; r<maxr; r*=1.2)
    {
      output << std::scientific << std::setprecision(9) << r << " " << N.N(r, xbj)  << endl;
    }
    
    return output.str();
}



int main(int argc, char* argv[])
{
    string datafile="data/proton/mve.dat";
    
    // Read data
    AmplitudeLib N(datafile);
   
    double xbj = 0.001;
    // Initialize interpolation for fast evaluation of dipole
    N.InitializeInterpolation(xbj);
    cout << "# " << N.GetString() << endl;
    
    const double Ns =  1.0-std::exp(-0.5);
    cout << "# Q_s^2(x=0.01) = " << N.SaturationScale(0.01, Ns) << " GeV^2" <<endl;
    cout << "# Q_s^2(x=0.0001) = " << N.SaturationScale(0.0001, Ns) << " GeV^2" << endl;


   
    cout << "# r [1/GeV]   N(r,x=" << xbj << ")" << endl;
    double minr = N.MinR()*1.01; double maxr=N.MaxR()*0.99;
    for (double r=minr; r<maxr; r*=1.2)
    {
        cout << std::scientific << std::setprecision(9) << r << " " << N.N(r, xbj)  << endl;
    }
   
   
    //double rs = N.SaturationScale(x, Ns);
    

  
    return 0;
}
    
// Define the module and export function
PYBIND11_MODULE(dipole_module, m) {
    m.def("SaturationScale", &AmplitudeLib::SaturationScale);
    m.def("dipole_amplitude", &DipoleAmplitude, "Execute dipole amplitude", py::arg("datafile") = "data/proton/mve.dat", py::arg("xbj") = 0.001);
}
