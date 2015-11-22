#include <iostream>
//~ #include <fstream>
//~ #include <cmath>

//~ #include "TH1.h"
//~ #include "TH2.h"
//~ #include "TH3.h"
//~ #include "TCanvas.h"
//~ #include "TRandom3.h"
//~ #include "TVector2.h"
//~ #include "TVector3.h"
//~ #include "TLorentzVector.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

/* Algorithms list
 * minName              algoName
 * ---------------------------------
 * Minuit/Minuit2       Migrad, Simplex, Scan, Combination  (default is Migrad)
 * Minuit2              Fumili2
 * Fumili
 * GSLMultiMin          ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent
 * GSLMultiFit
 * GSLSimAn
 * Genetic
 */

const int Ndim = 2;

using namespace ROOT::Math;
using namespace std;

double func1(const double* xx)
{
    double x = xx[0];
    double y = xx[1];
    return (x-3)*(x-3)+(y-2)*(y-2);
}

int minimizer()
{
    // Instantiate minimizer with a name/algorithm combo
    Minimizer* min = Factory::CreateMinimizer("Minuit", "");
    
    // The minimizer needs an IMultiGenFunction, which is easily provided
    // by a Functor, which is a generic wrapper class
    Functor fun = Functor(&func1, Ndim);
    
    // Give the function to the minimizer
    min->SetFunction(fun);
    
    // Give the function variables
    min->SetVariable(0, "x", 0, 0.01);
    min->SetVariable(1, "y", 0, 0.01);
    
    // Verbosity
    min->SetPrintLevel(0);
    
    // Minimize!
    min->Minimize();
    
    // Get out the values
    const double* minPoint = min->X();
    const double  minValue = min->MinValue();
    
    cout << "Value: " << minValue << endl;
    cout << "Point: ";
    for (int i=0; i<Ndim-1; ++i) {
        cout << minPoint[i] << " - ";
    }
    cout << minPoint[Ndim-1] << endl;
    
    return 0;
}

int main(int, char**)
{
    return minimizer();
}
