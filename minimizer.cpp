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
 * Minuit/Minuit2       Migrad, Simplex, Scan, Seek, Combined  (default is Migrad)
 * Minuit2              Fumili
 * Fumili
 * GSLMultiMin          ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent
 * GSLMultiFit
 * GSLSimAn
 * Linear
 * Genetic
 * 
 * (Fumili needs gradient)
 * 
 * https://root.cern.ch/fitting
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
    // Instantiate minimizer with a name/algorithm combo (case sensitive!)
    Minimizer* min = Factory::CreateMinimizer("Minuit", "Migrad");
    //Minimizer* min = Factory::CreateMinimizer("Minuit", "Simplex");
    //Minimizer* min = Factory::CreateMinimizer("Minuit", "Combined");
    //Minimizer* min = Factory::CreateMinimizer("Minuit", "Scan");
    //Minimizer* min = Factory::CreateMinimizer("Minuit", "Seek");
    
    //Minimizer* min = Factory::CreateMinimizer("Minuit2", "Migrad");
    //Minimizer* min = Factory::CreateMinimizer("Minuit2", "Simplex");
    //Minimizer* min = Factory::CreateMinimizer("Minuit2", "Combined");
    //Minimizer* min = Factory::CreateMinimizer("Minuit2", "Scan");
    //Minimizer* min = Factory::CreateMinimizer("Minuit", "Seek");
    //Minimizer* min = Factory::CreateMinimizer("Minuit2", "Fumili");
    
    //Minimizer* min = Factory::CreateMinimizer("GSLMultiMin", "ConjugateFR");
    //Minimizer* min = Factory::CreateMinimizer("GSLMultiMin", "ConjugatePR");
    //Minimizer* min = Factory::CreateMinimizer("GSLMultiMin", "BFGS");
    //Minimizer* min = Factory::CreateMinimizer("GSLMultiMin", "BFGS2");
    //Minimizer* min = Factory::CreateMinimizer("GSLMultiMin", "SteepestDescent");
    
    //Minimizer* min = Factory::CreateMinimizer("GSLMultiFit", "");
    
    //Minimizer* min = Factory::CreateMinimizer("GSLSimAn", "");

    if (min==nullptr) {
        cout << "Invalid algorithm name" << endl;
        return 1;
    }
    
    // The minimizer needs an IMultiGenFunction, which is easily provided
    // by a Functor, which is a generic wrapper class
    Functor fun = Functor(&func1, Ndim);
    
    // Give the function to the minimizer
    min->SetFunction(fun);
    
    // Give the function variables
    min->SetVariable(0, "x", 0, 0.01);
    min->SetVariable(1, "y", 0, 0.01);
    
    // Verbosity
    min->SetPrintLevel(5);
    
    // Algorithm strategy (higher=slower & more accurate)
    // min->SetStrategy(1);
    
    // Algorithms parameters
    min->SetTolerance(0.001); // For simplex
    
    // Minimize!
    min->Minimize();
    
    // Get out the values
    const double  minValue = min->MinValue();
    const double* minPoint = min->X();
    const double* minErrors= min->Errors();
    
    cout << "\nValue: " << minValue << endl;
    for (int i=0; i<Ndim; ++i) {
        cout << "x" << i+1 << ": " << minPoint[i];
        if (minErrors!=nullptr) cout << " ± " << minErrors[i];
        cout << endl;
    }
    cout << endl;
    
    return 0;
}

int main(int, char**)
{
    return minimizer();
}
