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
 * Minuit/Minuit2       Migrad, Simplex, Seek, Combined (default is Migrad)
 * Minuit/Minuit2       Scan (good for scans, not minimizations)
 * Minuit2              Fumili (*** needs different function type ***), Fumili2
 * Fumili               *** needs different function type ***
 * GSLMultiMin          ConjugateFR, ConjugatePR, BFGS, BFGS2, SteepestDescent
 * GSLMultiFit          *** needs different function type ***
 * GSLSimAn
 * Linear               *** needs linear function ***
 * Genetic
 * 
 * 
 * https://root.cern.ch/fitting
 * https://root.cern.ch/numerical-minimization
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
    // List of valid minimizers
    vector<pair<string,string>> minVector;
    //                                      name      algorithm
    minVector.push_back(pair<string,string>("Minuit", "Migrad"));
    minVector.push_back(pair<string,string>("Minuit", "Simplex"));
    minVector.push_back(pair<string,string>("Minuit", "Combined"));
    minVector.push_back(pair<string,string>("Minuit", "Seek"));
    // minVector.push_back(pair<string,string>("Minuit", "Scan"));
    
    minVector.push_back(pair<string,string>("Minuit2", "Migrad"));
    minVector.push_back(pair<string,string>("Minuit2", "Simplex"));
    minVector.push_back(pair<string,string>("Minuit2", "Combined"));
    minVector.push_back(pair<string,string>("Minuit2", "Seek"));
    // minVector.push_back(pair<string,string>("Minuit2", "Scan"));
    // minVector.push_back(pair<string,string>("Minuit2", "Fumili"));
    minVector.push_back(pair<string,string>("Minuit2", "Fumili2"));
    
    // minVector.push_back(pair<string,string>("Fumili",  "Fumili"));
    
    minVector.push_back(pair<string,string>("GSLMultiMin", "ConjugateFR"));
    minVector.push_back(pair<string,string>("GSLMultiMin", "ConjugatePR"));
    minVector.push_back(pair<string,string>("GSLMultiMin", "BFGS"));
    minVector.push_back(pair<string,string>("GSLMultiMin", "BFGS2"));
    minVector.push_back(pair<string,string>("GSLMultiMin", "SteepestDescent"));
    
    // minVector.push_back(pair<string,string>("GSLMultiFit", ""));
    
    minVector.push_back(pair<string,string>("GSLSimAn", ""));

    // minVector.push_back(pair<string,string>("Linear", "");
    
    minVector.push_back(pair<string,string>("Genetic", ""));

    for (unsigned int i=0; i < minVector.size(); ++i) {
    
        Minimizer* min = Factory::CreateMinimizer(minVector[i].first, minVector[i].second);
        if (min==nullptr) {
            cout << "Invalid algorithm: " << minVector[i].first << " - " << minVector[i].second << endl;
            continue;
        }
        else cout << "Using algorithm: " << minVector[i].first << " - " << minVector[i].second << endl;
        
        // The minimizer needs an IMultiGenFunction, which is easily provided
        // by a Functor, which is a generic wrapper class
        Functor fun = Functor(&func1, Ndim);
        
        // Give the function to the minimizer
        min->SetFunction(fun);
        
        // Give the function variables
        min->SetVariable(0, "x", 0, 0.1);
        min->SetVariable(1, "y", 0, 0.1);
        
        // Verbosity
        min->SetPrintLevel(1);
        
        // Algorithm strategy (higher=slower & more accurate)
        // min->SetStrategy(1);
        
        // Algorithms parameters
        min->SetMaxFunctionCalls(100000);
        min->SetMaxIterations(10000);
        min->SetTolerance(0.001);
        
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
    }
    
    return 0;
}

int main(int, char**)
{
    return minimizer();
}
