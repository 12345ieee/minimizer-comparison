#include <iostream>
//~ #include <fstream>
//~ #include <cmath>
#include <ctime>

#include "TH1.h"
//~ #include "TH2.h"
//~ #include "TH3.h"
//~ #include "TCanvas.h"
#include "TRandom3.h"
//~ #include "TVector2.h"
//~ #include "TVector3.h"
//~ #include "TLorentzVector.h"
#include "TString.h"

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
 * Genetic              *** sloooow and inaccurate ***
 * 
 * 
 * https://root.cern.ch/fitting
 * https://root.cern.ch/numerical-minimization
 */

const int Ndim = 8;

using namespace ROOT::Math;
using namespace std;

double parabola1(const double* xx)
{
    double x = xx[0];
    double y = xx[1];
    return (x-3)*(x-3)+(y-2)*(y-2);
}

double parabolaN(const double* xx)
{
    double accumulator=0;
    for (int i=0; i< Ndim; ++i) {
        double c = i+1;
        accumulator += (xx[i]-1/c)*(xx[i]-1/c);
    }
    return accumulator;
}

double rosenbrockN(const double* xx)
{
    double accumulator=0;
    for (int i=0; i< Ndim-1; ++i) {
        accumulator += 100*pow(xx[i+1]-xx[i]*xx[i], 2) + pow(xx[i]-1, 2);
    }
    return accumulator;
}

class MCMinimizer
{
protected:
    TRandom3 rng;
    
    int Ndim = 0;
    double (*function)(const double*); // = &rosenbrockN;
    
    vector<double> params;
    vector<double> params_min;
    vector<double> params_max;
    vector<string> names;
    double minValue;
    
    int max_function_calls = 0;
    int printLevel;
    
public:
    MCMinimizer()
    {
        this->rng = TRandom3(12345);
    }
    
    inline void SetMCFunction(double (*function)(const double*), int Ndim)
    {
        this->function = function;
        this->Ndim = Ndim;
    }
    
    inline void SetMCVariable(string name, double val, double min, double max)
    {
        this->params.push_back(val);
        this->params_min.push_back(min);
        this->params_max.push_back(max);
        this->names.push_back(name);
    }
    
    inline void SetMaxFunctionCalls(int calls)
    {
        this->max_function_calls = calls;
    }
    
    bool Minimize()
    {
        if (printLevel) cout << "Minimize using MCMinimizer" << endl;
        this->minValue = this->function(&params[0]); // evaluate in given initial point
        
        for (int i=0; i<this->max_function_calls; ++i) {
            double pars_array[Ndim];
            for (int par=0; par < Ndim; ++par) {     // get random params vector
                pars_array[par] = rng.Uniform(params_min[par], params_max[par]);
            }
            double nmin = (*function)(pars_array);
            if (nmin < this->minValue) {
                this->minValue = nmin;
                params = vector<double>(pars_array, pars_array + Ndim);
            }
        }
        if (printLevel) {
            cout << "FVAL         = " << this->minValue << endl;
            for (int i=0; i<Ndim; ++i) {
                cout << names[i] << "\t  = " << params[i] << endl;
            }
        }
        return true;
    }
    
    inline void SetPrintLevel(int printLevel)
    {
        this->printLevel = printLevel;
    }
    
    inline double MinValue()
    {
        return this->minValue;
    }
    
    inline double* X()
    {
        return &(this->params[0]);
    }
};

class MCZoomMinimizer
{
protected:
    TRandom3 rng;
    
    int Ndim = 0;
    double (*function)(const double*); // = &rosenbrockN;
    
    vector<double> params;
    vector<double> params_min;
    vector<double> params_max;
    vector<string> names;
    double minValue;
    
    int max_function_calls = 0;
    int splits = 0;
    int printLevel;
    
public:
    MCZoomMinimizer()
    {
        this->rng = TRandom3(12345);
    }
    
    inline void SetMCFunction(double (*function)(const double*), int Ndim)
    {
        this->function = function;
        this->Ndim = Ndim;
    }
    
    inline void SetMCVariable(string name, double val, double min, double max)
    {
        this->params.push_back(val);
        this->params_min.push_back(min);
        this->params_max.push_back(max);
        this->names.push_back(name);
    }
    
    inline void SetMaxFunctionCalls(int calls)
    {
        this->max_function_calls = calls;
    }
    
    inline void SetSplitNumber(int splits)
    {
        this->splits = splits;
    }
    
    bool Minimize()
    {
        if (printLevel) cout << "Minimize using MCZoomMinimizer" << endl;
        this->minValue = this->function(&params[0]); // evaluate in given initial point
        
        for (int ns=0; ns < splits; ++ns) {
            double average[2][Ndim];
            int calls[2][Ndim];
            int call_block = this->max_function_calls/this->splits;
            for (int i=0; i<call_block; ++i) {
                int call_over[Ndim];
                double pars_array[Ndim];
                for (int par=0; par < Ndim; ++par) {     // get random params vector
                    pars_array[par] = rng.Uniform(params_min[par], params_max[par]);
                    int ifup = pars_array[par] > ((params_min[par] + params_max[par])/2);
                    calls[ifup][par]++;
                    call_over[par] = ifup;
                }
                double nmin = (*function)(pars_array);
                for (int par=0; par < Ndim; ++par) {
                    average[call_over[par]][par] += nmin;
                }
                if (nmin < this->minValue) {
                    this->minValue = nmin;
                    params = vector<double>(pars_array, pars_array + Ndim);
                }
            }
            
            for (int par=0; par < Ndim; ++par) {
                if (average[0][par]/calls[0][par] > average[1][par]/calls[1][par]) {
                    params_max[par] = (params_min[par] + params_max[par])/2;
                }
                else {
                    params_min[par] = (params_min[par] + params_max[par])/2;
                }
            }
        }
        
        if (printLevel) {
            cout << "FVAL         = " << this->minValue << endl;
            for (int i=0; i<Ndim; ++i) {
                cout << names[i] << "\t  = " << params[i] << endl;
            }
        }
        return true;
    }
    
    inline void SetPrintLevel(int printLevel)
    {
        this->printLevel = printLevel;
    }
    
    inline double MinValue()
    {
        return this->minValue;
    }
    
    inline double* X()
    {
        return &(this->params[0]);
    }
};

int minimizer()
{
    // List of valid minimizers
    vector<pair<string,string>> minVector;
    //                                      name      algorithm
    // minVector.push_back(pair<string,string>("Minuit", "Migrad"));
    // minVector.push_back(pair<string,string>("Minuit", "Simplex"));
    // minVector.push_back(pair<string,string>("Minuit", "Combined"));
    // minVector.push_back(pair<string,string>("Minuit", "Seek"));
    // minVector.push_back(pair<string,string>("Minuit", "Scan"));
    
    minVector.push_back(pair<string,string>("Minuit2", "Migrad"));
    minVector.push_back(pair<string,string>("Minuit2", "Simplex"));
    // minVector.push_back(pair<string,string>("Minuit2", "Combined"));
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
    
    // minVector.push_back(pair<string,string>("Genetic", ""));

    for (unsigned int i=0; i < minVector.size(); ++i) {
    
        Minimizer* min = Factory::CreateMinimizer(minVector[i].first, minVector[i].second);
        if (min==nullptr) {
            cout << "Invalid algorithm: " << minVector[i].first << " - " << minVector[i].second << endl;
            continue;
        }
        else cout << "Using algorithm: " << minVector[i].first << " - " << minVector[i].second << endl;
        
        // The minimizer needs an IMultiGenFunction, which is easily provided
        // by a Functor, which is a generic wrapper class
        Functor fun = Functor(&rosenbrockN, Ndim);
        
        // Give the function to the minimizer
        min->SetFunction(fun);
        
        // Give the function variables
        if (minVector[i].first=="Genetic") {
            // Limited domain (genetic algorithm only)
            for (int i=0; i<Ndim; ++i) {
                min->SetLimitedVariable(i, Form("x%d", i), 0, 0.05, -0.5, +0.5);
            }
        }
        else {
            for (int i=0; i<Ndim; ++i) {
                min->SetVariable(i, Form("x%d", i), 0, 0.05);
            }
        }
        
        // Verbosity
        min->SetPrintLevel(1);
        
        // Algorithm strategy (higher=slower & more accurate)
        // min->SetStrategy(1);
        
        // Algorithms parameters
        min->SetMaxFunctionCalls(100000);
        min->SetMaxIterations(10000);
        min->SetTolerance(0.001);
        
        // Minimize!
        clock_t time = clock();
        bool success = min->Minimize();
        time = clock() - time;
        cout << "Success: " << success << endl;
        cout << "Time: " << 1000*(double)time/CLOCKS_PER_SEC << " ms " << endl;
        
        // Get out the values
        const double  minValue = min->MinValue();
        // const double* minPoint = min->X();
        // const double* minErrors= min->Errors();
        
        cout << "Value: " << minValue << endl;
        //~ for (int i=0; i<Ndim; ++i) {
            //~ cout << "x" << i+1 << ": " << minPoint[i];
            //~ if (minErrors!=nullptr) cout << " ± " << minErrors[i];
            //~ cout << endl;
        //~ }
        cout << endl;
    }
    
    // Now my minimizer
    {
        MCMinimizer* min = new MCMinimizer();
        cout << "Using algorithm: MCMinimizer" << endl;
        
        // Give the function to the minimizer
        min->SetMCFunction(rosenbrockN, Ndim);
        
        // Give the function variables
        for (int i=0; i<Ndim; ++i) {
            min->SetMCVariable(Form("x%d", i), 0, 0.9, 1.1);
        }
        
        // Verbosity
        min->SetPrintLevel(1);
        
        // Algorithms parameters
        min->SetMaxFunctionCalls(100000);
        
        // Minimize!
        clock_t time = clock();
        bool success = min->Minimize();
        time = clock() - time;
        cout << "Success: " << success << endl;
        cout << "Time: " << 1000*(double)time/CLOCKS_PER_SEC << " ms " << endl;
        
        // Get out the values
        const double  minValue = min->MinValue();
        // const double* minPoint = min->X();
        // const double* minErrors= min->Errors();
        
        cout << "Value: " << minValue << endl;
        //~ for (int i=0; i<Ndim; ++i) {
            //~ cout << "x" << i+1 << ": " << minPoint[i];
            //~ if (minErrors!=nullptr) cout << " ± " << minErrors[i];
            //~ cout << endl;
        //~ }
        cout << endl;
    }
    
    // Now my minimizer
    {
        MCZoomMinimizer* min = new MCZoomMinimizer();
        cout << "Using algorithm: MCMinimizer" << endl;
        
        // Give the function to the minimizer
        min->SetMCFunction(rosenbrockN, Ndim);
        
        // Give the function variables
        for (int i=0; i<Ndim; ++i) {
            min->SetMCVariable(Form("x%d", i), 0, 0.9, 1.1);
        }
        
        // Verbosity
        min->SetPrintLevel(1);
        
        // Algorithms parameters
        min->SetMaxFunctionCalls(100000);
        min->SetSplitNumber(10);
        
        // Minimize!
        clock_t time = clock();
        bool success = min->Minimize();
        time = clock() - time;
        cout << "Success: " << success << endl;
        cout << "Time: " << 1000*(double)time/CLOCKS_PER_SEC << " ms " << endl;
        
        // Get out the values
        const double  minValue = min->MinValue();
        // const double* minPoint = min->X();
        // const double* minErrors= min->Errors();
        
        cout << "Value: " << minValue << endl;
        //~ for (int i=0; i<Ndim; ++i) {
            //~ cout << "x" << i+1 << ": " << minPoint[i];
            //~ if (minErrors!=nullptr) cout << " ± " << minErrors[i];
            //~ cout << endl;
        //~ }
        cout << endl;
    }
    
    return 0;
}

int main(int, char**)
{
    return minimizer();
}
