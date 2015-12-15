//
//  simanneal.hpp
//  polymersim
//
//  Created by Bin Xu on 12/14/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef simanneal_hpp
#define simanneal_hpp

#include <iostream>
#include <fstream>
#include <ctime>
#include "app.hpp"

using std::ofstream;
using std::ifstream;
using std::string;
using std::to_string;

template <class S, class P>
class SimAnnealing
{
    App<S, P> app;
    string signature;
    string raw_filename;
    double beta;
public:
    SimAnnealing(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly, string run_signature) : app(nsim, nsumo, lsim, lsumo, lx, ly), signature(run_signature), beta(0.1)
    {
        raw_filename = "PolymerSim_";
        raw_filename += string("SimAnneal_");
        raw_filename += S::name + string("_");
        raw_filename += signature;
        raw_filename += string("_nsim") + to_string(nsim);
        raw_filename += string("_nsumo") + to_string(nsumo);
        raw_filename += string("_lsim") + to_string(lsim);
        raw_filename += string("_lsumo") + to_string(lsumo);
        raw_filename += string("_lx") + to_string(lx);
        raw_filename += string("_ly") + to_string(ly);
    }
    void Run(size_t rounds)
    {
        app.Initialize();
        ofstream r_out(raw_filename + "_running.txt");
        r_out<<"beta\t"<<"S_succ\t"<<"E_succ\t"<<"C_succ\t"<<"energy"<<endl;
        for (int i = 0; i < 100; i++)
        {
            for (int j = 1; j <= rounds/100; j++)
            {
                app.ResetMoveSucc('s');
                app.ResetMoveSucc('e');
                app.ResetMoveSucc('c');
                
                app.Proceed('s');
                app.Proceed('e');
                app.Proceed('c');
                
                if (j % (rounds/10000) == 0)
                {
                    r_out<<beta<<"\t";
                    r_out<<app.ResetMoveSucc('s')<<"\t";
                    r_out<<app.ResetMoveSucc('e')<<"\t";
                    r_out<<app.ResetMoveSucc('c')<<"\t";
                    r_out<<app.GetEnergy()<<endl;
                }
            }
            
            beta *= 1.05;
            app.SetBeta(beta);
            
            cout<<i+1<<" percent finished"<<endl;
            string filename = string(raw_filename) + string("_status_") + to_string(i+1) + string("_pct_beta_") + to_string(beta) + string(".txt");
            ofstream out(filename);
            app.Dump(out);
            out.close();
        }
        r_out.close();
    }
    
    void SelfTest()
    {
        app.TestPolymerAndSpace();
        app.TestPolymerConnection();
        app.TestReverseSpace();
        app.TestSpaceAndBond();
    }
    
    void Resume(int pct, string betas, size_t rounds)
    {
        app.Initialize();
        ofstream r_out(raw_filename + "_running"+to_string(pct)+".txt");
        r_out<<"beta\t"<<"S_succ\t"<<"E_succ\t"<<"C_succ\t"<<"energy"<<endl;
        string filename = string(raw_filename) + string("_status_") + to_string(pct) + string("_pct_beta_") + betas + string(".txt");
        ifstream in(filename);
        app.Resume(in);
        in.close();
        
        double beta = stof(betas);
        
        for (int i = pct; i < 100; i++)
        {
            for (int j = 1; j <= rounds/100; j++)
            {
                app.ResetMoveSucc('s');
                app.ResetMoveSucc('e');
                app.ResetMoveSucc('c');
                
                app.Proceed('s');
                app.Proceed('e');
                app.Proceed('c');
                
                if (j % (rounds/10000) == 0)
                {
                    r_out<<beta<<"\t";
                    r_out<<app.ResetMoveSucc('s')<<"\t";
                    r_out<<app.ResetMoveSucc('e')<<"\t";
                    r_out<<app.ResetMoveSucc('c')<<"\t";
                    r_out<<app.GetEnergy()<<endl;
                }
            }
            
            beta *= 1.05;
            app.SetBeta(beta);
            
            cout<<i+1<<" percent finished"<<endl;
            string filename = string(raw_filename) + string("_status_") + to_string(i+1) + string("_pct_beta_") + to_string(beta) + string(".txt");
            ofstream out(filename);
            app.Dump(out);
            out.close();
        }
        r_out.close();

    }
    
    
};



#endif /* simanneal_hpp */
