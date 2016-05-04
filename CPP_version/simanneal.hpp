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
    vector<tuple<int, double, double, int>> tasklist;
public:
    SimAnnealing(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly, string run_signature, vector<tuple<int, double, double, int>> thetasklist) : app(nsim, nsumo, lsim, lsumo, lx, ly), signature(run_signature), tasklist(thetasklist)
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
    
    void Run()
    {
        app.Initialize();
        ofstream r_out(raw_filename + "_running.txt");
        r_out<<"beta\t"<<"gamma\t"<<"S_succ\t"<<"E_succ\t"<<"C_succ\t"<<"T_succ\t"<<"R_succ\t"<<"CT_succ\t"<<"CR_succ\t"<<"energy"<<endl;
        for (auto task : tasklist)
        {
            int idx = std::get<0>(task);
            double beta = std::get<1>(task);
            double gamma = std::get<2>(task);
            int runs = std::get<3>(task);
            
            app.SetBeta(beta);
            app.SetGamma(gamma);
            
            cout<<"Running task #"<<idx<<" with beta = "<<beta<<" and gamma = "<<gamma<<", "<<runs<<" runs"<<endl;
            
            // reset the success statistics
            app.ResetMoveSucc('s');
            app.ResetMoveSucc('e');
            app.ResetMoveSucc('c');
            app.ResetMoveSucc('t');
            app.ResetMoveSucc('r');
            app.ResetMoveSucc('T');
            app.ResetMoveSucc('R');
            for(int i = 1; i < runs+1; i++)
            {
                app.Proceed('s');
                app.Proceed('e');
                app.Proceed('c');
                app.Proceed('t');
                app.Proceed('r');
                if (app.test_tube.LSim*2 <= app.test_tube.LSumo)
                {
              //      app.Proceed('T');
              //      app.Proceed('R');
                }
                
                
                if (i % (runs/10000) == 0)
                {
                    r_out<<beta<<"\t"<<gamma<<"\t";
                    r_out<<app.ResetMoveSucc('s')<<"\t";
                    r_out<<app.ResetMoveSucc('e')<<"\t";
                    r_out<<app.ResetMoveSucc('c')<<"\t";
                    r_out<<app.ResetMoveSucc('t')<<"\t";
                    r_out<<app.ResetMoveSucc('r')<<"\t";
                    r_out<<app.ResetMoveSucc('T')<<"\t";
                    r_out<<app.ResetMoveSucc('R')<<"\t";
                    r_out<<app.GetEnergy()<<endl;
                }
            }
            string filename = string(raw_filename) + string("_status_step_") + to_string(idx) + string("_beta_") + to_string(beta) + string("_gamma_") + to_string(gamma) + string(".txt");
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
    
};



#endif /* simanneal_hpp */
