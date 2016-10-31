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
    string initfile;
    vector<tuple<int, double, double, int, int, int>> tasklist;
    vector<char> move_list;
public:
    SimAnnealing(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly, string run_signature, vector<tuple<int, double, double, int, int, int>> thetasklist, string initfile) : app(nsim, nsumo, lsim, lsumo, lx, ly), signature(run_signature), tasklist(thetasklist), move_list({'s', 'e', 'c', 't', 'r', 'T', 'R', 'X', 'Y'}), initfile(initfile)
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
        if (initfile=="empty")
        {
            app.Initialize();
        }
        else
        {
            ifstream fin(initfile);
            app.Resume(fin);
            fin.close();
        }
        
        ofstream r_out(raw_filename + "_running.txt");
        r_out<<"beta\t"<<"gamma\t";
        for (auto move : move_list) r_out<<move<<"_succ\t";
        r_out<<"energy"<<endl;
        for (auto task : tasklist)
        {
            int idx = std::get<0>(task);
            double beta = std::get<1>(task);
            double gamma = std::get<2>(task);
            int runsim = std::get<3>(task);
            int runsumo = std::get<4>(task);
            int phos = std::get<5>(task);
            if (phos == 1)
            {
                app.PhosphorylateOneSiteEpyc();
            }
            else if(phos == -1)
            {
                app.DephosphorylateAllEpyc();
            }
            else if(phos == 5) // a code for changing the epyc length from 5 to 4
            {
                app.ReduceEpycLength();
            }
            app.SetBeta(beta);
            app.SetGamma(gamma);
            
            cout<<"Running task #"<<idx<<" with beta = "<<beta<<" and gamma = "<<gamma<<", "<<runsim<<" sim runs "<<runsumo<<" sumo runs"<<endl;
            
            // reset the success statistics
            for (auto move : move_list) app.ResetMoveSucc(move);
            if(runsim > runsumo)
            {
                for(int i = 1; i <= runsumo; i++)
                    for (auto move : move_list) app.Proceed(move);

                for (int i = 1; i <= runsim - runsumo; i++)
                    for (auto move : {'s', 'e', 'c'}) app.Proceed(move);
            }
            else
            {
                for(int i = 1; i <= runsim; i++)
                    for (auto move : move_list) app.Proceed(move);

                for (int i = 1; i <= runsumo - runsim; i++)
                    for (auto move : {'t', 'r'}) app.Proceed(move);
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
