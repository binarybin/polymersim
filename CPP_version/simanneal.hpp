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

typedef tuple<int, double, double, double, int, int, int> task_t;

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
    vector<task_t> tasklist;
    vector<char> move_list;
public:
    SimAnnealing(int nsim1, int nsim2, int nsumo, int lsim1, int lsim2, int lsumo, size_t lx, size_t ly, string run_signature, vector<tuple<int, double, double, double, int, int, int>> thetasklist, string initfile) : app(nsim1, nsim2, nsumo, lsim1, lsim2, lsumo, lx, ly), signature(run_signature), tasklist(thetasklist), move_list({'s', 'e', 'c', 't', 'r', 'T', 'R', 'X', 'Y'}), initfile(initfile)
    {
        raw_filename = "PolymerSim_";
#ifdef NO_TWO_END
        raw_filename += string("NoTwo_");
#endif
        raw_filename += string("SimAnneal_");
        raw_filename += S::name + string("_");
        raw_filename += signature;
        raw_filename += string("_nsim1_") + to_string(nsim1);
        raw_filename += string("_nsim2_") + to_string(nsim2);
        raw_filename += string("_nsumo_") + to_string(nsumo);
        raw_filename += string("_lsim1_") + to_string(lsim1);
        raw_filename += string("_lsim2_") + to_string(lsim2);
        raw_filename += string("_lsumo_") + to_string(lsumo);
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
            double gamma_intra = std::get<2>(task);
            double gamma_inter = std::get<3>(task);
            int runsim = std::get<4>(task);
            int runsumo = std::get<5>(task);
            int phos = std::get<6>(task);
            if (phos == 1)
            {
                app.PhosphorylateOneSiteEpyc();
            }
            else if(phos == -1)
            {
                app.DephosphorylateAllEpyc();
            }
            else if(phos == 5) // a code for reducing a phos number of epyc length by 1
            {
                app.ReduceEpycLength();
            }
            app.SetBeta(beta);
            app.SetGammaIntra(gamma_intra);
            app.SetGammaInter(gamma_inter);
            
            cout<<"Running task #"<<idx<<" with beta = "<<beta<<" and gamma_intra = "<<gamma_intra<< " and gamma_inter = "<<gamma_inter<<", "<<runsim<<" sim runs "<<runsumo<<" sumo runs"<<endl;
            
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
            
            string filename = string(raw_filename) + string("_status_step_") + to_string(idx) + string("_beta_") + to_string(beta) + string("_gammaintra_") + to_string(gamma_intra) + string("_gammainter_") + to_string(gamma_inter) + string(".txt");
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
