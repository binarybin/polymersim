//
//  main.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <ctime>
#include "app.hpp"

using std::ofstream;
using std::ifstream;

int main(int argc, const char * argv[])
{
    int nsim = 25;
    int nsumo = 25;
    int lsim = 10;
    int lsumo = 10;
    size_t lx = 100;
    size_t ly = 100;
    
    size_t nbr_run = 10000000;
    
    App<Space2D2L, Pos2d2l> app(nsim, nsumo, lsim, lsumo, lx, ly);
    app.Initialize();
    double beta = 0.1;
    ofstream r_out("runningdata.txt");
    r_out<<"beta\t"<<"S_succ\t"<<"E_succ\t"<<"C_succ\t"<<"energy"<<endl;
    for (int i = 0; i < 100; i++)
    {
        for (int j = 1; j <= nbr_run/100; j++)
        {
            app.Proceed('s');
            app.Proceed('e');
            app.Proceed('c');
            if (j % 10000 == 0)
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
        std::string filename = std::string("testout") + std::to_string(i+1) + std::string("pct.txt");
        ofstream out(filename);
        app.Dump(out);
        out.close();
    }
    r_out.close();
    
    app.ShowPolymer(std::cout);
    app.ShowSpace(std::cout);
    
    ofstream out("testout.txt");
    app.Dump(out);
    out.close();
    
/*    ifstream in("/Users/binxu/Desktop/testout.txt");
    app.Resume(in);
    in.close();*/
    
/*    for (int i = 0; i < nbr_run; i++)
    {
        app.Proceed('s');
        app.Proceed('e');
        app.Proceed('c');
    }*/
    
    app.TestPolymerAndSpace();
    app.TestPolymerConnection();
    app.TestReverseSpace();
    app.TestSpaceAndBond();
    
    return 0;
}

