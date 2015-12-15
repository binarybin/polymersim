//
//  main.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <map>
#include "simanneal.hpp"
using std::stoi;

int main(int argc, const char * argv[])
{
    cout<<"argc: "<<argc<<endl;
    int nsim = stoi(argv[1]);
    int nsumo = stoi(argv[2]);
    int lsim = stoi(argv[3]);
    int lsumo = stoi(argv[4]);
    size_t lx = stoi(argv[5]);
    size_t ly = stoi(argv[6]);
    
    size_t nbr_run = stoi(argv[7]);
    
    string signature = argv[8];
    
    SimAnnealing<Space2D1L, Pos2d1l> process(nsim, nsumo, lsim, lsumo, lx, ly, "test");
    
    if (argv[9][0] == 'r')
    {
        int pct = stoi(argv[10]);
        string beta = argv[11];
        process.Resume(pct, beta, nbr_run);
    }
    else
    {
        process.Run(nbr_run);
    }
    return 0;
}

