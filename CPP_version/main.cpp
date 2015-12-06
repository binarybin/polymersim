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
    
    size_t nbr_run = 100000;
    
    App<Space2D1L, Pos2d1l> app(nsim, nsumo, lsim, lsumo, lx, ly);
    app.Initialize();
    for (int i = 0; i < nbr_run; i++)
    {
        app.Proceed('s');
        app.Proceed('e');
        app.Proceed('c');
    }
    
    app.SetBeta(10);
    
    for (int i = 0; i < nbr_run; i++)
    {
        app.Proceed('s');
        app.Proceed('e');
        app.Proceed('c');
    }
    
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

