//
//  main.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <iostream>
#include "space2d1l.hpp"
#include "cornermove.hpp"
#include "endmove.hpp"


int main(int argc, const char * argv[])
{
    int nsim = 5;
    int nsumo = 5;
    int lsim = 5;
    int lsumo = 5;
    size_t lx = 20;
    size_t ly = 10;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    test_tube.Initialize();
    cout.flush();
    
    CornerMove cm(test_tube);
    EndMove em(test_tube);
    
    for (int i = 0; i < 10000; i++)
    {
        char typ_r = 0;
        int id_r = 0;
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        cout<<typ_r<<endl;
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        cout<<id_r<<endl;
        cm.ExecMove(id_r, typ_r);
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        cout<<typ_r<<endl;
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        cout<<id_r<<endl;
        em.ExecMove(id_r, typ_r);
    }
        
    
    PrintSpaceOccupation(std::cout, test_tube);
    
    
    return 0;
}
