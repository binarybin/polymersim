//
//  main.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <iostream>
#include "space2d1l.hpp"


int main(int argc, const char * argv[])
{
    int nsim = 50;
    int nsumo = 50;
    int lsim = 5;
    int lsumo = 5;
    size_t lx = 30;
    size_t ly = 30;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    test_tube.Initialize();
    cout.flush();
//    PrintSpaceOccupation(std::cout, test_tube);
    
    
    return 0;
}
