//
//  main.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <iostream>
#include <ctime>
#include "space2d1l.hpp"
#include "cornermove.hpp"
#include "endmove.hpp"
#include "snakemove.hpp"

void test_nomove();
void test_move();
void test_longmove();
void test_snakemove();
void test_longmovespeed();
int main(int argc, const char * argv[])
{
    test_longmovespeed();
    return 0;
}

void test_nomove()
{
    cout<<"Test of the polymer simulator"<<endl;
    int nsim = 5;
    int nsumo = 5;
    int lsim = 5;
    int lsumo = 5;
    size_t lx = 20;
    size_t ly = 10;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    
    cout<<"Test initialization"<<endl;
    test_tube.Initialize();
    cout.flush();
    
    CornerMove cm(test_tube);
    EndMove em(test_tube);
    
    PrintSpace(std::cout, test_tube);
    auto test_neighbor = test_tube.Neighbor(Position(1,2));
    cout<<"Test neighbor of (1,2)"<<endl;
    for (auto itr : test_neighbor)
        cout<<"x: "<<itr.x<<" y: "<<itr.y<<endl;
    cout<<"Test SafeRemove (0, 4)"<<endl;
    test_tube.SafeRemove(Position(0,4));
    
    
    cout<<"Test SafeCreate (0,4), value -1"<<endl;
    test_tube.SafeCreate(Position(0,4), -1);
    PrintSpace(std::cout, test_tube);
    
    cout<<"Test CreateBond (0,4), (0,3)"<<endl;
    test_tube.CreateBond(Position(0,4), Position(0,3));
    
    
    cout<<"Test SafeRemove (0,4), the case with bond"<<endl;
    test_tube.SafeRemove(Position(0,4));
    PrintSpace(std::cout, test_tube);
    
    cout<<"Test CanBuildBond (0,4), (0,3)"<<endl;
    cout<<test_tube.CanBuildBond(Position(0,3), Position(0,4))<<endl;
    cout<<"Create (0,4)"<<endl;
    test_tube.SafeCreate(Position(0,4), -1);
    cout<<test_tube.CanBuildBond(Position(0,3), Position(0,4))<<endl;
    
    cout<<"Test ExistBond (0,4), (0,3)"<<endl;
    cout<<test_tube.ExistBond(Position(0,4), Position(0,3))<<endl;
    cout<<"Build it and test"<<endl;
    test_tube.CreateBond(Position(0,3), Position(0,4));
    cout<<test_tube.ExistBond(Position(0,4), Position(0,3))<<endl;
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
}


void test_move()
{
    cout<<"Test of the polymer simulator"<<endl;
    int nsim = 5;
    int nsumo = 5;
    int lsim = 5;
    int lsumo = 5;
    size_t lx = 20;
    size_t ly = 10;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    
    cout<<"Test initialization"<<endl;
    test_tube.Initialize();
    cout.flush();
    
    CornerMove cm(test_tube);
    EndMove em(test_tube);
    cout<<"Redo initialization"<<endl;
    test_tube.Initialize();
    
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    cout<<"Test endmove"<<endl;
    em.ExecMove(1, 'i');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    cout<<"Test cornermove"<<endl;
    cm.ExecMove(1, 'i');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    cout<<"Test endmove"<<endl;
    em.ExecMove(5, 'i');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    cout<<"Test endmove"<<endl;
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'u');
    em.ExecMove(5, 'i');
    em.ExecMove(5, 'i');
    em.ExecMove(5, 'i');
    em.SetBeta(10);
    cm.SetBeta(10);
    cm.ExecMove(5, 'i');
    cm.ExecMove(5, 'u');
    cm.ExecMove(5, 'i');
    cm.ExecMove(5, 'u');
    cm.ExecMove(5, 'i');
    cm.ExecMove(5, 'u');
//    PrintSpace(std::cout, test_tube);
    cm.ExecMove(5, 'i');
    cm.ExecMove(5, 'u');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
}

void test_snakemove()
{
    cout<<"Test of the polymer simulator"<<endl;
    int nsim = 5;
    int nsumo = 5;
    int lsim = 5;
    int lsumo = 5;
    size_t lx = 20;
    size_t ly = 10;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    
    cout<<"Test initialization"<<endl;
    test_tube.Initialize();
    cout.flush();
    
    SnakeMove sm(test_tube);

    
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    sm.ExecMove(1, 'i');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    sm.ExecMove(1, 'i');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    sm.ExecMove(5, 'i');
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    
    sm.ExecMove(5, 'u');
    sm.ExecMove(5, 'u');
    sm.ExecMove(5, 'u');
    sm.ExecMove(5, 'u');
    sm.ExecMove(5, 'u');

    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    PrintPolymer(std::cout, test_tube);
}


void test_longmove()
{
    cout<<"Test of the polymer simulator"<<endl;
    int nsim = 5;
    int nsumo = 5;
    int lsim = 5;
    int lsumo = 5;
    size_t lx = 20;
    size_t ly = 10;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    
    cout<<"Test initialization"<<endl;
    test_tube.Initialize();
    cout.flush();
    
    CornerMove cm(test_tube);
    EndMove em(test_tube);
    
    for (int i = 0; i < 10000; i++)
    {
        char typ_r = 0;
        int id_r = 0;
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        cm.ExecMove(id_r, typ_r);
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        em.ExecMove(id_r, typ_r);
    }
    PrintSpace(std::cout, test_tube);
    PrintBond(std::cout, test_tube);
    PrintPolymer(std::cout, test_tube);
}

void test_longmovespeed()
{
    cout<<"Test of the polymer simulator"<<endl;
    int nsim = 150;
    int nsumo = 250;
    int lsim = 20;
    int lsumo = 20;
    size_t lx = 250;
    size_t ly = 250;
    Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
    
    cout<<"Test initialization"<<endl;
    test_tube.Initialize();
    
    SnakeMove sm(test_tube);
    EndMove em(test_tube);
    CornerMove cm(test_tube);
    
    for (int i = 0; i < 10000; i++)
    {
        char typ_r = 0;
        int id_r = 0;
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        sm.ExecMove(id_r, typ_r);
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        cm.ExecMove(id_r, typ_r);
        
        typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.LSim : test_tube.LSumo));
        em.ExecMove(id_r, typ_r);
    }
//    PrintSpace(std::cout, test_tube);
//    PrintBond(std::cout, test_tube);
//    PrintPolymer(std::cout, test_tube);
}

