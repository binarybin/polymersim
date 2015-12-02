//
//  app.hpp
//  polymersim
//
//  Created by Farzan Beroz on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef app_hpp
#define app_hpp

#include <iostream>
#include <ctime>
#include "space2d1l.hpp"
#include "space2d2l.hpp"
#include "move.hpp"

template <class S, class P>
class App
{
    S test_tube;
    Move<S, P, SnakeMove<S, P>> sm;
    Move<S, P, CornerMove<S, P>> cm;
    Move<S, P, EndMove<S, P>> em;
public:
    App(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly):
        test_tube(nsim, nsumo, lsim, lsumo, lx, ly), sm(test_tube), cm(test_tube), em(test_tube) {}
    App(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly, size_t lz):
        test_tube(nsim, nsumo, lsim, lsumo, lx, ly, lz), sm(test_tube), cm(test_tube), em(test_tube) {}
    void Initialize() {test_tube.Initialize();}
    
    void Proceed(char typ)
    {
        char typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
        char id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.NSim : test_tube.NSumo));
        switch (typ)
        {
            case 's':
                sm.ExecMove(id_r, typ_r);
                break;
                
            case 'e':
                em.ExecMove(id_r, typ_r);
                break;
            
            case 'c':
                cm.ExecMove(id_r, typ_r);
                break;
                
            default:
                break;
        }
    }
    
};


/*
 void test_longmovespeed()
 {
 cout<<"Test of the polymer simulator"<<endl;
 int nsim = 5;
 int nsumo = 5;
 int lsim = 10;
 int lsumo = 10;
 size_t lx = 30;
 size_t ly = 15;
 Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
 
 cout<<"Test initialization"<<endl;
 test_tube.Initialize();
 
 Move<Space2D1L, Pos2d1l, SnakeMove<Space2D1L, Pos2d1l>> sm(test_tube);
 Move<Space2D1L, Pos2d1l, CornerMove<Space2D1L, Pos2d1l>> cm(test_tube);
 Move<Space2D1L, Pos2d1l, EndMove<Space2D1L, Pos2d1l>> em(test_tube);
 clock_t begin = clock();
 for (int i = 0; i < 1000000; i++)
 {
 char typ_r = 0;
 int id_r = 0;
 
 typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
 id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.NSim : test_tube.NSumo));
 sm.ExecMove(id_r, typ_r);
 //        cout<<"sm"<<endl;
 //        PrintSpace(std::cout, test_tube);
 
 typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
 id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.NSim : test_tube.NSumo));
 cm.ExecMove(id_r, typ_r);
 //        cout<<"cm"<<endl;
 //        PrintSpace(std::cout, test_tube);
 
 typ_r = (double)rand()/RAND_MAX > 0.5 ? 'i' : 'u';
 id_r = 1 + ((double)rand()/RAND_MAX * (typ_r == 'i' ? test_tube.NSim : test_tube.NSumo));
 em.ExecMove(id_r, typ_r);
 //        cout<<"em"<<endl;
 //        PrintSpace(std::cout, test_tube);
 }
 clock_t end = clock();
 double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
 cout<<elapsed_secs<<" seconds"<<endl;
 PrintSpace(std::cout, test_tube);
 //    PrintBond(std::cout, test_tube);
 PrintPolymer(std::cout, test_tube);
 }
 */

/*
 void test_nomove()
 {
 cout<<"Test of the polymer simulator"<<endl;
 int nsim = 30;
 int nsumo = 20;
 int lsim = 30;
 int lsumo = 20;
 size_t lx = 120;
 size_t ly = 100;
 Space2D1L test_tube(nsim, nsumo, lsim, lsumo, lx, ly);
 
 cout<<"Test initialization"<<endl;
 test_tube.Initialize();
 cout.flush();
 
 CornerMove<Pos2d1l, Space2D1L> cm(test_tube);
 EndMove<Pos2d1l, Space2D1L> em(test_tube);
 
 PrintSpace(std::cout, test_tube);
 auto test_neighbor = test_tube.Neighbor(Pos2d1l(1,2));
 cout<<"Test neighbor of (1,2)"<<endl;
 for (auto itr : test_neighbor)
 cout<<"x: "<<itr.x<<" y: "<<itr.y<<endl;
 cout<<"Test SafeRemove (0, 4)"<<endl;
 test_tube.SafeRemove(Pos2d1l(0,4));
 
 
 cout<<"Test SafeCreate (0,4), value -1"<<endl;
 test_tube.SafeCreate(Pos2d1l(0,4), -1);
 PrintSpace(std::cout, test_tube);
 
 cout<<"Test CreateBond (0,4), (0,3)"<<endl;
 test_tube.CreateBond(Pos2d1l(0,4), Pos2d1l(0,3));
 
 
 cout<<"Test SafeRemove (0,4), the case with bond"<<endl;
 test_tube.SafeRemove(Pos2d1l(0,4));
 PrintSpace(std::cout, test_tube);
 
 cout<<"Test CanBuildBond (0,4), (0,3)"<<endl;
 cout<<test_tube.CanBuildBond(Pos2d1l(0,3), Pos2d1l(0,4))<<endl;
 cout<<"Create (0,4)"<<endl;
 test_tube.SafeCreate(Pos2d1l(0,4), -1);
 cout<<test_tube.CanBuildBond(Pos2d1l(0,3), Pos2d1l(0,4))<<endl;
 
 cout<<"Test ExistBond (0,4), (0,3)"<<endl;
 cout<<test_tube.ExistBond(Pos2d1l(0,4), Pos2d1l(0,3))<<endl;
 cout<<"Build it and test"<<endl;
 test_tube.CreateBond(Pos2d1l(0,3), Pos2d1l(0,4));
 cout<<test_tube.ExistBond(Pos2d1l(0,4), Pos2d1l(0,3))<<endl;
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
 
 CornerMove<Pos2d1l, Space2D1L> cm(test_tube);
 EndMove<Pos2d1l, Space2D1L> em(test_tube);
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
 
 SnakeMove<Pos2d1l, Space2D1L> sm(test_tube);
 
 
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
 
 CornerMove<Pos2d1l, Space2D1L> cm(test_tube);
 EndMove<Pos2d1l, Space2D1L> em(test_tube);
 
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
 }*/





#endif /* app_hpp */
