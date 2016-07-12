//
//  app.hpp
//  polymersim
//
//  Created by Farzan Beroz (Who is that?) on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef app_hpp
#define app_hpp

#include <iostream>
#include <sstream>
#include <ctime>
#include <tuple>
#include <string>
#include <stdexcept>
#include <random>
#include "space2d2l.hpp"
#include "move.hpp"
#include "rubimoves.hpp"


using std::make_tuple;
using std::tuple;
using std::string;
using std::istringstream;
using std::runtime_error;
using std::random_device;
using std::mt19937;
using std::uniform_int_distribution;
using std::generate_canonical;


template <class S, class P>
class App
{
    Move<S, P, SnakeMove<S, P>> sm;
    Move<S, P, CornerMove<S, P>> cm;
    Move<S, P, EndMove<S, P>> em;
    RubiMove<S, P, TranslationMove<S, P>> tm;
    RubiMove<S, P, RotationMove<S, P>> rm;
    RubiMove<S, P, TranslationMove<S, P>> Tm;
    RubiMove<S, P, RotationMove<S, P>> Rm;
    RubiMove<S, P, TranslationMove<S, P>> Xm;
    RubiMove<S, P, RotationMove<S, P>> Ym;
    random_device rd;
    mt19937 gen;
    int original_nbr_bond;
public:
    S test_tube;
    App(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly):
        test_tube(nsim, nsumo, lsim, lsumo, lx, ly),
        sm(test_tube), cm(test_tube), em(test_tube),
        tm(test_tube), rm(test_tube), Tm(test_tube), Rm(test_tube),
        Xm(test_tube), Ym(test_tube),
        gen(rd()), original_nbr_bond(0) {}
    
    void Initialize() {test_tube.Initialize();}
    
    void Proceed(char typ);
    
    void PhosphorylateOneSiteEpyc()
    {
        uniform_int_distribution<> dis(0, test_tube.LSim-1);
        for (int i = 0; i < test_tube.NSim; i++)
        {
            int id_r = dis(gen);
            test_tube.Phosphorylate(i, true, id_r);
            cout<<"Phosphorylate EPYC "<<i<<" site "<<id_r<<endl;
        }
    }
    void DephosphorylateAllEpyc()
    {
        for (int i = 0; i < test_tube.NSim; i++)
        {
            for (int j = 0; j < test_tube.LSim; j++)
            {
                if(test_tube.Phosphorylated(test_tube.Sims[i].locs[j]))
                    test_tube.Dephosphorylate(i, true, j);
            }
        }
    }
    void ReduceEpycLength()
    {
        test_tube.ReduceEpycLength();
    }
    
    int ResetMoveSucc(char typ)
    {
        int returnvalue;
        switch (typ)
        {
            case 's':
                returnvalue = sm.GetSucc();
                sm.ClearSucc();
                break;
                
            case 'e':
                returnvalue = em.GetSucc();
                em.ClearSucc();
                break;
                
            case 'c':
                returnvalue = cm.GetSucc();
                cm.ClearSucc();
                break;
                
            case 't':
                returnvalue = tm.GetSucc();
                tm.ClearSucc();
                break;
                
            case 'r':
                returnvalue = rm.GetSucc();
                rm.ClearSucc();
                break;
            
            case 'T':
                returnvalue = Tm.GetSucc();
                Tm.ClearSucc();
                break;
                
            case 'R':
                returnvalue = Rm.GetSucc();
                Rm.ClearSucc();
                break;
            case 'X':
                returnvalue = Xm.GetSucc();
                Xm.ClearSucc();
                break;
                
            case 'Y':
                returnvalue = Ym.GetSucc();
                Ym.ClearSucc();
                break;

            default:
                break;
        }
        return returnvalue;
    }
    
    double GetEnergy(){return -(original_nbr_bond + sm.bond_change + em.bond_change + cm.bond_change + tm.bond_change + rm.bond_change);}
    void ShowSpace(std::ostream &out) {PrintSpace(out, test_tube); }
    void ShowBond(std::ostream &out) {PrintBond(out, test_tube); }
    void ShowPolymer(std::ostream &out) {PrintPolymer(out, test_tube); }

    void SetBeta(double beta)
    {
        sm.SetBeta(beta);
        cm.SetBeta(beta);
        em.SetBeta(beta);
        tm.SetBeta(beta);
        rm.SetBeta(beta);
        Tm.SetBeta(beta);
        Rm.SetBeta(beta);
        Xm.SetBeta(beta);
        Ym.SetBeta(beta);
        
    }
    
    void SetGamma(double gamma)
    {
        sm.SetGamma(gamma);
        cm.SetGamma(gamma);
        em.SetGamma(gamma);
        tm.SetGamma(gamma);
        rm.SetGamma(gamma);
        Tm.SetGamma(gamma);
        Rm.SetGamma(gamma);
        Xm.SetGamma(gamma);
        Ym.SetGamma(gamma);
    }
    
    void Dump(std::ostream &out)
    {
        std::cout<<"Dumping current result to file"<<endl;
        out<<"[Spacelist]"<<endl;
        ShowSpace(out);
        out<<"[End]"<<endl;
        out<<"[Polymerlist]"<<endl;
        ShowPolymer(out);
        out<<"[End]"<<endl;
        out<<"[Bondlist]"<<endl;
        ShowBond(out);
        out<<"[End]"<<endl;
        std::cout<<"Finished dumping"<<endl;
    }
    void Resume(std::ifstream &in);
    vector<vector<string>> GenericSpaceAnal(vector<string>& spacebuf);
    vector<vector<string>> GenericPolyAnal(vector<string>& polybuf);
    vector<vector<string>> GenericBondAnal(vector<string>& bondbuf);
    
    void TestPolymerAndSpace();
    void TestPolymerConnection();
    void TestSpaceAndBond();
    void TestReverseSpace();
};

template <>
void App<Space2D2L, Pos2d2l>::TestPolymerAndSpace()
{
    auto tmpspace(test_tube.space);
    for (auto& poly : test_tube.Sims)
    {
        for (auto& loc : poly.locs)
        {
            int x = loc.x;
            int y = loc.y;
            int layer = loc.siml?0:1;
            assert(loc.siml);
            if (tmpspace[layer][x][y] != 1 && tmpspace[layer][x][y] != 2 && tmpspace[layer][x][y] != 3)
            {
                throw runtime_error("Polymer and space test failed for SIM");
            }
            tmpspace[layer][x][y] = 0;
        }
    }
    
    for (auto& poly : test_tube.Sumos)
    {
        for (auto& loc : poly.locs)
        {
            int x = loc.x;
            int y = loc.y;
            int layer = loc.siml?0:1;
            assert(!loc.siml);
            if (tmpspace[layer][x][y] != -1 && tmpspace[layer][x][y] != -2 && tmpspace[layer][x][y] != -3)
            {
                throw runtime_error("Polymer and space test failed for SUMO");
            }
            tmpspace[layer][x][y] = 0;
        }
    }
    
    for (int i = 0; i < tmpspace.size(); i++)
    {
        for (int j = 0; j < tmpspace[0].size(); j++)
        {
            if (tmpspace[0][i][j] != 0 || tmpspace[1][i][j] != 0)
            {
                throw runtime_error("Polymer and space test failed for space");
            }
        }
    }
    
    bool fail_deep_copy = true;
    for (int i = 0; i < tmpspace.size(); i++)
    {
        for (int j = 0; j < tmpspace[0].size(); j++)
        {
            if (tmpspace[0][i][j] != test_tube.space[0][i][j] || tmpspace[1][i][j] != test_tube.space[1][i][j])
            {
                fail_deep_copy = false;
            }
        }
    }
    
    if (fail_deep_copy)
    {
 //       throw runtime_error("Deep copy does not seem to work");
    }
    
    cout<<"Polymer and space test passed!"<<endl;
}

template <>
void App<Space2D2L, Pos2d2l>::TestSpaceAndBond()
{
    for (int i = 0; i < test_tube.space.size(); i++)
    {
        for (int j = 0;  j < test_tube.space[0].size(); j++)
        {
            if (test_tube.space[0][i][j] == 2 && test_tube.space[1][i][j] != -2)
            {
                throw runtime_error("Unmatched bond in the first layer");
            }
            else if (test_tube.space[0][i][j] != 2 && test_tube.space[1][i][j] == -2)
            {
                throw runtime_error("Unmatched bond in the second layer");
            }
        }
    }
    
    cout<<"Bond and space test passed!"<<endl;
}

template <>
void App<Space2D2L, Pos2d2l>::TestReverseSpace()
{
    auto tempspace(test_tube.space);
    for (int i = 0; i < test_tube.rspace.size(); i++)
    {
        for (int j = 0; j < test_tube.rspace[0].size(); j++)
        {
            if (test_tube.rspace[0][i][j][0] == NOBOND && test_tube.rspace[0][i][j][1] == NOBOND)
                continue;
            
            int id = test_tube.rspace[0][i][j][0];
            int pos = test_tube.rspace[0][i][j][1];
            Pos2d2l loc = test_tube.Sims[id].locs[pos];
            int rx = loc.x;
            int ry = loc.y;
            if (rx != i || ry != j)
            {
                throw runtime_error("Reverse space inconsistent for SIM");
            }
            tempspace[0][i][j] = 0;
            
        }
    }
    
    for (int i = 0; i < test_tube.rspace.size(); i++)
    {
        for (int j = 0; j < test_tube.rspace[0].size(); j++)
        {
            if (test_tube.rspace[1][i][j][0] == NOBOND && test_tube.rspace[1][i][j][1] == NOBOND)
                continue;
            
            int id = test_tube.rspace[1][i][j][0];
            int pos = test_tube.rspace[1][i][j][1];
            Pos2d2l loc = test_tube.Sumos[id].locs[pos];
            int rx = loc.x;
            int ry = loc.y;
            if (rx != i || ry != j)
            {
                throw runtime_error("Reverse space inconsistent for SUMO");
            }
            
            tempspace[1][i][j] = 0;
        }
    }
    
    for (int i = 0; i < tempspace.size(); i++)
    {
        for (int j = 0; j < tempspace[0].size(); j++)
        {
            if (tempspace[0][i][j] != 0 || tempspace[1][i][j] != 0)
            {
                throw runtime_error("Reverse space check failed, more monomers in space than in bond");
            }
        }
    }

    
   // cout<<"Reverse space test passed!"<<endl;
}

template <class S, class P>
void App<S,P>::TestPolymerConnection()
{
    for (auto& poly : test_tube.Sims)
        for (int pos = 0; pos < poly.locs.size()-1; pos++)
            if (! poly.locs[pos].IsNeighborOf(poly.locs[pos+1]))
                throw runtime_error("Polymer connection test failed for SIM");
    
    for (auto& poly : test_tube.Sumos)
        for (int pos = 0; pos < poly.locs.size()-1; pos++)
            if (! poly.locs[pos].IsNeighborOf(poly.locs[pos+1]))
                throw runtime_error("Polymer connection test failed for SUMO");
    cout<<"Polymer connection test passed!"<<endl;

}

template <class S, class P>
void App<S,P>::Proceed(char typ)
{
    int NPoly = ( (typ == 's' || typ == 'e' || typ == 'c') ? test_tube.NSim : test_tube.NSumo);
    uniform_int_distribution<> dis(0, NPoly-1);
    int id_r = dis(gen);
    switch (typ)
    {
        case 's':
            sm.ExecMove(id_r, 'i');
            break;
            
        case 'e':
            em.ExecMove(id_r, 'i');
            break;
            
        case 'c':
            cm.ExecMove(id_r, 'i');
            break;
            
        case 't':
            tm.ExecMove(id_r, 'u');
            break;
            
        case 'r':
            rm.ExecMove(id_r, 'u');
            break;
            
        case 'T':
            Tm.ExecDragMove(id_r, false);
            break;
        case 'R':
            Rm.ExecDragMove(id_r, false);
            break;
        case 'X':
            Xm.ExecDragMove(id_r, true);
        case 'Y':
            Ym.ExecDragMove(id_r, true);
            
        default:
            break;
    }
}

template <class S, class P>
void App<S,P>::Resume(std::ifstream &in)
{
    vector<string> spacebuf;
    vector<string> polybuf;
    vector<string> bondbuf;
    
    string tempstring;
    if (in.is_open())
    {
        vector<string> *working_buf = &spacebuf;
        while (std::getline(in, tempstring))
        {
            if (tempstring.compare("[Polymerlist]") == 0)
            {
                working_buf = &polybuf;
            }
            else if (tempstring.compare("[Bondlist]") == 0)
            {
                working_buf = &bondbuf;
            }
            
            (*working_buf).push_back(tempstring);
        }
    }
    
    vector<vector<string>> spaceanal = GenericSpaceAnal(spacebuf);
    vector<vector<string>> polyanal  = GenericPolyAnal(polybuf);
    vector<vector<string>> bondanal  = GenericBondAnal(bondbuf);
    cout<<"Finished!"<<endl;
    test_tube.Resume(spaceanal, polyanal, bondanal);
}

template <class S, class P>
vector<vector<string>> App<S, P>::GenericSpaceAnal(vector<string>& spacebuf)
{
    spacebuf.pop_back();
    spacebuf.erase(spacebuf.begin());
    assert(spacebuf[0].compare("Two dimensional single layer space") == 0 || spacebuf[0].compare("Two dimensional double layer space") == 0);
    assert(spacebuf[1].compare("Site occupation configuration: ") == 0);
    assert(spacebuf[spacebuf.size()-1].compare("END") == 0);

    spacebuf.erase(spacebuf.begin(), spacebuf.begin()+2);
    spacebuf.pop_back();
    char delim = '\t';
    vector<vector<string>> result;
    for (auto& ss : spacebuf)
    {
        istringstream split(ss);
        vector<string> tokens;
        for (std::string each; std::getline(split, each, delim); tokens.push_back(each));
        result.push_back(tokens);
    }
    assert(result.size() == test_tube.Ly);
    assert(result[0].size() == test_tube.Lx);
    assert(result[result.size()-1].size() == test_tube.Lx);
    
    return result;
}

template <class S, class P>
vector<vector<string>> App<S, P>::GenericPolyAnal(vector<string>& polybuf)
{
    polybuf.pop_back();
    polybuf.erase(polybuf.begin());
    assert(polybuf[0].compare("Two dimensional single layer space") == 0);
    assert(polybuf[polybuf.size()-1].compare("END") == 0);
    polybuf.erase(polybuf.begin());
    polybuf.pop_back();
    vector<vector<string>> result;
    vector<string> tokens;
    char delim = ' ';
    
    for (int i = 0; i < polybuf.size(); i++)
    {
        istringstream split(polybuf[i]);
        if (i % 2 == 0)
        {
            tokens.clear();
            for (std::string each; std::getline(split, each, delim); tokens.push_back(each));
        }
        else
        {
            for (std::string each; std::getline(split, each, delim); tokens.push_back(each));
            result.push_back(tokens);
        }
    }
    return result;
}

template <>
vector<vector<string>> App<Space2D2L, Pos2d2l>::GenericBondAnal(vector<string>& bondbuf)
{
    vector<vector<string>> result;
    result.push_back(bondbuf);
    return result;
}


#endif /* app_hpp */
