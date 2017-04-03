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
#include "space2d1l.hpp"
#include "space2d2l.hpp"
#include "move.hpp"


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
    S test_tube;
    Move<S, P, SnakeMove<S, P>> sm;
    Move<S, P, CornerMove<S, P>> cm;
    Move<S, P, EndMove<S, P>> em;
    Move<S, P, SnakeMove<S, P>> csm;
    Move<S, P, CornerMove<S, P>> ccm;
    Move<S, P, EndMove<S, P>> cem;
    random_device rd;
    mt19937 gen;
    int original_nbr_bond;
public:
    App(int nsim1, int nsim2, int nsumo, int lsim1, int lsim2, int lsumo, size_t lx, size_t ly):
        test_tube(nsim1, nsim2, nsumo, lsim1, lsim2, lsumo, lx, ly),
        sm(test_tube), cm(test_tube), em(test_tube),
        csm(test_tube), ccm(test_tube), cem(test_tube),
        gen(rd()), original_nbr_bond(0) {}
    
    App(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly, size_t lz):
        test_tube(nsim, nsumo, lsim, lsumo, lx, ly, lz),
        sm(test_tube), cm(test_tube), em(test_tube),
        csm(test_tube), ccm(test_tube), cem(test_tube),
        gen(rd()), original_nbr_bond(0) {}
    
    void Initialize() {test_tube.Initialize();}
    
    void Proceed(char typ);
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
                
            case 'S':
                returnvalue = csm.GetSucc();
                csm.ClearSucc();
                break;
                
            case 'E':
                returnvalue = cem.GetSucc();
                cem.ClearSucc();
                break;
                
            case 'C':
                returnvalue = ccm.GetSucc();
                ccm.ClearSucc();
                break;
                
            default:
                break;
        }
        return returnvalue;
    }
    
    double GetEnergy(){return -(original_nbr_bond + sm.bond_change + em.bond_change + cm.bond_change + csm.bond_change + ccm.bond_change + cem.bond_change);}
    void ShowSpace(std::ostream &out) {PrintSpace(out, test_tube); }
    void ShowBond(std::ostream &out) {PrintBond(out, test_tube); }
    void ShowPolymer(std::ostream &out) {PrintPolymer(out, test_tube); }

    void SetBeta(double beta)
    {
        sm.SetBeta(beta);
        cm.SetBeta(beta);
        em.SetBeta(beta);
        csm.SetBeta(beta);
        ccm.SetBeta(beta);
        cem.SetBeta(beta);
    }
    
    void SetGamma(double gamma)
    {
        sm.SetGamma(gamma);
        cm.SetGamma(gamma);
        em.SetGamma(gamma);
        csm.SetGamma(gamma);
        ccm.SetGamma(gamma);
        cem.SetGamma(gamma);
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
void App<Space2D1L, Pos2d1l>::TestPolymerAndSpace()
{
    auto tmpspace(test_tube.space);
    for (auto& poly : test_tube.Sims)
    {
        for (auto& loc : poly.locs)
        {
            int x = loc.x;
            int y = loc.y;
            if (tmpspace[x][y] != 1 && tmpspace[x][y] != 2)
            {
                throw runtime_error("Polymer and space test failed for SIM");
            }
            tmpspace[x][y] = 0;
        }
    }
    
    for (auto& poly : test_tube.Sumos)
    {
        for (auto& loc : poly.locs)
        {
            int x = loc.x;
            int y = loc.y;
            if (tmpspace[x][y] != -1 && tmpspace[x][y] != -2)
            {
                throw runtime_error("Polymer and space test failed for SUMO");
            }
            tmpspace[x][y] = 0;
        }
    }
    
    for (int i = 0; i < tmpspace.size(); i++)
    {
        for (int j = 0; j < tmpspace[0].size(); j++)
        {
            if (tmpspace[i][j] != 0)
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
            if (tmpspace[i][j] != test_tube.space[i][j])
            {
                fail_deep_copy = false;
            }
        }
    }
    
    if (fail_deep_copy)
    {
//        throw runtime_error("Deep copy does not seem to work");
    }
    cout<<"Polymer and space test passed!"<<endl;
}

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
            if (tmpspace[layer][x][y] != 1 && tmpspace[layer][x][y] != 2)
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
            if (tmpspace[layer][x][y] != -1 && tmpspace[layer][x][y] != -2)
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
void App<Space2D1L, Pos2d1l>::TestSpaceAndBond()
{
    auto tempspace(test_tube.space);
    for (int i = 0; i < test_tube.bond.size(); i++)
    {
        for (int j = 0;  j < test_tube.bond[0].size(); j++)
        {
            int x = test_tube.bond[i][j][0];
            int y = test_tube.bond[i][j][1];
            if (x != NOBOND && y != NOBOND)
            {
                if (abs(tempspace[i][j]) != 2)
                {
                    throw runtime_error("A bond suggested by bond is not fond on space");
                }
                tempspace[i][j] /= 2;
                int rx = test_tube.bond[x][y][0];
                int ry = test_tube.bond[x][y][1];
                if (rx != i || ry != j)
                {
                    throw runtime_error("Bond reversibility check failed");
                }
                
            }
        }
    }
    
    for (int i = 0; i < tempspace.size(); i++)
    {
        for (int j = 0; j < tempspace[0].size(); j++)
        {
            if (abs(tempspace[i][j]) == 2)
            {
                throw runtime_error("Bond check failed, more bond in space than in bond");
            }
        }
    }
    
    cout<<"Bond and space test passed!"<<endl;
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
void App<Space2D1L, Pos2d1l>::TestReverseSpace()
{
    auto tempspace(test_tube.space);
    for (int i = 0; i < test_tube.rspace.size(); i++)
    {
        for (int j = 0; j < test_tube.rspace[0].size(); j++)
        {
            if (test_tube.rspace[i][j][0] == NOBOND && test_tube.rspace[i][j][1] == NOBOND)
                continue;
            
            int id = test_tube.rspace[i][j][0];
            int pos = test_tube.rspace[i][j][1];
            Pos2d1l loc = test_tube.space[i][j] > 0 ? test_tube.Sims[id].locs[pos] : test_tube.Sumos[id].locs[pos];
            int rx = loc.x;
            int ry = loc.y;
            if (rx != i || ry != j)
            {
                throw runtime_error("Reverse space inconsistent");
            }
            tempspace[i][j] = 0;
        }
    }
    
    for (int i = 0; i < tempspace.size(); i++)
    {
        for (int j = 0; j < tempspace[0].size(); j++)
        {
            if (tempspace[i][j] != 0)
            {
                throw runtime_error("Reverse space check failed, more monomers in space than in bond");
            }
        }
    }

    cout<<"Reverse space test passed!"<<endl;
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

    
    cout<<"Reverse space test passed!"<<endl;
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
    char typ_r = generate_canonical<double, 3>(gen) >= 0.5 ? 'i' : 'u';
    int NPoly = (typ_r == 'i' ? test_tube.NSim : test_tube.NSumo);
    if(typ >= 'A' and typ <='Z') NPoly = test_tube.NSim;
    uniform_int_distribution<> dis(0, NPoly-1);
    int id_r = dis(gen);
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
            
        case 'S':
            csm.ExecCoMove(id_r);
            break;
            
        case 'E':
            cem.ExecCoMove(id_r);
            break;
            
        case 'C':
            ccm.ExecCoMove(id_r);
            break;
            
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

template <>
vector<vector<string>> App<Space2D1L, Pos2d1l>::GenericBondAnal(vector<string>& bondbuf)
{
    bondbuf.pop_back();
    bondbuf.erase(bondbuf.begin());
    assert(bondbuf[0].compare("Two dimensional single layer space") == 0);
    assert(bondbuf[1].compare("Bond configuration: ") == 0);
    assert(bondbuf[bondbuf.size()-1].compare("END") == 0);
    bondbuf.erase(bondbuf.begin());
    bondbuf.erase(bondbuf.begin());
    bondbuf.pop_back();
    vector<vector<string>> result;
    char delim = '\t';
    
    for (int i = 0; i < bondbuf.size(); i++)
    {
        vector<string> tokens;
        istringstream split(bondbuf[i]);
        for (std::string each; std::getline(split, each, delim); tokens.push_back(each));
        result.push_back(tokens);
    }
    
    return result;
}


#endif /* app_hpp */
