//
//  space2d2l.hpp
//  polymersim
//
//  Created by Bin Xu on 11/30/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef space2d2l_hpp
#define space2d2l_hpp

#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <string>
#include <tuple>
#include "position.h"
using std::vector;
using std::max;
using std::remove_if;
using std::cout;
using std::endl;
using std::invalid_argument;
using std::pair;
using std::string;
using std::tuple;
using std::make_tuple;


class Space2D2L
{
public:
    const static string name;
    int NSim, LSim, NSumo, LSumo;
    size_t Lx, Ly;
    int SimId, SumoId;
    vector<Polymer<Pos2d2l>> Sims;
    vector<Polymer<Pos2d2l>> Sumos;
    vector<vector<vector<int>>> space; //location is specified as space[layer][x][y]
    vector<vector<vector<vector<int>>>> rspace; // rspace[layer][x][y](polymerid, pos_on_poly)
    
    Space2D2L(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly) :
    NSim(nsim), NSumo(nsumo), LSim(lsim), LSumo(lsumo), Lx(lx), Ly(ly),
    SimId(0), SumoId(0),
    space(2, vector<vector<int>>(lx, vector<int>(ly))),
    rspace(2, vector<vector<vector<int>>>(lx, vector<vector<int>>(ly, vector<int>(2, NOBOND))))
    {}
    
    void Resume(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal, vector<vector<string>>& bondanal)
    {
        cout<<"Resuming from file"<<endl;
        ResumeBasicInfo(spaceanal, polyanal);
        ResumeSpace(spaceanal, polyanal);
        ResumePolymer(spaceanal, polyanal);
        ResumeBond(bondanal);
        ResumeReverse();
        cout<<"Finished resuming"<<endl;
    }
    
    void ResumeBasicInfo(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal);
    void ResumeSpace(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal);
    void ResumePolymer(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal);
    void ResumeBond(vector<vector<string>>& bondanal);
    void ResumeReverse();
    
    void Initialize()
    {
        if (max(LSim, LSumo) <= Ly && max(NSim, 2*NSumo) < Lx)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            DiluteInit('h');
        }
        else if (max(LSim, LSumo) <= Lx && max(NSim, 2*NSumo) < Ly)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            DiluteInit('v');
        }
        else
        {
            int rows_min_sim = ceil((double)NSim/Lx); // We need so many rows, otherwise it overflows in x direction
            int rows_max_sim = (int)Ly/LSim; // We cannot fill more than this many columns, otherwise it overflows in y direction
            int rows_sim = -1; // Initialize
            cout<<"min rows sim: "<<rows_min_sim<<"\t max rows sim: "<<rows_max_sim<<endl;
            if (rows_max_sim < rows_min_sim) // No valid range
            {
                throw invalid_argument("Dense Initialization failed for Sim, maybe too dense");
            }
            else // Take an intermediate number
                rows_sim = sqrt(rows_max_sim * rows_min_sim);
            cout << "Used dense initialization with " << rows_sim << " rows for Sim" << endl;
            DenseInit(rows_sim, 'i');
            
            
            int rows_min_sumo = ceil((double)NSumo/(Lx/2));
            int rows_max_sumo = (int) Ly / (LSumo/2);
            int rows_sumo = -1;
            if (rows_max_sumo < rows_min_sumo)
                throw invalid_argument("Dense Initialization failed for Sumo, maybe too dense");
            else
                rows_sumo = sqrt(rows_max_sumo * rows_min_sumo);
            cout << "Used dense initialization with " << rows_sumo << " rows for Sumo" << endl;
            DenseInit(rows_sumo, 'u');
        }
    }
    void CheckPhosphorylate()
    {
        for (auto poly : Sims)
        {
            cout<<"polymer"<<poly.poly_id<<endl;
            for (int i = 0; i < poly.locs.size(); i++)
            {
                if (Phosphorylated(poly.locs[i]))
                {
                    cout<<"One phos site"<<i<<endl;
                }
            }
        }
    }
    void Phosphorylate(int polyid, bool polytype, int siteid)
    {
        auto pt = polytype? Sims[polyid].locs[siteid] : Sumos[polyid].locs[siteid];
        UnsafeRemoveBond(pt);
        int val = GetSpacePoint(pt);
        assert(abs(val) == 1);
        SetSpacePoint(pt, val * 3);
#ifndef NDEBUG
        cout<<"Phosphorylate polymer "<<polyid<<" of type "<<polytype<<" at site "<<siteid<<endl;
        CheckPhosphorylate();
#endif
    }
    void Dephosphorylate(int polyid, bool polytype, int siteid)
    {
        auto pt = polytype? Sims[polyid].locs[siteid] : Sumos[polyid].locs[siteid];
        int val = GetSpacePoint(pt);
        assert(abs(val) == 3);
        SetSpacePoint(pt, val/3);
#ifndef NDEBUG
        cout<<"Dephosphorylate polymer "<<polyid<<" of type "<<polytype<<" at site "<<siteid<<endl;
#endif

    }
    bool Phosphorylated(Pos2d2l point)
    {
        int x = point.x;
        int y = point.y;
        int layer = point.siml? 0:1;
        return abs(space[layer][x][y]) == 3;
    }
    
    vector<Pos2d2l> Neighbor(Pos2d2l point)
    {
        int x = point.x;
        int y = point.y;
        bool siml = point.siml;
        return {Pos2d2l((x+1)%Lx, y, siml), Pos2d2l((x-1+(int)Lx)%Lx, y, siml), Pos2d2l(x, (y+1)%Ly, siml), Pos2d2l(x, (y-1+(int)Ly)%Ly, siml)};
    }
    
    vector<Pos2d2l> BondNeighbor(Pos2d2l point)
    {
        return {Pos2d2l(point.x, point.y, !point.siml)};
    }
    
    void UnsafeRemoveBond(Pos2d2l point)
    {
        int x = point.x;
        int y = point.y;
        int layer = point.siml? 0:1;
        if(abs(space[layer][x][y]) == 2)
        {
            assert(space[1-layer][x][y] * space[layer][x][y] == -4);
            space[layer][x][y] /= 2;
            space[1-layer][x][y] /= 2;
        }
    }
    
    void SafeRemove(Pos2d2l point)
    {
        int x = point.x;
        int y = point.y;
        int layer = point.siml? 0:1;
        assert(space[layer][x][y] != 0);
        
        if (abs(space[layer][x][y])== 2) // we have a bond to remove, so remove it
        {
            assert(abs(space[1-layer][x][y])==2); // assert bond consistency
            space[1-layer][x][y] = [](int &x) {return (x>0) - (x<0);}(space[1-layer][x][y]);
        }
        space[layer][x][y] = 0;
    }
    
    void SafeCreate(Pos2d2l point, int polyvalue)
    {
        int x = point.x;
        int y = point.y;
        int layer = point.siml? 0:1;
        assert(space[layer][x][y] == 0);
        space[layer][x][y] = polyvalue;
    }
    
    void CreateBond(Pos2d2l point1, Pos2d2l point2)
    {
        int x1 = point1.x;
        int y1 = point1.y;
        int layer1 = point1.siml? 0:1;
        int x2 = point2.x;
        int y2 = point2.y;
        int layer2 = point2.siml? 0:1;
        
        assert(layer1+layer2 == 1 && x1==x2 && y1==y2);
        assert(abs(space[layer1][x1][y1]) == 1 && abs(space[layer2][x2][y2]) == 1);
        
        space[layer1][x1][y1] *= 2;
        space[layer2][x2][y2] *= 2;
    }
    void CreateBond(Pos2d2l point)
    {
        CreateBond(point, BondNeighbor(point)[0]);
    }
    
    bool CanBuildBond(Pos2d2l point1, Pos2d2l point2)
    {
        return space[point1.siml? 0:1][point1.x][point1.y] * space[point1.siml? 0:1][point2.x][point2.y] == -1;
    }
    
    bool ExistBond(Pos2d2l point1, Pos2d2l point2)
    {
        return space[point1.siml? 0:1][point1.x][point1.y] * space[point1.siml? 0:1][point2.x][point2.y] == -4;
    }
    
    bool InABond(Pos2d2l point)
    {
        return abs(space[point.siml? 0:1][point.x][point.y]) == 2;
    }
    
    bool CanBuildBondTentative(Pos2d2l bpoint, int sitevalue)
    {
        return space[bpoint.siml? 0:1][bpoint.x][bpoint.y] * sitevalue == -1;
    }
        
    void DiluteInit(char direction);
    void DenseInit(int rows, char typ);
    void Place(char typ, int id, vector<Pos2d2l> locs) // i for sim, u for sumo
    {
        int spacetype = (typ == 'i')? 1 : -1;
        vector<Polymer<Pos2d2l>>& polymers = (typ=='i')? Sims : Sumos;
        
        Polymer<Pos2d2l> poly(id, locs.size());
        poly.locs = locs;
        polymers.push_back(poly);
        
        for (int idx = 0; idx < locs.size(); idx++)
        {
            int x = locs[idx].x;
            int y = locs[idx].y;
            int layer = locs[idx].siml? 0:1;
            if(space[layer][x][y] != 0)
                throw(std::invalid_argument("Initialization: trying to assign one point with two monomers."));
            space[layer][x][y] = spacetype;
            rspace[layer][x][y][0] = id;
            rspace[layer][x][y][1] = idx;
        }
    }
    
    bool EmptyPos(Pos2d2l& pos)
    {
        return space[pos.siml? 0:1][pos.x][pos.y] == 0;
    }
    
    void RSpacePointMove(Pos2d2l& oldpoint, Pos2d2l& newpoint)
    {
        int xnew = newpoint.x; int ynew = newpoint.y; int layernew = newpoint.siml? 0:1;
        int xold = oldpoint.x; int yold = oldpoint.y; int layerold = newpoint.siml? 0:1;
     
        assert(layernew == layerold);
        assert(rspace[layernew][xnew][ynew][0] == NOBOND && rspace[layernew][xnew][ynew][1] == NOBOND);
        rspace[layernew][xnew][ynew][0] = rspace[layernew][xold][yold][0];
        rspace[layernew][xnew][ynew][1] = rspace[layernew][xold][yold][1];
        rspace[layernew][xold][yold][0] = NOBOND;
        rspace[layernew][xold][yold][1] = NOBOND;
    }
    int GetSpacePoint(Pos2d2l& point)
    {
        return space[point.siml?0:1][point.x][point.y];
    }
    vector<int> GetRspacePoint(Pos2d2l& point)
    {
        return {rspace[point.siml? 0:1][point.x][point.y][0], rspace[point.siml? 0:1][point.x][point.y][1]};
    }
    void SetRspacePoint(Pos2d2l& point, int polyid, int locid)
    {
        rspace[point.siml? 0:1][point.x][point.y][0] = polyid;
        rspace[point.siml? 0:1][point.x][point.y][1] = locid;
    }
    
    void SetSpacePoint(Pos2d2l& point, int val)
    {
        space[point.siml?0:1][point.x][point.y] = val;
    }

    void MoveTogether(Pos2d2l& oldpoint, Pos2d2l& newpoint)
    {
        int xold = oldpoint.x; int yold = oldpoint.y;
        int xnew = newpoint.x; int ynew = newpoint.y;
        
        assert(space[0][xnew][ynew] == 0 && space[1][xnew][ynew] == 0);
        assert(space[0][xold][yold] != 0 && space[1][xold][yold] != 0);
        
        space[0][xnew][ynew] = space[0][xold][yold]; space[0][xold][yold] = 0;
        space[1][xnew][ynew] = space[1][xold][yold]; space[1][xold][yold] = 0;
    }
    
    tuple<int, int> GetFromReverseSpace(Pos2d2l& cooldpoint)
    {
        int polyid  = rspace[cooldpoint.siml? 0:1][cooldpoint.x][cooldpoint.y][0];
        int polypos = rspace[cooldpoint.siml? 0:1][cooldpoint.x][cooldpoint.y][0];
        return make_tuple(polyid, polypos);
    }
    
    void ReduceEpycLength()
    {
        LSim -= 1;
        for (auto& poly : Sims)
        {
            auto pt = poly.locs.back();
            SafeRemove(pt);
            int layer = pt.siml? 0:1;
            rspace[layer][pt.x][pt.y][0] = NOBOND;
            rspace[layer][pt.x][pt.y][1] = NOBOND;
            poly.locs.pop_back();
        }
    }
};

std::ostream& PrintSpace(std::ostream &out, Space2D2L &space);
std::ostream& PrintBond(std::ostream &out, Space2D2L &space);
std::ostream& PrintPolymer(std::ostream &out, Space2D2L &space);
std::ostream& PrintRSpace(std::ostream &out, Space2D2L &space);
#endif /* space2d2l_hpp */
