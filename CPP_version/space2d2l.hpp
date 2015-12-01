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
#include "position.h"
using std::vector;
using std::max;
using std::remove_if;
using std::cout;
using std::endl;
using std::invalid_argument;
using std::pair;

class Space2D2L
{
public:
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
    rspace(2, vector<vector<vector<int>>>(lx, vector<vector<int>>(ly, vector<int>(2))))
    {}
    
    void Initialize()
    {
        if (max(LSim, LSumo) <= Ly && max(NSim, NSumo) < Lx)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            DiluteInit('h');
        }
        else if (max(LSim, LSumo) <= Lx && max(NSim, NSumo) < Ly)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            DiluteInit('v');
        }
        else
        {
            int rows_min = ceil((NSim + NSumo)/Lx);
            int rows_max = floor(Ly / max(LSim, LSumo));
            int rows = -1;
            
            if (rows_max < rows_min)
                throw invalid_argument("Dense Initialization failed, maybe too dense");
            else
                rows = floor(sqrt(rows_max * rows_min));
            cout << "Used dense initialization with " << rows << " rows" << endl;
            DenseInit(rows);
        }
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
    
    void SafeRemove(Pos2d2l point)
    {
        int x = point.x;
        int y = point.y;
        int layer = point.siml? 0:1;
        assert(space[layer][x][y] != 0);
        space[layer][x][y] = 0;
        
        if (abs(space[layer][x][y])!= 1) // we have a bond to remove, so remove it
        {
            assert(abs(space[1-layer][x][y])==2); // assert bond consistency
            space[1-layer][x][y] = [](int &x) {return (x>0) - (x<0);}(space[1-layer][x][y]);
        }
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
        assert(abs(space[layer1][x1][y1]) == 1 && abs(space[layer2][x2][y2]));
        
        space[layer1][x1][y1] *= 2;
        space[layer2][x2][y2] *= 2;
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
        return abs(space[point.siml? 0:1][point.x][point.y]) != 1 && abs(space[point.siml? 0:1][point.x][point.y]) != 0;
    }
    
    bool CanBuildBondTentative(Pos2d2l bpoint, int sitevalue)
    {
        return space[bpoint.siml? 0:1][bpoint.x][bpoint.y] * sitevalue == -1;
    }
    
    void DiluteInit(char direction);
    void DenseInit(int rows);
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
        assert(rspace[layernew][xnew][ynew][0] == 0);
        rspace[layernew][xnew][ynew][0] = rspace[layernew][xold][yold][0];
        rspace[layernew][xnew][ynew][1] = rspace[layernew][xold][yold][1];
        rspace[layernew][xold][yold][0] = 0;
        rspace[layernew][xnew][ynew][1] = 0;
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
};

std::ostream& PrintSpace(std::ostream &out, Space2D2L &space);
std::ostream& PrintBond(std::ostream &out, Space2D2L &space);
std::ostream& PrintPolymer(std::ostream &out, Space2D2L &space);

#endif /* space2d2l_hpp */
