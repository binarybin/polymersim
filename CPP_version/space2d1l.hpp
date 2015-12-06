//
//  space2d1l.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef space2d1l_hpp
#define space2d1l_hpp
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <cassert>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "position.h"
using std::vector;
using std::max;
using std::remove_if;
using std::cout;
using std::endl;
using std::invalid_argument;
using std::string;

class Space2D1L
{
public:
    int NSim, LSim, NSumo, LSumo;
    size_t Lx, Ly;
    int SimId, SumoId;
    vector<Polymer<Pos2d1l>> Sims;
    vector<Polymer<Pos2d1l>> Sumos;
    vector<vector<int>> space;
    vector<vector<vector<int>>> bond;
    vector<vector<vector<int>>> rspace;
    

    Space2D1L(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly) :
    NSim(nsim), NSumo(nsumo), LSim(lsim), LSumo(lsumo), Lx(lx), Ly(ly),
    SimId(0), SumoId(0),
        space(lx, vector<int>(ly)),
        bond(lx, vector<vector<int>>(ly, vector<int>(2, NOBOND))),
        rspace(lx, vector<vector<int>>(ly, vector<int>(2, NOBOND)))
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
        if (max(LSim, LSumo) <= Ly && NSim + NSumo < Lx)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            DiluteInit('h');
        }
        else if (max(LSim, LSumo) <= Lx && NSim + NSumo < Ly)
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
    
    vector<Pos2d1l> Neighbor(Pos2d1l point)
    {
        int x = point.x;
        int y = point.y;
        return {Pos2d1l((x+1)%Lx, y), Pos2d1l((x-1+(int)Lx)%Lx, y), Pos2d1l(x, (y+1)%Ly), Pos2d1l(x, (y-1+(int)Ly)%Ly)};
    }
    
    vector<Pos2d1l> BondNeighbor(Pos2d1l point)
    {
        return Neighbor(point);
    }
    
    void SafeRemove(Pos2d1l point)
    {
        int x = point.x;
        int y = point.y;
        assert(space[x][y] != 0);
        space[x][y] = 0;
        
        int xb = bond[x][y][0];
        int yb = bond[x][y][1];
        
        if (xb != NOBOND && yb != NOBOND) // we have a bond to remove, so remove it
        {
            assert(bond[xb][yb][0] == x && bond[xb][yb][1] == y); // assert bond reversibility
            assert(abs(space[xb][yb])==2); // assert bond consistency
            bond[x][y][0] = NOBOND;
            bond[x][y][1] = NOBOND;
            
            bond[xb][yb][0] = NOBOND;
            bond[xb][yb][1] = NOBOND;
            
            space[xb][yb] = [](int &x) {return (x>0) - (x<0);}(space[xb][yb]);
            assert(space[xb][yb] != 0);
        }
    }
    
    void SafeCreate(Pos2d1l point, int polyvalue)
    {
        int x = point.x;
        int y = point.y;
        assert(space[x][y] == 0);
        space[x][y] = polyvalue;
    }
    
    void CreateBond(Pos2d1l point1, Pos2d1l point2)
    {
        int x1 = point1.x;
        int y1 = point1.y;
        int x2 = point2.x;
        int y2 = point2.y;
        
        assert(bond[x1][y1][0] == NOBOND && bond[x1][y1][1] == NOBOND);
        assert(bond[x2][y2][0] == NOBOND && bond[x2][y2][1] == NOBOND);
        assert(abs(space[x1][y1]) == 1 && abs(space[x2][y2]));
        
        bond[x1][y1][0] = x2, bond[x1][y1][1] = y2;
        bond[x2][y2][0] = x1, bond[x2][y2][1] = y1;
        
        space[x1][y1] *= 2;
        space[x2][y2] *= 2;
    }
    
    bool CanBuildBond(Pos2d1l point1, Pos2d1l point2)
    {
        return space[point1.x][point1.y] * space[point2.x][point2.y] == -1;
    }
    
    bool ExistBond(Pos2d1l point1, Pos2d1l point2)
    {
        return space[point1.x][point1.y] * space[point2.x][point2.y] == -4;
    }
    
    bool InABond(Pos2d1l point)
    {
        return abs(space[point.x][point.y]) != 1 && abs(space[point.x][point.y]) != 0;
    }
    
    bool CanBuildBondTentative(Pos2d1l bpoint, int sitevalue)
    {
        return space[bpoint.x][bpoint.y] * sitevalue == -1;
    }
    
    void DiluteInit(char direction);
    void DenseInit(int rows);
    void Place(char typ, int id, vector<Pos2d1l> locs) // i for sim, u for sumo
    {
        int spacetype = (typ == 'i')? 1 : -1;
        vector<Polymer<Pos2d1l>>& polymers = (typ=='i')? Sims : Sumos;
        
        Polymer<Pos2d1l> poly(id, locs.size());
        poly.locs = locs;
        polymers.push_back(poly);
        
        for (int idx = 0; idx < locs.size(); idx++)
        {
            int x = locs[idx].x;
            int y = locs[idx].y;
            space[x][y] = spacetype;
            rspace[x][y][0] = id;
            rspace[x][y][1] = idx;
        }
    }
    
    bool EmptyPos(Pos2d1l& pos)
    {
        return space[pos.x][pos.y] == 0;
    }
    
    void RSpacePointMove(Pos2d1l& oldpoint, Pos2d1l& newpoint)
    {
        int xnew = newpoint.x; int ynew = newpoint.y;
        int xold = oldpoint.x; int yold = oldpoint.y;
        
        assert(rspace[xnew][ynew][0] == NOBOND);
        rspace[xnew][ynew][0] = rspace[xold][yold][0];
        rspace[xnew][ynew][1] = rspace[xold][yold][1];
        rspace[xold][yold][0] = NOBOND;
        rspace[xnew][ynew][1] = NOBOND;
    }
    vector<int> GetRspacePoint(Pos2d1l& point)
    {
        return {rspace[point.x][point.y][0], rspace[point.x][point.y][1]};
    }
    int GetSpacePoint(Pos2d1l& point)
    {
        return space[point.x][point.y];
    }
    void SetSpacePoint(Pos2d1l& point, int val)
    {
        space[point.x][point.y] = val;
    }
    void SetRspacePoint(Pos2d1l& point, int polyid, int locid)
    {
        rspace[point.x][point.y][0] = polyid;
        rspace[point.x][point.y][1] = locid;
    }

};

std::ostream& PrintSpace(std::ostream &out, Space2D1L &space);
std::ostream& PrintBond(std::ostream &out, Space2D1L &space);
std::ostream& PrintPolymer(std::ostream &out, Space2D1L &space);
#endif /* space2d1l_hpp */
