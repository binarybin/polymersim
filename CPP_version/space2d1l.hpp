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
#include "position.h"
#include "space.hpp"
using std::vector;
using std::max;
using std::remove_if;
using std::cout;
using std::endl;
using std::invalid_argument;

class Space2D1L : protected Space
{
    friend std::ostream& PrintSpaceOccupation(std::ostream &out, Space2D1L &space);
    vector<vector<int>> space;
    vector<vector<vector<int>>> bond;
    vector<vector<vector<int>>> rspace;
    
public:
    Space2D1L(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly) :
        Space(nsim, nsumo, lsim, lsumo, lx, ly, 0),
        space(lx, vector<int>(ly)),
        bond(lx, vector<vector<int>>(ly, vector<int>(2))),
        rspace(lx, vector<vector<int>>(ly, vector<int>(2)))
    {}
    
    void Initialize()
    {
        if (max(this->LSim, this->LSumo) <= this->Ly
            && this->NSim + this->NSumo < this->Lx)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            this->DiluteInit('h');
        }
        else if (max(this->LSim, this->LSumo) <= this->Lx
            && this->NSim + this->NSumo < this->Ly)
        {
            cout << "Used dilute initialization, horizontal" <<endl;
            this->DiluteInit('v');
        }
        else
        {
            int rows_min = ceil((this->NSim + this->NSumo)/this->Lx);
            int rows_max = floor(this->Ly / max(this->LSim, this->LSumo));
            int rows = -1;
            
            if (rows_max < rows_min)
                throw invalid_argument("Dense Initialization failed, maybe too dense");
            else
                rows = floor(sqrt(rows_max * rows_min));
            cout << "Used dense initialization with " << rows << " rows" << endl;
            this -> DenseInit(rows);
        }
    }
    
    vector<Position> Neighbor(Position point)
    {
        int x = point.x;
        int y = point.y;
        return {Position((x+1)%Lx, y), Position((x-1)%Lx, y), Position(x, (y+1)%Ly), Position(x, (y-1)%Ly)};
    }
    
    void SafeRemove(Position point)
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
    
    void SafeCreate(Position point, int polyvalue)
    {
        int x = point.x;
        int y = point.y;
        assert(space[x][y] == 0);
        space[x][y] = polyvalue;
    }
    
    void CreateBond(Position point1, Position point2)
    {
        int x1 = point1.x;
        int y1 = point1.y;
        int x2 = point2.x;
        int y2 = point2.y;
        
        assert(bond[x1][y1][0] == NOBOND && bond[x1][y1][1] == NOBOND);
        assert(bond[x2][y2][0] == NOBOND && bond[x2][y2][1] == NOBOND);
        
        bond[x1][y1][0] = x2, bond[x1][y1][1] = y2;
        bond[x2][y2][0] = x1, bond[x2][y2][1] = y1;
        
        space[x1][y1] *= 2;
        space[x2][y2] *= 2;
    }
    
    bool CanBuildBond(Position point1, Position point2)
    {
        return space[point1.x][point1.y] * space[point2.x][point2.y] == -1;
    }
    
    bool ExistBond(Position point1, Position point2)
    {
        return space[point1.x][point1.y] * space[point2.x][point2.y] == -4;
    }
    
    bool InABond(Position point)
    {
        return abs(space[point.x][point.y]) != 1;
    }
    
    bool CanBuildBondTentative(Position bpoint, int sitevalue)
    {
        return space[bpoint.x][bpoint.y] * sitevalue == -1;
    }
    
    void DiluteInit(char direction);
    void DenseInit(int rows);
    void Place(char typ, int id, vector<Position> locs) // i for sim, u for sumo
    {
        int spacetype = 0;
        vector<Polymer>& polymers = Sims;
        if (typ=='i')
        {
            spacetype = 1;
            polymers = Sims;
            assert(locs.size() == LSim);
        }
        else if (typ == 'u')
        {
            spacetype = -1;
            polymers = Sumos;
            assert(locs.size() == LSumo);
        }
        else
            throw std::invalid_argument("Unrecognized polymer type");
        
        Polymer poly(id, locs.size());
        polymers.push_back(poly);
        
        for (int idx = 0; idx < locs.size(); idx++)
        {
            int x = locs[idx].x;
            int y = locs[idx].y;
            space[x][y] = spacetype;
            rspace[x][y][0] = id+1;
            rspace[x][y][1] = idx;
        }
    }
    

    
};

std::ostream& PrintSpaceOccupation(std::ostream &out, Space2D1L &space);

#endif /* space2d1l_hpp */
