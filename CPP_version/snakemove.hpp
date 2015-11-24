//
//  snakemove.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef snakemove_hpp
#define snakemove_hpp

#include "move.hpp"

class SnakeMove : public Move
{
public:
    SnakeMove(Space2D1L& space) : Move(space) {}
    
    void UpdatePolymer(Polymer& poly, int pointid, Position& newpoint)
    {
        if (pointid != 0) // move towards head
        {
            poly.locs.pop_back();
            poly.locs.insert(poly.locs.begin(), newpoint);
        }
        else
        {
            poly.locs.erase(poly.locs.begin());
            poly.locs.push_back(newpoint);
        }
    }
    
    void UpdateReverseCheckingSpace(Position& oldpoint, Position& newpoint, Polymer& poly)
    {
        int xnew = newpoint.x; int ynew = newpoint.y;
        int xold = oldpoint.x; int yold = oldpoint.y;
        
        assert(space.rspace[xnew][ynew][0] == 0);
        
        int polyid = space.rspace[xold][yold][0];
        space.rspace[xold][yold][0] = 0;
        space.rspace[xnew][ynew][1] = 0;
        
        for (int i = 0; i < poly.locs.size(); i++)
        {
            int xtemp = poly.locs[i].x;
            int ytemp = poly.locs[i].y;
            
            space.rspace[xtemp][ytemp][0] = polyid;
            space.rspace[xtemp][ytemp][1] = i;
        }
    }
    
    vector<tuple<int, Position>> GetPossibleMoves(Polymer& poly)
    {
        vector<tuple<int, Position>> possible_moves;
        int endloc_head = 0;
        int endloc_tail = (int)poly.locs.size() - 1;
        for (auto loc : space.Neighbor(poly.locs[endloc_head]))
            if (space.space[loc.x][loc.y] == 0)
                possible_moves.push_back(make_tuple(endloc_tail, loc)); // kill endloc_tail and create at loc
        
        for (auto loc : space.Neighbor(poly.locs[endloc_tail]))
            if (space.space[loc.x][loc.y] == 0)
                possible_moves.push_back(make_tuple(endloc_head, loc)); // kill endloc_head and create at loc
        return possible_moves;
    }
    
};


#endif /* snakemove_hpp */
