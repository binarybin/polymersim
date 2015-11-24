//
//  endmove.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef endmove_hpp
#define endmove_hpp

#include "move.hpp"

class EndMove : public Move
{
public:
    EndMove(Space2D1L& space) : Move(space) {}
    
    void UpdatePolymer(Polymer& poly, int pointid, Position& newpoint)
    {
        poly.locs[pointid] = newpoint;
    }
    
    void UpdateReverseCheckingSpace(Position& oldpoint, Position& newpoint, Polymer& poly)
    {
        int xnew = newpoint.x; int ynew = newpoint.y;
        int xold = oldpoint.x; int yold = oldpoint.y;
        
        assert(space.rspace[xnew][ynew][0] == 0);
        space.rspace[xnew][ynew][0] = space.rspace[xold][yold][0];
        space.rspace[xnew][ynew][1] = space.rspace[xold][yold][1];
        space.rspace[xold][yold][0] = 0;
        space.rspace[xnew][ynew][1] = 0;
    }
    
    vector<tuple<int, Position>> GetPossibleMoves(Polymer& poly)
    {
        int endloc_head = 0; int nextloc_head = 1;
        int endloc_tail = (int)poly.locs.size() - 1; int nextloc_tail = (int)poly.locs.size() - 2;
        vector<tuple<int, Position>> possible_moves;
        
        for (auto pos : space.Neighbor(poly.locs[nextloc_head]))
            if (space.space[pos.x][pos.y] == 0)
                possible_moves.push_back(make_tuple(endloc_head, pos));

        for (auto pos : space.Neighbor(poly.locs[nextloc_tail]))
            if (space.space[pos.x][pos.y] == 0)
                possible_moves.push_back(make_tuple(endloc_tail, pos));
        
        return possible_moves;
    }
};


#endif /* endmove_hpp */
