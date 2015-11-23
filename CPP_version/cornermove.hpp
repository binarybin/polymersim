//
//  cornermove.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef cornermove_hpp
#define cornermove_hpp

#include "move.hpp"

class CornerMove : public Move
{
public:
    CornerMove(Space2D1L& space) : Move(space) {}
    
    void UpdatePolymer(Polymer& poly, int pointid, Position& newpoint)
    {
        poly.locs[pointid] = newpoint;
    }
    
    void UpdateReverseCheckingSpace(Position& oldpoint, Position& newpoint)
    {
        int xnew = newpoint.x; int ynew = newpoint.y;
        int xold = oldpoint.x; int yold = oldpoint.y;
        
        assert(space.rspace[xnew][ynew][0] == 0 && space.rspace[xnew][ynew][1] == 0);
        space.rspace[xnew][ynew][0] = space.rspace[xold][yold][0];
        space.rspace[xnew][ynew][1] = space.rspace[xold][yold][1];
        space.rspace[xold][yold][0] = 0;
        space.rspace[xnew][ynew][1] = 0;
    }
    
    vector<tuple<int, Position>> GetPossibleMoves(Polymer& poly)
    {
        vector<tuple<int, Position>> possible_moves;
        for (int pos = 0; pos < poly.locs.size() - 2; pos++)
        {
            int x1 = poly.locs[pos].x; int y1 = poly.locs[pos].y;
            int x2 = poly.locs[pos+1].x; int y2 = poly.locs[pos+1].y;
            int x3 = poly.locs[pos+2].x; int y3 = poly.locs[pos+2].y;
            
            if (x1 != x3 && y1 != y3)
            {
                if (x1 == x2)
                {
                    assert(y1 != y2 && y2 == y3 && x2 != x3);
                    if (space.space[x3][y1] == 0)
                        possible_moves.push_back(make_tuple(pos+1, Position(x3, y1)));
                }
                else
                {
                    assert(y1 == y2 && x2 == x3 && y2 != y3);
                    if (space.space[x1][x3] == 0)
                        possible_moves.push_back(make_tuple(pos+1, Position(x1, y3)));
                }
            }
        }
        return possible_moves;
    }

};

#endif /* cornermove_hpp */
