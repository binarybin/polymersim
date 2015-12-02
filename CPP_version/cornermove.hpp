//
//  cornermove.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef cornermove_hpp
#define cornermove_hpp

#include <tuple>
#include <vector>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P> class CornerMove
{
public:
    S& space;
    CornerMove(S& thespace) : space(thespace) {}
    
    void UpdatePolymer(Polymer<P>& poly, int pointid, P& newpoint)
    {
        poly.locs[pointid] = newpoint;
    }
    
    void UpdateReverseCheckingSpace(P& oldpoint, P& newpoint, Polymer<P>& poly)
    {
        space.RSpacePointMove(oldpoint, newpoint);
    }
    vector<tuple<int, P>> GetPossibleMoves(Polymer<P>& poly);
    
};

template <> vector<tuple<int, Pos2d1l>> CornerMove<Space2D1L, Pos2d1l>::GetPossibleMoves(Polymer<Pos2d1l>& poly)
{
    vector<tuple<int, Pos2d1l>> possible_moves;
    for (int pos = 0; pos < poly.locs.size() - 2; pos++)
    {
        int x1 = poly.locs[pos].x; int y1 = poly.locs[pos].y;
        int x2 = poly.locs[pos+1].x; int y2 = poly.locs[pos+1].y;
        int x3 = poly.locs[pos+2].x; int y3 = poly.locs[pos+2].y;
        
        
        
        if (x1 != x3 && y1 != y3)
        {
      //      cout<<x1<<", "<<y1<<endl<<x2<<", "<<y2<<endl<<x3<<", "<<y3<<endl<<endl;
            if (x1 == x2)
            {
                assert(y1 != y2 && y2 == y3 && x2 != x3);
                if (space.space[x3][y1] == 0)
                    possible_moves.push_back(make_tuple(pos+1, Pos2d1l(x3, y1)));
            }
            else
            {
                assert(y1 == y2 && x2 == x3 && y2 != y3);
                if (space.space[x1][y3] == 0)
                    possible_moves.push_back(make_tuple(pos+1, Pos2d1l(x1, y3)));
            }
        }
    }
    return possible_moves;
}


#endif /* cornermove_hpp */
