//
//  endmove.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef endmove_hpp
#define endmove_hpp
#include <tuple>
#include <vector>
#include "space2d1l.hpp"
#include "space2d2l.cpp"

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P> class EndMove
{
private:
    S& space;
public:
    EndMove(S& thespace) : space(thespace) {}
    
    void UpdatePolymer(Polymer<P>& poly, int pointid, P& newpoint)
    {
        poly.locs[pointid] = newpoint;
    }
    
    vector<tuple<int, P>> GetPossibleMoves(Polymer<P>& poly)
    {
        int endloc_head = 0; int nextloc_head = 1;
        int endloc_tail = (int)poly.locs.size() - 1; int nextloc_tail = (int)poly.locs.size() - 2;
        vector<tuple<int, P>> possible_moves;
        
        for (auto pos : space.Neighbor(poly.locs[nextloc_head]))
            if (space.EmptyPos(pos))
                possible_moves.push_back(make_tuple(endloc_head, pos));
        
        for (auto pos : space.Neighbor(poly.locs[nextloc_tail]))
            if (space.EmptyPos(pos))
                possible_moves.push_back(make_tuple(endloc_tail, pos));
        
        return possible_moves;
    }

    void UpdateReverseCheckingSpace(P& oldpoint, P& newpoint, Polymer<P>& poly)
    {
        space.RSpacePointMove(oldpoint, newpoint);
    }
};


#endif /* endmove_hpp */
