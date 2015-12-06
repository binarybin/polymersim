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
#include "space2d1l.hpp"
#include "space2d2l.cpp"

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P> class CornerMove
{
private:
    S& space;
public:
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




#endif /* cornermove_hpp */
