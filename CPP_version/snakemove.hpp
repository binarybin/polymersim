//
//  snakemove.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef snakemove_hpp
#define snakemove_hpp

#include <tuple>
#include <vector>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P>
class SnakeMove
{
private:
    S& space;
public:
    SnakeMove(S& thespace): space(thespace) {}
    void UpdatePolymer(Polymer<P>& poly, int pointid, P& newpoint)
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
    
    void UpdateReverseCheckingSpace(P& oldpoint, P& newpoint, Polymer<P>& poly)
    {
        assert(space.GetRspacePoint(newpoint)[0] == NOBOND);
        
        int polyid = space.GetRspacePoint(oldpoint)[0];
        space.SetRspacePoint(oldpoint, NOBOND, NOBOND);
        for (int i = 0; i < poly.locs.size(); i++)
            space.SetRspacePoint(poly.locs[i], polyid, i);
    }

    vector<tuple<int, int, P>> GetPossibleMoves(Polymer<P>& poly)
    {
        int endloc_head = 0;
        int endloc_tail = (int)poly.locs.size() - 1;
        vector<tuple<int, int, P>> possible_moves;
        
        for (auto pos : space.Neighbor(poly.locs[endloc_head]))
            if (space.EmptyPos(pos))
                possible_moves.push_back(make_tuple(endloc_tail, endloc_head, pos));
        
        for (auto pos : space.Neighbor(poly.locs[endloc_tail]))
            if (space.EmptyPos(pos))
                possible_moves.push_back(make_tuple(endloc_head, endloc_tail, pos));
        
        return possible_moves;
    }

};



#endif /* snakemove_hpp */
