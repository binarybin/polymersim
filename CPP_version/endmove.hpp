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
    
    vector<tuple<int, P, int, int, P>> GetPossibleCoMoves(Polymer<P>& poly);
   
    void UpdateReverseCheckingSpace(P& oldpoint, P& newpoint, Polymer<P>& poly)
    {
        space.RSpacePointMove(oldpoint, newpoint);
    }
    
};

template<>
vector<tuple<int, Pos2d2l, int, int, Pos2d2l>> EndMove<Space2D2L, Pos2d2l>::GetPossibleCoMoves(Polymer<Pos2d2l>& polysim)
{
    vector<tuple<int, Pos2d2l, int, int, Pos2d2l>> result;
    vector<tuple<int, Pos2d2l>> possible_sim_moves = GetPossibleMoves(polysim);
    for (auto simmove : possible_sim_moves)
    {
        int simpointid=0; Pos2d2l point; std::tie(simpointid, point) = simmove;
        Pos2d2l otheroldpoint = polysim.locs[simpointid].OtherLayer();
        if (!space.EmptyPos(otheroldpoint))
        {
            int sumoid = space.rspace[1][otheroldpoint.x][otheroldpoint.y][0];
            Polymer<Pos2d2l>& polysumo = space.Sumos[sumoid];
            vector<tuple<int, Pos2d2l>> possible_sumo_moves = GetPossibleMoves(polysumo);
            for (auto sumomove : possible_sumo_moves)
            {
                int sumopointid = 0; Pos2d2l sumopoint; std::tie(sumopointid, sumopoint) = sumomove;
                Pos2d2l oldsumopoint = polysumo.locs[sumopointid];
                if(sumopoint.x == point.x && sumopoint.y == point.y && otheroldpoint.x == oldsumopoint.x&& otheroldpoint.y == oldsumopoint.y)
                    result.push_back(std::make_tuple(simpointid, point, sumoid, sumopointid, sumopoint));
            }
        }
    }
    
    return result;
}




#endif /* endmove_hpp */
