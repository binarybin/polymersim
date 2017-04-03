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
    vector<tuple<int, P, int, int, P>> GetPossibleCoMoves(Polymer<P>& poly);
    
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

template <> vector<tuple<int, Pos2d2l>> CornerMove<Space2D2L, Pos2d2l>::GetPossibleMoves(Polymer<Pos2d2l>& poly)
{
    vector<tuple<int, Pos2d2l>> possible_moves;
    for (int pos = 0; pos < poly.locs.size() - 2; pos++)
    {
        int x1 = poly.locs[pos].x; int y1 = poly.locs[pos].y;
        int x2 = poly.locs[pos+1].x; int y2 = poly.locs[pos+1].y;
        int x3 = poly.locs[pos+2].x; int y3 = poly.locs[pos+2].y;
        bool siml = poly.locs[pos].siml;
        int layer = siml? 0:1;
        
        if (x1 != x3 && y1 != y3)
        {
            if (x1 == x2)
            {
                assert(y1 != y2 && y2 == y3 && x2 != x3);
                if (space.space[layer][x3][y1] == 0)
                    possible_moves.push_back(make_tuple(pos+1, Pos2d2l(x3, y1, siml)));
            }
            else
            {
                assert(y1 == y2 && x2 == x3 && y2 != y3);
                if (space.space[layer][x1][y3] == 0)
                    possible_moves.push_back(make_tuple(pos+1, Pos2d2l(x1, y3, siml)));
            }
        }
    }
    return possible_moves;
}

template<>
vector<tuple<int, Pos2d2l, int, int, Pos2d2l>> CornerMove<Space2D2L, Pos2d2l>::GetPossibleCoMoves(Polymer<Pos2d2l>& polysim)
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




#endif /* cornermove_hpp */
