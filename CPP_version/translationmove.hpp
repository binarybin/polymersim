//
//  translationmove.hpp
//  polymersim
//
//  Created by Bin Xu on 4/8/16.
//  Copyright © 2016 Bin Xu. All rights reserved.
//

#ifndef translationmove_hpp
#define translationmove_hpp

#include <tuple>
#include <vector>
#include <cassert>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P>
class TranslationMove
{
private:
    S& space;
public:
    TranslationMove(S& thespace): space(thespace) {}
    vector<vector<P>> GetPossibleMoves(Polymer<P>& polysim, int polyid);
};


template<>
vector<vector<Pos2d2l>> TranslationMove<Space2D2L, Pos2d2l>::GetPossibleMoves(Polymer<Pos2d2l>& polysim, int polyid)
{
    vector<vector<Pos2d2l>> possible_moves;
    bool can_do = true;
    vector<Pos2d2l> vector_r;
    for (auto point : polysim.locs)
    {
        point.x = (point.x + 1) % space.Lx;
        int new_polyid = space.GetRspacePoint(point)[0];
        if (new_polyid != NOBOND && new_polyid != polyid)
        {
            can_do = false;
            break;
        }
        vector_r.push_back(point);
    }
    if (can_do) possible_moves.push_back(vector_r);
    
    vector<Pos2d2l> vector_l;
    can_do = true;
    for (auto point : polysim.locs)
    {
        point.x = ((int)(point.x) - 1 + (int)(space.Lx)) % (int)(space.Lx);
        int new_polyid = space.GetRspacePoint(point)[0];
        if (new_polyid != NOBOND && new_polyid != polyid)
        {
            can_do = false;
            break;
        }
        vector_l.push_back(point);
    }
    if (can_do) possible_moves.push_back(vector_l);
    
    
    can_do = true;
    vector<Pos2d2l> vector_d;
    for (auto point : polysim.locs)
    {
        point.y = (point.y + 1) % space.Ly;
        int new_polyid = space.GetRspacePoint(point)[0];
        if (new_polyid != NOBOND && new_polyid != polyid)
        {
            can_do = false;
            break;
        }
        vector_d.push_back(point);
    }
    if (can_do) possible_moves.push_back(vector_d);
    
    can_do = true;
    vector<Pos2d2l> vector_u;
    for (auto point : polysim.locs)
    {
        point.y = ((int)(point.y) - 1 + (int)(space.Ly)) % (int)(space.Ly);
        int new_polyid = space.GetRspacePoint(point)[0];
        if (new_polyid != NOBOND && new_polyid != polyid)
        {
            can_do = false;
            break;
        }
        vector_u.push_back(point);
    }
    if (can_do) possible_moves.push_back(vector_u);
    
    
    return possible_moves;
    
    
}



#endif /* translationmove_hpp */
