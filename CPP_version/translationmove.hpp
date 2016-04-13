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
    void UpdatePolymer(Polymer<P>& poly, vector<P> newpoints)
    {
        poly.locs = newpoints;
    }
    void UpdateReverseCheckingSpace(vector<P>& oldpoints, vector<P> newpoints)
    {
        int polyid = space.GetRspacePoint(oldpoints[0])[0];
        
        for (auto oldpt: oldpoints)
            space.SetRspacePoint(oldpt, NOBOND, NOBOND);
        
        for (int i = 0; i < newpoints.size(); i++)
            space.SetRspacePoint(newpoints[i], polyid, i);
    }
    vector<vector<P>> GetPossibleMoves(Polymer<P>& polysim, int polyid);
    tuple<int, vector<int>> ComputeBondInc(Polymer<P> poly, vector<P> newpoints);
    int ComputePSInc(Polymer<P> poly, vector<P> newpoints, int polyid)
    {
        int old_nbr_ps = 0, new_nbr_ps = 0;
        for (auto oldpoint : poly.locs)
        {
            for (auto pt : space.Neighbor(oldpoint))
            {
                if (space.EmptyPos(pt))
                {
                    old_nbr_ps ++;
                }
            }
        }
        
        for (auto newpoint : newpoints)
        {
            for (auto pt : space.Neighbor(newpoint))
            {
                if (space.EmptyPos(pt) || space.GetRspacePoint(pt)[0] == polyid)
                {
                    new_nbr_ps ++;
                }
            }
        }
        return new_nbr_ps - old_nbr_ps;
    }
    void BuildNewBonds(vector<P> newpoints, vector<int> bond_id_list)
    {
        for(int id : bond_id_list)
        {
            space.CreateBond(newpoints[id]);
        }
    }
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

template <class S, class P>
tuple<int, vector<int>> TranslationMove<S,P>::ComputeBondInc(Polymer<P> poly, vector<P> newpoints)
{
    int old_nbr_bond = 0, new_nbr_bond = 0;
    for (auto oldpoint : poly.locs)
    {
        if (space.InABond(oldpoint))
        {
            old_nbr_bond ++;
        }
    }
    
    std::map<int, vector<int>> rubisco_epyc; //map from epyc_id to point_in_rubisco
    for (int i = 0; i < space.LSumo; i++)
    {
        P bpoint = space.BondNeighbor(newpoints[i])[0];
        int epycid = space.GetRspacePoint(bpoint)[0];
        if (epycid != NOBOND)
        {
            rubisco_epyc[epycid].push_back(i);
        }
    }
    
    vector<int> result_pos;
    
    vector<pair<int, vector<int>>> intersect;
    
    for(auto epyc_link : rubisco_epyc)
    {
        if (epyc_link.second.size() > 1)
            intersect.push_back(epyc_link);
        else
            result_pos.push_back(epyc_link.second[0]);
    }
    
    for (auto epyc_link : intersect)
        // each epyc_link is a pair with an epyc_id and a vector of rubisco points connected to it
    {
        int nbr_first_half = 0;
        int nbr_second_half = 0;
        
        for (auto pointid : epyc_link.second)
        {
            if (pointid >= 0 && pointid < 4) nbr_first_half ++;
            else if (pointid < 8 && pointid >=4) nbr_second_half ++;
            else
                throw(std::invalid_argument("ComputeBondInc"));
        }
        
        
        if (nbr_first_half == nbr_second_half) // The equal case, randomly decide which half to pick
        {
            if (rand()%2 == 0)
                nbr_first_half ++;
            else
                nbr_second_half ++;
        }
        
        
        if (nbr_first_half > nbr_second_half)
        {
            for (auto pointid : epyc_link.second)
                if (pointid >= 0 && pointid < 4)
                    result_pos.push_back(pointid);
        }
        else
        {
            for (auto pointid : epyc_link.second)
                if (pointid < 8 && pointid >= 4)
                    result_pos.push_back(pointid);
        }
    }
    
    new_nbr_bond = result_pos.size();
    return std::make_tuple(new_nbr_bond - old_nbr_bond, result_pos);
}

#endif /* translationmove_hpp */
