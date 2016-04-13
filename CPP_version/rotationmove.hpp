//
//  rotationmove.hpp
//  polymersim
//
//  Created by Bin Xu on 4/11/16.
//  Copyright © 2016 Bin Xu. All rights reserved.
//

#ifndef rotationmove_hpp
#define rotationmove_hpp


#include <tuple>
#include <vector>
#include <cassert>
#include <map>
#include <cstdlib>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P>
class RotationMove
{
private:
    S& space;
public:
    RotationMove(S& thespace): space(thespace) {}
    
    void UpdatePolymer(Polymer<P>& poly, vector<P> newpoints) // "l" for cLockwise, "n" for couNterclockwise
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
    
    vector<vector<Pos2d2l>> GetPossibleMoves(Polymer<P>& polysim, int polyid);
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
vector<vector<Pos2d2l>> RotationMove<Space2D2L, Pos2d2l>::GetPossibleMoves(Polymer<Pos2d2l>& polysim, int polyid)
{
    vector<vector<Pos2d2l>> possible_moves;
    vector<Pos2d2l> pure_new_points;
    
    int Lx = (int)(space.Lx);
    int Ly = (int)(space.Ly);
    
    auto Tentative = [Lx, Ly, polysim] (int first, int second)
    {
        Pos2d2l tentative_point(polysim.locs[second]);
        tentative_point.x += (Lx + polysim.locs[second].x - polysim.locs[first].x);
        tentative_point.y += (Ly + polysim.locs[second].y - polysim.locs[first].y);
        tentative_point.x %= Lx;
        tentative_point.y %= Ly;
        return tentative_point;
    };
    
    auto tent7 = Tentative(2,3);
    if (space.space[1][tent7.x][tent7.y] != 0)
    {
        return possible_moves;
    }
    auto tent6 = Tentative(5,4);
    if (space.space[1][tent6.x][tent6.y] != 0)
    {
        return possible_moves;
    }
    auto tent0 = Tentative(3,2);
    if (space.space[1][tent0.x][tent0.y] != 0)
    {
        return possible_moves;
    }
    auto tent1 = Tentative(4,5);
    if (space.space[1][tent1.x][tent1.y] != 0)
    {
        return possible_moves;
    }
    
    vector<Pos2d2l> cclk;
    cclk.push_back(tent0);
    cclk.push_back(tent1);
    cclk.push_back(polysim.locs[5]);
    cclk.push_back(polysim.locs[2]);
    cclk.push_back(polysim.locs[3]);
    cclk.push_back(polysim.locs[4]);
    cclk.push_back(tent6);
    cclk.push_back(tent7);
    
    vector<Pos2d2l> clk;
    clk.push_back(cclk[6]);
    clk.push_back(cclk[7]);
    clk.push_back(cclk[4]);
    clk.push_back(cclk[5]);
    clk.push_back(cclk[2]);
    clk.push_back(cclk[3]);
    clk.push_back(cclk[0]);
    clk.push_back(cclk[1]);
    
    possible_moves.push_back(cclk);
    possible_moves.push_back(clk);
    
    return possible_moves;
    
    
}

template <class S, class P>
tuple<int, vector<int>> RotationMove<S,P>::ComputeBondInc(Polymer<P> poly, vector<P> newpoints)
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



#endif /* rotationmove_hpp */
