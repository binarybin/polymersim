//
//  move.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef move_hpp
#define move_hpp

#include <cstdlib>
#include <stdexcept>
#include <tuple>
#include <random>
#include "space2d2l.hpp"
#include "endmove.hpp"
#include "snakemove.hpp"
#include "cornermove.hpp"

using std::invalid_argument;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::random_device;
using std::mt19937;
using std::uniform_int_distribution;
using std::generate_canonical;

template <class S, class P, class M> class Move
{
private:
    S& space;
    M move;
    double beta;
    double gamma;
    
    int succ;
    
    random_device rd;
    mt19937 gen;
    
public:
    int bond_change;
    Move(S &thespace): move(thespace), space(thespace), beta(0), gamma(0), bond_change(0), succ(0), gen(rd()){}
    
    void ClearSucc(){succ=0;}
    int GetSucc(){return succ;}
    
    tuple<bool, int> ExecMove(int polyid, char polytyp);
    void SetBeta(double setbeta)  {beta = setbeta;}
    void SetGamma(double setgamma) {gamma = setgamma;}
    
    P ChooseBond(vector<P> bond_choice)
    {
        uniform_int_distribution<> dis(0, (int)bond_choice.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < bond_choice.size());
        return bond_choice[id];
    }
    
    tuple<int, int, P> ChooseMove(vector<tuple<int, int, P>> possible_moves)
    {
        uniform_int_distribution<> dis(0, (int)possible_moves.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    
    double Weight(int nbr_bond_inc, int nbr_ps_inc) {return exp(nbr_bond_inc * beta - nbr_ps_inc * gamma);}
    tuple<int, vector<int>> ComputeBondInc(Polymer<P>& poly, vector<P> newpoints);
    int ComputePSIncCrossLayer(P oldpoint, P newpoint, int polyid)
    {
        int nbr_ps_inc = 0;
        for (auto& pt : space.Neighbor(space.BondNeighbor(oldpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                nbr_ps_inc --;
            }
        }
        
        for (auto& pt : space.Neighbor(space.BondNeighbor(newpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                nbr_ps_inc ++;
            }
        }
        return nbr_ps_inc;
    }
    int ComputePSIncSameLayer(P oldpoint, P newpoint, int polyid)
    {
        int nbr_ps_inc = 0;
        for (auto& pt : space.Neighbor(oldpoint))
        {
            if (space.EmptyPos(pt))
            {
                nbr_ps_inc --;
            }
        }
        
        for (auto& pt : space.Neighbor(newpoint))
        {
            if (space.EmptyPos(pt) || space.GetRspacePoint(pt)[0] == polyid)
            {
                nbr_ps_inc ++;
            }
        }
        return nbr_ps_inc;
    }
};


// TODO: modify this
/*
template <class S, class P, class M>
tuple<int, vector<int>> Move<S, P, M>::ComputeBondInc(Polymer<P>& poly, vector<P> newpoints)
{
    int old_nbr_bond = 0, new_nbr_bond = 0;
    for (auto oldpoint : poly.locs)
    {
        if (space.InABond(oldpoint))
        {
            old_nbr_bond ++;
        }
    }
    
    std::map<int, vector<pair<int, int>>> epyc_rubisco; //map from rubisco_id to (point_in_epic, point_in_rubisco)
    for (int i = 0; i < space.LSim; i++)
        if (!space.Phosphorylated(poly.locs[i]))
        {
            P bpoint = space.BondNeighbor(newpoints[i])[0];
            int rubiscoid = space.GetRspacePoint(bpoint)[0];
            int rubiscopos = space.GetRspacePoint(bpoint)[1];
            if (rubiscoid != NOBOND)
                if (!space.Phosphorylated(bpoint))
                {
                    epyc_rubisco[rubiscoid].push_back(std::make_pair(i, rubiscopos));
                }
        }
    
    vector<int> result_pos;
    
    vector<pair<int, vector<pair<int, int>>>> intersect;
    
    for(auto epyc_link : epyc_rubisco)
    {
        if (epyc_link.second.size() > 1)
            intersect.push_back(epyc_link);
        else
            result_pos.push_back(epyc_link.second[0].first);
    }
    
    for (auto epyc_link : intersect)
        // each epyc_link is a pair with an epyc_id and a vector of rubisco points connected to it
    {
        int nbr_first_half = 0;
        int nbr_second_half = 0;
        
        for (auto pt_epyc_rubi : epyc_link.second)
        {
            if (pt_epyc_rubi.second >= 0 && pt_epyc_rubi.second < 4) nbr_first_half ++;
            else if (pt_epyc_rubi.second < 8 && pt_epyc_rubi.second >=4) nbr_second_half ++;
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
            for (auto pt_epyc_rubi : epyc_link.second)
                if (pt_epyc_rubi.second >= 0 && pt_epyc_rubi.second < 4)
                    result_pos.push_back(pt_epyc_rubi.first);
        }
        else
        {
            for (auto pt_epyc_rubi : epyc_link.second)
                if (pt_epyc_rubi.second < 8 && pt_epyc_rubi.second >= 4)
                    result_pos.push_back(pt_epyc_rubi.first);
        }
    }
    
    for (auto id : result_pos)
    {
        assert(!space.Phosphorylated(poly.locs[id]));
    }
    
    new_nbr_bond = (int)(result_pos.size());
    return std::make_tuple(new_nbr_bond - old_nbr_bond, result_pos);
}
*/
// The version that does not impose the non-two-end-interacting rule
template <class S, class P, class M>
tuple<int, vector<int>> Move<S, P, M>::ComputeBondInc(Polymer<P>& poly, vector<P> newpoints)
{
    int old_nbr_bond = 0, new_nbr_bond = 0;
    for (auto oldpoint : poly.locs)
    {
        if (space.InABond(oldpoint))
        {
            old_nbr_bond ++;
        }
    }
    vector<int> result_pos;
    for (int i = 0; i < space.LSim; i++)
        if (!space.Phosphorylated(poly.locs[i]))
        {
            P bpoint = space.BondNeighbor(newpoints[i])[0];
            int rubiscoid = space.GetRspacePoint(bpoint)[0];
            if (rubiscoid != NOBOND)
                if (!space.Phosphorylated(bpoint))
                {
                    result_pos.push_back(i);
                }
        }
    
    for (auto id : result_pos)
    {
        assert(!space.Phosphorylated(poly.locs[id]));
    }
    
    new_nbr_bond = (int)(result_pos.size());
    return std::make_tuple(new_nbr_bond - old_nbr_bond, result_pos);
}


template <class S, class P, class M>
tuple<bool, int> Move<S, P, M>::ExecMove(int polyid, char polytyp)
{
    Polymer<P>& poly = polytyp=='i' ? space.Sims[polyid] : space.Sumos[polyid];
    
    vector<tuple<int, int, P>> possible_moves = move.GetPossibleMoves(poly);
    
    if (possible_moves.empty()) return make_tuple(false, 0);
    
    int kill_pointid = 0; int create_pointid; P newpoint;
    tie(kill_pointid, create_pointid, newpoint) = ChooseMove(possible_moves);
    
    P oldpoint = poly.locs[kill_pointid];
    
    int nbr_bond_inc = 0; vector<int> bond_id_list;
    Polymer<P> temp_new_polymer = poly;
    move.UpdatePolymer(temp_new_polymer, kill_pointid, newpoint);
#ifndef NDEBUG
    cout<<"dealing with EPYC "<<polyid<<endl;
#endif
    tie(nbr_bond_inc, bond_id_list) = ComputeBondInc(poly, temp_new_polymer.locs);
    
    
    int nbr_ps_inc = ComputePSIncSameLayer(oldpoint, newpoint, polyid) + ComputePSIncCrossLayer(oldpoint, newpoint, polyid);
    
    if (generate_canonical<double, 10>(gen) < Weight(nbr_bond_inc, nbr_ps_inc))
    {
        int sitevalue = space.GetSpacePoint(oldpoint);
        sitevalue = sitevalue/abs(sitevalue);
        
        for (auto oldpoint : poly.locs)
        {
            space.UnsafeRemoveBond(oldpoint);
        }
        
        vector<int> phosphorylated_pts;
        for (int i = 0; i < poly.locs.size(); i++)
        {
            if (space.Phosphorylated(poly.locs[i]))
            {
                phosphorylated_pts.push_back(i);
                space.Dephosphorylate(polyid, polytyp, i);
            }
        }
        
        assert(!space.EmptyPos(oldpoint));
        space.SafeRemove(oldpoint);
        space.SafeCreate(newpoint, sitevalue);
        move.UpdatePolymer(poly, kill_pointid, newpoint);
        move.UpdateReverseCheckingSpace(oldpoint, newpoint, poly);
        
        
        for (auto pointid : bond_id_list)
        {
            space.CreateBond(poly.locs[pointid]);
        }
        for (auto phosid : phosphorylated_pts)
        {
            space.Phosphorylate(polyid, polytyp, phosid);
        }
        
        succ += 1;
        bond_change += nbr_bond_inc;
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);
}

#endif /* move_hpp */
