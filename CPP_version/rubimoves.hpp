//
//  rubimoves.hpp
//  polymersim
//
//  Created by Bin Xu on 4/12/16.
//  Copyright © 2016 Bin Xu. All rights reserved.
//

#ifndef rubimoves_hpp
#define rubimoves_hpp

#include <cstdlib>
#include <stdexcept>
#include <tuple>
#include <random>
#include <set>
#include <algorithm>
#include <unordered_map>
#include "space2d1l.hpp"
#include "space2d2l.hpp"
#include "translationmove.hpp"
#include "rotationmove.hpp"

using std::invalid_argument;
using std::tuple;
using std::make_tuple;
using std::tie;
using std::random_device;
using std::mt19937;
using std::uniform_int_distribution;
using std::generate_canonical;

template <class S, class P, class M> class RubiMove
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
    RubiMove(S &thespace): move(thespace), space(thespace), beta(0), gamma(0), bond_change(0), succ(0), gen(rd()){}
    
    void ClearSucc(){succ=0;}
    int GetSucc(){return succ;}
    
    void SetBeta(double setbeta)  {beta = setbeta;}
    void SetGamma(double setgamma) {gamma = setgamma;}
    
    tuple<bool, int> ExecMove(int polyid, char polytyp);
    vector<P> ChooseMove(vector<vector<P>> possible_moves)
    {
        uniform_int_distribution<> dis(0, (int)possible_moves.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    
    double Weight(int nbr_bond_inc, int nbr_ps_inc) {return exp(nbr_bond_inc * beta - nbr_ps_inc * gamma);}
    
    void UpdateReverseCheckingSpace(vector<P>& oldpoints, vector<P> newpoints)
    {
        int polyid = space.GetRspacePoint(oldpoints[0])[0];
        
        for (auto oldpt: oldpoints)
            space.SetRspacePoint(oldpt, NOBOND, NOBOND);
        
        for (int i = 0; i < newpoints.size(); i++)
            space.SetRspacePoint(newpoints[i], polyid, i);
    }
    
    void BuildNewBonds(vector<P> newpoints, vector<int> bond_id_list)
    {
        for(int id : bond_id_list)
        {
            space.CreateBond(newpoints[id]);
        }
    }
    tuple<int, vector<int>> ComputeBondInc(Polymer<P> poly, vector<P> newpoints);
    
    bool PointInVector(P point, vector<P>& vec)
    {
        for (auto pt : vec)
        {
            assert(pt.siml == point.siml);
            if (pt.x == point.x && pt.y == point.y)
            {
                return true;
            }
        }
        return false;
    }
    int ComputePSIncCrossLayer(Polymer<P> poly, vector<P> newpoints, int polyid);
    int ComputePSIncSameLayer(Polymer<P> poly, vector<P> newpoints, int polyid);
    
    
    
    //TO BE MODIFIED
    
    bool ExecDragMove(int polyid);
    void UpdateReverseCheckingSpaceDragMove(vector<P>& oldpoints1, vector<P> newpoints1, vector<P>& oldpoints2, vector<P> newpoints2)
    {
        int polyid1 = space.GetRspacePoint(oldpoints1[0])[0];
        int polyid2 = space.GetRspacePoint(oldpoints2[0])[0];
        
        for (auto oldpt: oldpoints1)
            space.SetRspacePoint(oldpt, NOBOND, NOBOND);
        
        for (auto oldpt: oldpoints2)
            space.SetRspacePoint(oldpt, NOBOND, NOBOND);
        
        for (int i = 0; i < newpoints1.size(); i++)
            space.SetRspacePoint(newpoints1[i], polyid1, i);
        for (int i = 0; i < newpoints2.size(); i++)
            space.SetRspacePoint(newpoints2[i], polyid2, i);
    }

    int ComputePSIncDragSameLayer(TriMoveInfo<P>& the_trimove);
    int ComputePSIncDragCrossLayer(TriMoveInfo<P>& the_trimove);
};


template <class S, class P, class M>
tuple<bool, int> RubiMove<S, P, M>::ExecMove(int polyid, char polytyp)
{
    Polymer<P>& poly = polytyp=='i' ? space.Sims[polyid] : space.Sumos[polyid];
    int sitevalue = polytyp=='i' ? 1 : -1;
    
    auto possible_moves = move.GetPossibleMoves(poly, polyid);
    
    if (possible_moves.empty()) return make_tuple(false, 0);
    
    auto newpoints = ChooseMove(possible_moves);
    
    int nbr_bond_inc = 0; vector<int> bond_id_list;
    tie(nbr_bond_inc, bond_id_list) = ComputeBondInc(poly, newpoints);
    
    int nbr_ps_inc = ComputePSIncCrossLayer(poly, newpoints, polyid) + ComputePSIncSameLayer(poly, newpoints, polyid);
    
    if (generate_canonical<double, 10>(gen) < Weight(nbr_bond_inc, nbr_ps_inc))
    {
        for (auto oldpoint : poly.locs)
        {
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        
        for (auto newpoint : newpoints)
        {
            assert(space.EmptyPos(newpoint));
            space.SafeCreate(newpoint, sitevalue);
        }
        
        
        UpdateReverseCheckingSpace(poly.locs, newpoints);
        poly.locs = newpoints;
        
        BuildNewBonds(newpoints, bond_id_list);
        
        succ ++;
        bond_change += nbr_bond_inc;
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);
}

template <class S, class P, class M>
tuple<int, vector<int>> RubiMove<S,P,M>::ComputeBondInc(Polymer<P> poly, vector<P> newpoints)
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


template <class S, class P, class M>
int RubiMove<S,P,M>::ComputePSIncCrossLayer(Polymer<P> poly, vector<P> newpoints, int polyid)
{
    int old_nbr_ps = 0, new_nbr_ps = 0;
    for (auto oldpoint : poly.locs)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(oldpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                old_nbr_ps ++;
            }
        }
    }
    
    for (auto newpoint : newpoints)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(newpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                new_nbr_ps ++;
            }
        }
    }
    return new_nbr_ps - old_nbr_ps;
}

template <class S, class P, class M>
int RubiMove<S,P,M>::ComputePSIncSameLayer(Polymer<P> poly, vector<P> newpoints, int polyid)
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
            bool will_be_killed = space.GetRspacePoint(pt)[0] == polyid;
            bool currently_empty = space.EmptyPos(pt);
            
            bool will_be_created = PointInVector(pt, newpoints); // pt is in newpoints
            
            bool will_be_PS = (will_be_killed || currently_empty) && (! will_be_created);
            
            if (will_be_PS)
            {
                new_nbr_ps ++;
            }
        }
    }
    return new_nbr_ps - old_nbr_ps;
}



// To be modified

template <class S, class P, class M>
bool RubiMove<S, P, M>::ExecDragMove(int polyid)
{
    Polymer<P>& poly = space.Sumos[polyid];
    
    bool is_trimer; int epyc1id; int epyc2id;
    tie(is_trimer, epyc1id, epyc2id) = IsATrimerState(polyid);
    if (!is_trimer) return false; // The state is not a trimer, no move
    
    auto possible_moves = FilterForTriMoves(move.GetPossibleMoves(poly, polyid), epyc1id, epyc2id);
    
    if (possible_moves.empty()) return false; // No available move for the trimer, no move
    
    auto newpoints = ChooseMove(possible_moves); // Choose the target points of the rubisco
    
    TriMoveInfo<P> tri_move_info = CleanUpTheTriMove(polyid, epyc1id, epyc2id, newpoints); // Based on those info, determine the move details
    
    int nbr_ps_inc = ComputePSIncTriSameLayer(tri_move_info) + ComputePSIncTriCrossLayer(tri_move_info);
    
    if (generate_canonical<double, 10>(gen) < Weight(0, nbr_ps_inc))
    {
        for (auto oldpoint : poly.locs)
        {
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        
        for (auto newpoint : tri_move_info.rubisconewpoints)
        {
            assert(space.EmptyPos(newpoint));
            space.SafeCreate(newpoint, -1);
        }
        
        
        for (auto oldpoint : space.Sims[epyc1id].locs)
        {
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        
        for (auto oldpoint : space.Sims[epyc2id].locs)
        {
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        
        for (auto newpoint : tri_move_info.epyc1newpoints)
        {
            assert(space.EmptyPos(newpoint));
            space.SafeCreate(newpoint, 1);
        }
        
        for (auto newpoint : tri_move_info.epyc2newpoints)
        {
            assert(space.EmptyPos(newpoint));
            space.SafeCreate(newpoint, 1);
        }
        
        for (auto pt : space.Sims[epyc1id].locs)
        {
            assert(space.GetRspacePoint(pt)[0] == epyc1id);
        }
        
        for (auto pt : space.Sims[epyc2id].locs)
        {
            assert(space.GetRspacePoint(pt)[0] == epyc2id);
        }
        
        for (auto pt : poly.locs)
        {
            assert(space.GetRspacePoint(pt)[0] == polyid);
        }
        UpdateReverseCheckingSpace(poly.locs, tri_move_info.rubisconewpoints);
        UpdateReverseCheckingSpaceEpycCoMove(space.Sims[epyc1id].locs, tri_move_info.epyc1newpoints, space.Sims[epyc2id].locs, tri_move_info.epyc2newpoints);
        
        poly.locs = tri_move_info.rubisconewpoints;
        space.Sims[epyc1id].locs = tri_move_info.epyc1newpoints;
        space.Sims[epyc2id].locs = tri_move_info.epyc2newpoints;
        
        
        
        BuildNewBonds(tri_move_info.rubisconewpoints, tri_move_info.rubiscoinbond);
        succ ++;
        return true;
    }
    else
        return false;
}

template <class S, class P, class M>
int RubiMove<S,P,M>::ComputePSIncTriSameLayer(TriMoveInfo<P>& the_trimove)
{
    int ps_inc_rubisco = ComputePSIncSameLayer(space.Sumos[the_trimove.rubiscoid], the_trimove.rubisconewpoints, the_trimove.rubiscoid);
    int ps_inc_onlyepyc1 = ComputePSIncSameLayer(space.Sims[the_trimove.epyc1id], the_trimove.epyc1newpoints, the_trimove.epyc1id);
    int ps_inc_onlyepyc2 = ComputePSIncSameLayer(space.Sims[the_trimove.epyc2id], the_trimove.epyc2newpoints, the_trimove.epyc2id);
    int ps_epycs_interface = 0;
    for (auto newpoint : space.Sims[the_trimove.epyc1id].locs)
    {
        for (auto pt : space.Neighbor(newpoint))
        {
            if (space.GetRspacePoint(pt)[0] == the_trimove.epyc2id)
            {
                ps_epycs_interface ++;
            }
        }
    }
    
    return ps_inc_rubisco + ps_inc_onlyepyc1 + ps_inc_onlyepyc2 - ps_epycs_interface;
}


template <class S, class P, class M>
int RubiMove<S,P,M>::ComputePSIncTriCrossLayer(TriMoveInfo<P>& the_trimove)
{
    int old_nbr_ps = 0;
    
    for (auto oldpoint : space.Sumos[the_trimove.rubiscoid].locs)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(oldpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                old_nbr_ps ++;
            }
        }
    }
    
    for (auto oldpoint : space.Sims[the_trimove.epyc1id].locs)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(oldpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                old_nbr_ps ++;
            }
        }
    }
    
    for (auto oldpoint : space.Sims[the_trimove.epyc2id].locs)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(oldpoint)[0]))
        {
            if (space.EmptyPos(pt))
            {
                old_nbr_ps ++;
            }
        }
    }
    
    // Count new number of PS bonds
    int new_nbr_ps = 0;
    auto epycnewpoints = the_trimove.epyc1newpoints;
    epycnewpoints.insert(epycnewpoints.end(), the_trimove.epyc2newpoints.begin(), the_trimove.epyc2newpoints.end());
    for (auto newpoint : epycnewpoints)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(newpoint)[0]))
        {
            bool will_be_killed = space.GetRspacePoint(pt)[0] == the_trimove.rubiscoid;
            bool currently_empty = space.EmptyPos(pt);
            
            bool will_be_created = PointInVector(pt, the_trimove.rubisconewpoints); // pt is in newpoints
            
            bool will_be_PS = (will_be_killed || currently_empty) && (!will_be_created);
            if (will_be_PS)
            {
                new_nbr_ps ++;
            }
        }
    }
    
    for (auto newpoint : the_trimove.rubisconewpoints)
    {
        for (auto pt : space.Neighbor(space.BondNeighbor(newpoint)[0]))
        {
            bool will_be_killed = (space.GetRspacePoint(pt)[0] == the_trimove.epyc1id) || (space.GetRspacePoint(pt)[0] == the_trimove.epyc2id);
            bool currently_empty = space.EmptyPos(pt);
            
            bool will_be_created = PointInVector(pt, epycnewpoints); // pt is in newpoints
            
            bool will_be_PS = (will_be_killed || currently_empty) && (!will_be_created);
            if (will_be_PS)
            {
                new_nbr_ps ++;
            }
        }
    }
    
    return new_nbr_ps - old_nbr_ps;
}



#endif /* rubimoves_hpp */
