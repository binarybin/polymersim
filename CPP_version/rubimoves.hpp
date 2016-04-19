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

template<class P> struct TriMoveInfo
{
    int rubiscoid;
    int epyc1id;
    int epyc2id;
    vector<P> rubisconewpoints;
    vector<P> epyc1newpoints;
    vector<P> epyc2newpoints;
    vector<int> rubiscoinbond;
};

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
    
    tuple<bool, int> ExecMove(int polyid, char polytyp);
    bool ExecTriMove(int polyid);
    
    void SetBeta(double setbeta)  {beta = setbeta;}
    void SetGamma(double setgamma) {gamma = setgamma;}
    
    vector<P> ChooseMove(vector<vector<P>> possible_moves)
    {
        uniform_int_distribution<> dis(0, (int)possible_moves.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    
    double Weight(int nbr_bond_inc, int nbr_ps_inc) {return exp(nbr_bond_inc * beta - nbr_ps_inc * gamma);}
    
    tuple<bool, int, int> IsATrimerState(int rubiscoid);
    TriMoveInfo<P> CleanUpTheTriMove(int rubiscoid, int epyc1id, int epyc2id, vector<P> new_rubisco_points);
    vector<vector<P>> FilterForTriMoves(vector<vector<P>> possible_moves, int epyc1id, int epyc2id);
    
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
    
    int ComputePSIncTri(TriMoveInfo<P>& the_trimove)
    {
        int ps_inc_rubisco = ComputePSInc(space.Sumos[the_trimove.rubiscoid], the_trimove.rubisconewpoints, the_trimove.rubiscoid);
        int ps_inc_onlyepyc1 = ComputePSInc(space.Sims[the_trimove.epyc1id], the_trimove.epyc1newpoints, the_trimove.epyc1id);
        int ps_inc_onlyepyc2 = ComputePSInc(space.Sims[the_trimove.epyc2id], the_trimove.epyc2newpoints, the_trimove.epyc2id);
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
    
    int nbr_ps_inc = ComputePSInc(poly, newpoints, polyid);
    
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
bool RubiMove<S, P, M>::ExecTriMove(int polyid)
{
    Polymer<P>& poly = space.Sumos[polyid];
    
    bool is_trimer; int epyc1id; int epyc2id;
    tie(is_trimer, epyc1id, epyc2id) = IsATrimerState(polyid);
    if (!is_trimer) return false; // The state is not a trimer, no move
    
    auto possible_moves = FilterForTriMoves(move.GetPossibleMoves(poly, polyid), epyc1id, epyc2id);
    
    if (possible_moves.empty()) return false; // No available move for the trimer, no move
    
    auto newpoints = ChooseMove(possible_moves); // Choose the target points of the rubisco
    
    TriMoveInfo<P> tri_move_info = CleanUpTheTriMove(polyid, epyc1id, epyc2id, newpoints); // Based on those info, determine the move details
    
    int nbr_ps_inc = ComputePSIncTri(tri_move_info);
    
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
        
        
        for (auto oldpoint : space.Sims[tri_move_info.epyc1id].locs)
        {
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        
        for (auto newpoint : tri_move_info.epyc1newpoints)
        {
            assert(space.EmptyPos(newpoint));
            space.SafeCreate(newpoint, 1);
        }
        
        for (auto oldpoint : space.Sims[tri_move_info.epyc2id].locs)
        {
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        
        for (auto newpoint : tri_move_info.epyc2newpoints)
        {
            assert(space.EmptyPos(newpoint));
            space.SafeCreate(newpoint, 1);
        }
        
        UpdateReverseCheckingSpace(poly.locs, tri_move_info.rubisconewpoints);
        UpdateReverseCheckingSpace(space.Sims[tri_move_info.epyc1id].locs, tri_move_info.epyc1newpoints);
        UpdateReverseCheckingSpace(space.Sims[tri_move_info.epyc2id].locs, tri_move_info.epyc2newpoints);
        poly.locs = tri_move_info.rubisconewpoints;
        space.Sims[tri_move_info.epyc1id].locs = tri_move_info.epyc1newpoints;
        space.Sims[tri_move_info.epyc2id].locs = tri_move_info.epyc2newpoints;
        
        BuildNewBonds(tri_move_info.rubisconewpoints, tri_move_info.rubiscoinbond);
        
        succ ++;
        
        return true;
    }
    else
        return false;
}


template <class S, class P, class M>
tuple<bool, int, int> RubiMove<S,P,M>::IsATrimerState(int rubiscoid)//use the simplist creteria: two epycs, full overlap
{
    std::set<int> epycids;
    for (auto point : space.Sumos[rubiscoid].locs)
    {
        int epycid = space.GetRspacePoint(point)[0];
        if (epycid != NOBOND)
            epycids.insert(epycid);
    }
    if (epycids.size() != 2)
        return false;
    
    for (int epycid : epycids)
    {
        for(auto point : space.Sims[epycid].locs)
        {
            if (space.GetRspacePoint(point)[0] != rubiscoid)    return make_tuple(false, 0, 0);
        }
    }
    vector<int> epids(epycids.begin(), epycids.end());
    
    return make_tuple(true, epids[0], epids[1]);
}

template <class S, class P, class M>
vector<vector<P>> RubiMove<S,P,M>::FilterForTriMoves(vector<vector<P>> possible_moves, int epyc1id, int epyc2id)
{
    vector<vector<P>> result;
    for (vector<P> newpoints : possible_moves)
    {
        bool good = true;
        for (P newpoint : newpoints)
        {
            int newepycid = space.GetRspacePoint(space.BondNeighbor(newpoint)[0])[0];
            if (newepycid != NOBOND && newepycid != epyc1id && newepycid != epyc2id) // The new point is occupied by an epyc other than epyc1 and epyc2
                good = false;
        }
        if (good)
            result.push_back(newpoints);
    }
    return result;
}

template <class S, class P, class M>
TriMoveInfo<P> RubiMove<S,P,M>::CleanUpTheTriMove(int rubiscoid, int epyc1id, int epyc2id, vector<P> new_rubisco_points)
{
    TriMoveInfo<P> result;
    
    result.epyc1id = epyc1id;
    result.epyc2id = epyc2id;
    result.rubisconewpoints = new_rubisco_points;
    
    
    for (auto pointepyc : space.Sims[result.epyc1id].locs)
    {
        auto pos_rubi = space.GetRspacePoint(space.BondNeighbor(pointepyc)[0])[1];
        assert(pos_rubi >= 0 && pos_rubi<space.LSumo);
        auto new_epyc_point = space.BondNeighbor(new_rubisco_points[pos_rubi])[0];
        result.epyc1newpoints.push_back(new_epyc_point);
    }
    
    for (auto pointepyc : space.Sims[result.epyc2id].locs)
    {
        auto pos_rubi = space.GetRspacePoint(space.BondNeighbor(pointepyc)[0])[1];
        assert(pos_rubi >= 0 && pos_rubi<space.LSumo);
        auto new_epyc_point = space.BondNeighbor(new_rubisco_points[pos_rubi])[0];
        result.epyc2newpoints.push_back(new_epyc_point);
    }
    
    for (int i = 0; i < space.LSumo; i++)
    {
        if (space.InABond(space.Sumos[rubiscoid].locs[i]))
        {
            result.rubiscoinbond.push_back(i);
        }
    }
    return result;
    
    
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













#endif /* rubimoves_hpp */