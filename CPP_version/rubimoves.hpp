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
#include <utility>
#include <random>
#include <set>
#include <algorithm>
#include <unordered_map>
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
    
    DragMoveInfo<P> ChooseMove(vector<DragMoveInfo<P>> possible_moves)
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
    
    
    bool ExecDragMove(int polyid, bool blob);
    void UpdateRealSpaceDrag(DragMoveInfo<P>& the_move);
    void UpdateRSpaceDrag(DragMoveInfo<P>& the_move);
    void UpdatePolymersDrag(DragMoveInfo<P>& the_move);
    void BuildNewBondsDrag(DragMoveInfo<P>& the_move);
    pair<int, int> ComputePSIncDrag(DragMoveInfo<P>& the_move);
};


template <class S, class P, class M>
tuple<bool, int> RubiMove<S, P, M>::ExecMove(int polyid, char polytyp)
{
    Polymer<P>& poly = polytyp=='i' ? space.Sims[polyid] : space.Sumos[polyid];
    
    auto possible_moves = move.GetPossibleMoves(poly, polyid);
    
    if (possible_moves.empty()) return make_tuple(false, 0);
    
    auto newpoints = ChooseMove(possible_moves);
    
    int nbr_bond_inc = 0; vector<int> bond_id_list;
    tie(nbr_bond_inc, bond_id_list) = ComputeBondInc(poly, newpoints);
    
    int nbr_ps_inc = ComputePSIncCrossLayer(poly, newpoints, polyid) + ComputePSIncSameLayer(poly, newpoints, polyid);
    
    if (generate_canonical<double, 10>(gen) < Weight(nbr_bond_inc, nbr_ps_inc))
    {
        // NOTE: for this part, a better way is to dephosphorylate all of them, finish the move, and phosphorylate them. So sitevalues will still be +/- 1
        vector<int> sitevalues;
        
        for (auto oldpoint : poly.locs)
        {
            assert(!space.EmptyPos(oldpoint));
            int sitevalue = space.GetSpacePoint(oldpoint);
            if (abs(sitevalue) == 2)
            {
                sitevalue /= 2;
            }
            sitevalues.push_back(sitevalue);
            space.SafeRemove(oldpoint);
        }
        
        int idx = 0;
        for (auto newpoint : newpoints)
        {
            int sitevalue = sitevalues[idx];
            idx ++;
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


// TODO: modify this

/*
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
        if (epycid != NOBOND )
            if (!space.Phosphorylated(bpoint))
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
*/
// The version that does not impose the non-two-end-interacting rule
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
    
    vector<int> result_pos;
    
    for (int i = 0; i < space.LSumo; i++)
    {
        P bpoint = space.BondNeighbor(newpoints[i])[0];
        int epycid = space.GetRspacePoint(bpoint)[0];
        if (epycid != NOBOND )
            if (!space.Phosphorylated(bpoint))
            {
                result_pos.push_back(i);
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

template <class S, class P, class M>
bool RubiMove<S, P, M>::ExecDragMove(int polyid, bool blob)
{

    vector<DragMoveInfo<P>> possible_moves;
    if (blob)
    {
        possible_moves = move.GetPossibleBlobMoves(polyid);
    }
    else
    {
        possible_moves = move.GetPossibleDragMoves(polyid);
    }
    
    
    
    if (possible_moves.empty()) return false; // No available move for the trimer, no move
    
    auto the_move = ChooseMove(possible_moves); // Choose the target points of the rubisco
#ifndef NDEBUG
    PrintSpace(cout, space);
#endif
    auto nbr_ps_inc_pair = ComputePSIncDrag(the_move);
#ifndef NDEBUG
    PrintSpace(cout, space);
#endif
    int nbr_ps_inc = nbr_ps_inc_pair.first + nbr_ps_inc_pair.second;
    
    
    if (generate_canonical<double, 10>(gen) < Weight(0, nbr_ps_inc))
    {
        UpdateRealSpaceDrag(the_move);
        UpdateRSpaceDrag(the_move);
        UpdatePolymersDrag(the_move);
        BuildNewBondsDrag(the_move);
        succ ++;
        return true;
    }
    else
        return false;
}
template <class S, class P, class M>
void RubiMove<S, P, M>::UpdateRealSpaceDrag(DragMoveInfo<P>& the_move)
{
    vector<vector<bool>> phosepyc;

    for (auto epycid : the_move.epycIDs)
    {
        vector<bool> phos;
        for (auto oldpoint : space.Sims[epycid].locs)
        {
            if (abs(space.GetSpacePoint(oldpoint)) == 3)
                phos.push_back(true);
            else
                phos.push_back(false);
#ifndef NDEBUG
            cout<<epycid<<" "<<oldpoint.x<<" "<<oldpoint.y<<" "<<oldpoint.siml<<" "<<space.GetSpacePoint(oldpoint)<<endl;
#endif
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        phosepyc.push_back(phos);
    }
#ifndef NDEBUG
    cout<<"EPYC phos"<<endl;
    for (auto line : phosepyc)
    {
        for (auto p : line)
        {
            cout<<p<<'\t';
        }
        cout<<endl;
    }
#endif
    
    int epycidx = 0;
    for (const auto& epycnewpoints : the_move.epycNewPoints)
    {
        int idx = 0;
        for (auto newpoint : epycnewpoints)
        {
            assert(space.EmptyPos(newpoint));
            if (phosepyc[epycidx][idx])
                space.SafeCreate(newpoint, 3);
            else
                space.SafeCreate(newpoint, 1);
            
            idx ++;
        }
        epycidx ++;
    }
    
    vector<vector<bool>> phosrubi;
    for (auto rubiid : the_move.rubiscoIDs)
    {
        vector<bool> phos;
        for (auto oldpoint : space.Sumos[rubiid].locs)
        {
            if (abs(space.GetSpacePoint(oldpoint)) == 3)
                phos.push_back(true);
            else
                phos.push_back(false);
            assert(!space.EmptyPos(oldpoint));
            space.SafeRemove(oldpoint);
        }
        phosepyc.push_back(phos);
    }
    
    int rubiidx = 0;
    for (const auto& rubinewpoints : the_move.rubiscoNewPoints)
    {
        int idx = 0;
        for (auto newpoint : rubinewpoints)
        {
            assert(space.EmptyPos(newpoint));
            if (phosepyc[epycidx][idx])
                space.SafeCreate(newpoint, -3);
            else
                space.SafeCreate(newpoint, -1);
            idx ++;
        }
        rubiidx ++;
    }
}


template <class S, class P, class M>
void RubiMove<S, P, M>::UpdateRSpaceDrag(DragMoveInfo<P>& the_move)
{
    for (auto rubiid : the_move.rubiscoIDs)
        for (auto pt : space.Sumos[rubiid].locs) assert(space.GetRspacePoint(pt)[0] == rubiid);

    for (int i = 0; i < the_move.rubiscoIDs.size(); i++)
    {
        int polyid = the_move.rubiscoIDs[i];
        const auto& oldpoints = space.Sumos[polyid].locs;
        for (auto oldpt: oldpoints)
            space.SetRspacePoint(oldpt, NOBOND, NOBOND);
    }
    for (int i = 0; i < the_move.rubiscoIDs.size(); i++)
    {
        int polyid = the_move.rubiscoIDs[i];
        auto& newpoints = the_move.rubiscoNewPoints[i];
        for (int j = 0; j < newpoints.size(); j++)
            space.SetRspacePoint(newpoints[j], polyid, j);
    }
    
    for (auto epycid : the_move.epycIDs)
        for (auto pt : space.Sims[epycid].locs) assert(space.GetRspacePoint(pt)[0] == epycid);
    for (int i = 0; i < the_move.epycIDs.size(); i++)
    {
        int polyid = the_move.epycIDs[i];
        const auto& oldpoints = space.Sims[polyid].locs;
        for (auto oldpt: oldpoints)
            space.SetRspacePoint(oldpt, NOBOND, NOBOND);
    }
    for (int i = 0; i < the_move.epycIDs.size(); i++)
    {
        int polyid = the_move.epycIDs[i];
        auto& newpoints = the_move.epycNewPoints[i];
        for (int j = 0; j < newpoints.size(); j++)
            space.SetRspacePoint(newpoints[j], polyid, j);
    }
}

template <class S, class P, class M>
void RubiMove<S, P, M>::UpdatePolymersDrag(DragMoveInfo<P>& the_move)
{
    for (int i = 0; i < the_move.epycIDs.size(); i++)
        space.Sims[the_move.epycIDs[i]].locs = the_move.epycNewPoints[i];
    for (int i = 0; i < the_move.rubiscoIDs.size(); i++)
        space.Sumos[the_move.rubiscoIDs[i]].locs = the_move.rubiscoNewPoints[i];
}

template <class S, class P, class M>
void RubiMove<S, P, M>::BuildNewBondsDrag(DragMoveInfo<P>& the_move)
{
    for (int i = 0; i < the_move.rubiscoIDs.size(); i++)
    {
        for(int id : the_move.rubiscoInBondIDs[i])
        {
            space.CreateBond(the_move.rubiscoNewPoints[i][id]);
        }
    }
}


template <class S, class P, class M>
pair<int, int> RubiMove<S,P,M>::ComputePSIncDrag(DragMoveInfo<P>& the_move)
{
    // step1 : compute the number of old PS bonds
    int ps_rubisco_old = 0;
    int ps_rubi_epyc_old = 0;
    for (auto rubiID: the_move.rubiscoIDs)
    {
        const auto& OldPts = space.Sumos[rubiID].locs;
        for (auto pt : OldPts)
        {
            for (auto nbrpt : space.Neighbor(pt))
                ps_rubisco_old += (space.EmptyPos(nbrpt)? 1 : 0);
            for (auto nbrpt : space.Neighbor(space.BondNeighbor(pt)[0]))
                ps_rubi_epyc_old += (space.EmptyPos(nbrpt)? 1 : 0);
        }
    }
    
    int ps_epyc_old = 0;
    int ps_epyc_rubi_old = 0;
    for (auto epycID : the_move.epycIDs)
    {
        const auto& OldPts = space.Sims[epycID].locs;
        for (auto pt : OldPts)
        {
            for (auto nbrpt : space.Neighbor(pt))
                ps_epyc_old += (space.EmptyPos(nbrpt)? 1 : 0);
            for (auto nbrpt : space.Neighbor(space.BondNeighbor(pt)[0]))
                ps_epyc_rubi_old += (space.EmptyPos(nbrpt)? 1 : 0);
        }
    }
    
    // step2 : make tentative moves in space.space
    vector<int> rubiVal;
    for (auto rubiID: the_move.rubiscoIDs)
    {
        const auto& OldPts = space.Sumos[rubiID].locs;
        for (auto pt : OldPts)
        {
            rubiVal.push_back(space.GetSpacePoint(pt));
            space.SetSpacePoint(pt, 0);
        }
    }
    int rubicount = 0;
    for (const auto& NewPts : the_move.rubiscoNewPoints)
    {
        for (auto pt : NewPts)
        {
            space.SetSpacePoint(pt, rubiVal[rubicount]);
            rubicount++;
        }
    }
    vector<int> epycVal;
    for (auto epycID : the_move.epycIDs)
    {
        const auto& OldPts = space.Sims[epycID].locs;
        for (auto pt : OldPts)
        {
            epycVal.push_back(space.GetSpacePoint(pt));
            space.SetSpacePoint(pt, 0);
        }
    }
    int epyccount = 0;
    for (const auto& NewPts : the_move.epycNewPoints)
    {
        for (auto pt : NewPts)
        {
            space.SetSpacePoint(pt, epycVal[epyccount]);
            epyccount++;
        }
    }
    
    // step3 : Compute new number of bonds
    int ps_rubisco_new = 0;
    int ps_rubi_epyc_new = 0;
    for (const auto& NewPts : the_move.rubiscoNewPoints)
    {
        for (auto pt : NewPts)
        {
            for (auto nbrpt : space.Neighbor(pt))
                ps_rubisco_new += (space.EmptyPos(nbrpt)? 1 : 0);
            for (auto nbrpt : space.Neighbor(space.BondNeighbor(pt)[0]))
                ps_rubi_epyc_new += (space.EmptyPos(nbrpt)? 1 : 0);
        }
    }
    
    int ps_epyc_new = 0;
    int ps_epyc_rubi_new = 0;
    for (const auto& NewPts : the_move.epycNewPoints)
    {
        for (auto pt : NewPts)
        {
            for (auto nbrpt : space.Neighbor(pt))
                ps_epyc_new += (space.EmptyPos(nbrpt)? 1 : 0);
            for (auto nbrpt : space.Neighbor(space.BondNeighbor(pt)[0]))
                ps_epyc_rubi_new += (space.EmptyPos(nbrpt)? 1 : 0);
        }
    }

    // step4 : restore the original space
    rubicount = 0;
    for (const auto& NewPts : the_move.rubiscoNewPoints)
        for (auto pt : NewPts)
            space.SetSpacePoint(pt, 0);

    for (auto rubiID: the_move.rubiscoIDs)
    {
        const auto& OldPts = space.Sumos[rubiID].locs;
        for (auto pt : OldPts)
        {
            space.SetSpacePoint(pt, rubiVal[rubicount]);
            rubicount++;
        }
    }
    
    epyccount = 0;
    for (const auto& NewPts : the_move.epycNewPoints)
        for (auto pt : NewPts)
            space.SetSpacePoint(pt, 0);

    for (auto epycID : the_move.epycIDs)
    {
        const auto& OldPts = space.Sims[epycID].locs;
        for (auto pt : OldPts)
        {
            space.SetSpacePoint(pt, epycVal[epyccount]);
            epyccount++;
        }
    }
    return std::make_pair(ps_epyc_new + ps_rubisco_new - ps_epyc_old - ps_rubisco_old,
                          ps_epyc_rubi_new + ps_rubi_epyc_new - ps_epyc_rubi_old - ps_rubi_epyc_old);
    
}


#endif /* rubimoves_hpp */
