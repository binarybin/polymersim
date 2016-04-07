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
#include "space2d1l.hpp"
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
    tuple<bool, int> ExecCoMove(int simid);
    void SetBeta(double setbeta)  {beta = setbeta;}
    void SetGamma(double setgamma) {gamma = setgamma;}
    
    P ChooseBond(vector<P> bond_choice)
    {
        uniform_int_distribution<> dis(0, (int)bond_choice.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < bond_choice.size());
        return bond_choice[id];
    }
    
    tuple<int, P> ChooseMove(vector<tuple<int, P>> possible_moves)
    {
        uniform_int_distribution<> dis(0, (int)possible_moves.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    tuple<int, P, int, int, P> ChooseCoMove(vector<tuple<int, Pos2d2l, int, int, Pos2d2l>> possible_moves)
    {
        uniform_int_distribution<> dis(0, (int)possible_moves.size()-1);
        size_t id = dis(gen);
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    
    double Weight(int nbr_bond_inc, int nbr_ps_inc) {return exp(nbr_bond_inc * beta - nbr_ps_inc * gamma);}
   
};


template <class S, class P, class M>
tuple<bool, int> Move<S, P, M>::ExecMove(int polyid, char polytyp)
{
    Polymer<P>& poly = polytyp=='i' ? space.Sims[polyid] : space.Sumos[polyid];
    int sitevalue = polytyp=='i' ? 1 : -1;
    
    int nbr_bond_inc = 0;
    vector<tuple<int, P>> possible_moves = move.GetPossibleMoves(poly);
    
    if (possible_moves.empty()) return make_tuple(false, nbr_bond_inc);
    
    int pointid = 0; P newpoint; tie(pointid, newpoint) = ChooseMove(possible_moves);
    
    P oldpoint = poly.locs[pointid];
    
    if (space.InABond(oldpoint)) nbr_bond_inc --;
    
    vector<P> bond_choice;
    
    // The version for SIM-SUMO system
    /*for (auto bpoint : space.BondNeighbor(newpoint))
        if (space.CanBuildBondTentative(bpoint, sitevalue))
            bond_choice.push_back(bpoint);*/
    
    // The version for Rubisco system
    for (auto bpoint : space.BondNeighbor(newpoint))
    {
        if (space.CanBuildRubiscoBondTentative(poly, polytyp, polyid, pointid, bpoint, sitevalue))
            bond_choice.push_back(bpoint);
    }
    
    if (!bond_choice.empty()) nbr_bond_inc ++;
    
    int nbr_ps_inc = 0;
    for (auto& pt : space.Neighbor(oldpoint))
    {
        if (space.EmptyPos(pt))
        {
            nbr_ps_inc --;
        }
        else
        {
            nbr_ps_inc ++;
        }
    }
    
    for (auto& pt : space.Neighbor(newpoint))
    {
        if (space.EmptyPos(pt))
        {
            nbr_ps_inc ++;
        }
        else
        {
            nbr_ps_inc --;
        }
    }
    
    if (generate_canonical<double, 10>(gen) < Weight(nbr_bond_inc, nbr_ps_inc))
    {
        assert(!space.EmptyPos(oldpoint));
        space.SafeRemove(oldpoint);
        space.SafeCreate(newpoint, sitevalue);
        move.UpdatePolymer(poly, pointid, newpoint);
        move.UpdateReverseCheckingSpace(oldpoint, newpoint, poly);
        
        if (!bond_choice.empty())
        {
            P bpoint = ChooseBond(bond_choice);
            space.CreateBond(newpoint, bpoint);
        }
        succ += 1;
        bond_change += nbr_bond_inc;
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);
}

template <class S, class P, class M>
tuple<bool, int> Move<S, P, M>::ExecCoMove(int simid)
{
    Polymer<P>& poly = space.Sims[simid];
    
    vector<tuple<int, P, int, int, P>> possible_comoves = move.GetPossibleCoMoves(poly);
    //move sims[simid]$1 to $2, move sumos[$3]$4 to $5
    
    if (possible_comoves.empty()) return make_tuple(false, 0);
    
    int simpointid = 0; int sumoid = 0; int sumopointid = 0; P newsimpoint; P newsumopoint;
    tie(simpointid, newsimpoint, sumoid, sumopointid, newsumopoint) = ChooseCoMove(possible_comoves);
    
    P oldsimpoint = poly.locs[simpointid];
    P oldsumopoint = oldsimpoint.OtherLayer();
    
    int nbr_bond_inc = 0;
    
    if(!space.InABond(oldsimpoint))
        nbr_bond_inc ++;
    
    int nbr_ps_inc = 0;
    for (auto& pt : space.Neighbor(oldsimpoint))
    {
        if (space.EmptyPos(pt))
        {
            nbr_ps_inc --;
        }
        else
        {
            nbr_ps_inc ++;
        }
    }
    
    for (auto& pt : space.Neighbor(newsimpoint))
    {
        if (space.EmptyPos(pt))
        {
            nbr_ps_inc ++;
        }
        else
        {
            nbr_ps_inc --;
        }
    }
    for (auto& pt : space.Neighbor(oldsumopoint))
    {
        if (space.EmptyPos(pt))
        {
            nbr_ps_inc --;
        }
        else
        {
            nbr_ps_inc ++;
        }
    }
    
    for (auto& pt : space.Neighbor(newsumopoint))
    {
        if (space.EmptyPos(pt))
        {
            nbr_ps_inc ++;
        }
        else
        {
            nbr_ps_inc --;
        }
    }


    
    if (generate_canonical<double, 10>(gen) < Weight(nbr_bond_inc, nbr_ps_inc))
    {
        space.SafeRemove(oldsimpoint);
        space.SafeRemove(oldsumopoint);
        space.SafeCreate(newsimpoint, 1);
        space.SafeCreate(newsumopoint, -1);
        
        move.UpdatePolymer(poly, simpointid, newsimpoint);
        move.UpdatePolymer(space.Sumos[sumoid], sumopointid, newsumopoint);
        
        move.UpdateReverseCheckingSpace(oldsimpoint, newsimpoint, poly);
        move.UpdateReverseCheckingSpace(oldsumopoint, newsumopoint, space.Sumos[sumoid]);
        
        space.CreateBond(newsimpoint, newsumopoint);
        succ += 1;
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);

    
    
}


#endif /* move_hpp */
