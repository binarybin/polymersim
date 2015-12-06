//
//  move.cpp
//  polymersim
//
//  Created by Bin Xu on 12/6/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include "move.hpp"


template <class S, class P, class M>
tuple<bool, int> Move<S, P, M>::ExecMove(int polyid, char polytyp)
{
    Polymer<P>& poly = polytyp=='i' ? space.Sims[polyid-1] : space.Sumos[polyid-1];
    int sitevalue = polytyp=='i' ? 1 : -1;
    
    int nbr_bond_inc = 0;
    vector<tuple<int, P>> possible_moves = move.GetPossibleMoves(poly);
    
    if (possible_moves.empty()) return make_tuple(false, nbr_bond_inc);
    
    int pointid = 0; P newpoint; tie(pointid, newpoint) = ChooseMove(possible_moves);
    
    P oldpoint = poly.locs[pointid];
    
    if (space.InABond(oldpoint)) nbr_bond_inc -= 1;
    
    vector<P> bond_choice;
    
    for (auto bpoint : space.BondNeighbor(newpoint))
        if (space.CanBuildBondTentative(bpoint, sitevalue))
            bond_choice.push_back(bpoint);
    
    if (!bond_choice.empty()) nbr_bond_inc += 1;
    
    if ((double)rand()/(RAND_MAX) < Weight(nbr_bond_inc))
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
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);
}
