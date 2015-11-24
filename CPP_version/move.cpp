//
//  move.cpp
//  polymersim
//
//  Created by Bin Xu on 11/19/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <stdio.h>
#include "move.hpp"

tuple<bool, int> Move::ExecMove(int polyid, char polytyp)
{
    Polymer& poly = polytyp=='i' ? space.Sims[polyid-1] : space.Sumos[polyid-1];
    int sitevalue = polytyp=='i' ? 1 : -1;
    
    int nbr_bond_inc = 0;
    vector<tuple<int, Position>> possible_moves = GetPossibleMoves(poly);

    if (possible_moves.empty()) return make_tuple(false, nbr_bond_inc);
    
    int pointid = 0; Position newpoint; tie(pointid, newpoint) = ChooseMove(possible_moves);
    
    Position oldpoint = poly.locs[pointid];
    
    if (space.InABond(oldpoint)) nbr_bond_inc -= 1;
    
    vector<Position> bond_choice;
    
    for (auto bpoint : space.Neighbor(newpoint))
        if (space.CanBuildBondTentative(bpoint, sitevalue))
            bond_choice.push_back(bpoint);
    
    if (!bond_choice.empty()) nbr_bond_inc += 1;
    
    if ((double)rand()/(RAND_MAX) < Weight(nbr_bond_inc))
    {
        space.SafeRemove(oldpoint);
        space.SafeCreate(newpoint, sitevalue);
        UpdatePolymer(poly, pointid, newpoint);
        UpdateReverseCheckingSpace(oldpoint, newpoint, poly);
        
        if (!bond_choice.empty())
        {
            Position bpoint = ChooseBond(bond_choice);
            space.CreateBond(newpoint, bpoint);
        }
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);
}
