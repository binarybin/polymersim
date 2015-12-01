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
#include "space2d1l.hpp"
#include "endmove.hpp"
#include "snakemove.hpp"
#include "cornermove.hpp"

using std::invalid_argument;
using std::tuple;
using std::make_tuple;
using std::tie;


template <class S, class P, class M> class Move
{
public:
    S& space;
    M move;
    double beta;

    Move(S &thespace): move(thespace), space(thespace), beta(0){}
    
    tuple<bool, int> ExecMove(int polyid, char polytyp);
    void SetBeta(double setbeta)  {this->beta = setbeta;}
    
    P ChooseBond(vector<P> bond_choice)
    {
        double ran = (double)rand()/(RAND_MAX);
        size_t id = floor(ran * bond_choice.size());
        assert(id >= 0 && id < bond_choice.size());
        return bond_choice[id];
    }
    
    tuple<int, P> ChooseMove(vector<tuple<int, P>> possible_moves)
    {
        double ran = (double)rand()/(RAND_MAX);
        size_t id = floor(ran * possible_moves.size());
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    
    double Weight(int nbr_bond_inc) {return exp(nbr_bond_inc * beta);}
   
};


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


#endif /* move_hpp */
