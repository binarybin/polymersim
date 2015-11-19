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
#include "space.hpp"
#include "space2d1l.hpp"

using std::invalid_argument;
using std::tuple;
using std::make_tuple;
using std::tie;


class PossibleMoves
{
public:
    int i;
    bool isEmpty();
    tuple<int, Position> RandSelect();
};


class Move
{
    Space& space;
    double beta;
public:
    Move(Space &thespace): space(thespace), beta(0){}
    SetBeta(double setbeta)  {this->beta = setbeta;}
    tuple<bool, int> ExecMove(int polyid, char polytyp)
    {
        Polymer& poly = space.Sims[1];
        int sitevalue = 0;
        
        if (polytyp == 'i')
        {
            poly = space.Sims[polyid];
            sitevalue = 1;
        }
        else if(polytyp == 'u')
        {
            poly = space.Sumos[polyid];
            sitevalue = -1;
        }
        else
            throw invalid_argument("Polymer type undefined.");
        
        int nbr_bond_inc = 0;
        PossibleMoves &possible_moves = GetPossibleMoves(Polymer& poly);
        if (possible_moves.isEmpty())
            return make_tuple(false, nbr_bond_inc);
        
        int pointid = 0;
        Position newpoint;
        tie(pointid, newpoint) = possible_moves.RandSelect();
        Position oldpoint = poly.locs[pointid];
        
        if (space.InABond(oldpoint))
        {
            nbr_bond_inc -= 1;
        }
        
        vector<Position> bond_choice;
        
        for (auto bpoint : space.Neighbor(newpoint))
        {
            if (space.CanBuildBondTentative(sitevalue, bpoint))
            {
                bond_choice.push_back(bpoint);
            }
        }
        
        if (!bond_choice.empty())
        {
            nbr_bond_inc += 1;
        }
        
        if ((double)rand()/(RAND_MAX) < this->Weight(nbr_bond_inc))
        {
            space.SafeRemove(oldpoint);
            space.SafeCreate(newpoint, sitevalue);
            this->UpdatePolymer(poly, pointid, newpoint);
            this->UpdateReverseCheckingSpace(oldpoint, newpoint);
            
            if (!bond_choice.empty())
            {
                Position bpoint = this->ChooseBond(bond_choice);
                space.CreateBond(newpoint, bpoint);
            }
            
            return make_tuple(true, nbr_bond_inc);
            
        }
        else
        {
            return make_tuple(false, nbr_bond_inc);
        }
        
    }
    
    Position ChooseBond(vector<Position> bond_choice);
    virtual void UpdatePolymer(Polymer& poly, int pointid, Position& newpoint);
    virtual void UpdateReverseCheckingSpace(Position& oldpoint, Position& newpoint);
    virtual PossibleMoves& GetPossibleMoves = 0;
    
    double Weight(int nbr_bond_inc) {return exp(nbr_bond_inc * beta);}
    
    
};
#endif /* move_hpp */
