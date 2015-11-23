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

using std::invalid_argument;
using std::tuple;
using std::make_tuple;
using std::tie;


class Move
{
public:
    Space2D1L& space;
    double beta;

    Move(Space2D1L &thespace): space(thespace), beta(0){}
    
    tuple<bool, int> ExecMove(int polyid, char polytyp);
    void SetBeta(double setbeta)  {this->beta = setbeta;}
    
    Position ChooseBond(vector<Position> bond_choice)
    {
        double ran = (double)rand()/(RAND_MAX);
        size_t id = floor(ran * bond_choice.size());
        assert(id >= 0 && id < bond_choice.size());
        return bond_choice[id];
    }
    
    tuple<int, Position> ChooseMove(vector<tuple<int, Position>> possible_moves)
    {
        double ran = (double)rand()/(RAND_MAX);
        size_t id = floor(ran * possible_moves.size());
        assert(id >= 0 && id < possible_moves.size());
        return possible_moves[id];
    }
    
    virtual void UpdatePolymer(Polymer& poly, int pointid, Position& newpoint) = 0;
    virtual void UpdateReverseCheckingSpace(Position& oldpoint, Position& newpoint) = 0;
    virtual vector<tuple<int, Position>> GetPossibleMoves(Polymer& poly) = 0;
    
    double Weight(int nbr_bond_inc) {return exp(nbr_bond_inc * beta);}
};


#endif /* move_hpp */
