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
#include "space2d2l.cpp"
#include "endmove.hpp"
#include "snakemove.hpp"
#include "cornermove.hpp"

using std::invalid_argument;
using std::tuple;
using std::make_tuple;
using std::tie;


template <class S, class P, class M> class Move
{
private:
    S& space;
    M move;
    double beta;
public:
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




#endif /* move_hpp */
