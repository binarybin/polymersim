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
    
    tuple<bool, int> ExecMove(int polyid, char polytyp);
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
    tie(nbr_bond_inc, bond_id_list) = move.ComputeBondInc(poly, newpoints);
    
    int nbr_ps_inc = move.ComputePSInc(poly, newpoints, polyid);
    
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
        
        
        move.UpdateReverseCheckingSpace(poly.locs, newpoints);
        move.UpdatePolymer(poly, newpoints);
        
        move.BuildNewBonds(newpoints, bond_id_list);
        
        succ ++;
        bond_change += nbr_bond_inc;
        
        return make_tuple(true, nbr_bond_inc);
    }
    else
        return make_tuple(false, nbr_bond_inc);
}


#endif /* rubimoves_hpp */
