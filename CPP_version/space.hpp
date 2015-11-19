//
//  space.hpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef space_hpp
#define space_hpp

#include <vector>
#include "position.h"
using std::vector;

int static NOBOND = -1;

class Polymer {
public:
    int poly_id;
    vector<Position> locs;
    
    Polymer(int id, size_t size) : poly_id(id)
    {
        this->locs.resize(size);
    }
    Polymer(){}
};

class Space {
protected:
    int NSim, LSim, NSumo, LSumo;
    size_t Lx, Ly, Lz;
    int SimId, SumoId;
    
    vector<Polymer> Sims;
    vector<Polymer> Sumos;
    
public:
    Space(int nsim, int nsumo, int lsim, int lsumo, size_t lx, size_t ly, size_t lz)
    : NSim(nsim), NSumo(nsumo), LSim(lsim), LSumo(lsumo), Lx(lx), Ly(ly), Lz(lz),
    SimId(0), SumoId(0), Sims(nsim), Sumos(nsumo) {}
    
    virtual void Initialize() = 0;
};


#endif /* space_hpp */
