//
//  position.h
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef position_h
#define position_h

#include <vector>
using std::vector;

class Position
{
public:
    int x, y;
    Position(int input_x, int input_y): x(input_x), y(input_y) {}
    Position(): x(0), y(0) {}
};

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


#endif /* position_h */
