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

int static NOBOND = -1;


class Pos2d1l
{
public:
    int x, y;
    Pos2d1l(int input_x, int input_y): x(input_x), y(input_y) {}
    Pos2d1l(): x(0), y(0) {}
    bool IsNeighborOf(Pos2d1l& other)
    {
        return (x == other.x) != (y == other.y);
    }
};

class Pos2d2l
{
public:
    int x, y;
    bool siml; // in the sim layer
    Pos2d2l(int input_x, int input_y, bool input_layer): x(input_x), y(input_y), siml(input_layer) {}
    Pos2d2l(): x(0), y(0), siml(true){}
    bool IsNeighborOf(Pos2d2l& other)
    {
        return (siml == other.siml) && (x == other.x) != (y == other.y);
    }
    
    Pos2d2l OtherLayer()
    {
        return Pos2d2l(x, y, !siml);
    }
};

class Pos3d1l
{
public:
    int x, y, z;
    bool siml; // in the sim layer
    Pos3d1l(int input_x, int input_y, int input_z, bool input_layer):
    x(input_x), y(input_y), z(input_z){}
    Pos3d1l(): x(0), y(0), z(0) {}
};

class Pos3d2l
{
public:
    int x, y, z;
    bool siml; // in the sim layer
    Pos3d2l(int input_x, int input_y, int input_z, bool input_layer):
        x(input_x), y(input_y),
        z(input_z), siml(input_layer) {}
    Pos3d2l(): x(0), y(0), z(0), siml(true){}
};



template <class P>
class Polymer {
public:
    int poly_id;
    vector<P> locs;
    
    Polymer(int id, size_t size) : poly_id(id)
    {
        this->locs.resize(size);
    }
    Polymer(){}
};

template <class P>
struct DragMoveInfo
{
    vector<int> rubiscoIDs;
    vector<int> epycIDs;
    vector<vector<P>> rubiscoNewPoints;
    vector<vector<P>> epycNewPoints;
    vector<vector<int>> rubiscoInBondIDs;
};



#endif /* position_h */
