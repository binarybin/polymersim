//
//  position.h
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#ifndef position_h
#define position_h

class Position
{
public:
    int x, y, z;
    char layer;
    Position(int input_x, int input_y): x(input_x), y(input_y), z(0), layer(0) {}
    Position(): x(0), y(0), z(0), layer(0) {}
};

#endif /* position_h */
