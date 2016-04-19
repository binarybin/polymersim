//
//  rotationmove.hpp
//  polymersim
//
//  Created by Bin Xu on 4/11/16.
//  Copyright © 2016 Bin Xu. All rights reserved.
//

#ifndef rotationmove_hpp
#define rotationmove_hpp


#include <tuple>
#include <vector>
#include <cassert>
#include <map>
#include <cstdlib>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P>
class RotationMove
{
private:
    S& space;
public:
    RotationMove(S& thespace): space(thespace) {}
    vector<vector<Pos2d2l>> GetPossibleMoves(Polymer<P>& polysim, int polyid);
};


template<>
vector<vector<Pos2d2l>> RotationMove<Space2D2L, Pos2d2l>::GetPossibleMoves(Polymer<Pos2d2l>& polysim, int polyid)
{
    vector<vector<Pos2d2l>> possible_moves;
    vector<Pos2d2l> pure_new_points;
    
    int Lx = (int)(space.Lx);
    int Ly = (int)(space.Ly);
    
    auto Tentative = [Lx, Ly, polysim] (int first, int second)
    {
        Pos2d2l tentative_point(polysim.locs[second]);
        tentative_point.x += (Lx + polysim.locs[second].x - polysim.locs[first].x);
        tentative_point.y += (Ly + polysim.locs[second].y - polysim.locs[first].y);
        tentative_point.x %= Lx;
        tentative_point.y %= Ly;
        return tentative_point;
    };
    
    auto tent7 = Tentative(2,3);
    if (space.space[1][tent7.x][tent7.y] != 0)
    {
        return possible_moves;
    }
    auto tent6 = Tentative(5,4);
    if (space.space[1][tent6.x][tent6.y] != 0)
    {
        return possible_moves;
    }
    auto tent0 = Tentative(3,2);
    if (space.space[1][tent0.x][tent0.y] != 0)
    {
        return possible_moves;
    }
    auto tent1 = Tentative(4,5);
    if (space.space[1][tent1.x][tent1.y] != 0)
    {
        return possible_moves;
    }
    
    vector<Pos2d2l> cclk;
    cclk.push_back(tent0);
    cclk.push_back(tent1);
    cclk.push_back(polysim.locs[5]);
    cclk.push_back(polysim.locs[2]);
    cclk.push_back(polysim.locs[3]);
    cclk.push_back(polysim.locs[4]);
    cclk.push_back(tent6);
    cclk.push_back(tent7);
    
    vector<Pos2d2l> clk;
    clk.push_back(cclk[6]);
    clk.push_back(cclk[7]);
    clk.push_back(cclk[4]);
    clk.push_back(cclk[5]);
    clk.push_back(cclk[2]);
    clk.push_back(cclk[3]);
    clk.push_back(cclk[0]);
    clk.push_back(cclk[1]);
    
    possible_moves.push_back(cclk);
    possible_moves.push_back(clk);
    
    return possible_moves;
    
    
}



#endif /* rotationmove_hpp */
