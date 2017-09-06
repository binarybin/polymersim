//
//  flipmove.hpp
//  
//
//  Created by Guanhua He on 7/18/17.
//
//

#ifndef flipmove_hpp
#define flipmove_hpp

#include <tuple>
#include <vector>
#include <cassert>
#include <deque>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <algorithm>
#include <stdexcept>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P>
class FlipMove
{
private:
    S& space;
public:
    FlipMove(S& thespace): space(thespace) {}
    void UpdatePolymer(Polymer<P>& poly, vector<int> new_wiring)
    {
        vector<P> temp_new_polymer;
        for (int i = 0; i < poly.locs.size(); i++)
            temp_new_polymer.push_back(poly.locs[new_wiring[i]]);
        poly.locs = temp_new_polymer;
    }
    
    void UpdateReverseCheckingSpace(Polymer<P>& poly)
    {
        int polyid = space.GetRspacePoint(poly.locs[0])[0];
        for (int i = 0; i < poly.locs.size(); i++)
            space.SetRspacePoint(poly.locs[i], polyid, i);
    }
    
    vector<vector<int>> GetPossibleMoves(Polymer<P>& poly)
    {
        vector<tuple<int, int, P>> possible_moves;
        vector<int> xs;
        vector<int> ys;
        
        for (int i = 0; i < poly.locs.size(); i++)
        {
            xs.push_back(poly.locs[i].x);
            ys.push_back(poly.locs[i].y);
        }
        
        int maxx = *max_element(std::begin(xs), std::end(xs));
        int maxy = *max_element(std::begin(ys), std::end(ys));
        int minx = *min_element(std::begin(xs), std::end(xs));
        int miny = *min_element(std::begin(ys), std::end(ys));
        
        vector<vector<int>> moves;
        
        vector<int> new_wiring;
        int cross = 0;

        if (maxx - minx == 1  && maxy - miny == 3 )
        {
            for (int i = 0; i < poly.locs.size()-1; i++)
            {
                
                if (poly.locs[i].y == poly.locs[i+1].y && poly.locs[i].y != maxy && poly.locs[i].y != miny)
                {

                    cross += i;

                }
            }
            
            if (cross == 0)
            {
                if (poly.locs[0].y == poly.locs[5].y)
                    cross = -4;
                else if(poly.locs[0].x == poly.locs[5].x)
                    cross = -2;
                else if (poly.locs[4].y == poly.locs[7].y)
                    cross = -6;
                else
                    cross = -1;
            }
            
        }
        else if (maxx - minx == 3 && maxy - miny == 1)
        {
            for (int i = 0; i < poly.locs.size()-1; i++)
            {
                if (poly.locs[i].x == poly.locs[i+1].x && poly.locs[i].x != maxx && poly.locs[i].x != minx)
                {

                    cross += i;
                }
            }
            if (cross == 0)
            {
                if (poly.locs[0].x == poly.locs[5].x)
                    cross = -4;
                else if(poly.locs[0].y == poly.locs[5].y)
                    cross = -2;
                else if (poly.locs[4].x == poly.locs[7].x)
                    cross = -6;
                else
                    cross = -1;
            }
        }
        
        if (cross == 0 || cross == -1)
            return moves;
        else if (abs(cross)==2)
            new_wiring = {0,1,2,7,6,5,4,3};
        else if (abs(cross)==4)
            new_wiring = {4,3,2,1,0,5,6,7};
        else if (abs(cross)==6)
            new_wiring = {2,1,0,3,4,7,6,5};
        else
        {
            cout<<"cross: "<<cross<<endl;
            throw(std::invalid_argument("Flip move unknown type."));
        }
        moves.push_back(new_wiring);
        return moves;
        
    }
};




#endif /* flipmove_hpp */
