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
#include <set>
#include <algorithm>
#include "position.h"
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
    P GetRefPoint(int rubiscoid)
    {
        auto locs = space.Sumos[rubiscoid].locs;
        int max_x = locs[2].x, max_y = locs[2].y;
        for (int i = 3; i <= 5; i++)
        {
            if (locs[i].x == (max_x + 1) % space.Lx)
                max_x = locs[i].x;
            if (locs[i].y == (max_y + 1) % space.Ly)
                max_y = locs[i].y;
        }
        for (int i = 2; i <= 5; i++)
        {
            if (locs[i].x == max_x && locs[i].y == max_y)
            {
                return locs[i];
            }
        }
        throw std::invalid_argument("Get ref point error.");
    }
    
    P RotationPoint(P oldpt, const P& refpt, int order)
    {
        if(oldpt.x < refpt.x - space.LSim - 1) oldpt.x += space.Lx;
        if(oldpt.y < refpt.y - space.LSim - 1) oldpt.y += space.Ly;
        if(oldpt.x > refpt.x + space.LSim + 1) oldpt.x -= space.Lx;
        if(oldpt.y > refpt.y + space.LSim + 1) oldpt.y -= space.Ly;
        assert(abs(oldpt.x - refpt.x) < space.LSim + 1);
        assert(abs(oldpt.y - refpt.y) < space.LSim + 1);
        
        P newpt;
        newpt.siml = oldpt.siml; //I don't care about the refpt's layer
        if (order == 1)
        {
            newpt.x = (-oldpt.y+refpt.x+refpt.y-1+(int)(space.Lx))%space.Lx;
            newpt.y = (oldpt.x+refpt.y-refpt.x+(int)(space.Ly))%space.Ly;
        }
        else if (order == 2)
        {
            newpt.x = (-oldpt.x - 1 + 2*refpt.x + (int)(space.Lx))%space.Lx;
            newpt.y = (-oldpt.y - 1 + 2*refpt.y + (int)(space.Ly))%space.Ly;
        }
        else if (order == 3)
        {
            newpt.x = (oldpt.y + refpt.x - refpt.y + (int)(space.Lx))%space.Lx;
            newpt.y = (-oldpt.x - 1 + refpt.x + refpt.y + (int)(space.Ly)) % space.Ly;
        }
        else
            throw(std::invalid_argument("Wrong rotation order."));
        assert(newpt.x>=0 && newpt.y>=0);
        
        return newpt;
    }
    vector<vector<P>> GetPossibleMoves(const Polymer<P>& polysim, int polyid);
    bool PassSingleMoveFilter(const vector<P>& newpoints, int polyid)
    {
        for (auto pt : newpoints)
        {
            int ptid = space.GetRspacePoint(pt)[0];
            if (ptid != NOBOND && ptid != polyid)
            {
                return false;
            }
        }
        return true;
    }
    
    vector<DragMoveInfo<P>> GetPossibleMultipleMoves(const vector<int>& rubiscoIDs, const vector<int>& epycIDs, const P& refpt)
    {
        vector<DragMoveInfo<P>> possible_moves;
        
        vector<vector<int>> rubiscoInBondIDs;
        for (int rubiscoID : rubiscoIDs)
        {
            vector<int> rubisco_sites_in_a_bond;
            for (int l = 0; l < space.LSumo; l++)
            {
                if(space.InABond(space.Sumos[rubiscoID].locs[l]))
                    rubisco_sites_in_a_bond.push_back(l);
            }
            rubiscoInBondIDs.push_back(rubisco_sites_in_a_bond);
        }
        
        
        for (int i = 1; i < 3; i++)
        {
            DragMoveInfo<P> the_move;
            
            for (int epycID : epycIDs)
            {
                vector<P> epycNewPoints;
                for (auto pt : space.Sims[epycID].locs)
                {
                    epycNewPoints.push_back(RotationPoint(pt, refpt, i));
                }
                the_move.epycNewPoints.push_back(epycNewPoints);
            }
            if (!PassDragMoveLayerFilter(the_move, 'e')) continue;
            
            for (int rubiscoID : rubiscoIDs)
            {
                vector<P> rubiscoNewPoints;
                for (auto pt : space.Sumos[rubiscoID].locs)
                {
                    rubiscoNewPoints.push_back(RotationPoint(pt, refpt, i));
                }
                the_move.rubiscoNewPoints.push_back(rubiscoNewPoints);
            }
            if (!PassDragMoveLayerFilter(the_move, 'r'))  continue;
            
            the_move.rubiscoIDs = rubiscoIDs;
            the_move.epycIDs = epycIDs;
            the_move.rubiscoInBondIDs = rubiscoInBondIDs;
            possible_moves.push_back(the_move);
        }
        return possible_moves;
    }
    bool PassDragMoveLayerFilter(const DragMoveInfo<P>& dragmove, char layer)
    {
        auto NewPoints = (layer == 'e')? dragmove.epycNewPoints : dragmove.rubiscoNewPoints;
        auto IDs = (layer == 'e')? dragmove.epycIDs : dragmove.rubiscoIDs;
        for (auto new_pts : NewPoints)
        {
            for (auto newpt : new_pts)
            {
                int ptid = space.GetRspacePoint(newpt)[0];
                if (ptid != NOBOND)
                {
                    // This site is not intrisically empty, we should check if it's gonna be removed
                    bool not_in_the_list = true;
                    for (auto id : IDs) if (id == ptid) not_in_the_list = false;
                    
                    
                    //I see a point that I am trying to move to
                    //that is neither empty nor in the list to be cleared
                    //so that is not a feasible move
                    if (not_in_the_list) return false;
                }
            }
        }
        return true;
    }
    
    vector<DragMoveInfo<P>> GetPossibleDragMoves(int rubiscoID)
    {
        vector<int> rubiscoIDs;
        rubiscoIDs.push_back(rubiscoID);
        
        auto rubiscoOldPoints = space.Sumos[rubiscoID].locs;
        
        std::set<int> connected_epycs;
        for(auto pt : rubiscoOldPoints)
        {
            int epycID = space.GetRspacePoint(space.BondNeighbor(pt)[0])[0];
            if (epycID != NOBOND) connected_epycs.insert(epycID);
        }
        vector<int> epycIDs;
        std::copy(connected_epycs.begin(), connected_epycs.end(), std::back_inserter(epycIDs));
        
        P refpt = GetRefPoint(rubiscoID);
        
        return GetPossibleMultipleMoves(rubiscoIDs, epycIDs, refpt);
    }
    
};


template<>
vector<vector<Pos2d2l>> RotationMove<Space2D2L, Pos2d2l>::GetPossibleMoves(const Polymer<Pos2d2l>& polysim, int polyid)
{
    vector<vector<Pos2d2l>> possible_moves;
    
    auto refpoint = GetRefPoint(polyid);
    
    for (int i = 1; i <= 3; i++)
    {
        vector<Pos2d2l> new_points;
        for (int l=0; l<space.LSumo; l++)
        {
            new_points.push_back(RotationPoint(polysim.locs[l], refpoint, i));
        }
        if (PassSingleMoveFilter(new_points, polyid))
        {
            possible_moves.push_back(new_points);
        }
        for(int idx = 1; idx < space.LSumo; idx++)
        {
            assert(new_points[idx].x == new_points[idx-1].x || (new_points[idx].x - new_points[idx-1].x + space.Lx)%space.Lx == 1 || (new_points[idx-1].x - new_points[idx].x + space.Lx)%space.Lx == 1);
            assert(new_points[idx].y == new_points[idx-1].y || (new_points[idx].y - new_points[idx-1].y + space.Ly)%space.Ly == 1 || (new_points[idx-1].y - new_points[idx].y + space.Ly)%space.Ly == 1);
        }
    }
    
    return possible_moves;
}



#endif /* rotationmove_hpp */
