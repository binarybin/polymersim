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
#include <deque>
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
    
    bool TouchedOtherPolymers(const vector<int>& rubiscoIDs, const vector<int>& epycIDs)
    {
        for (auto rubiID : rubiscoIDs)
        {
            for (auto rubiPT : space.Sumos[rubiID].locs)
            {
                int correspondingID = space.GetRspacePoint(space.BondNeighbor(rubiPT)[0])[0];
                if (correspondingID == NOBOND) continue;
                bool danger = true;
                for (auto epycID : epycIDs)
                    if (epycID == correspondingID) danger = false;
                if (danger)
                    return true;
            }
        }
        
        for (auto epycID : epycIDs)
        {
            for (auto epycPT : space.Sims[epycID].locs)
            {
                int correspondingID = space.GetRspacePoint(space.BondNeighbor(epycPT)[0])[0];
                if (correspondingID == NOBOND) continue;
                bool danger = true;
                for (auto rubiID : rubiscoIDs)
                    if (rubiID == correspondingID) danger = false;
                if (danger)
                    return true;
            }
        }
        return false;
    }
    
    vector<DragMoveInfo<P>> GetPossibleMultipleMoves(const vector<int>& rubiscoIDs, const vector<int>& epycIDs, const P& refpt)
    {
        vector<DragMoveInfo<P>> possible_moves;
        
        if (TouchedOtherPolymers(rubiscoIDs, epycIDs))
            return possible_moves;
        
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
        
        
        for (int i = 1; i <= 3; i++)
        {
            DragMoveInfo<P> the_move;
            the_move.rubiscoIDs = rubiscoIDs;
            the_move.epycIDs = epycIDs;
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
            
            
            
            the_move.rubiscoInBondIDs = rubiscoInBondIDs;
            possible_moves.push_back(the_move);
#ifndef NDEBUG
    /*        if (rand() % 100 == 0)
            {
                cout<<"Rotational blob move type: "<<i<<endl;
                cout<<"Rubiscos"<<endl;
                for (auto idx : rubiscoIDs)
                {
                    cout<<idx<<'\t';
                }
                cout<<endl;
                cout<<"EPYCs"<<endl;
                for (auto idx : epycIDs)
                {
                    cout<<idx<<'\t';
                }
                cout<<endl;
            }
*/
#endif
        }
        return possible_moves;
    }
    bool PassDragMoveLayerFilter(const DragMoveInfo<P>& dragmove, char layer)
    {
        auto NewPoints = (layer == 'e')? dragmove.epycNewPoints : dragmove.rubiscoNewPoints;
        auto IDs = (layer == 'e')? dragmove.epycIDs : dragmove.rubiscoIDs;
        auto CoIDs = (layer != 'e')? dragmove.epycIDs : dragmove.rubiscoIDs;
        
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
                int coptid = space.GetRspacePoint(space.BondNeighbor(newpt)[0])[0];
                if (coptid != NOBOND)
                {
                    // The corresponding site on the other layer is not intrinsically empty, we should check if it's gonna be removed
                    bool not_in_the_list = true;
                    for (auto coid : CoIDs) if (coid == coptid) not_in_the_list = false;
                    // I see a point that I am trying to move to of which the other layer point
                    //is neither empty nor in the list to be cleared so that it is not a feasible move
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
    
    vector<DragMoveInfo<P>> GetPossibleBlobMoves(int rubiscoID)
    {
        std::map<int, int> rubiTable;
        std::map<int, int> epycTable;
        std::deque<pair<int, bool>> tasks;
        tasks.push_back(std::make_pair(rubiscoID, false));
        rubiTable[rubiscoID] = 1;
        // Get the lists of rubiscos and epycs involved
        while (!tasks.empty())
        {
            auto task = tasks.front();
            tasks.pop_front();
            const vector<P>& pts = (task.second ? space.Sims[task.first].locs : space.Sumos[task.first].locs);
            for (const auto & pt : pts)
            {
                bool newType = !task.second;
                int linked_to_id = space.GetRspacePoint(space.BondNeighbor(pt)[0])[0];
                if (linked_to_id != NOBOND)
                {
                    auto& Table = (newType ? epycTable : rubiTable);
                    if (Table[linked_to_id] == 0)
                    {
                        Table[linked_to_id] = 1;
                        tasks.push_back(std::make_pair(linked_to_id, newType));
                    }
                }
            }
        }
        
        // Copy the involved rubiscos and epycs to their lists
        vector<int> rubis;
        for (auto a_pair : rubiTable)
            if (a_pair.second == 1) rubis.push_back(a_pair.first);
        
        vector<int> epycs;
        for (auto a_pair : epycTable)
            if (a_pair.second == 1) epycs.push_back(a_pair.first);
        
        vector<DragMoveInfo<P>> the_result;
        
        
        // Check if there is a cycle
        // A cycle will make rotation ambiguous
        // So we don't like it
        std::map<int, int> x_occ, y_occ;
        for (auto rubiid : rubis)
        {
            for (auto pt : space.Sumos[rubiid].locs)
            {
                x_occ[pt.x] = 1;
                y_occ[pt.y] = 1;
            }
        }
        
        for (auto epycid : epycs)
        {
            for (auto pt : space.Sims[epycid].locs)
            {
                x_occ[pt.x] = 1;
                y_occ[pt.y] = 1;
            }
        }
        bool x_cycle = true, y_cycle = true;
        for (int x = 0; x<space.Lx; x++)
        {
            if (x_occ[x] != 1) x_cycle = false;
        }
        for (int y = 0; y<space.Ly; y++)
        {
            if (y_occ[y] != 1) y_cycle = false;
        }
        if (x_cycle || y_cycle)
        {
            return the_result;
        }

        // Now we just push them into the GetPossibleMultipleMoves function
        // to get rotations around different points
        for (auto rubiid : rubis)
        {
            P refpt = GetRefPoint(rubiid);
            auto moves_one_point = GetPossibleMultipleMoves(rubis, epycs, refpt);
            the_result.insert( the_result.end(), moves_one_point.begin(), moves_one_point.end() );
        }
        return the_result;
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
