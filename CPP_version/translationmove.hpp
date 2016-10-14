//
//  translationmove.hpp
//  polymersim
//
//  Created by Bin Xu on 4/8/16.
//  Copyright © 2016 Bin Xu. All rights reserved.
//

#ifndef translationmove_hpp
#define translationmove_hpp

#include <tuple>
#include <vector>
#include <cassert>
#include <deque>
#include <iostream>
#include <fstream>
#include <cstdlib>

using std::tuple;
using std::vector;
using std::make_tuple;

template <class S, class P>
class TranslationMove
{
private:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    S& space;
public:
    TranslationMove(S& thespace): space(thespace) {}
    vector<vector<P>> GetPossibleMoves(const Polymer<P>& polysim, int polyid);
    
    P TranslationPoint(P oldpt, int direction)
    {
        P newpt;
        newpt.siml = oldpt.siml;
        if (direction == 1)
        {
            newpt.x = (oldpt.x+1)%space.Lx;
            newpt.y = oldpt.y;
        }
        else if (direction == 2)
        {
            newpt.x = oldpt.x;
            newpt.y = (oldpt.y + 1)%space.Ly;
        }
        else if (direction == 3)
        {
            newpt.x = (oldpt.x - 1 + (int)(space.Lx))%space.Lx;
            newpt.y = oldpt.y;
        }
        else if (direction == 4)
        {
            newpt.x = oldpt.x;
            newpt.y = (oldpt.y - 1 + (int)(space.Ly))%space.Ly;
        }
        else
            throw(std::invalid_argument("Wrong translation order."));
        return newpt;
    }
    
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
        
        
        for (int i = 1; i <= 4; i++)
        {
            DragMoveInfo<P> the_move;
            the_move.rubiscoIDs = rubiscoIDs;
            the_move.epycIDs = epycIDs;
            for (int epycID : epycIDs)
            {
                vector<P> epycNewPoints;
                for (auto pt : space.Sims[epycID].locs)
                {
                    epycNewPoints.push_back(TranslationPoint(pt, i));
                }
                the_move.epycNewPoints.push_back(epycNewPoints);
            }
            if (!PassDragMoveLayerFilter(the_move, 'e')) continue;
            
            for (int rubiscoID : rubiscoIDs)
            {
                vector<P> rubiscoNewPoints;
                for (auto pt : space.Sumos[rubiscoID].locs)
                {
                    rubiscoNewPoints.push_back(TranslationPoint(pt, i));
                }
                the_move.rubiscoNewPoints.push_back(rubiscoNewPoints);
            }
            if (!PassDragMoveLayerFilter(the_move, 'r'))  continue;
            
            
            
            the_move.rubiscoInBondIDs = rubiscoInBondIDs;
            possible_moves.push_back(the_move);
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
        
        P refpt;
        
        return GetPossibleMultipleMoves(rubiscoIDs, epycIDs, refpt);
    }
    vector<DragMoveInfo<P>> GetPossibleBlobMoves(int rubiscoID)
    {
        std::map<int, int> rubiTable;
        std::map<int, int> epycTable;
        std::deque<pair<int, bool>> tasks;
        tasks.push_back(std::make_pair(rubiscoID, false));
        rubiTable[rubiscoID] = 1;
        while (!tasks.empty())
        {
            auto task = tasks.front();
            tasks.pop_front();
            const vector<P>& pts = (task.second ? space.Sims[task.first].locs : space.Sumos[task.first].locs);
            for (const auto & pt : pts)
            {
                int linked_to_id = space.GetRspacePoint(space.BondNeighbor(pt)[0])[0];
                bool newType = !task.second;
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
        vector<int> rubis;
        for (auto a_pair : rubiTable)
            if (a_pair.second == 1) rubis.push_back(a_pair.first);
        
        vector<int> epycs;
        for (auto a_pair : epycTable)
            if (a_pair.second == 1) epycs.push_back(a_pair.first);
        
        P refpt;
        
        auto the_result = GetPossibleMultipleMoves(rubis, epycs, refpt);
        
        return the_result;
    }

};
template<class S, class P>
vector<vector<P>> TranslationMove<S, P>::GetPossibleMoves(const Polymer<P>& polysim, int polyid)
{
    vector<vector<P>> possible_moves;
    
    for (int i = 1; i <= 4; i++)//Moves
    {
        vector<P> new_points;
        for (int l=0; l<space.LSumo; l++) //Points
        {
            new_points.push_back(TranslationPoint(polysim.locs[l], i));
        }
        if (PassSingleMoveFilter(new_points, polyid))
        {
            possible_moves.push_back(new_points);
        }
    }
    
    return possible_moves;
}


#endif /* translationmove_hpp */
