//
//  space2d1l.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include "space2d1l.hpp"

std::ostream& PrintSpaceOccupation(std::ostream &out, Space2D1L &space)
{
    out<<"Two dimensional single layer space"<<endl;
    out<<"Site occupation configuration: "<<endl;
    for (int i = 0; i < space.Lx; i++)
    {
        for (int j = 0; j < space.Ly; j++)
        {
            out<<std::setw(3)<<space.space[i][j];
        }
        out<<endl;
    }
    out<<"END"<<endl;
    return out;
}

void Space2D1L::DiluteInit(char direction)
{
    if (direction == 'h')
    {
        int dl = (int)Lx/(NSim + NSumo);
        for (int x = 0; x < NSim; x++)
        {
            vector<Position> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].x = x * dl;
                locs[l].y = l;
            }
            this->Place('i', x+1, locs);
        }
        
        for (int x = 0; x < NSumo; x++)
        {
            vector<Position> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].x = (int)Lx - 1 - x * dl;
                locs[l].y = l;
            }
            this->Place('u', x+1, locs);
        }
    }
    else if (direction == 'v')
    {
        int dl = (int)Ly/(NSim + NSumo);
        for (int y = 0; y < NSim; y++)
        {
            vector<Position> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].x = l;
                locs[l].y = y*dl;
            }
            this->Place('i', y+1, locs);
        }
        
        for (int y = 0; y < NSumo; y++)
        {
            vector<Position> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].x = l;
                locs[l].y = (int)Ly - 1 - y * dl;
            }
            this->Place('u', y+1, locs);
        }
    }
    else
        throw std::invalid_argument("");
}

void Space2D1L::DenseInit(int rows)
{
    if (LSim * rows > Ly || LSumo * rows > Ly)
    {
        throw std::invalid_argument("Too many rows to initialize polymers.");
    }
    
    int lines_sim = ceil((float)NSim/rows);
    int lines_sumo = ceil((float)NSumo/rows);
    
    int dx = (int) Lx / (lines_sim + lines_sumo );
    cout << "Dense initialization: dx = "<<dx<<endl;
    cout << "dx = " << Lx << " /("<<lines_sim<<" + "<<lines_sumo<<" )"<<endl;
    if(dx==0)
        throw std::invalid_argument("Just a little bit too dense. This can be settled but the current implementation does not support");
    int dy = (int) Ly/rows;
    int count_sim = 0;
    
    for (int xidx = 0; xidx < lines_sim; xidx++)
    {
        if (count_sim >= NSim)
            break;
        
        int x = xidx*dx;
        
        for (int yidx = 0; yidx < rows; yidx++)
        {
            if (count_sim >= NSim)
                break;
            int y = yidx*dy;
            count_sim += 1;
            cout<<"Initializing SIM "<<count_sim<<" row "<<x<<" column "<<y<<endl;
            

            vector<Position> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].x = x;
                locs[l].y = y+l;
            }
            
            this->Place('i', count_sim, locs);
        }
    }
    
    int count_sumo = 0;
    
    for (int xidx = 0 + lines_sim; xidx < lines_sim + lines_sumo; xidx++)
    {
        if (count_sumo >= NSumo)
            break;
        
        int x = xidx*dx;
        
        for (int yidx = 0; yidx < rows; yidx++)
        {
            if (count_sumo >= NSumo)
                break;
            int y = yidx*dy;
            count_sumo += 1;
            cout<<"Initializing SUMO "<<count_sumo<<" row "<<x<<" column "<<y<<endl;
            
            
            vector<Position> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].x = x;
                locs[l].y = y+l;
            }
            
            this->Place('u', count_sumo, locs);
        }
    }
}








