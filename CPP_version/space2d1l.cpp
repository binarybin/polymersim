//
//  space2d1l.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//
#include <sstream>
#include "space2d1l.hpp"

std::ostream& PrintSpace(std::ostream &out, Space2D1L &space)
{
    out<<"Two dimensional single layer space"<<endl;
    out<<"Site occupation configuration: "<<endl;
    for (int i = 0; i < space.Lx; i++)
    {
        for (int j = 0; j < space.Ly; j++)
        {
            out<<space.space[i][j]<<'\t';
        }
        out<<endl;
    }
    out<<"END"<<endl;
    return out;
}

std::ostream& PrintBond(std::ostream &out, Space2D1L &space)
{
    out<<"Two dimensional single layer space"<<endl;
    out<<"Bond configuration: "<<endl;
    for (int i = 0; i < space.Lx; i++)
    {
        for (int j = 0; j < space.Ly; j++)
        {
            out<<"("<<space.bond[i][j][0]<<","<<space.bond[i][j][1]<<")\t";
        }
        out<<endl;
    }
    out<<"END"<<endl;
    return out;
}

std::ostream& PrintPolymer(std::ostream &out, Space2D1L &space)
{
    out<<"Two dimensional single layer space"<<endl;
    for (int i = 0; i < space.Sims.size(); i++)
    {
        out<<"Sim Nbr "<<i<<endl;
        for (auto loc : space.Sims[i].locs)
            out<<"("<<loc.x<<","<<loc.y<<")"<<" ";
        out<<endl;
    }
    for (int i = 0; i < space.Sumos.size(); i++)
    {
        out<<"Sumo Nbr "<<i<<endl;
        for (auto loc : space.Sumos[i].locs)
            out<<"("<<loc.x<<","<<loc.y<<")"<<" ";
        out<<endl;
    }
    out<<"END"<<endl;
    return out;
}



void Space2D1L::ResumeBasicInfo(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal)
{
    int SimCount = 0;
    while (polyanal[SimCount][0].compare("Sim") == 0)
        SimCount++;

    NSim = SimCount;
    NSumo = (int)polyanal.size() - SimCount;
    LSim = (int)polyanal[0].size() - 3;
    LSumo = (int)polyanal[SimCount].size() - 3;
    Lx = spaceanal.size();
    Ly = spaceanal[0].size();
    
    space = vector<vector<int>>(Lx, vector<int>(Ly));
    bond = vector<vector<vector<int>>>(Lx, vector<vector<int>>(Ly, vector<int>(2, NOBOND)));
    rspace = vector<vector<vector<int>>>(Lx, vector<vector<int>>(Ly, vector<int>(2, NOBOND)));

    Sims.clear();
    Sumos.clear();
    
}

void Space2D1L::ResumeSpace(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal)
{
    for (int i = 0; i < Lx; i++)
    {
        for (int j = 0; j < Ly; j++)
        {
            space[i][j] = std::stoi(spaceanal[i][j]);
        }
    }
    
}

void Space2D1L::ResumePolymer(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal)
{
    for (int id = 0; id < NSim; id++)
    {
        Polymer<Pos2d1l> poly;
        poly.poly_id = id;
        for (int pos = 0; pos < polyanal[pos].size(); pos++)
        {
            if (pos>2)
            {
                string ss(polyanal[id][pos]);
                std::replace(ss.begin(), ss.end(), '(', ' ');
                std::replace(ss.begin(), ss.end(), ')', ' ');
                vector<string> tokens;
                std::istringstream split(ss);
                for (std::string each; std::getline(split, each, ','); tokens.push_back(each));
                int x = std::stoi(tokens[0]);
                int y = std::stoi(tokens[1]);
                poly.locs.push_back(Pos2d1l(x, y));
            }
        }
        Sims.push_back(poly);
    }
    
    for (int id = NSim; id < NSim+NSumo; id++)
    {
        Polymer<Pos2d1l> poly;
        poly.poly_id = id;
        for (int pos = 0; pos < polyanal[id].size(); pos++)
        {
            if (pos>2)
            {
                string ss(polyanal[id][pos]);
                std::replace(ss.begin(), ss.end(), '(', ' ');
                std::replace(ss.begin(), ss.end(), ')', ' ');
                vector<string> tokens;
                std::istringstream split(ss);
                for (std::string each; std::getline(split, each, ','); tokens.push_back(each));
                int x = std::stoi(tokens[0]);
                int y = std::stoi(tokens[1]);
                poly.locs.push_back(Pos2d1l(x, y));
            }
        }
        Sumos.push_back(poly);
    }
}

void Space2D1L::ResumeBond(vector<vector<string>>& bondanal)
{
    for (int i = 0; i < Lx; i++)
    {
        for (int j = 0; j < Ly; j++)
        {
            string ss(bondanal[i][j]);
            std::replace(ss.begin(), ss.end(), '(', ' ');
            std::replace(ss.begin(), ss.end(), ')', ' ');
            vector<string> tokens;
            std::istringstream split(ss);
            for (std::string each; std::getline(split, each, ','); tokens.push_back(each));
            int x = std::stoi(tokens[0]);
            int y = std::stoi(tokens[1]);
            bond[i][j][0] = x;
            bond[i][j][1] = y;
        }
    }
}

void Space2D1L::ResumeReverse()
{
    for (int id = 0; id < NSim; id++)
    {
        for (int pos = 0; pos < LSim; pos++)
        {
            int x = Sims[id].locs[pos].x;
            int y = Sims[id].locs[pos].y;
            rspace[x][y][0] = id;
            rspace[x][y][1] = pos;
        }
    }
    
    for (int id = 0; id < NSumo; id++)
    {
        for (int pos = 0; pos < LSumo; pos++)
        {
            int x = Sumos[id].locs[pos].x;
            int y = Sumos[id].locs[pos].y;
            rspace[x][y][0] = id;
            rspace[x][y][1] = pos;
        }
    }
    cout<<"Finished"<<endl;
}

void Space2D1L::DiluteInit(char direction)
{
    if (direction == 'h')
    {
        int dl = (int)Lx/(NSim + NSumo);
        for (int x = 0; x < NSim; x++)
        {
            vector<Pos2d1l> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].x = x * dl;
                locs[l].y = l;
            }
            this->Place('i', x, locs);
        }
        
        for (int x = 0; x < NSumo; x++)
        {
            vector<Pos2d1l> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].x = (int)Lx - 1 - x * dl;
                locs[l].y = l;
            }
            this->Place('u', x, locs);
        }
    }
    else if (direction == 'v')
    {
        int dl = (int)Ly/(NSim + NSumo);
        for (int y = 0; y < NSim; y++)
        {
            vector<Pos2d1l> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].x = l;
                locs[l].y = y*dl;
            }
            this->Place('i', y, locs);
        }
        
        for (int y = 0; y < NSumo; y++)
        {
            vector<Pos2d1l> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].x = l;
                locs[l].y = (int)Ly - 1 - y * dl;
            }
            this->Place('u', y, locs);
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
            
            cout<<"Initializing SIM "<<count_sim<<" row "<<x<<" column "<<y<<endl;
            

            vector<Pos2d1l> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].x = x;
                locs[l].y = y+l;
            }
            
            this->Place('i', count_sim, locs);
            count_sim += 1;
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
            
            cout<<"Initializing SUMO "<<count_sumo<<" row "<<x<<" column "<<y<<endl;
            
            
            vector<Pos2d1l> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].x = x;
                locs[l].y = y+l;
            }
            
            this->Place('u', count_sumo, locs);
            count_sumo += 1;
        }
    }
}

const string Space2D1L::name = "2d1l";
