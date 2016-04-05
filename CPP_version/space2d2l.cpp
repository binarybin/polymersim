//
//  space2d2l.cpp
//  polymersim
//
//  Created by Bin Xu on 11/30/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//
#include <sstream>
#include "space2d2l.hpp"

std::ostream& PrintSpace(std::ostream &out, Space2D2L &space)
{
    out<<"Two dimensional double layer space"<<endl;
    out<<"Site occupation configuration: "<<endl;
    for (int i = 0; i < space.Lx; i++)
    {
        for (int j = 0; j < space.Ly; j++)
        {
            out<<space.space[0][i][j]<<","<<space.space[1][i][j]<<'\t';
        }
        out<<endl;
    }
    out<<"END"<<endl;
    return out;
}

std::ostream& PrintRSpace(std::ostream &out, Space2D2L &space)
{
    out<<"Two dimensional double layer space"<<endl;
    out<<"Reverse space configuration: "<<endl;
    for (int i = 0; i < space.Lx; i++)
    {
        for (int j = 0; j < space.Ly; j++)
        {
            if (space.rspace[0][i][j][0] < 0)
                out<<"*,";
            else
                out<<space.rspace[0][i][j][0]<<",";
            if (space.rspace[0][i][j][1] < 0)
                out<<"*--";
            else
                out<<space.rspace[0][i][j][1]<<"--";
            if (space.rspace[1][i][j][0] < 0)
                out<<"*,";
            else
                out<<space.rspace[1][i][j][0]<<",";
            if (space.rspace[1][i][j][1] < 0)
                out<<"*\t";
            else
                out<<space.rspace[1][i][j][1]<<"\t";
        }
        out<<endl;
    }
    out<<"END"<<endl;
    return out;
}

std::ostream& PrintPolymer(std::ostream &out, Space2D2L &space)
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


std::ostream& PrintBond(std::ostream &out, Space2D2L &space)
{
    out<<"No bond info for this geometry"<<endl;
    return out;
}

void Space2D2L::DiluteInit(char direction)
{
    if (direction == 'h')
    {
        int dl = (int)Lx/max(NSim, NSumo);
        for (int x = 0; x < NSim; x++)
        {
            vector<Pos2d2l> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].siml = true;
                locs[l].x = x * dl;
                locs[l].y = l;
            }
            this->Place('i', x, locs);
        }
        
        for (int x = 0; x < NSumo; x++)
        {
            vector<Pos2d2l> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].siml = false;
                locs[l].x = x * dl;
                locs[l].y = l;
            }
            this->Place('u', x, locs);
        }
    }
    else if (direction == 'v')
    {
        int dl = (int)Ly/max(NSim, NSumo);
        for (int y = 0; y < NSim; y++)
        {
            vector<Pos2d2l> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].siml = true;
                locs[l].x = l;
                locs[l].y = y*dl;
            }
            this->Place('i', y, locs);
        }
        
        for (int y = 0; y < NSumo; y++)
        {
            vector<Pos2d2l> locs(LSumo);
            for (int l = 0; l < LSumo; l++)
            {
                locs[l].siml = false;
                locs[l].x = l;
                locs[l].y = y * dl;
            }
            this->Place('u', y, locs);
        }
    }
    else
        throw std::invalid_argument("");
}

void Space2D2L::DenseInit(int rows, char typ)
{
    int Length = (typ == 'i' ? LSim : LSumo);
    
    int lines = (typ == 'i' ? ceil((float)NSim/rows) : ceil((float)NSumo/rows));
    int dx = (int) Lx / lines;
    cout << "Dense initialization: dx = "<<dx<<endl;
    cout << "dx = " << Lx << " /"<<lines<<endl;
    if(dx==0)
        throw std::invalid_argument("Just a little bit too dense. This can be settled but the current implementation does not support");
    
    int dy = (int) Ly/rows;
    if (dy < Length)
        throw std::invalid_argument("Too many rows to initialize polymers.");
    
    int count = 0;
    int N = (typ == 'i' ? NSim : NSumo);
    
    for (int xidx = 0; xidx < lines; xidx++)
    {
        if (count >= N)
            break;
        
        int x = xidx*dx;
        
        for (int yidx = 0; yidx < rows; yidx++)
        {
            if (count >= N)
                break;
            int y = yidx*dy;
            
            cout<<"Initializing "<< (typ == 'i' ? "SIM " : "SUMO ")<<count<<" row "<<x<<" column "<<y<<endl;
            
            
            vector<Pos2d2l> locs(LSim);
            for (int l = 0; l < LSim; l++)
            {
                locs[l].siml = (typ=='i');
                locs[l].x = x;
                locs[l].y = y+l;
            }
            
            Place(typ, count, locs);
            count += 1;
        }
    }
}

void Space2D2L::ResumeBasicInfo(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal)
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
    
    space = vector<vector<vector<int>>>(2, vector<vector<int>>(Lx, vector<int>(Ly)));
    rspace = vector<vector<vector<vector<int>>>>(2, vector<vector<vector<int>>>(Lx, vector<vector<int>>(Ly, vector<int>(2, NOBOND))));
    
    Sims.clear();
    Sumos.clear();
    
}

void Space2D2L::ResumeSpace(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal)
{
    for (int i = 0; i < Lx; i++)
    {
        for (int j = 0; j < Ly; j++)
        {
            string ss(spaceanal[i][j]);
            vector<string> tokens;
            std::istringstream split(ss);
            for (std::string each; std::getline(split, each, ','); tokens.push_back(each));
            int simspace = std::stoi(tokens[0]);
            int sumospace = std::stoi(tokens[1]);
            space[0][i][j] = simspace;
            space[1][i][j] = sumospace;
        }
    }
}

void Space2D2L::ResumePolymer(vector<vector<string>>& spaceanal, vector<vector<string>>& polyanal)
{
    for (int id = 0; id < NSim; id++)
    {
        Polymer<Pos2d2l> poly;
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
                poly.locs.push_back(Pos2d2l(x, y, true));
            }
        }
        Sims.push_back(poly);
    }
    
    for (int id = NSim; id < NSim+NSumo; id++)
    {
        Polymer<Pos2d2l> poly;
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
                poly.locs.push_back(Pos2d2l(x, y, false));
            }
        }
        Sumos.push_back(poly);
    }
}

void Space2D2L::ResumeBond(vector<vector<string>>& bondanal)
{
}

void Space2D2L::ResumeReverse()
{
    for (int id = 0; id < NSim; id++)
    {
        for (int pos = 0; pos < LSim; pos++)
        {
            int x = Sims[id].locs[pos].x;
            int y = Sims[id].locs[pos].y;
            rspace[0][x][y][0] = id;
            rspace[0][x][y][1] = pos;
        }
    }
    
    for (int id = 0; id < NSumo; id++)
    {
        for (int pos = 0; pos < LSumo; pos++)
        {
            int x = Sumos[id].locs[pos].x;
            int y = Sumos[id].locs[pos].y;
            rspace[1][x][y][0] = id;
            rspace[1][x][y][1] = pos;
        }
    }
}

const string Space2D2L::name = "2d2l";
