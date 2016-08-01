//
//  main.cpp
//  polymersim
//
//  Created by Bin Xu on 11/18/15.
//  Copyright © 2015 Bin Xu. All rights reserved.
//

#include <map>
#include "simanneal.hpp"
using std::stoi;

vector<tuple<int, double, double, int>> GetTaskList(std::string filename)
{
    vector<tuple<int, double, double, int>> result;
    
    ifstream in(filename);
    if(in.fail())
    {
        throw runtime_error(std::string("Failed to open file ") + filename);
    }
    else
    {
        cout<<"File open successful: "<<filename<<endl;
    }
    string tempstring;
    vector<string> spacebuf;
    while (getline(in, tempstring)) spacebuf.push_back(tempstring);
    
    char delim = '\t';
    for (auto& ss : spacebuf)
    {
        istringstream split(ss);
        vector<string> tokens;
        for (std::string each; std::getline(split, each, delim); tokens.push_back(each));
        int idx = stoi(tokens[0]);
        double beta = stof(tokens[1]);
        double gamma = stof(tokens[2]);
        int runs = stoi(tokens[3]);
        result.push_back(make_tuple(idx, beta, gamma, runs));
    }

    
    in.close();
    
    return result;
}

int main(int argc, const char * argv[])
{
    cout<<"Polymer Simulation Started!"<<endl;
    
    cout<<"argc: "<<argc<<endl;
    int nsim = stoi(argv[1]);
    int nsumo = stoi(argv[2]);
    int lsim = stoi(argv[3]);
    int lsumo = stoi(argv[4]);
    size_t lx = stoi(argv[5]);
    size_t ly = stoi(argv[6]);
    
    vector<tuple<int, double, double, int>> tasklist = GetTaskList(argv[7]);
    string signature = argv[8];
    
    SimAnnealing<Space2D2L, Pos2d2l> process(nsim, nsumo, lsim, lsumo, lx, ly, signature, tasklist);
    
    process.Run();
    return 0;
}