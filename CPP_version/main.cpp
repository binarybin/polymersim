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

vector<tuple<int, double, double, int, int, int>> GetTaskList(std::string filename)
{
    vector<tuple<int, double, double, int, int, int>> result;
    
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
        int runsim = stoi(tokens[3]);
        int runsumo = stoi(tokens[4]);
        int phos = stoi(tokens[5]);
        result.push_back(make_tuple(idx, beta, gamma, runsim, runsumo, phos));
    }

    in.close();
    
    return result;
}

int main(int argc, const char * argv[])
{
    cout<<"Polymer Simulation Started!"<<endl;
    
    cout<<"argc: "<<argc<<endl;
    int nsim1 = stoi(argv[1]);
    int nsim2 = stoi(argv[2]);
    int nsumo = stoi(argv[3]);
    int lsim1 = stoi(argv[4]);
    int lsim2 = stoi(argv[5]);
    int lsumo = stoi(argv[6]);
    size_t lx = stoi(argv[7]);
    size_t ly = stoi(argv[8]);
    
    vector<tuple<int, double, double, int, int, int>> tasklist = GetTaskList(argv[9]);
    string signature = argv[10];
    string initfile = argv[11];
    SimAnnealing<Space2D2L, Pos2d2l> process(nsim1, nsim2, nsumo, lsim1, lsim2, lsumo, lx, ly, signature, tasklist, initfile);
    
    process.Run();
    return 0;
}
