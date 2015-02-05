//
//  main.cpp
//  Project1_555
//
//  Created by Bing Hong Fu on 1/31/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#include <iostream>

#include "fsm.h"
#include "base.h"
#include "myUtil.h"
#include "utils.h"
#include "channels.h"
#include <cmath>

using namespace std;
using namespace gr::trellis;

enum module{
    PSK,
    PAM,
    ORTH
};

vector<int> convEncoder(fsm f, vector<int> infoBits);

vector<int> convDecoder(fsm f, vector<int> input);

vector<int> grayEncoder(vector<int> &input, int M);

vector<int> grayDecoder(vector<int> &input, int M);

vector<vector<double>> modulation(module mod, vector<int> bits, int M);

vector<int> minimumDistanceDetect(vector<vector<double>> map, vector<vector<double>> receivedSignal);

int additionalDist(int from, int to, fsm f, vector<int> observation);

int findInput(int i, int s, fsm f);

int main(int argc, const char * argv[]) {
    
    // Get input of fsm from text file
    
    const char * filename = "input.txt";
    fsm f(filename);
    vector<int> info({0, 0, 0, 1, 1, 0, 1, 1});
    
    /*
    
    decoded = convEncoder(f, info);
    decoded = grayEncoder(info, 16);
    decoded = grayDecoder(decoded, 16);
     
    for(auto iter = decoded.begin(); iter < decoded.end(); iter ++)
        cout << *iter << " ";

    */
    
    
    /*---------------modulation test------------------
     
    vector<vector<double>> moded = modulation(module::ORTH, info, 4);
    
    //printMatrix(moded);
    
    AWGN(moded, 1);
    //printMatrix(moded);
    
    //--------------------------------------------------*/
    
    /*---------------guassian and mean/var test---------
    int seed = -123;
    int *ptr = &seed;
    
    vector<double> ran;
    for(int i = 1; i < 10000; i++)
        ran.push_back(gauss1(ptr));
    
    cout << mean(ran) << "   " << variance(ran) << endl;
     
    //-------------------------------------------------*/
    
    
    return 0;
}

vector<int> convEncoder(fsm f, vector<int> infoBits){
    vector<int>encoded;
    int state = 0;
    int n = log2(f.I());
    int k = log2(f.O());
    const vector<int> NS = f.NS();
    const vector<int> OS = f.OS();
    vector<int> decEncode;
    
    for(auto iter = infoBits.begin(); iter < infoBits.end(); iter += n){
        vector<int> temp(iter, iter+n);
        int input = base2dec(temp, 2);
        int output = OS[state*f.I() + input];
        
        state = NS[state*f.I() + input];
        decEncode.push_back(output);
    }
    
    de2bi(k, decEncode, encoded);
    return encoded;
}

vector<int> grayEncoder(vector<int> &input, int M){
    vector<int> grayEncoded;
    if(4 == M){
        vector<vector<int>> grayMap({{0, 0}, {0, 1}, {1, 1}, {1, 0}});
        for(auto iter = input.begin(); iter < input.end(); iter += 2){
            vector<int> temp(iter, iter+2);
            int decimal = base2dec(temp, 2);
            grayEncoded.insert(grayEncoded.end(), grayMap[decimal].begin(), grayMap[decimal].end());
        }
    }
    
    else if(8 == M){
        vector<vector<int>> grayMap({{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 0}, {1, 1, 1}, {1, 0, 1}, {1, 0, 0}});
        for(auto iter = input.begin(); iter < input.end(); iter += 3){
            vector<int> temp(iter, iter+3);
            int decimal = base2dec(temp, 2);
            grayEncoded.insert(grayEncoded.end(), grayMap[decimal].begin(), grayMap[decimal].end());
        }
    }
    
    else if(16 == M){
        vector<vector<int>> grayMap({{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 0, 1, 0}, {0, 1, 1, 0}, {0, 1, 1, 1}, {0, 1, 0, 1}, {0, 1, 0, 0}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 1}, {1, 1, 1, 0}, {1, 0, 1, 0}, {1, 0, 1, 1}, {1, 0, 0, 1}, {1, 0, 0, 0}});
        for(auto iter = input.begin(); iter < input.end(); iter += 4){
            vector<int> temp(iter, iter+4);
            int decimal = base2dec(temp, 2);
            grayEncoded.insert(grayEncoded.end(), grayMap[decimal].begin(), grayMap[decimal].end());
        }
    }

    return grayEncoded;
}

vector<int> grayDecoder(vector<int> &input, int M){
    vector<int> grayDecoded;
    if(4 == M){
        vector<vector<int>> grayMap({{0, 0}, {0, 1}, {1, 1}, {1, 0}});
        for(auto iter = input.begin(); iter < input.end(); iter += 2){
            vector<int> temp(iter, iter+2);
            int decimal = base2dec(temp, 2);
            grayDecoded.insert(grayDecoded.end(), grayMap[decimal].begin(), grayMap[decimal].end());
        }
    }
    
    else if(8 == M){
        vector<vector<int>> grayMap({{0, 0, 0}, {0, 0, 1}, {0, 1, 1}, {0, 1, 0}, {1, 1, 1}, {1, 1, 0}, {1, 0,  0}, {1, 0, 1}});
        for(auto iter = input.begin(); iter < input.end(); iter += 3){
            vector<int> temp(iter, iter+3);
            int decimal = base2dec(temp, 2);
            grayDecoded.insert(grayDecoded.end(), grayMap[decimal].begin(), grayMap[decimal].end());
        }
    }
    
    else if(16 == M){
        vector<vector<int>> grayMap({{0, 0, 0, 0}, {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 0, 1, 0}, {0, 1, 1, 1}, {0, 1, 1, 0}, {0, 1, 0, 0}, {0, 1, 0, 1}, {1, 1, 1, 1}, {1, 1, 1, 0}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 0, 0, 0}, {1, 0, 0, 1}, {1, 0, 1, 1}, {1, 0, 1, 0}});
        for(auto iter = input.begin(); iter < input.end(); iter += 4){
            vector<int> temp(iter, iter+4);
            int decimal = base2dec(temp, 2);
            grayDecoded.insert(grayDecoded.end(), grayMap[decimal].begin(), grayMap[decimal].end());
        }
    }
    return grayDecoded;
}

vector<vector<double>> modulation(module mod, vector<int> bits, int M){
    //  Modulation constellation of PAM, PSK, Orthogonal
    //  E_b = 1

    
    vector<vector<double>> symbol;
    int n = log2(M);
    if(mod == PAM){

    }
    
    else if(mod == PSK){
        for(auto iter = bits.begin(); iter < bits.end(); iter += n){
            vector<int> temp(iter, iter + n);
            int decimal = base2dec(temp, 2);
            vector<double> currentSymbol({log2(M) * cos(decimal * 2 * M_PI / M), log2(M) * sin(decimal * 2 * M_PI / M)});
            symbol.push_back(currentSymbol);
        }
    }
    
    else if(mod == ORTH){
        vector<vector<double>> map;
        for(int i = 0; i < M; i++){
            vector<double> temp(M, 0);
            temp[i] = log2(M);
            map.push_back(temp);
        }
        
        for(auto iter = bits.begin(); iter < bits.end(); iter += n){
            vector<int> temp(iter, iter + n);
            int decimal = base2dec(temp, 2);
            symbol.push_back(map[decimal]);
        }
    }
    
    return symbol;
}

vector<int> minimumDistanceDetect(vector<vector<double>> map, vector<vector<double>> receivedSignal){
    vector<int>temp;
    vector<unsigned long> index = nearestNeighbors(map, receivedSignal);
    for(auto iter = index.begin(); iter < index.end(); iter++)
        temp.push_back((int) *iter);
    return temp;
}




vector<int> convDecoder(fsm f, vector<int> input){
    const int s = f.S();
    const int k = log2(f.O());
    const int n = log2(f.I());
    const vector<int> NS = f.NS();
    const vector<int> OS = f.OS();
    int si;
    
    vector<double> td(s, 999);
    td[0] = 0;
    vector<vector<int>> surviver(s, {});
    vector<vector<int>> surviver_next(s, {});
    vector<int> x_star(s, 0);
    vector<int> s_star(s, 0);
    vector<double> td_next(s, 999);
    
    for(auto iter = input.begin(); iter < input.end(); iter += k){
        vector<int> t(iter, iter+k);
        //int input = base2dec(t, 2);
        
        for(auto siter = surviver.begin(); siter < surviver.end(); siter ++){
            si = (int)(siter - surviver.begin());
            
            td_next[si] = td[0] + additionalDist(0, si, f, t);
            s_star[si] = 0;
            for(int i = 0; i < s; i++){
                if(td[i] + additionalDist(i, si, f, t) < td_next[si]){
                    td_next[si] = td[i] + additionalDist(i, si, f, t);
                    s_star[si] = i;
                }
            }
            x_star[si] = findInput(s_star[si], s, f);
            surviver_next[si] = surviver[si];
            surviver_next[si].push_back(x_star[si]);
        }
        surviver = surviver_next;
        td = td_next;
    }
    
    int index = minIndex(td);
    vector<int>decodedSymbol = surviver[index];
    vector<int> temp;
    de2bi(n, decodedSymbol, temp);

    return temp;
}

int findInput(int from, int to, fsm f){
    const vector<int> NS = f.NS();
    for(int i = 0; i < f.I(); i++){
        if(NS[from * f.I() + i] == to)
            return i;
    }
    
    cerr << " cannot find link" << endl;
    return -1;
}



int additionalDist(int from, int to, fsm f, vector<int> observation){
    const int I = f.I();
    const int O = f.O();
    const vector<int> OS = f.OS();
    const vector<int> NS = f.NS();
    
    if(from > f.S()-1 || to > f.S()-1){
        cerr << "State out of range" << endl;
        return -1;
    }
    
    int n = log2(O);
    
    int symbol = -1;
    
    for(int i = from * I; i < (from+1) * I; i++){
        if(NS[i] == to){
            symbol = OS[i];
        }
    }
    
    if(symbol == -1)
        return 999;
    
    vector<int> output(n, 0);
    dec2base(symbol, 2, output);
    return hammingDist(output, observation);
}



