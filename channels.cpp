//
//  channels.cpp
//  Project1_555
//
//  Created by Bing Hong Fu on 2/2/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#include "channels.h"
#include "utils.h"
#include <cmath>

using namespace std;

void AWGN(vector<vector<double>> &input, double N_0){
    // input should be a matrix
    
/*    vector<vector<double>> temp;


    
    auto iter12 = input.begin();
    auto iter22 = iter12->begin();
    
    for(auto iter11 = temp.begin(), iter12 = input.begin(); iter11 < temp.end(); iter11 ++, iter12 ++){
        for(auto iter21 = iter11->begin(), iter22 = iter12->begin(); iter21 < iter11->end(); iter21++, iter22 ++){
            *iter21 = *iter22 + sqrt(N_0/2) * gauss1(ptr);
        }
    }
*/
    int seed = -123;
    int *ptr = &seed;
    
    for(auto iter = input.begin(); iter < input.end(); iter++){
        for(auto iter1 = iter->begin(); iter1 < iter->end(); iter1 ++){
            *iter1 += (sqrt(N_0/2) * gauss1(ptr));
        }
    }

}