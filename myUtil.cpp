//
//  myUtil.cpp
//  Project1_555
//
//  Created by Bing Hong Fu on 2/3/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#include "myUtil.h"
#include "base.h"
#include <iostream>

using namespace std;


bool de2bi(int num, vector<int> &de, vector<int> &bi){
// num: number of bit, de: decimal vector, bi: binary vector
    for(auto iter = de.begin(); iter < de.end(); iter++){
        vector<int> temp(num, 0);
        if(!gr::trellis::dec2base(*iter, 2, temp))
            return false;
        bi.insert(bi.end(), temp.begin(), temp.end());
    }
    return true;
}

bool bi2de(int num, vector<int> &de, vector<int> &bi){
    if(! bi.size() % ((unsigned int) num)){
        std::cout << "binary bits not fit" << endl;
        return false;
    }
    for(auto iter = bi.begin(); iter < bi.end(); iter += num){
        vector<int> temp(iter, iter+num);
        de.push_back(gr::trellis::base2dec(temp, 2));
    }
    return true;
}

double mean(vector<double> input){
    double sum = 0;
    for(auto i = input.begin(); i < input.end(); i++)
        sum += *i;
    return sum/(double)input.size();
}

double mean(vector<int> input){
    int sum = 0;
    for(auto i = input.begin(); i < input.end(); i++)
        sum += *i;
    return sum/(double)input.size();
}


double variance(vector<double> input){
    auto s = input.size();
    double square_sum = 0;
    for(auto i = input.begin(); i < input.end(); i++){
        square_sum += (*i * *i);
    }
    double m = mean(input);
    return square_sum/s - m * m;
}

double variance(vector<int> input){
    auto s = input.size();
    int square_sum = 0;
    for(auto i = input.begin(); i < input.end(); i++){
        square_sum += (*i * *i);
    }
    double m = mean(input);
    return square_sum/s - m * m;
}

void printMatrix(vector<vector<double>> input){
    for(auto iter = input.begin(); iter < input.end(); iter++){
        for(auto iter1 = iter->begin(); iter1 < iter->end(); iter1++)
            cout << *iter1 << "\t";
        cout << endl;
    }
}

double dist2(vector<double> vec1, vector<double> vec2){
    if(vec1.size() != vec2.size()){
        cerr << "Distance not available: Vector length don't match!" << endl;
        exit(0);
    }
    
    double dist = 0;
    vector<int>::iterator iter2;
    for(auto iter1 = vec1.begin(), iter2 = vec2.begin(); iter1 < vec1.end(); iter1++, iter2++)
        dist += (*iter1 - *iter2) * (*iter1 - *iter2);
    return dist;
}

void printMatrix(vector<vector<int>> input){
    for(auto iter = input.begin(); iter < input.end(); iter++){
        for(auto iter1 = iter->begin(); iter1 < iter->end(); iter1++)
            cout << *iter1 << "\t";
        cout << endl;
    }
}

unsigned long nearestNeighbor(vector<vector<double>>map, vector<double> input){
    if(map.size() == 0){
        cerr << "Neibour Map Null" << endl;
        exit(0);
    }
    
    auto min_iter = map.begin();
    double min_dist = dist2(*min_iter, input);
    for(auto iter = map.begin(); iter < map.end(); iter++){
        if(dist2(*iter, input) < min_dist){
            min_iter = iter;
            min_dist = dist2(*min_iter, input);
        }
    }
    return min_iter - map.begin();
}

int hammingDist(vector<int> vec1, vector<int> vec2){
    if(vec1.size() != vec2.size()){
        cerr << "Distance not available: Vector length don't match!" << endl;
        exit(0);
    }
    int count = 0;
    vector<int>::iterator iter2;
    for(auto iter1 = vec1.begin(), iter2 = vec2.begin(); iter1 < vec1.end(); iter1++, iter2++){
        if(*iter1 != *iter2)
            count  ++;
    }
    return count;
}

int minIndex(vector<double> vec){
    double min = *vec.begin();
    int index = 0;
    for(auto iter = vec.begin(); iter < vec.end(); iter ++){
        if(*iter < min){
            min = *iter;
            index = int(iter - vec.begin());
        }
    }
    return index;
}

vector<unsigned long> nearestNeighbors(vector<vector<double>> map, vector<vector<double>> input){
    vector<unsigned long> temp;
    for(auto iter = input.begin(); iter < input.end(); iter++){
        temp.push_back(nearestNeighbor(map, *iter));
    }
    return temp;
}

