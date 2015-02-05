//
//  myUtil.h
//  Project1_555
//
//  Created by Bing Hong Fu on 2/3/15.
//  Copyright (c) 2015 Bing Hong Fu. All rights reserved.
//

#ifndef __Project1_555__myUtil__
#define __Project1_555__myUtil__

#include <stdio.h>
#include <vector>
bool de2bi(int num, std::vector<int> &de, std::vector<int> &bi);
//  Decimal vector to binary vector. Each decimal number is converted to num bits

bool bi2de(int num, std::vector<int> &de, std::vector<int> &bi);
//  Binary vector to decimal vector. Each num bits are converted to one number

double variance(std::vector<double> input);
double variance(std::vector<int> input);

double mean(std::vector<double> input);
double mean(std::vector<int> input);

void printMatrix(std::vector<std::vector<double>> input);
void printMatrix(std::vector<std::vector<int>> input);

double dist2(std::vector<double> vec1, std::vector<double> vec2);

unsigned long nearestNeighbor(std::vector<std::vector<double>>map, std::vector<double> input);
std::vector<unsigned long> nearestNeighbors(std::vector<std::vector<double>> map, std::vector<std::vector<double>> input);

int hammingDist(std::vector<int> vec1, std::vector<int> vec2);

int minIndex(std::vector<double> vec)
#endif /* defined(__Project1_555__myUtil__) */
