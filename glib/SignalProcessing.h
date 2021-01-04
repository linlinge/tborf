#pragma once
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
void DaubechiesWavelet(vector<double>& dat,vector<double>& cA, vector<double>& cD);
void DaubechiesWavelet(vector<float>& dat,vector<double>& cA, vector<double>& cD);