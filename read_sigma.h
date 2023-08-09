#pragma once
#include "parameters.h"
#include <iostream>
#include <vector>
#include <complex> //this contains complex numbers and trig functions
#include <fstream>
#include <cmath>
#include <string>


void get_sigma(Parameters &parameters, std::vector<std::vector<dcomp>> &sigma_mb_r_up, std::vector<std::vector<dcomp>> &sigma_mb_r_down,
    std::vector<std::vector<dcomp>> &sigma_mb_l_up, std::vector<std::vector<dcomp>> &sigma_mb_l_down);

//void read_self_energy(Parameters &parameters, std::vector<std::vector<dcomp>> &sigma_mb_up, std::vector<std::vector<dcomp>> &sigma_mb_down, string filename);