//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/17/19.
//
#pragma once

#include <map>
#include <string>
#include <vector>
#include "Config.h"
#include <numeric>

namespace MCMPS {
    class Class_DataMeasurement {
    private:
        std::vector<MCMPS::type_RealVal> Data;

    public:
        Class_DataMeasurement() = default;
        ~Class_DataMeasurement() = default;

        // Add a measurement value to a specific quantity(key)
        void AppendValue(const MCMPS::type_RealVal _which_value){
            Data.emplace_back(_which_value);
        }
        MCMPS::type_RealVal Get_AveValue(){
            return std::accumulate(Data.begin(), Data.end(), 0.) / Data.size();
        }

        void ClearValue(){
            Data.clear();
        }

    };
}



