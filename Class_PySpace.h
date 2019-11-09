//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/16/19.
//
#pragma once

#include "Config.h"
#include <vector>

namespace MCMPS {
    class Class_PySpace {
    private:
        // Array to store the values of the physical degrees of freedom
        std::vector<MCMPS::type_PyDegree> Array_PyDegree;

    public:
        Class_PySpace() = default;
        explicit Class_PySpace(type_NumSite _NumSite);
        ~Class_PySpace() = default;

        const MCMPS::type_NumSite size() const;

        const MCMPS::type_PyDegree Get_Spin_0(type_NumSite _Which_Site) const;
        const MCMPS::type_PyDegree Get_OppoSpin_0(type_NumSite _Which_Site) const;

        const MCMPS::type_PyDegree Get_Spin_Z(type_NumSite _Which_Site) const;
        const MCMPS::type_PyDegree Get_OppoSpin_Z(type_NumSite _Which_Site) const;

        const MCMPS::type_PyDegree Get_Spin_R(type_NumSite _Which_Site) const;
        const MCMPS::type_PyDegree Get_OppoSpin_R(type_NumSite _Which_Site) const;

        const MCMPS::type_PyDegree Get_Spin_RZ(type_NumSite _Which_Site) const;
        const MCMPS::type_PyDegree Get_OppoSpin_RZ(type_NumSite _Which_Site) const;

        const type_NumSite Get_RSite(type_NumSite _Which_Site) const;

        void Flip_Spin(type_NumSite _Which_Site);
    };
}




