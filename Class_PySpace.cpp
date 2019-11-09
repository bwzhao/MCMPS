//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/16/19.
//

#include "Class_PySpace.h"
#include <algorithm>
#include <random>

MCMPS::Class_PySpace::Class_PySpace(MCMPS::type_NumSite _NumSite):
Array_PyDegree(_NumSite, -1)
{
    //initialize the spin and VB
    std::vector<type_NumSite> SitaA_Array;
    std::vector<type_NumSite> SitaB_Array;
    for (type_NumSite index_site = 0; index_site != Array_PyDegree.size(); ++index_site) {
        if (index_site % 2 == 0) {
            SitaA_Array.push_back(index_site);
        }
        else {
            SitaB_Array.push_back(index_site);
        }
    }

    std::random_device ran_dev;
    std::mt19937 g(ran_dev());
    static std::uniform_int_distribution<MCMPS::type_NumSite> Ran_Int(0, 1);

    std::shuffle(SitaA_Array.begin(), SitaA_Array.end(), g);
    std::shuffle(SitaB_Array.begin(), SitaB_Array.end(), g);

    for(type_NumSite index = 0; index != SitaA_Array.size(); ++index){
        auto siteA = SitaA_Array[index];
        auto siteB = SitaB_Array[index];

        //Site
        auto temp_color = Ran_Int(ran_dev) * 2 - 1;
        Array_PyDegree[siteA] = temp_color;
        Array_PyDegree[siteB] = -temp_color;
    }
}

const MCMPS::type_NumSite MCMPS::Class_PySpace::size() const{
    return Array_PyDegree.size();
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_Spin_0(type_NumSite _Which_Site) const {
    return Array_PyDegree[_Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_OppoSpin_0(type_NumSite _Which_Site) const {
    return -Array_PyDegree[_Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_Spin_Z(type_NumSite _Which_Site) const {
    return -Array_PyDegree[_Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_OppoSpin_Z(type_NumSite _Which_Site) const {
    return Array_PyDegree[_Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_Spin_R(type_NumSite _Which_Site) const {
    return Array_PyDegree[Array_PyDegree.size() - 1 - _Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_OppoSpin_R(type_NumSite _Which_Site) const {
    return -Array_PyDegree[Array_PyDegree.size() - 1 - _Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_Spin_RZ(type_NumSite _Which_Site) const {
    return -Array_PyDegree[Array_PyDegree.size() - 1 - _Which_Site];
}

const MCMPS::type_PyDegree MCMPS::Class_PySpace::Get_OppoSpin_RZ(type_NumSite _Which_Site) const {
    return Array_PyDegree[Array_PyDegree.size() - 1 - _Which_Site];
}

void MCMPS::Class_PySpace::Flip_Spin(MCMPS::type_NumSite _Which_Site) {
    Array_PyDegree[_Which_Site] = -Array_PyDegree[_Which_Site];
}

const MCMPS::type_NumSite MCMPS::Class_PySpace::Get_RSite(MCMPS::type_NumSite _Which_Site) const {
    return (Array_PyDegree.size() - 1 - _Which_Site);
}


