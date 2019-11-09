//
// Created by Bowen Zhao on 10/16/19.
//
#pragma once
// Define some library-related parameters
//#define ARMA_NO_DEBUG
//#define ARMA_EXTRA_DEBUG
#include <armadillo>

#include <cstdint>
#include <string>

namespace MCMPS {
    using type_RealVal = double;
    using type_IntVal = int;
    using type_Mat = arma::Mat<MCMPS::type_RealVal>;

    using type_BondDim = arma::uword;
    using type_PyDegree = std::int_fast32_t;
    using type_NumSite = std::int_fast32_t;


    constexpr type_PyDegree SPINUP = 1;
    constexpr type_PyDegree SPINDOWN = -1;

    constexpr type_RealVal MATELE_Heisenberg_DIA = 0.25;
    constexpr type_RealVal MATELE_Heisenbergy_OFF = 0.5;

    constexpr type_RealVal MATELE_Ising_DIA = 1;
    constexpr type_RealVal MATELE_Ising_OFF = 1;

    constexpr type_RealVal DELTA0 = 1;
    constexpr type_RealVal DELTA_EXPONENT = 0.75;
    constexpr type_RealVal DELTA_Q= 0.9;

    constexpr type_IntVal P_0= 10;
    constexpr type_IntVal F_0= 100;

    constexpr type_IntVal K_MAX = 8;

    const std::vector<std::string> VEC_NAME = {{'E'}};

    inline MCMPS::type_IntVal P_k(int _Which_Step) {
        return P_0 * _Which_Step;
    }

    inline MCMPS::type_IntVal F_k(int _Which_Step) {
        return F_0 * _Which_Step;
    }
}
