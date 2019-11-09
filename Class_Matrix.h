//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/17/19.
//

#pragma once

#include "Config.h"
#include <armadillo>

namespace MCMPS{
    using Class_Matrix = MCMPS::type_Mat;
//class Class_Matrix {
//    private:
//        MCMPS::type_Mat Matrix;
//    public:
//        Class_Matrix() = default;
//        explicit Class_Matrix(MCMPS::type_BondDim _BondDim);
//        explicit Class_Matrix(const MCMPS::type_Mat &_Which_Matrix);
//        ~Class_Matrix() = default;
//
//        // Set the value for a specific element
//        void Set_Ele(MCMPS::type_BondDim Which_Row, MCMPS::type_BondDim Which_Column, MCMPS::type_RealVal Which_Val);
//        // Get the value of a specific element
//        const MCMPS::type_RealVal & Get_Ele(MCMPS::type_BondDim Which_Row, MCMPS::type_BondDim Which_Column) const;
//        // Matrix Multiplication
//        inline MCMPS::Class_Matrix operator * (const MCMPS::Class_Matrix & _Other_Matrix) const;
//
//        inline MCMPS::Class_Matrix operator * (double _Other_Number) const;
//
//        inline void operator += (const MCMPS::Class_Matrix & _Other_Matrix);
//        inline void operator *= (const MCMPS::type_RealVal _Which_number);
//
//        inline MCMPS::Class_Matrix Get_sign () const;
//
//        // Get the trace of this matrix
//        MCMPS::type_RealVal Get_trace() const;
//
//        // Transpose
//        MCMPS::Class_Matrix t() const;
//
//        // Set to zero
//        void zeros();
//    };
}



