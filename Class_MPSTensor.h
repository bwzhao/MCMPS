//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/16/19.
//

#pragma once

#include "Config.h"
#include "Class_Matrix.h"

namespace MCMPS {
    class Class_MPSTensor {
    private:
        // Store the elements in this tensor, in this case, the rank-2 tensor is just a matrix.
        // The matrix correspond to s=+1
        MCMPS::Class_Matrix Matrix_p;
        // The matrix correspond to s=-1
        MCMPS::Class_Matrix Matrix_m;

    public:
        Class_MPSTensor() = default;
        explicit Class_MPSTensor(type_BondDim _BondDim);
        ~Class_MPSTensor() = default;

        MCMPS::Class_Matrix & Get_Matrix(type_PyDegree _Which_Spin);
//        MCMPS::Class_Matrix & Get_Matrix_Z(type_PyDegree _Which_Spin);
        void Set_Zeros();

        void operator *= (type_RealVal _Which_number);
        void operator += (const Class_MPSTensor & _Other_Tensor);

        void Update_Tensor(const Class_MPSTensor &_D1, const Class_MPSTensor &_D2, type_RealVal _Val_ES,
                           type_RealVal _Val_delta);
    };
}



