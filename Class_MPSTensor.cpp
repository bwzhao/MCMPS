//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/16/19.
//

#include "Class_MPSTensor.h"

MCMPS::Class_MPSTensor::Class_MPSTensor(MCMPS::type_BondDim _BondDim):
        Matrix_p(_BondDim, _BondDim),
        Matrix_m(_BondDim, _BondDim)
{
    //initialize
    Matrix_p.randu();
    Matrix_m.randu();

    Matrix_p %= arma::sign(Matrix_p);
    Matrix_m %= arma::sign(Matrix_m);
}

MCMPS::Class_Matrix &MCMPS::Class_MPSTensor::Get_Matrix(type_PyDegree _Which_Spin) {
    if (_Which_Spin == SPINUP) {
        return Matrix_p;
    }
    else{
        return Matrix_m;
    }
}

//MCMPS::Class_Matrix &MCMPS::Class_MPSTensor::Get_Matrix_Z(type_PyDegree _Which_Spin) {
//    if (_Which_Spin == SPINUP) {
//        return Matrix_m;
//    }
//    else{
//        return Matrix_p;
//    }
//}

void MCMPS::Class_MPSTensor::Set_Zeros() {
    Matrix_p.zeros();
    Matrix_m.zeros();
}

void MCMPS::Class_MPSTensor::operator*=(MCMPS::type_RealVal _Which_number) {
    Matrix_p *= _Which_number;
    Matrix_m *= _Which_number;
}

void MCMPS::Class_MPSTensor::operator+=(const Class_MPSTensor & _Other_Tensor) {
    Matrix_m += _Other_Tensor.Matrix_m;
    Matrix_p += _Other_Tensor.Matrix_p;

}

void MCMPS::Class_MPSTensor::Update_Tensor(const Class_MPSTensor &_D1, const Class_MPSTensor &_D2, type_RealVal _Val_ES,
                                           type_RealVal _Val_delta) {
    // Update p matrix
    decltype(Matrix_p) dMatrix_p = arma::sign(_D2.Matrix_p * 2 - _D1.Matrix_p * (2 * _Val_ES));
    decltype(Matrix_p) random_p = arma::randu(arma::size(dMatrix_p)) * _Val_delta;
    Matrix_p -= (random_p) % (dMatrix_p);

    // Update m matrix
    decltype(Matrix_m) dMatrix_m = arma::sign(_D2.Matrix_m * 2 - _D1.Matrix_m * (2 * _Val_ES));
    decltype(Matrix_m) random_m = arma::randu(arma::size(dMatrix_m)) * _Val_delta;
    Matrix_m -= (random_m) % (dMatrix_m);

//    Matrix_p.print();
//    std::cout << std::endl;
//    Matrix_m.print();
//    std::cout << std::endl;

}
