//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/17/19.
//

//#include "Class_Matrix.h"
//
//MCMPS::Class_Matrix::Class_Matrix(MCMPS::type_BondDim _BondDim):
//Matrix(_BondDim, _BondDim, arma::fill::randu)
//{
//
//}
//
//void MCMPS::Class_Matrix::Set_Ele(MCMPS::type_BondDim Which_Row, MCMPS::type_BondDim Which_Column,
//                                  MCMPS::type_RealVal Which_Val) {
//    this->Matrix(Which_Row, Which_Column) = Which_Val;
//}
//
//const MCMPS::type_RealVal & MCMPS::Class_Matrix::Get_Ele(MCMPS::type_BondDim Which_Row, MCMPS::type_BondDim Which_Column) const{
//    return this->Matrix(Which_Row, Which_Column);
//}
//
//MCMPS::Class_Matrix MCMPS::Class_Matrix::operator * (const MCMPS::Class_Matrix &_Other_Matrix) const{
//    return MCMPS::Class_Matrix(this->Matrix * _Other_Matrix.Matrix);
//}
//
//
//
//MCMPS::Class_Matrix::Class_Matrix(const MCMPS::type_Mat &_Which_Matrix) :
//Matrix(_Which_Matrix)
//{
//}
//
//MCMPS::type_RealVal MCMPS::Class_Matrix::Get_trace() const{
//    return arma::trace(this->Matrix);
//}
//
//MCMPS::Class_Matrix MCMPS::Class_Matrix::t() const {
//    return MCMPS::Class_Matrix(this->Matrix.t());
//}
//
//MCMPS::Class_Matrix MCMPS::Class_Matrix::operator*(double _Other_Number) const {
//    return MCMPS::Class_Matrix(this->Matrix * _Other_Number);
//}
//
//void MCMPS::Class_Matrix::operator+=(const MCMPS::Class_Matrix & _Other_Matrix) {
//    this->Matrix += _Other_Matrix.Matrix;
//}
//
//void MCMPS::Class_Matrix::zeros() {
//    this->Matrix.zeros();
//}
//
//void MCMPS::Class_Matrix::operator*= (const MCMPS::type_RealVal _Which_number) {
//    this->Matrix *= _Which_number;
//}
//
//MCMPS::Class_Matrix MCMPS::Class_Matrix::Get_sign() const {
//    return MCMPS::Class_Matrix(arma::sign(Matrix));
//}

