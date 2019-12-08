//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/16/19.
//

#include "Class_Network.h"

MCMPS::Class_Network::Class_Network(MCMPS::type_NumSite _NumSite, MCMPS::type_BondDim _BondDim, int _Which_CPU, std::string _TypeUpdate):
        TypeUpdate(_TypeUpdate),
        PySpace(_NumSite),
        Tensor_A(_BondDim),
        Tensor_B(_BondDim),
// Initialize the pre-calculated matrices as the random matrices
        Array_L_0(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_R_0(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_G_0(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_L_R(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_R_R(_NumSite,   Class_Matrix(_BondDim, _BondDim)),
        Array_G_R(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_L_Z(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_R_Z(_NumSite,   Class_Matrix(_BondDim, _BondDim)),
        Array_G_Z(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_L_RZ(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Array_R_RZ(_NumSite,   Class_Matrix(_BondDim, _BondDim)),
        Array_G_RZ(_NumSite, Class_Matrix(_BondDim, _BondDim)),
        Matrix_Id(MCMPS::type_Mat(_BondDim, _BondDim, arma::fill::eye)),
// Initialize the weight as 1
        Val_WS_0(1.),
        Val_WS_Z(1.),
        Val_WS_R(1.),
        Val_WS_RZ(1.),
        Which_CPU(_Which_CPU),
        Num_Measure(0),
        Tensor_BinD1_A(_BondDim),
        Tensor_BinD1_B(_BondDim),
        Tensor_BinD2_A(_BondDim),
        Tensor_BinD2_B(_BondDim),
        Val_BinEnergy(0.),
        Val_BinEnergy2(0.)
{
    Tensor_BinD1_A.Set_Zeros();
    Tensor_BinD1_B.Set_Zeros();
    Tensor_BinD2_A.Set_Zeros();
    Tensor_BinD2_B.Set_Zeros();

    Cal_R_0();
    Cal_L_0();
    Cal_R_Z();
    Cal_L_Z();
    Cal_R_R();
    Cal_L_R();
    Cal_R_RZ();
    Cal_L_RZ();

    MCMPS::Class_DataMeasurement temp_Class_DataMeasurement;
    for (const auto &which_key: MCMPS::VEC_NAME) {
        MeasureData.insert(std::make_pair(which_key, temp_Class_DataMeasurement));
    }

    char temp_chars_nameclass[100];
    snprintf(temp_chars_nameclass, 100, "Heisenberg_%d_%d_%d_",
             _NumSite,
             _BondDim,
             _Which_CPU);
    Str_NameFile = std::string(temp_chars_nameclass) + _TypeUpdate;
}

void MCMPS::Class_Network::Cal_L_0() {
    // Firstly set the matrix at the first site
    Array_L_0[0] = this->Get_Tensor_0(0).Get_Matrix(PySpace.Get_Spin_0(0));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = index_count;
        Array_L_0[index_site] =
                Array_L_0[index_site - 1] * this->Get_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site));
    }
}
void MCMPS::Class_Network::Cal_L_Z() {
    // Firstly set the matrix at the first site
    Array_L_Z[0] = this->Get_Tensor_Z(0).Get_Matrix(PySpace.Get_Spin_Z(0));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = index_count;
        Array_L_Z[index_site] =
                Array_L_Z[index_site - 1] * this->Get_Tensor_Z(index_site).Get_Matrix(PySpace.Get_Spin_Z(index_site));
    }
}
void MCMPS::Class_Network::Cal_L_R() {
    // Firstly set the matrix at the first site
    Array_L_R[0] = this->Get_Tensor_R(0).Get_Matrix(PySpace.Get_Spin_R(0));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = index_count;
        Array_L_R[index_site] =
                Array_L_R[index_site - 1] * this->Get_Tensor_R(index_site).Get_Matrix(PySpace.Get_Spin_R(index_site));
    }
}
void MCMPS::Class_Network::Cal_L_RZ() {
    // Firstly set the matrix at the first site
    Array_L_RZ[0] = this->Get_Tensor_RZ(0).Get_Matrix(PySpace.Get_Spin_RZ(0));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = index_count;
        Array_L_RZ[index_site] =
                Array_L_RZ[index_site - 1] * this->Get_Tensor_RZ(index_site).Get_Matrix(
                        PySpace.Get_Spin_RZ(index_site));
    }
}

void MCMPS::Class_Network::Cal_R_0() {
    // Firstly set the matrix at the last site
    Array_R_0[PySpace.size() - 1] = this->Get_Tensor_0(PySpace.size() - 1).Get_Matrix(
            PySpace.Get_Spin_0(PySpace.size() - 1));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = PySpace.size() - 1 - index_count;
        Array_R_0[index_site] =
                this->Get_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site)) * Array_R_0[index_site + 1];
    }
}
void MCMPS::Class_Network::Cal_R_Z() {
    // Firstly set the matrix at the last site
    Array_R_Z[PySpace.size() - 1] = this->Get_Tensor_Z(PySpace.size() - 1).Get_Matrix(
            PySpace.Get_Spin_Z(PySpace.size() - 1));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = PySpace.size() - 1 - index_count;
        Array_R_Z[index_site] =
                this->Get_Tensor_Z(index_site).Get_Matrix(PySpace.Get_Spin_Z(index_site)) * Array_R_Z[index_site + 1];
    }
}
void MCMPS::Class_Network::Cal_R_R() {
    // Firstly set the matrix at the last site
    Array_R_R[PySpace.size() - 1] = this->Get_Tensor_R(PySpace.size() - 1).Get_Matrix(
            PySpace.Get_Spin_R(PySpace.size() - 1));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = PySpace.size() - 1 - index_count;
        Array_R_R[index_site] =
                this->Get_Tensor_R(index_site).Get_Matrix(PySpace.Get_Spin_R(index_site)) * Array_R_R[index_site + 1];
    }
}
void MCMPS::Class_Network::Cal_R_RZ() {
    // Firstly set the matrix at the last site
    Array_R_RZ[PySpace.size() - 1] = this->Get_Tensor_RZ(PySpace.size() - 1).Get_Matrix(
            PySpace.Get_Spin_RZ(PySpace.size() - 1));
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size(); ++index_count) {
        auto index_site = PySpace.size() - 1 - index_count;
        Array_R_RZ[index_site] =
                this->Get_Tensor_RZ(index_site).Get_Matrix(PySpace.Get_Spin_RZ(index_site)) * Array_R_RZ[index_site + 1];
    }
}

void MCMPS::Class_Network::Cal_G_0() {
    // Firstly set the matrix at the first site
    Array_G_0[0] = Array_R_0[1];
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site = index_count;
        Array_G_0[index_site] = (Array_R_0[index_site + 1] * Array_L_0[index_site - 1]);
    }
    Array_G_0[PySpace.size() - 1] = Array_L_0[PySpace.size() - 2];
}
void MCMPS::Class_Network::Cal_G_Z() {
    // Firstly set the matrix at the first site
    Array_G_Z[0] = Array_R_Z[1];
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site = index_count;
        Array_G_Z[index_site] = (Array_R_Z[index_site + 1] * Array_L_Z[index_site - 1]);
    }
    Array_G_Z[PySpace.size() - 1] = Array_L_Z[PySpace.size() - 2];
}
void MCMPS::Class_Network::Cal_G_R() {
    // Firstly set the matrix at the first site
    Array_G_R[0] = Array_R_R[1];
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site = index_count;
        Array_G_R[index_site] = (Array_R_R[index_site + 1] * Array_L_R[index_site - 1]);
    }
    Array_G_R[PySpace.size() - 1] = Array_L_R[PySpace.size() - 2];
}
void MCMPS::Class_Network::Cal_G_RZ() {
    // Firstly set the matrix at the first site
    Array_G_RZ[0] = Array_R_RZ[1];
    for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site = index_count;
        Array_G_RZ[index_site] = (Array_R_RZ[index_site + 1] * Array_L_RZ[index_site - 1]);
    }
    Array_G_RZ[PySpace.size() - 1] = Array_L_RZ[PySpace.size() - 2];
}

MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_Tensor_0(MCMPS::type_NumSite Which_Site) {
    if (Which_Site % 2 == 0) {
        return Tensor_A;
    }
    else{
        return Tensor_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_Tensor_Z(MCMPS::type_NumSite Which_Site) {
    if (Which_Site % 2 == 0) {
        return Tensor_A;
    }
    else{
        return Tensor_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_Tensor_R(MCMPS::type_NumSite Which_Site) {
    if (PySpace.Get_RSite(Which_Site) % 2 == 0) {
        return Tensor_A;
    }
    else{
        return Tensor_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_Tensor_RZ(MCMPS::type_NumSite Which_Site) {
    if (PySpace.Get_RSite(Which_Site) % 2 == 0){
        return Tensor_A;
    }
    else{
        return Tensor_B;
    }
}

MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D1_Tensor_0(MCMPS::type_NumSite Which_Site) {
    if (Which_Site % 2 == 0) {
        return Tensor_BinD1_A;
    }
    else{
        return Tensor_BinD1_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D1_Tensor_Z(MCMPS::type_NumSite Which_Site) {
    if (Which_Site % 2 == 0) {
        return Tensor_BinD1_A;
    }
    else{
        return Tensor_BinD1_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D1_Tensor_R(MCMPS::type_NumSite Which_Site) {
    if (PySpace.Get_RSite(Which_Site) % 2 == 0) {
        return Tensor_BinD1_A;
    }
    else{
        return Tensor_BinD1_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D1_Tensor_RZ(MCMPS::type_NumSite Which_Site) {
    if (PySpace.Get_RSite(Which_Site) % 2 == 0) {
        return Tensor_BinD1_A;
    }
    else{
        return Tensor_BinD1_B;
    }
}

MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D2_Tensor_0(MCMPS::type_NumSite Which_Site) {
    if (Which_Site % 2 == 0) {
        return Tensor_BinD2_A;
    }
    else{
        return Tensor_BinD2_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D2_Tensor_Z(MCMPS::type_NumSite Which_Site) {
    if (Which_Site % 2 == 0) {
        return Tensor_BinD2_A;
    }
    else{
        return Tensor_BinD2_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D2_Tensor_R(MCMPS::type_NumSite Which_Site) {
    if (PySpace.Get_RSite(Which_Site) % 2 == 0) {
        return Tensor_BinD2_A;
    }
    else{
        return Tensor_BinD2_B;
    }
}
MCMPS::Class_MPSTensor &MCMPS::Class_Network::Get_D2_Tensor_RZ(MCMPS::type_NumSite Which_Site) {
    if (PySpace.Get_RSite(Which_Site) % 2 == 0) {
        return Tensor_BinD2_A;
    }
    else{
        return Tensor_BinD2_B;
    }
}

void MCMPS::Class_Network::Update_Spin_Pair() {
    // Random number generator
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr) + Which_CPU));
    static std::uniform_real_distribution<double> Ran_Double(0., 1.);

    // Update from left to right, so calculate array_R first
    Cal_R_0();
    // Get current weight
    Val_WS_0 =  arma::trace(Array_R_0[0]);

    // Update pairs of spin
    // First for all the pairs except for the last
    const Class_Matrix* ptr_Lmm1 = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp2 = &this->Matrix_Id;

    for (MCMPS::type_NumSite index_count = 0; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site_0 = index_count;
        auto index_site_1 = index_count + 1;
        if (index_site_0 == 0){
            ptr_Lmm1 = &this->Matrix_Id;
        }
        else{
            auto index_site_m1 = index_site_0 - 1;
            ptr_Lmm1 = &Array_L_0[index_site_m1];
        }

        if (PySpace.Get_Spin_0(index_site_0) != PySpace.Get_Spin_0(index_site_1)){
            if (index_count == PySpace.size() - 2){
                ptr_Rmp2 = &this->Matrix_Id;
            }
            else{
                auto index_site_2 = index_count + 2;
                ptr_Rmp2 = &Array_R_0[index_site_2];
            }

            auto Val_WSp = arma::trace((*ptr_Lmm1)
                    * this->Get_Tensor_0(index_site_0).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_0))
                    * this->Get_Tensor_0(index_site_1).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_1))
                    * (*ptr_Rmp2));
            auto temp_ratio = Val_WSp * Val_WSp / Val_WS_0 / Val_WS_0;

            // If we accept the change
            if (Ran_Double(engine) < temp_ratio){
                Val_WS_0 = Val_WSp;
                PySpace.Flip_Spin(index_site_0);
                PySpace.Flip_Spin(index_site_1);
                Array_L_0[index_site_0] = (*ptr_Lmm1)
                        * this->Get_Tensor_0(index_site_0).Get_Matrix(PySpace.Get_Spin_0(index_site_0));
            }
            else{
                Array_L_0[index_site_0] = (*ptr_Lmm1)
                        * this->Get_Tensor_0(index_site_0).Get_Matrix(PySpace.Get_Spin_0(index_site_0));
            }
        }
        else{
            Array_L_0[index_site_0] = (*ptr_Lmm1)
                    * this->Get_Tensor_0(index_site_0).Get_Matrix(PySpace.Get_Spin_0(index_site_0));
        }
    }
    // Then is the pair across the boundary
    {
        const MCMPS::type_NumSite index_site_left = 0;
        const MCMPS::type_NumSite index_site_right = PySpace.size() - 1;

        if (PySpace.Get_Spin_0(index_site_left) != PySpace.Get_Spin_0(index_site_right)){
            // Calculate current Array_R
            Array_R_0[index_site_right] = this->Get_Tensor_0(index_site_right).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_right));
            for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
                auto index_site = index_site_right - index_count;
                Array_R_0[index_site] =
                        this->Get_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site)) * Array_R_0[index_site + 1];
            }

            Array_R_0[index_site_left] = this->Get_Tensor_0(index_site_left).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_left))
                    * Array_R_0[index_site_left + 1];
            auto Val_WSp = arma::trace(Array_R_0[index_site_left]);
            auto temp_ratio = Val_WSp * Val_WSp / Val_WS_0 / Val_WS_0;

            // If we accept the change
            if (Ran_Double(engine) < temp_ratio){
                Val_WS_0 = Val_WSp;
                PySpace.Flip_Spin(index_site_left);
                PySpace.Flip_Spin(index_site_right);
                Cal_L_0();
//                Pre_Cal_R();
                assert(arma::approx_equal(Array_L_0[index_site_right], Array_R_0[0], "reldiff", 0.0002));
            }
            else{
                Array_L_0[index_site_right] = Array_L_0[index_site_right - 1] *
                        this->Get_Tensor_0(index_site_right).Get_Matrix(
                                PySpace.Get_Spin_0(index_site_right));
                Cal_R_0();
//                Pre_Cal_L();
                assert(arma::approx_equal(Array_L_0[index_site_right], Array_R_0[0], "reldiff", 0.0002));
            }
        }
        else{
            Array_L_0[index_site_right] = Array_L_0[index_site_right - 1] *
                    this->Get_Tensor_0(index_site_right).Get_Matrix(PySpace.Get_Spin_0(index_site_right));
            Cal_R_0();
            assert(arma::approx_equal(Array_L_0[index_site_right], Array_R_0[0], "reldiff", 0.0002));
        }
    }
}
void MCMPS::Class_Network::Update_Spin_Pair_Sym() {
    // Random number generator
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr) + Which_CPU));
    static std::uniform_real_distribution<double> Ran_Double(0., 1.);

    // Update from left to right, so calculate array_R first
    Cal_R_0();
    Cal_R_Z();
    Cal_L_R();
    Cal_L_RZ();

    // Get current weight
    Val_WS_0 =  arma::trace(Array_R_0[0]);
    Val_WS_Z =  arma::trace(Array_R_Z[0]);
    Val_WS_R =  arma::trace(Array_L_R[PySpace.size() - 1]);
    Val_WS_RZ =  arma::trace(Array_L_RZ[PySpace.size() - 1]);

    // Update pairs of spin
    // First for all the pairs except for the last

    const Class_Matrix* ptr_Lmm1_0 = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp2_0 = &this->Matrix_Id;

    const Class_Matrix* ptr_Lmm1_Z = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp2_Z = &this->Matrix_Id;

    const Class_Matrix* ptr_Rmp1_R = &this->Matrix_Id;
    const Class_Matrix* ptr_Lmm2_R = &this->Matrix_Id;

    const Class_Matrix* ptr_Rmp1_RZ = &this->Matrix_Id;
    const Class_Matrix* ptr_Lmm2_RZ = &this->Matrix_Id;

    for (MCMPS::type_NumSite index_count = 0; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site0_0 = index_count;
        auto index_site1_0 = index_count + 1;
        auto index_site0_R = PySpace.Get_RSite(index_site0_0);
        auto index_site1_R = PySpace.Get_RSite(index_site1_0);

        if (index_site0_0 == 0){
            ptr_Lmm1_0 = &this->Matrix_Id;
            ptr_Lmm1_Z = &this->Matrix_Id;
            ptr_Rmp1_R = &this->Matrix_Id;
            ptr_Rmp1_RZ = &this->Matrix_Id;
        }
        else{
            auto index_sitem1_0 = index_site0_0 - 1;
            auto index_sitem1_R = PySpace.Get_RSite(index_sitem1_0);
            ptr_Lmm1_0 = &Array_L_0[index_sitem1_0];
            ptr_Lmm1_Z = &Array_L_Z[index_sitem1_0];
            ptr_Rmp1_R = &Array_R_R[index_sitem1_R];
            ptr_Rmp1_RZ = &Array_R_RZ[index_sitem1_R];
        }

        if (PySpace.Get_Spin_0(index_site0_0) != PySpace.Get_Spin_0(index_site1_0)){
            if (index_count == PySpace.size() - 2){
                ptr_Rmp2_0 = &this->Matrix_Id;
                ptr_Rmp2_Z = &this->Matrix_Id;
                ptr_Lmm2_R = &this->Matrix_Id;
                ptr_Lmm2_RZ = &this->Matrix_Id;
            }
            else{
                auto index_site2_0 = index_count + 2;
                auto index_site2_R = PySpace.Get_RSite(index_site2_0);
                ptr_Rmp2_0 = &Array_R_0[index_site2_0];
                ptr_Rmp2_Z = &Array_R_Z[index_site2_0];
                ptr_Lmm2_R = &Array_L_R[index_site2_R];
                ptr_Lmm2_RZ = &Array_L_RZ[index_site2_R];
            }

            auto Val_WSp_0 = arma::trace((*ptr_Lmm1_0) *
                                         this->Get_Tensor_0(index_site0_0).Get_Matrix(
                                                         PySpace.Get_OppoSpin_0(index_site0_0)) *
                                         this->Get_Tensor_0(index_site1_0).Get_Matrix(
                                                         PySpace.Get_OppoSpin_0(index_site1_0)) *
                                         (*ptr_Rmp2_0));

            auto Val_WSp_Z = arma::trace((*ptr_Lmm1_Z) *
                                         this->Get_Tensor_Z(index_site0_0).Get_Matrix(
                                                         PySpace.Get_OppoSpin_Z(index_site0_0)) *
                                         this->Get_Tensor_Z(index_site1_0).Get_Matrix(
                                                         PySpace.Get_OppoSpin_Z(index_site1_0)) *
                                         (*ptr_Rmp2_Z));

            auto Val_WSp_R = arma::trace((*ptr_Lmm2_R) *
                                         this->Get_Tensor_R(index_site1_R).Get_Matrix(
                                                 PySpace.Get_OppoSpin_R(index_site1_R)) *
                                         this->Get_Tensor_R(index_site0_R).Get_Matrix(
                                                         PySpace.Get_OppoSpin_R(index_site0_R)) *
                                         (*ptr_Rmp1_R));

            auto Val_WSp_RZ = arma::trace((*ptr_Lmm2_RZ) *
                                          this->Get_Tensor_RZ(index_site1_R).Get_Matrix(
                                                          PySpace.Get_OppoSpin_RZ(index_site1_R)) *
                                          this->Get_Tensor_RZ(index_site0_R).Get_Matrix(
                                                          PySpace.Get_OppoSpin_RZ(index_site0_R)) *
                                          (*ptr_Rmp1_RZ));

            auto Val_WSp = Val_WSp_0 + Val_WSp_Z + Val_WSp_R + Val_WSp_RZ;
            auto Val_WS = Val_WS_0 + Val_WS_Z + Val_WS_R + Val_WS_RZ;
            auto temp_ratio = Val_WSp * Val_WSp / Val_WS / Val_WS;

            // If we accept the change
            if (Ran_Double(engine) < temp_ratio){
                Val_WS_0 = Val_WSp_0;
                Val_WS_R = Val_WSp_R;
                Val_WS_Z = Val_WSp_Z;
                Val_WS_RZ = Val_WSp_RZ;

                PySpace.Flip_Spin(index_site0_0);
                PySpace.Flip_Spin(index_site1_0);
            }
        }
        Array_L_0[index_site0_0] = (*ptr_Lmm1_0) * this->Get_Tensor_0(index_site0_0).Get_Matrix(
                PySpace.Get_Spin_0(index_site0_0));

        Array_L_Z[index_site0_0] = (*ptr_Lmm1_Z) * this->Get_Tensor_Z(index_site0_0).Get_Matrix(
                PySpace.Get_Spin_Z(index_site0_0));

        Array_R_R[index_site0_R] = this->Get_Tensor_R(index_site0_R).Get_Matrix(
                PySpace.Get_Spin_R(index_site0_R)) * (*ptr_Rmp1_R);

        Array_R_RZ[index_site0_R] = this->Get_Tensor_RZ(index_site0_R).Get_Matrix(
                PySpace.Get_Spin_RZ(index_site0_R)) * (*ptr_Rmp1_RZ);
    }
    // Then is the pair across the boundary
    {
        const MCMPS::type_NumSite index_site_left = 0;
        const MCMPS::type_NumSite index_site_right = PySpace.size() - 1;

        if (PySpace.Get_Spin_0(index_site_left) != PySpace.Get_Spin_0(index_site_right)){
            // Assume we can flip the spins
            PySpace.Flip_Spin(index_site_left);
            PySpace.Flip_Spin(index_site_right);

            // Calculate current Array_R
            Cal_R_0();
            Cal_R_Z();
            Cal_L_R();
            Cal_L_RZ();

            auto Val_WSp_0 = arma::trace(Array_R_0[index_site_left]);
            auto Val_WSp_Z = arma::trace(Array_R_Z[index_site_left]);
            auto Val_WSp_R = arma::trace(Array_L_R[index_site_right]);
            auto Val_WSp_RZ = arma::trace(Array_L_RZ[index_site_right]);

            auto Val_WSp = Val_WSp_0 + Val_WSp_Z + Val_WSp_R + Val_WSp_RZ;
            auto Val_WS = Val_WS_0 +  Val_WS_Z + Val_WS_R + Val_WS_RZ;
            auto temp_ratio = Val_WSp * Val_WSp / Val_WS / Val_WS;

            // If we accept the change
            if (Ran_Double(engine) < temp_ratio){
                Val_WS_0 = Val_WSp_0;
                Val_WS_R = Val_WSp_R;
                Val_WS_Z = Val_WSp_Z;
                Val_WS_RZ = Val_WSp_RZ;

                Cal_L_0();
                Cal_L_Z();
                Cal_R_R();
                Cal_R_RZ();

                assert(arma::approx_equal(Array_L_0[index_site_right], Array_R_0[0], "reldiff", 0.0002));
                assert(arma::approx_equal(Array_L_Z[index_site_right], Array_R_Z[0], "reldiff", 0.0002));
                assert(arma::approx_equal(Array_L_R[index_site_right], Array_R_R[0], "reldiff", 0.0002));
                assert(arma::approx_equal(Array_L_RZ[index_site_right], Array_R_RZ[0], "reldiff", 0.0002));

                assert(abs((Val_WS_0 - arma::trace(Array_R_0[index_site_left]))/ Val_WS_0) <= 0.00001);
                assert(abs((Val_WS_Z - arma::trace(Array_R_Z[index_site_left])) / Val_WS_Z) <= 0.00001);
                assert(abs((Val_WS_R - arma::trace(Array_R_R[index_site_left])) / Val_WS_Z) <= 0.00001);
                assert(abs((Val_WS_RZ - arma::trace(Array_R_RZ[index_site_left])) / Val_WS_RZ) <= 0.00001);
            }
            else{
                // Flip back
                PySpace.Flip_Spin(index_site_left);
                PySpace.Flip_Spin(index_site_right);

                Array_L_0[index_site_right] = Array_L_0[index_site_right - 1] *
                        this->Get_Tensor_0(index_site_right).Get_Matrix(
                                PySpace.Get_Spin_0(index_site_right));

                Array_L_Z[index_site_right] = Array_L_Z[index_site_right - 1] *
                        this->Get_Tensor_Z(index_site_right).Get_Matrix(
                                PySpace.Get_Spin_Z(index_site_right));

                Array_R_R[index_site_left] =
                        this->Get_Tensor_R(index_site_left).Get_Matrix(
                                PySpace.Get_Spin_R(index_site_left)) * Array_R_R[index_site_left + 1];

                Array_R_RZ[index_site_left] =
                        this->Get_Tensor_RZ(index_site_left).Get_Matrix(
                                PySpace.Get_Spin_RZ(index_site_left)) * Array_R_RZ[index_site_left + 1];

                Cal_R_0();
                Cal_R_Z();
                Cal_L_R();
                Cal_L_RZ();

                assert(arma::approx_equal(Array_L_0[index_site_right], Array_R_0[0], "reldiff", 0.0002));
                assert(arma::approx_equal(Array_L_Z[index_site_right], Array_R_Z[0], "reldiff", 0.0002));
                assert(arma::approx_equal(Array_L_R[index_site_right], Array_R_R[0], "reldiff", 0.0002));
                assert(arma::approx_equal(Array_L_RZ[index_site_right], Array_R_RZ[0], "reldiff", 0.0002));

                assert(abs((Val_WS_0 - arma::trace(Array_R_0[index_site_left]))/ Val_WS_0) <= 0.00001);
                assert(abs((Val_WS_Z - arma::trace(Array_R_Z[index_site_left])) / Val_WS_Z) <= 0.00001);
                assert(abs((Val_WS_R - arma::trace(Array_R_R[index_site_left])) / Val_WS_Z) <= 0.00001);
                assert(abs((Val_WS_RZ - arma::trace(Array_R_RZ[index_site_left])) / Val_WS_RZ) <= 0.00001);
            }
        }
        else{
            Array_L_0[index_site_right] = Array_L_0[index_site_right - 1] *
                                          this->Get_Tensor_0(index_site_right).Get_Matrix(
                                                  PySpace.Get_Spin_0(index_site_right));
            Array_L_Z[index_site_right] =  Array_L_Z[index_site_right - 1] *
                                           this->Get_Tensor_Z(index_site_right).Get_Matrix(
                                                   PySpace.Get_Spin_Z(index_site_right));
            Array_R_R[index_site_left] =
                    this->Get_Tensor_R(index_site_left).Get_Matrix(
                            PySpace.Get_Spin_R(index_site_left)) * Array_R_R[index_site_left + 1];
            Array_R_RZ[index_site_left] =
                    this->Get_Tensor_RZ(index_site_left).Get_Matrix(
                            PySpace.Get_Spin_RZ(index_site_left)) * Array_R_RZ[index_site_left + 1];

            Cal_R_0();
            Cal_R_Z();
            Cal_L_R();
            Cal_L_RZ();

            assert(arma::approx_equal(Array_L_0[index_site_right], Array_R_0[0], "reldiff", 0.0002));
            assert(arma::approx_equal(Array_L_Z[index_site_right], Array_R_Z[0], "reldiff", 0.0002));
            assert(arma::approx_equal(Array_L_R[index_site_right], Array_R_R[0], "reldiff", 0.0002));
            assert(arma::approx_equal(Array_L_RZ[index_site_right], Array_R_RZ[0], "reldiff", 0.0002));

            assert(abs((Val_WS_0 - arma::trace(Array_R_0[index_site_left]))/ Val_WS_0) <= 0.00001);
            assert(abs((Val_WS_Z - arma::trace(Array_R_Z[index_site_left])) / Val_WS_Z) <= 0.00001);
            assert(abs((Val_WS_R - arma::trace(Array_R_R[index_site_left])) / Val_WS_Z) <= 0.00001);
            assert(abs((Val_WS_RZ - arma::trace(Array_R_RZ[index_site_left])) / Val_WS_RZ) <= 0.00001);
        }
    }



    // At the end of the update, both Array_L and Array_R have been calculated.
    // But it seems that is a bug that I have to explicitly calculate them
//    Pre_Cal_L();
//    Pre_Cal_R();
}
void MCMPS::Class_Network::Update_Spin_Single() {
    // Random number generator
    static std::mt19937_64 engine(static_cast<unsigned>(std::time(nullptr) + Which_CPU));
    static std::uniform_real_distribution<double> Ran_Double(0., 1.);

    // Update from left to right, so calculate array_R first
    Cal_R_0();
    // Get current weight
    Val_WS_0 =  arma::trace(Array_R_0[0]);

    // Update pairs of spin
    // First for all the pairs except for the last
    const Class_Matrix* ptr_Lmm1 = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp1 = &this->Matrix_Id;
    for (MCMPS::type_NumSite index_count = 0; index_count != PySpace.size(); ++index_count) {
        auto index_site_0 = index_count;
        auto index_site_1 = index_count + 1;

        if (index_count == 0){
            ptr_Lmm1 = &this->Matrix_Id;
        }
        else{
            auto index_site_m1 = index_count - 1;
            ptr_Lmm1 = &Array_L_0[index_site_m1];
        }

        if (index_count == PySpace.size() - 1){
            ptr_Rmp1 = &this->Matrix_Id;
        }
        else{
            ptr_Rmp1 = &Array_R_0[index_site_1];
        }

        auto Val_WSp = arma::trace(( (*ptr_Lmm1) * this->Get_Tensor_0(index_site_0).Get_Matrix(
                PySpace.Get_OppoSpin_0(index_site_0)) *
                                    (*ptr_Rmp1)));
        auto temp_ratio = Val_WSp * Val_WSp / Val_WS_0 / Val_WS_0;

        // If we accept the change
        if (Ran_Double(engine) < temp_ratio){
            Val_WS_0 = Val_WSp;
            PySpace.Flip_Spin(index_site_0);
            PySpace.Flip_Spin(index_site_1);
            Array_L_0[index_site_0] = (*ptr_Lmm1) * this->Get_Tensor_0(index_site_0).Get_Matrix(
                    PySpace.Get_Spin_0(index_site_0));
        }
        else{
            Array_L_0[index_site_0] = (*ptr_Lmm1) * this->Get_Tensor_0(index_site_0).Get_Matrix(
                    PySpace.Get_Spin_0(index_site_0));
        }
    }
    Cal_R_0();

    assert(arma::approx_equal(Array_L_0[PySpace.size() - 1], Array_R_0[0], "reldiff", 0.0002));
    // At the end of the update, both Array_L and Array_R have been calculated.
}

void MCMPS::Class_Network::Measure_Heisenberg() {
    // Measure the energy at current spin configuration
    MCMPS::type_RealVal temp_ES_Dia = 0.;
    MCMPS::type_RealVal temp_ES_Off = 0.;

    // Update pairs of spin
    // First for all the pairs except for the last
    const Class_Matrix* ptr_Lmm1 = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp2 = &this->Matrix_Id;
    for (MCMPS::type_NumSite index_count = 0; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site_0 = index_count;
        auto index_site_1 = index_count + 1;

        if (PySpace.Get_Spin_0(index_site_0) != PySpace.Get_Spin_0(index_site_1)){
            // Diagonal energy
            temp_ES_Dia += -MATELE_Heisenberg_DIA;

            // Off diagonal energy
            if (index_count == 0){
                ptr_Lmm1 = &this->Matrix_Id;
            }
            else{
                auto index_site_m1 = index_count - 1;
                ptr_Lmm1 = &Array_L_0[index_site_m1];
            }
            if (index_count == PySpace.size() - 2){
                ptr_Rmp2 = &this->Matrix_Id;
            }
            else{
                auto index_site_2 = index_count + 2;
                ptr_Rmp2 = &Array_R_0[index_site_2];
            }
            auto Val_WSp = arma::trace(((*ptr_Lmm1)
                    * this->Get_Tensor_0(index_site_0).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_0))
                    * this->Get_Tensor_0(index_site_1).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_1))
                    * (*ptr_Rmp2)));
            temp_ES_Off += Val_WSp / Val_WS_0 * MATELE_Heisenbergy_OFF;
        }
        else{
            // Diagonal energy
            temp_ES_Dia += MATELE_Heisenberg_DIA;
        }
    }
    // Then is the pair across the boundary
    {
        MCMPS::type_NumSite index_site_left = 0;
        MCMPS::type_NumSite index_site_right = PySpace.size() - 1;

        if (PySpace.Get_Spin_0(index_site_left) != PySpace.Get_Spin_0(index_site_right)){
            // Diagonal energy
            temp_ES_Dia += -MATELE_Heisenberg_DIA;

            // Off-Diagonal energy
            auto temp_Matrix = this->Get_Tensor_0(index_site_right).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_right));
            for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
                auto index_site = index_site_right - index_count;
                temp_Matrix = this->Get_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site)) * temp_Matrix;
            }
            temp_Matrix = this->Get_Tensor_0(index_site_left).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_left)) * temp_Matrix;

            auto Val_WSp = arma::trace(temp_Matrix);
            temp_ES_Off += Val_WSp / Val_WS_0 * MATELE_Heisenbergy_OFF;
        }
        else{
            // Diagonal energy
            temp_ES_Dia += MATELE_Heisenberg_DIA;
        }
    }

    // Total energy at this configuration
    auto temp_ES = temp_ES_Dia + temp_ES_Off;

    // Increase the number of measurements
    ++Num_Measure;
    // Measure energy
    Val_BinEnergy += temp_ES;
    Val_BinEnergy2 += temp_ES*temp_ES;

    // Measure derivative
    Cal_G_0();
    for (MCMPS::type_NumSite index_site = 0; index_site != PySpace.size() ; ++index_site) {
        auto & Which_D1_Matrix = Get_D1_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site));
        auto & Which_D2_Matrix = Get_D2_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site));
        Which_D1_Matrix += Array_G_0[index_site].t() * (1. / Val_WS_0);
        Which_D2_Matrix += Array_G_0[index_site].t() * (1. / Val_WS_0 * temp_ES);
    }
}
void MCMPS::Class_Network::Measure_Heisenberg_Sym() {
    // Measure the energy at current spin configuration
    MCMPS::type_RealVal temp_ES_Dia = 0.;
    MCMPS::type_RealVal temp_ES_Off = 0.;

    auto Val_WS = Val_WS_0 + Val_WS_Z + Val_WS_R + Val_WS_RZ;

    // Update pairs of spin
    // First for all the pairs except for the last
    // Update pairs of spin
    // First for all the pairs except for the last
    const Class_Matrix* ptr_Lmm1_0 = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp2_0 = &this->Matrix_Id;

    const Class_Matrix* ptr_Lmm1_Z = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp2_Z = &this->Matrix_Id;

    const Class_Matrix* ptr_Rmp1_R = &this->Matrix_Id;
    const Class_Matrix* ptr_Lmm2_R = &this->Matrix_Id;

    const Class_Matrix* ptr_Rmp1_RZ = &this->Matrix_Id;
    const Class_Matrix* ptr_Lmm2_RZ = &this->Matrix_Id;

    for (MCMPS::type_NumSite index_count = 0; index_count != PySpace.size() - 1; ++index_count) {
        auto index_site_0 = index_count;
        auto index_site_1 = index_count + 1;
        if (index_site_0 == 0){
            ptr_Lmm1_0 = &this->Matrix_Id;
            ptr_Lmm1_Z = &this->Matrix_Id;
            ptr_Rmp1_R = &this->Matrix_Id;
            ptr_Rmp1_RZ = &this->Matrix_Id;
        }
        else{
            auto index_site_m1 = index_site_0 - 1;
            ptr_Lmm1_0 = &Array_L_0[index_site_m1];
            ptr_Lmm1_Z = &Array_L_Z[index_site_m1];
            ptr_Rmp1_R = &Array_R_R[PySpace.Get_RSite(index_site_m1)];
            ptr_Rmp1_RZ = &Array_R_RZ[PySpace.Get_RSite(index_site_m1)];
        }

        if (PySpace.Get_Spin_0(index_site_0) != PySpace.Get_Spin_0(index_site_1)){
            // Diagonal energy
            temp_ES_Dia += -MATELE_Heisenberg_DIA;

            if (index_count == PySpace.size() - 2){
                ptr_Rmp2_0 = &this->Matrix_Id;
                ptr_Rmp2_Z = &this->Matrix_Id;
                ptr_Lmm2_R = &this->Matrix_Id;
                ptr_Lmm2_RZ = &this->Matrix_Id;
            }
            else{
                auto index_site_2 = index_count + 2;
                ptr_Rmp2_0 = &Array_R_0[index_site_2];
                ptr_Rmp2_Z = &Array_R_Z[index_site_2];
                ptr_Lmm2_R = &Array_L_R[PySpace.Get_RSite(index_site_2)];
                ptr_Lmm2_RZ = &Array_L_RZ[PySpace.Get_RSite(index_site_2)];
            }

            auto Val_WSp_0 = arma::trace((*ptr_Lmm1_0) *
                                         this->Get_Tensor_0(index_site_0).Get_Matrix(
                                                 PySpace.Get_OppoSpin_0(index_site_0)) *
                                         this->Get_Tensor_0(index_site_1).Get_Matrix(
                                                 PySpace.Get_OppoSpin_0(index_site_1)) *
                                         (*ptr_Rmp2_0));

            auto Val_WSp_Z = arma::trace((*ptr_Lmm1_Z) *
                                         this->Get_Tensor_Z(index_site_0).Get_Matrix(
                                                 PySpace.Get_OppoSpin_Z(index_site_0)) *
                                         this->Get_Tensor_Z(index_site_1).Get_Matrix(
                                                 PySpace.Get_OppoSpin_Z(index_site_1)) *
                                         (*ptr_Rmp2_Z));

            auto Val_WSp_R = arma::trace((*ptr_Lmm2_R) *
                                         this->Get_Tensor_R(PySpace.Get_RSite(index_site_1)).Get_Matrix(
                                                 PySpace.Get_OppoSpin_R(PySpace.Get_RSite(index_site_1))) *
                                         this->Get_Tensor_R(PySpace.Get_RSite(index_site_0)).Get_Matrix(
                                                 PySpace.Get_OppoSpin_R(PySpace.Get_RSite(index_site_0))) *
                                         (*ptr_Rmp1_R));

            auto Val_WSp_RZ = arma::trace((*ptr_Lmm2_RZ) *
                                          this->Get_Tensor_RZ(PySpace.Get_RSite(index_site_1)).Get_Matrix(
                                                  PySpace.Get_OppoSpin_RZ(PySpace.Get_RSite(index_site_1))) *
                                          this->Get_Tensor_RZ(PySpace.Get_RSite(index_site_0)).Get_Matrix(
                                                  PySpace.Get_OppoSpin_RZ(PySpace.Get_RSite(index_site_0))) *
                                          (*ptr_Rmp1_RZ));

            auto Val_WSp = Val_WSp_0 + Val_WSp_Z + Val_WSp_R + Val_WSp_RZ;

            temp_ES_Off += Val_WSp / Val_WS * MATELE_Heisenbergy_OFF;
        }
        else{
            // Diagonal energy
            temp_ES_Dia += MATELE_Heisenberg_DIA;
        }
    }
    // Then is the pair across the boundary
    {
        MCMPS::type_NumSite index_site_left = 0;
        MCMPS::type_NumSite index_site_right = PySpace.size() - 1;

        if (PySpace.Get_Spin_0(index_site_left) != PySpace.Get_Spin_0(index_site_right)){
            // Diagonal energy
            temp_ES_Dia += -MATELE_Heisenberg_DIA;

            // Off-Diagonal energy
            auto temp_Matrix_0 = this->Get_Tensor_0(index_site_right).Get_Matrix(
                    PySpace.Get_OppoSpin_0(index_site_right));
            auto temp_Matrix_Z = this->Get_Tensor_Z(index_site_right).Get_Matrix(
                    PySpace.Get_OppoSpin_Z(index_site_right));
            auto temp_Matrix_R = this->Get_Tensor_R(PySpace.Get_RSite(index_site_right)).Get_Matrix(
                    PySpace.Get_OppoSpin_R(PySpace.Get_RSite(index_site_right)));
            auto temp_Matrix_RZ = this->Get_Tensor_RZ(PySpace.Get_RSite(index_site_right)).Get_Matrix(
                    PySpace.Get_OppoSpin_RZ(PySpace.Get_RSite(index_site_right)));

            for (MCMPS::type_NumSite index_count = 1; index_count != PySpace.size() - 1; ++index_count) {
                auto index_site = index_site_right - index_count;
                temp_Matrix_0 =
                        this->Get_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site)) * temp_Matrix_0;
                temp_Matrix_Z =
                        this->Get_Tensor_Z(index_site).Get_Matrix(PySpace.Get_Spin_Z(index_site)) * temp_Matrix_Z;
                temp_Matrix_R = temp_Matrix_R *
                        this->Get_Tensor_R(PySpace.Get_RSite(index_site)).Get_Matrix(PySpace.Get_Spin_R(PySpace.Get_RSite(index_site)));
                temp_Matrix_RZ = temp_Matrix_RZ * this->Get_Tensor_RZ(PySpace.Get_RSite(index_site)).Get_Matrix(
                        PySpace.Get_Spin_RZ(PySpace.Get_RSite(index_site)));
            }
            temp_Matrix_0 = this->Get_Tensor_0(index_site_left).Get_Matrix(PySpace.Get_OppoSpin_0(index_site_left)) * temp_Matrix_0;
            temp_Matrix_Z = this->Get_Tensor_Z(index_site_left).Get_Matrix(PySpace.Get_OppoSpin_Z(index_site_left)) * temp_Matrix_Z;
            temp_Matrix_R = temp_Matrix_R * this->Get_Tensor_R(PySpace.Get_RSite(index_site_left)).Get_Matrix(
                    PySpace.Get_OppoSpin_R(PySpace.Get_RSite(index_site_left)));
            temp_Matrix_RZ = temp_Matrix_RZ * this->Get_Tensor_RZ(PySpace.Get_RSite(index_site_left)).Get_Matrix(
                    PySpace.Get_OppoSpin_RZ(PySpace.Get_RSite(index_site_left)));

            auto Val_WSp = arma::trace(temp_Matrix_0 + temp_Matrix_Z + temp_Matrix_R + temp_Matrix_RZ);
            temp_ES_Off += Val_WSp / Val_WS * MATELE_Heisenbergy_OFF;
        }
        else{
            // Diagonal energy
            temp_ES_Dia += MATELE_Heisenberg_DIA;
        }
    }

    // Total energy at this configuration
    auto temp_ES = temp_ES_Dia + temp_ES_Off;

    // Increase the number of measurements
    ++Num_Measure;
    // Measure energy
    Val_BinEnergy += temp_ES;
    Val_BinEnergy2 += temp_ES*temp_ES;

    // Measure derivative
    Cal_G_0();
    Cal_G_Z();
    Cal_G_R();
    Cal_G_RZ();
    for (MCMPS::type_NumSite index_site = 0; index_site != PySpace.size() ; ++index_site) {
        Get_D1_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site)) += Array_G_0[index_site].t() * (1. / Val_WS);
        Get_D2_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site)) += Array_G_0[index_site].t() * (1. / Val_WS * temp_ES);

        Get_D1_Tensor_Z(index_site).Get_Matrix(PySpace.Get_Spin_Z(index_site)) += Array_G_Z[index_site].t() * (1. / Val_WS);
        Get_D2_Tensor_Z(index_site).Get_Matrix(PySpace.Get_Spin_Z(index_site)) += Array_G_Z[index_site].t() * (1. / Val_WS * temp_ES);

        Get_D1_Tensor_R(index_site).Get_Matrix(PySpace.Get_Spin_R(index_site)) += Array_G_R[index_site].t() * (1. / Val_WS);
        Get_D2_Tensor_R(index_site).Get_Matrix(PySpace.Get_Spin_R(index_site)) += Array_G_R[index_site].t() * (1. / Val_WS * temp_ES);

        Get_D1_Tensor_RZ(index_site).Get_Matrix(PySpace.Get_Spin_RZ(index_site)) += Array_G_RZ[index_site].t() * (1. / Val_WS);
        Get_D2_Tensor_RZ(index_site).Get_Matrix(PySpace.Get_Spin_RZ(index_site)) += Array_G_RZ[index_site].t() * (1. / Val_WS * temp_ES);
    }
}

void MCMPS::Class_Network::Measure_Ising() {
    // Measure the energy at current spin configuration
    MCMPS::type_RealVal temp_ES_Dia = 0.;
    MCMPS::type_RealVal temp_ES_Off = 0.;

    // Update pairs of spin
    // First for all the pairs except for the last
    const Class_Matrix* ptr_Lmm1 = &this->Matrix_Id;
    const Class_Matrix* ptr_Rmp1 = &this->Matrix_Id;
    for (MCMPS::type_NumSite index_count = 0; index_count != PySpace.size(); ++index_count) {
        auto index_site_0 = index_count;
        auto index_site_1 = (index_count + 1) %  PySpace.size();

        if (PySpace.Get_Spin_0(index_site_0) != PySpace.Get_Spin_0(index_site_1)){
            // Diagonal energy
            temp_ES_Dia += -MATELE_Ising_DIA;
        }
        else{
            // Diagonal energy
            temp_ES_Dia += MATELE_Ising_DIA;
        }

        // Off diagonal energy
        if (index_count == 0){
            ptr_Lmm1 = &this->Matrix_Id;
        }
        else{
            auto index_site_m1 = index_count - 1;
            ptr_Lmm1 = &Array_L_0[index_site_m1];
        }
        if (index_count == PySpace.size() - 1){
            ptr_Rmp1 = &this->Matrix_Id;
        }
        else{
            ptr_Rmp1 = &Array_R_0[index_site_1];
        }
        auto Val_WSp = arma::trace(((*ptr_Lmm1) * this->Get_Tensor_0(index_site_0).Get_Matrix(
                PySpace.Get_OppoSpin_0(index_site_0)) *
                                    (*ptr_Rmp1)));
        temp_ES_Off += Val_WSp / Val_WS_0 * MATELE_Ising_OFF;
    }
    // Total energy at this configuration
    auto temp_ES = temp_ES_Dia + temp_ES_Off;

    // Increase the number of measurements
    ++Num_Measure;
    // Measure energy
    Val_BinEnergy += temp_ES;
    Val_BinEnergy2 += temp_ES*temp_ES;

    // Measure derivative
    Cal_G_0();
    for (MCMPS::type_NumSite index_site = 0; index_site != PySpace.size() ; ++index_site) {
        auto & Which_D1_Matrix = Get_D1_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site));
        auto & Which_D2_Matrix = Get_D2_Tensor_0(index_site).Get_Matrix(PySpace.Get_Spin_0(index_site));
        Which_D1_Matrix += Array_G_0[index_site].t() * (1. / Val_WS_0);
        Which_D2_Matrix += Array_G_0[index_site].t() * (1. / Val_WS_0 * temp_ES);
    }
}
void MCMPS::Class_Network::Update_Tensor(int _Which_Step) {
    // Update Each tensor
    auto Val_delta = delta_k(_Which_Step + 1);
    Tensor_A.Update_Tensor(Tensor_BinD1_A, Tensor_BinD2_A, Val_BinEnergy, Val_delta);
    Tensor_B.Update_Tensor(Tensor_BinD1_B, Tensor_BinD2_B, Val_BinEnergy, Val_delta);


}

MCMPS::type_RealVal MCMPS::Class_Network::delta_k(int _Which_Step) {
//    return DELTA0 * std::pow(_Which_Step, -DELTA_EXPONENT);
    return DELTA0 * std::pow(DELTA_Q, _Which_Step);
}

void MCMPS::Class_Network::CleanMeasure() {
    // Clear the Bin values
    Num_Measure = 0;
    Val_BinEnergy = 0;
    Val_BinEnergy2 = 0;
    Tensor_BinD1_A.Set_Zeros();
    Tensor_BinD1_B.Set_Zeros();
    Tensor_BinD2_A.Set_Zeros();
    Tensor_BinD2_B.Set_Zeros();
}

void MCMPS::Class_Network::AveMeasure() {
    // Firstly calculate the Bin mean value
    Val_BinEnergy *= 1./Num_Measure;
    Val_BinEnergy2 *= 1./Num_Measure;
    Tensor_BinD1_A *= 1./Num_Measure;
    Tensor_BinD1_B *= 1./Num_Measure;
    Tensor_BinD2_A *= 1./Num_Measure;
    Tensor_BinD2_B *= 1./Num_Measure;
}

void MCMPS::Class_Network::Measure_Energy(int _Which_Step) {
    MeasureData.at("E").AppendValue(Val_BinEnergy);
    MeasureData.at("E^2").AppendValue(Val_BinEnergy2);
}

void MCMPS::Class_Network::Write_Energy(int _Which_Step) {
    auto ave_E = MeasureData.at("E").Get_AveValue() / this->PySpace.size();
    auto ave_E2 = MeasureData.at("E^2").Get_AveValue() / this->PySpace.size() / this->PySpace.size();
    auto size = MeasureData.at("E").Get_Size();

    std::cout << _Which_Step << "\t"
              << ave_E << "\t"
              << (ave_E2  - ave_E * ave_E) / size << "\t"
              << std::endl;

    // Outfile
    std::ofstream outfile;
    outfile.open(Str_NameFile + ".txt", std::ofstream::out | std::ofstream::app);
    outfile.setf(std::ios::fixed);
    outfile.precision(16);

    outfile << _Which_Step << "\t"
              << ave_E << "\t"
              << (ave_E2 - ave_E * ave_E) / size << "\t"
              << std::endl;

    outfile.close();

    MeasureData.at("E").ClearValue();
    MeasureData.at("E^2").ClearValue();
}

std::string MCMPS::Class_Network::Get_TypeUpdate() {
    return this->TypeUpdate;
}


