//
// Created by Bowen Zhao (bwzhao@bu.edu) on 10/16/19.
//

#pragma once

#include "Config.h"
#include "Class_MPSTensor.h"
#include "Class_PySpace.h"
#include <vector>
#include "Class_DataMeasurement.h"
#include <cmath>

namespace MCMPS{
    class Class_Network {
    private:
        // Store the physical degrees of freedom
        MCMPS::Class_PySpace PySpace;

        // Store the elements for each identical tensor
        // In this case, we use different tensors for different sublattices
        MCMPS::Class_MPSTensor Tensor_A;
        MCMPS::Class_MPSTensor Tensor_B;

        // Store the pre-calculated matrices
        // Left multiplication
        std::vector<MCMPS::Class_Matrix> Array_L_0;
        // Right multiplication
        std::vector<MCMPS::Class_Matrix> Array_R_0;
        // matrix G(m): omit the m-matrix, and take transpose
        std::vector<MCMPS::Class_Matrix> Array_G_0;

        // Pre-calculated matrices for different sectors:
        // 1. Reflection
        std::vector<MCMPS::Class_Matrix> Array_L_R;
        std::vector<MCMPS::Class_Matrix> Array_R_R;
        std::vector<MCMPS::Class_Matrix> Array_G_R;
        // 2. Z-inversion
        std::vector<MCMPS::Class_Matrix> Array_L_Z;
        std::vector<MCMPS::Class_Matrix> Array_R_Z;
        std::vector<MCMPS::Class_Matrix> Array_G_Z;
        // 3. Z-inversion + Reflection
        std::vector<MCMPS::Class_Matrix> Array_L_RZ;
        std::vector<MCMPS::Class_Matrix> Array_R_RZ;
        std::vector<MCMPS::Class_Matrix> Array_G_RZ;

        // Identity Matrix for convenience
        MCMPS::Class_Matrix Matrix_Id;


        // Store the weight W(S) of the current spin configuration
        MCMPS::type_RealVal Val_WS_0;
        MCMPS::type_RealVal Val_WS_Z;
        MCMPS::type_RealVal Val_WS_R;
        MCMPS::type_RealVal Val_WS_RZ;

        // Distinguish with indiviual run
        int Which_CPU;

        //// For Measurement
        // Count the number of current "local" measurement for updating the tensors
        int Num_Measure;
        MCMPS::Class_MPSTensor Tensor_BinD1_A;
        MCMPS::Class_MPSTensor Tensor_BinD1_B;
        MCMPS::Class_MPSTensor Tensor_BinD2_A;
        MCMPS::Class_MPSTensor Tensor_BinD2_B;

        MCMPS::type_RealVal Val_BinEnergy;
        MCMPS::type_RealVal Val_BinEnergy2;

        // For measurements:
        std::map<std::string, MCMPS::Class_DataMeasurement> MeasureData;

    public:
        Class_Network() = default;
        explicit Class_Network(type_NumSite _NumSite, type_BondDim _BondDim, int _Which_CPU, std::string _TypeUpdate);
        ~Class_Network() = default;

        // Determine if we need symmetric wavefunction
        std::string TypeUpdate;
        std::string Str_NameFile;

        // Calculate the pre-calculated matrices
        void Cal_L_0();
        void Cal_R_0();
        void Cal_G_0();
        void Cal_L_R();
        void Cal_R_R();
        void Cal_G_R();
        void Cal_L_Z();
        void Cal_R_Z();
        void Cal_G_Z();
        void Cal_L_RZ();
        void Cal_R_RZ();
        void Cal_G_RZ();

        // Get the tensor at certain state
        MCMPS::Class_MPSTensor &Get_Tensor_0(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D1_Tensor_0(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D2_Tensor_0(type_NumSite Which_Site);
        // Get the tensor of certain quantum number:
        MCMPS::Class_MPSTensor &Get_Tensor_R(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D1_Tensor_R(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D2_Tensor_R(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_Tensor_Z(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D1_Tensor_Z(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D2_Tensor_Z(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_Tensor_RZ(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D1_Tensor_RZ(type_NumSite Which_Site);
        MCMPS::Class_MPSTensor &Get_D2_Tensor_RZ(type_NumSite Which_Site);

        // Single update step for the spin configurations
        // This update will update two spins at a time. And it will conserve the total spin
        void Update_Spin_Pair();
        void Update_Spin_Pair_Sym();
        // This update will just update one single spin
        void Update_Spin_Single();

        // Single Measurement step
        void Measure_Heisenberg();
        void Measure_Heisenberg_Sym();
        void Measure_Ising();

        void AveMeasure();
        void CleanMeasure();

        // After several measurement, we can have one update of the tensor
        void Update_Tensor(int _Which_Step);

        void Measure_Energy(int _Which_Step);
        void Write_Energy(int _Which_Step);
        inline type_RealVal delta_k(int _Which_Step);

        std::string Get_TypeUpdate();
    };
}



