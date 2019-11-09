#include <iostream>
#include "Config.h"
#include "Class_Network.h"

int main() {
    MCMPS::type_NumSite Numsite = 12;
    MCMPS::type_BondDim BondDim = 12;
    int Which_CPU = 1;

    // Construct the MPS
    auto MPS = MCMPS::Class_Network(Numsite, BondDim, Which_CPU);

    for (int index_loop = 0; index_loop != 1; ++index_loop){
        for (int index_Step_k = 0; index_Step_k != MCMPS::K_MAX; ++index_Step_k) {
            for (int index_Bin = 0; index_Bin != MCMPS::P_k(index_Step_k); ++index_Bin){
                for (int index_Measure = 0; index_Measure != MCMPS::F_k(index_Step_k); ++index_Measure) {
//                    MPS.Update_Spin_Pair();
//                    MPS.Measure_Heisenberg();
                    MPS.Update_Spin_Pair_Sym();
                    MPS.Measure_Heisenberg_Sym();
//                    MPS.Update_Spin_Single();
//                    MPS.Measure_Ising();
                }
                MPS.AveMeasure();
                MPS.Update_Tensor(index_Step_k);
                MPS.CleanMeasure();
            }
        }
    }
//
//    for (int index_loop = 0; index_loop != 3; ++index_loop){
//        for (int index_Step_k = 0; index_Step_k != MCMPS::K_MAX; ++index_Step_k) {
//            for (int index_Bin = 0; index_Bin != MCMPS::P_k(index_Step_k); ++index_Bin){
//                for (int index_Measure = 0; index_Measure != MCMPS::F_k(index_Step_k); ++index_Measure) {
//                    MPS.Update_Spin_Pair_Sym();
//                    MPS.Measure_Heisenberg_Sym();
//                }
//                MPS.AveMeasure();
//                MPS.Update_Tensor(index_Step_k);
//                MPS.CleanMeasure();
//            }
//        }
//    }





    return 0;
}