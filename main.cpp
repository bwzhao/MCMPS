#include <iostream>
#include "Config.h"
#include "Class_Network.h"

int main(int argc, char *argv[]) {
    std::vector<std::string> para_main;
    for (decltype(argc) nl_arg = 0; nl_arg != argc; ++nl_arg) {
        std::string temp_string(argv[nl_arg]);
        para_main.push_back(temp_string);
    }

    MCMPS::type_NumSite Numsite = std::stoi(para_main[1]);
    MCMPS::type_BondDim BondDim = std::stoi(para_main[2]);
    int Which_CPU = std::stoi(para_main[3]);
    // "Simple", "Neel", "Sym"
    std::string TypeUpdate = para_main[4];
    MCMPS::type_RealVal Val_J1 = std::stod(para_main[5]);
    MCMPS::type_RealVal Val_J2 = std::stod(para_main[6]);

    MCMPS::type_RealVal Val_Sign_Z = std::stod(para_main[7]);
    MCMPS::type_RealVal Val_Sign_R = std::stod(para_main[8]);

    // Construct the MPS
    auto MPS = MCMPS::Class_Network(Numsite, BondDim, Which_CPU, TypeUpdate, Val_J1, Val_J2, Val_Sign_Z, Val_Sign_R);

    for (int index_loop = 0; index_loop != 1; ++index_loop){
        for (int index_Step_k = 1; index_Step_k != MCMPS::K_MAX; ++index_Step_k) {
            for (int index_Bin = 0; index_Bin != MCMPS::P_k(index_Step_k); ++index_Bin){
                for (int index_Measure = 0; index_Measure != MCMPS::F_k(index_Step_k); ++index_Measure) {
                    if (MPS.Get_TypeUpdate() == "Simple"){
                        MPS.Update_Spin_Pair();
                        MPS.Measure_Heisenberg();
                    }
                    else if (MPS.Get_TypeUpdate() == "Sym"){
                        MPS.Update_Spin_Pair_Sym(Val_Sign_Z, Val_Sign_R);
                        MPS.Measure_Heisenberg_Sym();
                    }
                    else if (MPS.Get_TypeUpdate() == "J1J2"){
                        MPS.Update_Spin_Pair_Sym(Val_Sign_Z, Val_Sign_R);
                        MPS.Measure_J1J2_Sym(Val_Sign_Z, Val_Sign_R);
                    }
                    else if (MPS.Get_TypeUpdate() == "Ising") {
                        MPS.Update_Spin_Single();
                        MPS.Measure_Ising();
                    }
                }
                MPS.AveMeasure();
                MPS.Update_Tensor(index_Step_k);
                MPS.CleanMeasure();
            }

            for (int index_Bin = 0; index_Bin != MCMPS::P_k(index_Step_k); ++index_Bin){
                for (int index_Measure = 0; index_Measure != MCMPS::F_k(index_Step_k); ++index_Measure) {
                    if (MPS.Get_TypeUpdate() == "Simple"){
                        MPS.Update_Spin_Pair();
                        MPS.Measure_Heisenberg();
                    }
                    else if (MPS.Get_TypeUpdate() == "Sym"){
                        MPS.Update_Spin_Pair_Sym(Val_Sign_Z, Val_Sign_R);
                        MPS.Measure_Heisenberg_Sym();
                    }
                    else if (MPS.Get_TypeUpdate() == "J1J2"){
                        MPS.Update_Spin_Pair_Sym(Val_Sign_Z, Val_Sign_R);
                        MPS.Measure_J1J2_Sym(Val_Sign_Z, Val_Sign_R);
                    }
                    else if (MPS.Get_TypeUpdate() == "Ising") {
                        MPS.Update_Spin_Single();
                        MPS.Measure_Ising();
                    }
                }
                MPS.AveMeasure();
                MPS.Measure_Energy(index_Step_k);
                MPS.CleanMeasure();
            }
            MPS.Write_Energy(index_Step_k);
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