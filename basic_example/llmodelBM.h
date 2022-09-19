#pragma once
#include <random>

class LL_model_BM{
    public:
        int n, cycle;
        double ***spins, **P2;
        bool *broken;
        int* neighbors_list;
        double E;
        double beta, rotation_angle, acceptance_rate;
        std::default_random_engine generator;
        std::uniform_int_distribution<int> random_I;
        std::uniform_real_distribution<double> p;

        LL_model_BM(int n);

        void neighbors(int I);
        void set_beta(double beta);
        double dot(int I, int i, int J, int j);
        double E_IJ(int I, int J);
        double E_neighbors(int I);
        void H();
        void calculate_P2();

        bool BarkerWatts_move();
        void adjust_rotation_angle();

        bool switch_move();
        void switch_sites(int I, int J);

        void initialize_lattice_parallel();
        void initialize_lattice_random();

        void thermalize_BarkerWatts(int N);

        void break_molecule(int I);
        void align_molecule(int I);

        void break_molecules(int count);

        

};