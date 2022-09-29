#pragma once
#include <random>

class LL_model_BM{
    public:
        int n, cycle;
        int* shuffled_I;
        double ***spins, **P2;
        bool *broken;
        int* neighbors_list;
        double E;
        double polar_order;
        double beta, rotation_angle, acceptance_rate;
        std::default_random_engine generator;
        std::uniform_int_distribution<int> random_I;
        std::uniform_real_distribution<double> p;

        LL_model_BM(int n);
        
        void calculate_polar_order();
        void neighbors(int I);
        void set_beta(double beta);
        double dot(int I, int i, int J, int j);
        double E_IJ(int I, int J);
        double E_neighbors(int I);
        void H();
        void calculate_P2();

        void shuffle_I();

        bool BarkerWatts_move(int I);
        void BarkerWatts_cycle();
        void adjust_rotation_angle();

        bool switch_move();
        void switch_sites(int I, int J);

        void initialize_lattice_parallel();
        void initialize_lattice_random();

        void thermalize_BarkerWatts(int N);

        void break_molecule(int I);
        void align_molecule(int I);

        void break_molecules(int count);

        int count_broken_molecules();

};