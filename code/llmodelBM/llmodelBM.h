#pragma once
#include <random>
#define pi 3.141592653589793238462


class LL_model_BM{
    public:
        
        int *shuffled_I, *neighbors_list;
        double ***spins, **P2;
        bool *broken;
        int n, cycle, broken_neighbors;
        double E, polar_order, beta, rotation_angle, acceptance_rate, swap_acceptance_rate;
        
        std::default_random_engine generator;
        std::uniform_int_distribution<int> random_I;
        std::uniform_real_distribution<double> p;
        std::uniform_int_distribution<int> random_neighbor;
        
        LL_model_BM(int n);
        
        //output
        std::string metadata();
        std::string names();
        friend std::ostream& operator <<(std::ostream& os, const LL_model_BM& model);
        
        //basic functions
        void set_beta(double beta);
        void neighbors(int I);
        double dot(int I, int i, int J, int j);
        void rotation3D(double* vector, double n0, double n1, double n2, double gama);

        //energy, P2 calculation
        double E_IJ(int I, int J);
        double E_neighbors(int I);
        void H();
        void calculate_P2();
        void count_neighbors_of_same_kind();
        

        //orientational MC
        void shuffle_I();
        bool BarkerWatts_move(int I);
        void BarkerWatts_cycle();
        void thermalize_BarkerWatts(int N);
        void adjust_rotation_angle();

        //swap MC
        bool swap_move();
        void swap_sites(int I, int J);
        void swap_cycle();

        //lattice initialization 
        void initialize_lattice_parallel();
        void initialize_lattice_random();
        void break_molecule(int I);
        void align_molecule(int I);
        void break_molecules(int count);

        //additional order parameters
        void calculate_polar_order();
        int count_broken_molecules();
        
};