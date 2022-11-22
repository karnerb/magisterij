#pragma once
#include <random>
#define pi 3.141592653589793238462



class LC_Elastomer{
    public:
        
        int **neighbors_list;
        int *shuffled_I;
        double ***spins, **P2;
        bool *broken;
        int n, cycle, broken_neighbors, broken_cluster_count;
        double Energy_LL_interaction, Energy_elastic, Energy_coupling, polar_order, beta, rotation_angle, resize_step, acceptance_rate, resize_acceptance_rate;
        double alpha;
        double lambda;
        double sigma;
        
        std::default_random_engine generator;
        std::uniform_int_distribution<int> random_I;
        std::uniform_real_distribution<double> p;
        std::uniform_int_distribution<int> random_neighbor;
        
        LC_Elastomer(int n);
        
        //output
        std::string metadata();
        std::string names();
        friend std::ostream& operator <<(std::ostream& os, const LC_Elastomer& model);
        
        //basic functions
        void set_beta(double beta);
        void set_lambda(double lambda);
        double dot(int I, int i, int J, int j);
        void rotation3D(double* vector, double n0, double n1, double n2, double gama);

        //energy, P2 calculation
        double E_IJ(int I, int J);
        double E_neighbors(int I);
        double E_coupling(int I);
        double E_elastic();

        void H_neighbors();
        void H_coupling();
        void H_elastic();
        void calculate_P2();
        void count_neighbors_of_same_kind();
        
        //coupling parameter
        double Q(double lambda);

        //orientational MC
        void shuffle_I();
        bool BarkerWatts_move(int I);
        void BarkerWatts_cycle();
        void thermalize_BarkerWatts(int N);
        void adjust_rotation_angle();
        void BarkerWatts_cycle2();
        void BarkerWatts_cycle3();

        //swap MC
        bool swap_move();
        void swap_sites(int I, int J);
        void swap_cycle();

        //resize MC
        bool resize_move();
        void adjust_resize_step();

        //lattice initialization 
        void initialize_lattice_parallel();
        void initialize_lattice_random();
        void break_molecule(int I);
        void align_molecule(int I);
        void break_molecules(int count);
        void make_one_cluster(int a);
        void make_two_clusters(int a);
        void checkerboard_pattern();
        void make_empty_cube_cluster(int a);
        void make_L_shaped_cluster(int a, int b);

        //additional order parameters
        void calculate_polar_order();
        void cluster_count_and_size();
        int count_broken_molecules();
        
};
