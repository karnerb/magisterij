#include <random>

class LL_model{
    public:
        int n, cycle;
        int* neighbors_list;
        double **spins, **P2;
        double E;
        double polar_order;
        double beta, rotation_angle, acceptance_rate;
        std::default_random_engine generator;
        std::uniform_int_distribution<int> random_I;
        std::uniform_real_distribution<double> p;
        
        
        LL_model(int n);

        void neighbors(int I);

        void set_beta(double beta);
        double dot(int I, int J);
        double E_IJ(int I, int J);
        double E_neighbors(int I);
        void H();
        void calculate_P2();

        bool BarkerWatts_move();
        void adjust_rotation_angle();

        void cluster_move();
        
        void initialize_lattice_parallel();
        void initialize_lattice_random();

        void thermalize_BarkerWatts(int N);
        void thermalize_cluster(int N);

        

};