//
// Created by monika on 15.10.19.
//

#include <iostream>
#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include <iomanip>

using namespace std;
using namespace arma;


inline double Pott2DEnergy(mat conf, double J) {
    int L = conf.n_rows;
    int neig0, neig1, neig2, neig3, main;
    double energy = 0;
    for(int i=0; i<L; i++) {
        for(int j=0; j<L; j++) {
            main = conf(i, j);
            neig0 = conf(i, (j+1)%L);
            neig1 = conf((i+L-1)%L, j);
            neig2 = conf(i, (j+L-1)%L);
            neig3 = conf((i+1)%L, j);
            energy += bool(main != neig0);
            energy += bool(main != neig1);
            energy += bool(main != neig2);
            energy += bool(main != neig3);


        }
    }
    return J*energy/2.0;
}

inline double magn(mat conf, int q) {
    int L = conf.n_rows;
    int size = L*L;
    double magn = 0;
    double angle;
    double angle_x = 0;
    double angle_y = 0;
    for(int i=0; i<L; i++){
        for(int j=0; j<L; j++) {
            angle = (double) conf(i, j)*2*M_PI/q;
            angle_x += cos(angle);
            angle_y += sin(angle);

        }
    }
    return sqrt(angle_x*angle_x + angle_y*angle_y)/size;

}

inline mat initial_conf(int L, int q){

    std::random_device rd; // obtain a random number from hardware
    std::mt19937 eng(rd()); // seed the generator
    std::uniform_int_distribution<> distr(0, q-1); // define the range

    mat initial(L, L);
    initial.zeros();

    for (int i=0; i<L; i++) {
        for(int j=0; j<L; j++) {
            initial(i, j) = distr(eng);
        }
    }

    return initial;

}

inline double beta_critical(int q) {
    return 0.5*log(1+sqrt(q));
}

inline mat WolffUpdate(mat conf, int q, double beta, double J, int N) {
    int L = conf.n_rows;
    random_device rd;
    mt19937 eng(rd());
    uniform_int_distribution<> distr(0, L-1);
    uniform_int_distribution<> distr2(0, q-1);

    vec sites_x(L*L);
    vec sites_y(L*L);

    vec neig_sites_x(4);
    vec neig_sites_y(4);

    cube neigb (4, L, L);
    mat already_cluster(L, L);


    int ss;
    int x_site, y_site;
    double r, p;





    // Number of iterations
    for(int it=0; it<N; it++) {

        for(int i=0; i<L; i++) {
            for(int j=0; j<L; j++) {
                neigb(0, i, j) = conf(i, (j+1)%L);
                neigb(1, i, j) = conf((i+L-1)%L, j);
                neigb(2, i, j) = conf(i, (j+L-1)%L);
                neigb(3, i, j) = conf((i+1)%L, j);
            }
        }

//        cout << it << endl;

        sites_x.zeros();
        sites_y.zeros();

        sites_x(0) = distr(eng);
        sites_y(0) = distr(eng);

        ss = conf(sites_x(0), sites_y(0));


        already_cluster.zeros();
        already_cluster(sites_x(0), sites_y(0)) = 1;


        int k = 1; // number of sites in the cluster
        int limit0 = 0;
        int limit1 = k;

        while (limit0 < limit1) {

            for (int i = limit0; i < limit1; i++) {
                x_site = (int) sites_x(i);
                y_site = (int) sites_y(i);

                neig_sites_x(0) = x_site;
                neig_sites_x(1) = (x_site + L - 1) % L;
                neig_sites_x(2) = x_site;
                neig_sites_x(3) = (x_site + 1) % L;

                neig_sites_y(0) = (y_site + 1) % L;
                neig_sites_y(1) = y_site;
                neig_sites_y(2) = (y_site - 1 + L) % L;
                neig_sites_y(3) = y_site;



                for (int j = 0; j < 4; j++) {
                    if (neigb(j, sites_x(i), sites_y(i)) == ss &&
                        already_cluster(neig_sites_x(j), neig_sites_y(j)) != 1.0) {
                        r = ((double) rand() / (RAND_MAX));
                        p = 1-exp(-2 * beta);
                        if (r < p) {
                            sites_x(k) = neig_sites_x(j);
                            sites_y(k) = neig_sites_y(j);
                            already_cluster(neig_sites_x(j), neig_sites_y(j)) = 1;
                            k++;

                        }

                    }
                }

            }


            limit0 = limit1;
            limit1 = k;

        }





//        sites_x = sites_x.head(k);
//        sites_y = sites_y.head(k);


        int new_s = ss;


        while (ss == new_s) {
            new_s = distr2(eng);
        }

        for (int n = 0; n < k; n++) {
            conf(sites_x(n), sites_y(n)) = new_s;
        }

//        conf.print();

    }
    return conf;

}

inline double heat_capacity(double beta, vec energy) {
    int size = energy.n_rows;
    double suma = accu(energy);
    double suma_kw = accu(energy%energy);

    return beta*beta*(suma_kw - suma*suma);
}

inline double susceptibility(double beta, vec magn) {
    int size = magn.n_rows;
    double suma = accu(magn);
    double suma_kw = accu(magn%magn);

    return beta*(suma_kw - suma*suma);
}


inline void saveintofile(mat config, ofstream &file) {
    int L = config.n_rows;

    for (int i = 0; i < L; i++) {
        for (int j = 0; j < L; j++) {
            file << setprecision(1) << config(i, j) << " ";
        }
    }
    file << "\n";
}

