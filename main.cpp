#include <iostream>
#include <armadillo>
#include <math.h>
#include <cstdlib>
#include <ctime>
#include "some_functions.cpp"

using namespace std;
using namespace arma;


int main() {

    // 1.  MONTE CARLO

    // a) variables

    int L = 50;
    int t = 1;
    int n = 1; // the multiple of space between the independent configurations
    double J = 1.0;
    double T = 0.1;
    double beta = 1.0 / T;
    double E_new;
    int q = 4;

    // files

    ostringstream ss1;

    ss1 << L;
    string Size = ss1.str();


    string conf_name, hist_energy_name, hist_magn_name, capacity_name, suscpet_name;

    // b) initial random configuration

//    mat configurations;
//    mat new_conf;
//    configurations.load("L=20q=5_T=1.61.txt");
//    new_conf = configurations.row(configurations.n_rows-1);
//    new_conf.reshape(20, 20);


    mat new_conf = initial_conf(L, q);
//    new_conf.print();

    // c) MC

    int N = 10000; // n_configurations
    int step = 50; // n_updates
    int therm_step = 1000;
    double Temp;
    int counter;
    vec energies(N);
    vec magnets(N);
//    vec heat_capacities(140);
//    vec suscpet(140);

    suscpet_name = "L=" + Size + "_suscp.txt";
    capacity_name = "L=" + Size  + "_heat.txt";



    int n_temp = 0;
    for (int i=300; i>0; i-=1) {
        Temp = 0.01*i;
        beta = 1.0/Temp;

        cout << "Temp " << Temp << endl;

        ostringstream ss2;

        ss2 << Temp;
        string Temp2 = ss2.str();

        ostringstream ss3;

        ss3 << q;
        string q_val = ss3.str();

        conf_name = "L=" + Size + "q=" + q_val + "_T=" + Temp2 + ".txt";
        hist_energy_name = "L=" + Size +  "q=" + q_val + "_T=" + Temp2 + "_energy.txt";
        hist_magn_name = "L=" + Size + "q=" + q_val + "_T=" + Temp2 + "_magnet.txt";

        ofstream myfile;
        myfile.open(conf_name.c_str(), ios_base::app);


        //thermalisation
//        cout << "thermalisation " << endl;
        new_conf  = WolffUpdate(new_conf, q, beta, J, therm_step);
        counter = 0;


        while(counter<N) {
            new_conf = WolffUpdate(new_conf, q, beta, J, step);
            magnets(counter) = magn(new_conf, q);
            energies(counter) = Pott2DEnergy(new_conf, J);
            saveintofile(new_conf, myfile);
            counter++;
        }

//        new_conf.print();


        energies.save(hist_energy_name, csv_ascii);
        magnets.save(hist_magn_name, csv_ascii);

//        heat_capacities(n_temp) = heat_capacity(beta, energies);
//        suscpet(n_temp) =  susceptibility(beta, magnets);

        n_temp++;

    }

//    heat_capacities.save(capacity_name, csv_ascii);
//    suscpet.save(suscpet_name, csv_ascii);


return 0;


}