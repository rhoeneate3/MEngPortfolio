//-------------------------------------------------------------------//
//         eduPIC : educational 1d3v PIC/MCC simulation code         //
//           version 1.0, release date: March 16, 2021               //
//                       :) Share & enjoy :)                         //
//-------------------------------------------------------------------//
// When you use this code, you are required to acknowledge the       //
// authors by citing the paper:                                      //
// Z. Donko, A. Derzsi, M. Vass, B. Horvath, S. Wilczek              //
// B. Hartmann, P. Hartmann:                                         //
// "eduPIC: an introductory particle based  code for radio-frequency //
// plasma simulation"                                                //
// Plasma Sources Science and Technology, vol 30, pp. 095017 (2021)  //
//-------------------------------------------------------------------//
// Disclaimer: The eduPIC (educational Particle-in-Cell/Monte Carlo  //
// Collisions simulation code), Copyright (C) 2021                   //
// Zoltan Donko et al. is free software: you can redistribute it     //
// and/or modify it under the terms of the GNU General Public License//
// as published by the Free Software Foundation, version 3.          //
// This program is distributed in the hope that it will be useful,   //
// but WITHOUT ANY WARRANTY; without even the implied warranty of    //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU  //
// General Public License for more details at                        //
// https://www.gnu.org/licenses/gpl-3.0.html.                        //
//-------------------------------------------------------------------//

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdbool>
#include <cmath>
#include <ctime>
#include <random>
#include <time.h>
#include <omp.h>
#include <string.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>     // for access() and chdir()
#include <sys/stat.h>   // for mkdir()
#include <omp.h>
#include <time.h>
#include <windows.h>  // for __declspec(dllexport)
#include <direct.h>   // for _mkdir on Windows
#include <io.h>       // for access()



using namespace::std;

#define MAX_THREADS 128
// constants

double time1, timedif1,timedif2, timedifHR, updted_time;

const double     PI             = 3.141592653589793;          // mathematical constant Pi
const double     TWO_PI         = 2.0 * PI;                   // two times Pi
const double     E_CHARGE       = 1.60217662e-19;             // electron charge [C]
const double     EV_TO_J        = E_CHARGE;                   // eV <-> Joule conversion factor
const double     E_MASS         = 9.10938356e-31;             // mass of electron [kg]
const double     AR_MASS        = 6.63352090e-26;             // mass of argon atom [kg]
const double     MU_ARAR        = AR_MASS / 2.0;              // reduced mass of two argon atoms [kg]
const double     K_BOLTZMANN    = 1.38064852e-23;             // Boltzmann's constant [J/K]
const double     EPSILON0       = 8.85418781e-12;             // permittivity of free space [F/m]
const double     MU0            = 4.0 * PI *1e-7;                // premeabiliy of free space [H/m]
const int     c_light              = 3e8 ;
// simulation parameters
const int        N_G            = 1000;                        // number of grid points
const int        N_T            = 4000;                       // time steps within an RF period
const double     FREQUENCY      = 13.56e6;                    // driving frequency [Hz]
const double     WEIGHT         = 7e4;                      // weight of superparticles
const int N_INIT = 1000; 
const double electrode_area = 1.0e-4;  // JUST A SCALING FACTOR!

double L = 0.025; // electrode gap [m]
double VOLTAGE = 250; // voltage amplitude [V]
double PRESSURE = 10;  // gas pressure [Pa]
double TEMPERATURE = 350.0; // background gas temperature [K]
double     N_TURNS = 5.0;                       //no. turns of solenoid
double POWER = 0.0;
double TARGET_EFFICIENCY = 60;

// Global variable to store heating efficiency
double heating_efficiency = 0.0;  // Will hold the value of heating efficiency
double R_coil = 10 ;
double B_STRENGTH  = MU0 * N_TURNS * (VOLTAGE/(R_coil)) ;                   // magnetic field strength T
    
double MIN_POWER = 0.0, MAX_POWER = 0.0;
double MIN_PRESSURE = 0.0, MAX_PRESSURE = 0.0;
double MIN_L = 0.0, MAX_L = 0.0;
double MIN_TEMPERATURE = 0.0, MAX_TEMPERATURE = 0.0;

double MIN_VOLTAGE;
double MAX_VOLTAGE;
double MIN_GAP;
double MAX_GAP;

// Define a struct to hold all the information returned from check_and_save_info()
typedef struct {
    double ecoll_freq;
    double icoll_freq;
    double density;
    double plas_freq;
    double meane;
    double debye_length;
    double power_e; 
    double power_i;
    double sim_time; 
    double R_plasma;
    double heating_efficiency;
    char param_name[32];  
} Info;

// for genetic algorithm bit:
// #define POP_SIZE 30  // Population size
// #define GENERATIONS 30  // Number of generations
// #define MUTATION_RATE 0.35  // 10% mutation chance

// for genetic algorithm bit:
#define POP_SIZE 10  // Population size
#define GENERATIONS 20  // Number of generations
#define MUTATION_RATE 0.1  // 10% mutation chance



// additional (derived) constants

const double     PERIOD         = 1.0 / FREQUENCY;                           // RF period length [s]
const double     DT_E           = PERIOD / (double)(N_T);                    // electron time step [s]
const int        N_SUB          = 20;                                        // ions move only in these cycles (subcycling)
const double     DT_I           = N_SUB * DT_E;                              // ion time step [s]
const double     DX             = L / (double)(N_G - 1);                     // spatial grid division [m]
const double     INV_DX         = 1.0 / DX;                                  // inverse of spatial grid size [1/m]
const double     GAS_DENSITY    = PRESSURE / (K_BOLTZMANN * TEMPERATURE);    // background gas density [1/m^3]
const double     OMEGA          = TWO_PI * FREQUENCY;                        // angular frequency [rad/s]

// electron and ion cross sections

const int        N_CS           = 5;                          // total number of processes / cross sections
const int        E_ELA          = 0;                          // process identifier: electron/elastic
const int        E_EXC          = 1;                          // process identifier: electron/excitation
const int        E_ION          = 2;                          // process identifier: electron/ionization
const int        I_ISO          = 3;                          // process identifier: ion/elastic/isotropic
const int        I_BACK         = 4;                          // process identifier: ion/elastic/backscattering
const double     E_EXC_TH       = 11.5;                       // electron impact excitation threshold [eV]
const double     E_ION_TH       = 15.8;                       // electron impact ionization threshold [eV]
const int        CS_RANGES      = 1000000;                    // number of entries in cross section arrays
const double     DE_CS          = 0.001;                      // energy division in cross section arrays [eV]
typedef float    cross_section[CS_RANGES];                    // cross section array
cross_section    sigma[N_CS];                                 // set of cross section arrays
cross_section    sigma_tot_e;                                 // total macroscopic cross section of electrons
cross_section    sigma_tot_i;                                 // total macroscopic cross section of ions

// particle coordinates

const int        MAX_N_P = 1000000;                           // maximum number of particles (electrons / ions)
typedef double   particle_vector[MAX_N_P];                    // array for particle properties
int              N_e = 0;                                     // number of electrons
int              N_i = 0;                                     // number of ions
particle_vector  x_e, vx_e, vy_e, vz_e;                       // coordinates of electrons (one spatial, three velocity components)
particle_vector  x_i, vx_i, vy_i, vz_i;                       // coordinates of ions (one spatial, three velocity components)

typedef double   xvector[N_G];                                // array for quantities defined at gird points
xvector          efield,pot;                                  // electric field and potential
xvector          e_density,i_density;                         // electron and ion densities
xvector          cumul_e_density,cumul_i_density;             // cumulative densities

typedef unsigned long long int Ullong;                        // compact name for 64 bit unsigned integer
Ullong       N_e_abs_pow  = 0;                                // counter for electrons absorbed at the powered electrode
Ullong       N_e_abs_gnd  = 0;                                // counter for electrons absorbed at the grounded electrode
Ullong       N_i_abs_pow  = 0;                                // counter for ions absorbed at the powered electrode
Ullong       N_i_abs_gnd  = 0;                                // counter for ions absorbed at the grounded electrode

// electron energy probability function

const int    N_EEPF  = 2000;                                 // number of energy bins in Electron Energy Probability Function (EEPF)
const double DE_EEPF = 0.05;                                 // resolution of EEPF [eV]
typedef double eepf_vector[N_EEPF];                          // array for EEPF
eepf_vector eepf     = {0.0};                                // time integrated EEPF in the center of the plasma

// ion flux-energy distributions

const int    N_IFED   = 200;                                 // number of energy bins in Ion Flux-Energy Distributions (IFEDs)
const double DE_IFED  = 1.0;                                 // resolution of IFEDs [eV]
typedef int  ifed_vector[N_IFED];                            // array for IFEDs
ifed_vector  ifed_pow = {0};                                 // IFED at the powered electrode
ifed_vector  ifed_gnd = {0};                                 // IFED at the grounded electrode
double       mean_i_energy_pow;                              // mean ion energy at the powered electrode
double       mean_i_energy_gnd;                              // mean ion energy at the grounded electrode

// spatio-temporal (XT) distributions

const int N_BIN                     = 20;                    // number of time steps binned for the XT distributions
const int N_XT                      = N_T / N_BIN;           // number of spatial bins for the XT distributions
typedef double xt_distr[N_G][N_XT];                          // array for XT distributions (decimal numbers)
xt_distr pot_xt                     = {0.0};                 // XT distribution of the potential
xt_distr efield_xt                  = {0.0};                 // XT distribution of the electric field
xt_distr ne_xt                      = {0.0};                 // XT distribution of the electron density
xt_distr ni_xt                      = {0.0};                 // XT distribution of the ion density
xt_distr ue_xt                      = {0.0};                 // XT distribution of the mean electron velocity
xt_distr ui_xt                      = {0.0};                 // XT distribution of the mean ion velocity
xt_distr je_xt                      = {0.0};                 // XT distribution of the electron current density
xt_distr ji_xt                      = {0.0};                 // XT distribution of the ion current density
xt_distr powere_xt                  = {0.0};                 // XT distribution of the electron powering (power absorption) rate
xt_distr poweri_xt                  = {0.0};                 // XT distribution of the ion powering (power absorption) rate
xt_distr meanee_xt                  = {0.0};                 // XT distribution of the mean electron energy
xt_distr meanei_xt                  = {0.0};                 // XT distribution of the mean ion energy
xt_distr counter_e_xt               = {0.0};                 // XT counter for electron properties
xt_distr counter_i_xt               = {0.0};                 // XT counter for ion properties
xt_distr ioniz_rate_xt              = {0.0};                 // XT distribution of the ionisation rate

double   mean_energy_accu_center    = 0;                     // mean electron energy accumulator in the center of the gap
Ullong   mean_energy_counter_center = 0;                     // mean electron energy counter in the center of the gap
Ullong   N_e_coll                   = 0;                     // counter for electron collisions
Ullong   N_i_coll                   = 0;                     // counter for ion collisions
double   Time;                                               // total simulated time (from the beginning of the simulation)
int      cycle,no_of_cycles,cycles_done;                     // current cycle and total cycles in the run, cycles completed
int      arg1;                                               // used for reading command line arguments
char     st0[80];                                            // used for reading command line arguments
FILE     *datafile;                                          // used for saving data
bool     measurement_mode;                                   // flag that controls measurements and data saving

//---------------------------------------------------------------------------//
// C++ Mersenne Twister 19937 generator                                      //
// R01(MTgen) will genarate uniform distribution over [0,1) interval         //
// RMB(MTgen) will generate Maxwell-Boltzmann distribution (of gas atoms)    //
//---------------------------------------------------------------------------//

std::random_device rd{}; 
std::mt19937 MTgen(rd());
std::uniform_real_distribution<> R01(0.0, 1.0);
std::normal_distribution<> RMB(0.0,sqrt(K_BOLTZMANN * TEMPERATURE / AR_MASS));

//----------------------------------------------------------------------------//
//  electron cross sections: A V Phelps & Z Lj Petrovic, PSST 8 R21 (1999)    //
//----------------------------------------------------------------------------//

void set_electron_cross_sections_ar(void){
    int    i;
    double en,qmel,qexc,qion;
    
    printf(">> eduPIC: Setting e- / Ar cross sections\n");
    for(i=0; i<CS_RANGES; i++){
        if (i == 0) {en = DE_CS;} else {en = DE_CS * i;}                            // electron energy
        qmel = fabs(6.0 / pow(1.0 + (en/0.1) + pow(en/0.6,2.0), 3.3)
                    - 1.1 * pow(en, 1.4) / (1.0 + pow(en/15.0, 1.2)) / sqrt(1.0 + pow(en/5.5, 2.5) + pow(en/60.0, 4.1)))
        + 0.05 / pow(1.0 + en/10.0, 2.0) + 0.01 * pow(en, 3.0) / (1.0 + pow(en/12.0, 6.0));
        if (en > E_EXC_TH)
            qexc = 0.034 * pow(en-11.5, 1.1) * (1.0 + pow(en/15.0, 2.8)) / (1.0 + pow(en/23.0, 5.5))
            + 0.023 * (en-11.5) / pow(1.0 + en/80.0, 1.9);
        else
            qexc = 0;
        if (en > E_ION_TH)
            qion = 970.0 * (en-15.8) / pow(70.0 + en, 2.0) + 0.06 * pow(en-15.8, 2.0) * exp(-en/9);
        else
            qion = 0;
        sigma[E_ELA][i] = qmel * 1.0e-20;       // cross section for e- / Ar elastic collision
        sigma[E_EXC][i] = qexc * 1.0e-20;       // cross section for e- / Ar excitation
        sigma[E_ION][i] = qion * 1.0e-20;       // cross section for e- / Ar ionization
    }
}

//------------------------------------------------------------------------------//
//  ion cross sections: A. V. Phelps, J. Appl. Phys. 76, 747 (1994)             //
//------------------------------------------------------------------------------//

void set_ion_cross_sections_ar(void){
    int    i;
    double e_com,e_lab,qmom,qback,qiso;
    
    printf(">> eduPIC: Setting Ar+ / Ar cross sections\n");
    for(i=0; i<CS_RANGES; i++){
        if (i == 0) {e_com = DE_CS;} else {e_com = DE_CS * i;}             // ion energy in the center of mass frame of reference
        e_lab = 2.0 * e_com;                                               // ion energy in the laboratory frame of reference
        qmom  = 1.15e-18 * pow(e_lab,-0.1) * pow(1.0 + 0.015 / e_lab, 0.6);
        qiso  = 2e-19 * pow(e_lab,-0.5) / (1.0 + e_lab) + 3e-19 * e_lab / pow(1.0 + e_lab / 3.0, 2.0);
        qback = (qmom-qiso) / 2.0;
        sigma[I_ISO][i]  = qiso;             // cross section for Ar+ / Ar isotropic part of elastic scattering
        sigma[I_BACK][i] = qback;            // cross section for Ar+ / Ar backward elastic scattering
    }
}

//----------------------------------------------------------------------//
//  calculation of total cross sections for electrons and ions          //
//----------------------------------------------------------------------//

void calc_total_cross_sections(void){
    int i;
    
    for(i=0; i<CS_RANGES; i++){
        sigma_tot_e[i] = (sigma[E_ELA][i] + sigma[E_EXC][i] + sigma[E_ION][i]) * GAS_DENSITY;   // total macroscopic cross section of electrons
        sigma_tot_i[i] = (sigma[I_ISO][i] + sigma[I_BACK][i]) * GAS_DENSITY;                    // total macroscopic cross section of ions
    }
}

//----------------------------------------------------------------------//
//  test of cross sections for electrons and ions                       //
//----------------------------------------------------------------------//

void test_cross_sections(void){
    FILE  * f;
    int   i,j;
    
    f = fopen("cross_sections.dat","w");        // cross sections saved in data file: cross_sections.dat
    for(i=0; i<CS_RANGES; i++){
        fprintf(f,"%12.4f ",i*DE_CS);
        for(j=0; j<N_CS; j++) fprintf(f,"%14e ",sigma[j][i]);
        fprintf(f,"\n");
    }
    fclose(f);
}

//---------------------------------------------------------------------//
// find upper limit of collision frequencies                           //
//---------------------------------------------------------------------//

double max_electron_coll_freq (void){
    int i;
    double e,v,nu,nu_max;
    nu_max = 0;
    for(i=0; i<CS_RANGES; i++){
        e  = i * DE_CS;
        v  = sqrt(2.0 * e * EV_TO_J / E_MASS);
        nu = v * sigma_tot_e[i];
        if (nu > nu_max) {nu_max = nu;}
    }
    return nu_max;
}

double max_ion_coll_freq (void){
    int i;
    double e,g,nu,nu_max;
    nu_max = 0;
    for(i=0; i<CS_RANGES; i++){
        e  = i * DE_CS;
        g  = sqrt(2.0 * e * EV_TO_J / MU_ARAR);
        nu = g * sigma_tot_i[i];
        if (nu > nu_max) nu_max = nu;
    }
    return nu_max;
}

//----------------------------------------------------------------------//
// initialization of the simulation by placing a given number of        //
// electrons and ions at random positions between the electrodes        //
//----------------------------------------------------------------------//

void init(int nseed){
    int i;
    
    for (i=0; i<nseed; i++){
        x_e[i]  = L * R01(MTgen);               // initial random position of the electron
        vx_e[i] = 0; vy_e[i] = 0; vz_e[i] = 0;  // initial velocity components of the electron
        x_i[i]  = L * R01(MTgen);               // initial random position of the ion
        vx_i[i] = 0; vy_i[i] = 0; vz_i[i] = 0;  // initial velocity components of the ion
    }
    N_e = nseed;    // initial number of electrons
    N_i = nseed;    // initial number of ions
}

//----------------------------------------------------------------------//
// e / Ar collision  (cold gas approximation)                           //
//----------------------------------------------------------------------//

void collision_electron (double xe, double *vxe, double *vye, double *vze, int eindex){
    const double F1 = E_MASS  / (E_MASS + AR_MASS);
    const double F2 = AR_MASS / (E_MASS + AR_MASS);
    double t0,t1,t2,rnd;
    double g,g2,gx,gy,gz,wx,wy,wz,theta,phi;
    double chi,eta,chi2,eta2,sc,cc,se,ce,st,ct,sp,cp,energy,e_sc,e_ej;
    
    // calculate relative velocity before collision & velocity of the centre of mass
    
    gx = (*vxe);
    gy = (*vye);
    gz = (*vze);
    g  = sqrt(gx * gx + gy * gy + gz * gz);
    wx = F1 * (*vxe);
    wy = F1 * (*vye);
    wz = F1 * (*vze);
    
    // find Euler angles
    
    if (gx == 0) {theta = 0.5 * PI;}
    else {theta = atan2(sqrt(gy * gy + gz * gz),gx);}
    if (gy == 0) {
        if (gz > 0){phi = 0.5 * PI;} else {phi = - 0.5 * PI;}
    } else {phi = atan2(gz, gy);}
    st  = sin(theta);
    ct  = cos(theta);
    sp  = sin(phi);
    cp  = cos(phi);
    
    // choose the type of collision based on the cross sections
    // take into account energy loss in inelastic collisions
    // generate scattering and azimuth angles
    // in case of ionization handle the 'new' electron
    
    t0   =     sigma[E_ELA][eindex];
    t1   = t0 +sigma[E_EXC][eindex];
    t2   = t1 +sigma[E_ION][eindex];
    rnd  = R01(MTgen);
    if (rnd < (t0/t2)){                              // elastic scattering
        chi = acos(1.0 - 2.0 * R01(MTgen));          // isotropic scattering
        eta = TWO_PI * R01(MTgen);                   // azimuthal angle
    } else if (rnd < (t1/t2)){                       // excitation
        energy = 0.5 * E_MASS * g * g;               // electron energy
        energy = fabs(energy - E_EXC_TH * EV_TO_J);  // subtract energy loss for excitation
        g   = sqrt(2.0 * energy / E_MASS);           // relative velocity after energy loss
        chi = acos(1.0 - 2.0 * R01(MTgen));          // isotropic scattering
        eta = TWO_PI * R01(MTgen);                   // azimuthal angle
    } else {                                         // ionization
        energy = 0.5 * E_MASS * g * g;               // electron energy
        energy = fabs(energy - E_ION_TH * EV_TO_J);  // subtract energy loss of ionization
        e_ej  = 10.0 * tan(R01(MTgen) * atan(energy/EV_TO_J / 20.0)) * EV_TO_J; // energy of the ejected electron
        e_sc = fabs(energy - e_ej);                  // energy of scattered electron after the collision
        g    = sqrt(2.0 * e_sc / E_MASS);            // relative velocity of scattered electron
        g2   = sqrt(2.0 * e_ej / E_MASS);            // relative velocity of ejected electron
        chi  = acos(sqrt(e_sc / energy));            // scattering angle for scattered electron
        chi2 = acos(sqrt(e_ej / energy));            // scattering angle for ejected electrons
        eta  = TWO_PI * R01(MTgen);                  // azimuthal angle for scattered electron
        eta2 = eta + PI;                             // azimuthal angle for ejected electron
        sc  = sin(chi2);
        cc  = cos(chi2);
        se  = sin(eta2);
        ce  = cos(eta2);
        gx  = g2 * (ct * cc - st * sc * ce);
        gy  = g2 * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
        gz  = g2 * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
        x_e[N_e]  = xe;                              // add new electron
        vx_e[N_e] = wx + F2 * gx;
        vy_e[N_e] = wy + F2 * gy;
        vz_e[N_e] = wz + F2 * gz;
        N_e++;
        x_i[N_i]  = xe;                              // add new ion
        vx_i[N_i] = RMB(MTgen);                      // velocity is sampled from background thermal distribution
        vy_i[N_i] = RMB(MTgen);
        vz_i[N_i] = RMB(MTgen);
        N_i++;
    }
    
    // scatter the primary electron
    
    sc = sin(chi);
    cc = cos(chi);
    se = sin(eta);
    ce = cos(eta);
    
    // compute new relative velocity:
    
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    
    // post-collision velocity of the colliding electron
    
    (*vxe) = wx + F2 * gx;
    (*vye) = wy + F2 * gy;
    (*vze) = wz + F2 * gz;
}

//----------------------------------------------------------------------//
// Ar+ / Ar collision                                                   //
//----------------------------------------------------------------------//

void collision_ion (double *vx_1, double *vy_1, double *vz_1,
                    double *vx_2, double *vy_2, double *vz_2, int e_index){
    double   g,gx,gy,gz,wx,wy,wz,rnd;
    double   theta,phi,chi,eta,st,ct,sp,cp,sc,cc,se,ce,t1,t2;
    
    // calculate relative velocity before collision
    // random Maxwellian target atom already selected (vx_2,vy_2,vz_2 velocity components of target atom come with the call)
    
    gx = (*vx_1)-(*vx_2);
    gy = (*vy_1)-(*vy_2);
    gz = (*vz_1)-(*vz_2);
    g  = sqrt(gx * gx + gy * gy + gz * gz);
    wx = 0.5 * ((*vx_1) + (*vx_2));
    wy = 0.5 * ((*vy_1) + (*vy_2));
    wz = 0.5 * ((*vz_1) + (*vz_2));
    
    // find Euler angles
    
    if (gx == 0) {theta = 0.5 * PI;} else {theta = atan2(sqrt(gy * gy + gz * gz),gx);}
    if (gy == 0) {
        if (gz > 0){phi = 0.5 * PI;} else {phi = - 0.5 * PI;}
    } else {phi = atan2(gz, gy);}
    
    // determine the type of collision based on cross sections and generate scattering angle
    
    t1  =      sigma[I_ISO][e_index];
    t2  = t1 + sigma[I_BACK][e_index];
    rnd = R01(MTgen);
    if  (rnd < (t1 /t2)){                        // isotropic scattering
        chi = acos(1.0 - 2.0 * R01(MTgen));      // scattering angle
    } else {                                     // backward scattering
        chi = PI;                                // scattering angle
    }
    eta = TWO_PI * R01(MTgen);                   // azimuthal angle
    sc  = sin(chi);
    cc  = cos(chi);
    se  = sin(eta);
    ce  = cos(eta);
    st  = sin(theta);
    ct  = cos(theta);
    sp  = sin(phi);
    cp  = cos(phi);
    
    // compute new relative velocity
    
    gx = g * (ct * cc - st * sc * ce);
    gy = g * (st * cp * cc + ct * cp * sc * ce - sp * sc * se);
    gz = g * (st * sp * cc + ct * sp * sc * ce + cp * sc * se);
    
    // post-collision velocity of the ion
    
    (*vx_1) = wx + 0.5 * gx;
    (*vy_1) = wy + 0.5 * gy;
    (*vz_1) = wz + 0.5 * gz;
}

//-----------------------------------------------------------------//
// solve Poisson equation (Thomas algorithm)                       //
//-----------------------------------------------------------------//

void solve_Poisson_CCP (xvector rho1, double tt){
    const double A =  1.0;
    const double B = -2.0;
    const double C =  1.0;
    const double S = 1.0 / (2.0 * DX);
    const double ALPHA = -DX * DX / EPSILON0;
    xvector      g, w, f;
    int          i;
    
    // apply potential to the electrodes - boundary conditions
    
    pot[0]     = VOLTAGE * cos(OMEGA * tt);         // potential at the powered electrode
    pot[N_G-1] = 0.0;                               // potential at the grounded electrode
    
    // solve Poisson equation
    
    for(i=1; i<=N_G-2; i++) f[i] = ALPHA * rho1[i];
    f[1] -= pot[0];
    f[N_G-2] -= pot[N_G-1];
    w[1] = C/B;
    g[1] = f[1]/B;
    for(i=2; i<=N_G-2; i++){
        w[i] = C / (B - A * w[i-1]);
        g[i] = (f[i] - A * g[i-1]) / (B - A * w[i-1]);
    }
    pot[N_G-2] = g[N_G-2];
    for (i=N_G-3; i>0; i--) pot[i] = g[i] - w[i] * pot[i+1];            // potential at the grid points between the electrodes
    
    // compute electric field
    
    for(i=1; i<=N_G-2; i++) efield[i] = (pot[i-1] - pot[i+1]) * S;      // electric field at the grid points between the electrodes
    efield[0]     = (pot[0]     - pot[1])     * INV_DX - rho1[0]     * DX / (2.0 * EPSILON0);   // powered electrode
    efield[N_G-1] = (pot[N_G-2] - pot[N_G-1]) * INV_DX + rho1[N_G-1] * DX / (2.0 * EPSILON0);   // grounded electrode
}
//---------------------------------------------------------------------//
// simulation of one radiofrequency cycle (CCP)                        //
//---------------------------------------------------------------------//

void do_one_cycle_CCP (void){
    const double DV       =  electrode_area * DX;
    const double FACTOR_W = WEIGHT / DV;
    const double FACTOR_E = DT_E / E_MASS * E_CHARGE;
    const double FACTOR_I = DT_I / AR_MASS * E_CHARGE;
    const double MIN_X    = 0.45 * L;                       // min. position for EEPF collection
    const double MAX_X    = 0.55 * L;                       // max. position for EEPF collection
    int      k, t, p, energy_index;
    double   g, g_sqr, gx, gy, gz, vx_a, vy_a, vz_a, e_x, energy, nu, p_coll, v_sqr, velocity;
    double   mean_v, c0, c1, c2, rate;
    bool     out;
    xvector  rho;
    int      t_index;
    
    for (t=0; t<N_T; t++){          // the RF period is divided into N_T equal time intervals (time step DT_E)
        Time += DT_E;               // update of the total simulated time
        t_index = t / N_BIN;        // index for XT distributions
        
        // step 1: compute densities at grid points
 // step 1: compute densities at grid points (electron density every time step)

 #pragma omp parallel
 {
     int tid = omp_get_thread_num();
     static double e_density_private[MAX_THREADS][N_G]; // Still risky, better to allocate outside!
     memset(e_density_private[tid], 0, sizeof(double) * N_G);
 
     #pragma omp for
     for (int k = 0; k < N_e; k++) {
         double c0 = x_e[k] * INV_DX;
         int p = (int)c0;
         if (p >= 0 && p < N_G - 1) {
             double w1 = (p + 1 - c0) * FACTOR_W;
             double w2 = (c0 - p) * FACTOR_W;
             e_density_private[tid][p]   += w1;
             e_density_private[tid][p+1] += w2;
         }
     }
 
     #pragma omp barrier
 
     #pragma omp single
     {
         for (int p = 0; p < N_G; p++) e_density[p] = 0;
 
         for (int t = 0; t < omp_get_num_threads(); t++) {
             for (int p = 0; p < N_G; p++) {
                 e_density[p] += e_density_private[t][p];
             }
         }
 
         e_density[0]     *= 2.0;
         e_density[N_G-1] *= 2.0;
 
         for (int p = 0; p < N_G; p++) cumul_e_density[p] += e_density[p];
     }
 }
 
 if ((t % N_SUB) == 0) {                                            // ion density - computed in every N_SUB-th time steps (subcycling) - CAN'T PARALLELISE!
    for(p=0; p<N_G; p++) i_density[p] = 0;
    for(k=0; k<N_i; k++){
        c0 = x_i[k] * INV_DX;
        p  = int(c0);
        i_density[p]   += (p + 1 - c0) * FACTOR_W;  
        i_density[p+1] += (c0 - p) * FACTOR_W;
    }
    i_density[0]     *= 2.0;
    i_density[N_G-1] *= 2.0;
}
for(p=0; p<N_G; p++) cumul_i_density[p] += i_density[p];
        
        // step 2: solve Poisson equation
        #pragma omp parallel for
        for(p=0; p<N_G; p++) rho[p] = E_CHARGE * (i_density[p] - e_density[p]);  // get charge density
        solve_Poisson_CCP(rho,Time);                                                 // compute potential and electric field
        
        // steps 3 & 4: move particles according to electric field interpolated to particle positions
        
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
        
            // Private thread-local arrays
            double counter_e_xt_private[N_G] = {0.0};
            double ue_xt_private[N_G] = {0.0};
            double meanee_xt_private[N_G] = {0.0};
            double ioniz_rate_xt_private[N_G] = {0.0};
            double eepf_private[N_EEPF] = {0.0};
            double mean_energy_accu_center_private = 0.0;
            double mean_energy_counter_center_private = 0.0;
        
            #pragma omp for private(c0, p, c1, c2, e_x, mean_v, v_sqr, energy, energy_index, velocity, rate)
            for (int k = 0; k < N_e; k++) {  // move all electrons
                c0  = x_e[k] * INV_DX;
                p   = (int)c0;
                if (p < 0 || p >= N_G - 1) continue;  // guard against out-of-bounds
                c1  = p + 1.0 - c0;
                c2  = c0 - p;
                e_x = c1 * efield[p] + c2 * efield[p+1];
        
                if (measurement_mode) {
                    mean_v = vx_e[k] - 0.5 * e_x * FACTOR_E;
        
                    counter_e_xt_private[p]   += c1;
                    counter_e_xt_private[p+1] += c2;
        
                    ue_xt_private[p]   += c1 * mean_v;
                    ue_xt_private[p+1] += c2 * mean_v;
        
                    v_sqr  = mean_v * mean_v + vy_e[k] * vy_e[k] + vz_e[k] * vz_e[k];
                    energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
        
                    meanee_xt_private[p]   += c1 * energy;
                    meanee_xt_private[p+1] += c2 * energy;
        
                    energy_index = (int)(energy / DE_CS + 0.5);
                    if (energy_index > CS_RANGES - 1) energy_index = CS_RANGES - 1;
                    velocity = sqrt(v_sqr);
                    rate = sigma[E_ION][energy_index] * velocity * DT_E * GAS_DENSITY;
        
                    ioniz_rate_xt_private[p]   += c1 * rate;
                    ioniz_rate_xt_private[p+1] += c2 * rate;
        
                    if ((MIN_X < x_e[k]) && (x_e[k] < MAX_X)) {
                        energy_index = (int)(energy / DE_EEPF);
                        if (energy_index < N_EEPF) {
                            eepf_private[energy_index] += 1.0;
                        }
                        mean_energy_accu_center_private += energy;
                        mean_energy_counter_center_private += 1.0;
                    }
                }
        
                // update velocity and position
                vx_e[k] -= e_x * FACTOR_E;
                x_e[k]  += vx_e[k] * DT_E;
            }
        
            // Reduction: sum private arrays into global arrays (only one thread at a time here)
            #pragma omp critical
            {
                for (int i = 0; i < N_G; i++) {
                    counter_e_xt[i][t_index] += counter_e_xt_private[i];
                    ue_xt[i][t_index] += ue_xt_private[i];
                    meanee_xt[i][t_index] += meanee_xt_private[i];
                    ioniz_rate_xt[i][t_index] += ioniz_rate_xt_private[i];
                }
                for (int i = 0; i < N_EEPF; i++) {
                    eepf[i] += eepf_private[i];
                }
                mean_energy_accu_center += mean_energy_accu_center_private;
                mean_energy_counter_center += mean_energy_counter_center_private;
            }
        }
        
        
        if ((t % N_SUB) == 0) {                       // move all ions in every N_SUB-th time steps (subcycling)
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int nthreads = omp_get_num_threads();
            
                // Private thread-local arrays
                double counter_i_xt_private[N_G] = {0.0};
                double ui_xt_private[N_G] = {0.0};
                double meanei_xt_private[N_G] = {0.0};
            
                #pragma omp for private(k, c0, p, c1, c2, e_x, mean_v, v_sqr, energy)
                for (k = 0; k < N_i; k++) {
                    c0  = x_i[k] * INV_DX;
                    p   = (int)c0;
                    if (p < 0 || p >= N_G - 1) continue;  // guard against out-of-bounds
                    c1  = p + 1 - c0;
                    c2  = c0 - p;
                    e_x = c1 * efield[p] + c2 * efield[p+1];
            
                    if (measurement_mode) {
                        mean_v = vx_i[k] + 0.5 * e_x * FACTOR_I;
            
                        counter_i_xt_private[p]   += c1;
                        counter_i_xt_private[p+1] += c2;
            
                        ui_xt_private[p]   += c1 * mean_v;
                        ui_xt_private[p+1] += c2 * mean_v;
            
                        v_sqr  = mean_v * mean_v + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                        energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
            
                        meanei_xt_private[p]   += c1 * energy;
                        meanei_xt_private[p+1] += c2 * energy;
                    }
            
                    // update velocity and position
                    vx_i[k] += e_x * FACTOR_I;
                    x_i[k]  += vx_i[k] * DT_I;
                }
            
                // Reduction: sum private arrays into global arrays
                #pragma omp critical
                {
                    for (int i = 0; i < N_G; i++) {
                        counter_i_xt[i][t_index] += counter_i_xt_private[i];
                        ui_xt[i][t_index] += ui_xt_private[i];
                        meanei_xt[i][t_index] += meanei_xt_private[i];
                    }
                }
            }
            
        }
        
        // step 5: check boundaries
        
        k = 0;
        while(k < N_e) {    // check boundaries for all electrons in every time step
            out = false;
            if (x_e[k] < 0) {N_e_abs_pow++; out = true;}    // the electron is out at the powered electrode
            if (x_e[k] > L) {N_e_abs_gnd++; out = true;}    // the electron is out at the grounded electrode
            if (out) {                                      // remove the electron, if out
                x_e [k] = x_e [N_e-1];
                vx_e[k] = vx_e[N_e-1];
                vy_e[k] = vy_e[N_e-1];
                vz_e[k] = vz_e[N_e-1];
                N_e--;
            } else k++;
        }
        
        if ((t % N_SUB) == 0) {   // check boundaries for all ions in every N_SUB-th time steps (subcycling)
            k = 0;
            while(k < N_i) {
                out = false;
                if (x_i[k] < 0) {       // the ion is out at the powered electrode
                    N_i_abs_pow++;
                    out    = true;
                    v_sqr  = vx_i[k] * vx_i[k] + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS *  v_sqr/ EV_TO_J;
                    energy_index = (int)(energy / DE_IFED);
                    if (energy_index < N_IFED) {ifed_pow[energy_index]++;}       // save IFED at the powered electrode
                }
                if (x_i[k] > L) {       // the ion is out at the grounded electrode
                    N_i_abs_gnd++;
                    out    = true;
                    v_sqr  = vx_i[k] * vx_i[k] + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
                    energy_index = (int)(energy / DE_IFED);
                    if (energy_index < N_IFED) {ifed_gnd[energy_index]++;}        // save IFED at the grounded electrode
                }
                if (out) {  // delete the ion, if out
                    x_i [k] = x_i [N_i-1];
                    vx_i[k] = vx_i[N_i-1];
                    vy_i[k] = vy_i[N_i-1];
                    vz_i[k] = vz_i[N_i-1];
                    N_i--;
                } else k++;
            }
        }
        
        // step 6: collisions
        #pragma omp parallel for private(v_sqr, velocity, energy, energy_index, nu, p_coll)
        for (k=0; k<N_e; k++){                              // checking for occurrence of a collision for all electrons in every time step
            v_sqr = vx_e[k] * vx_e[k] + vy_e[k] * vy_e[k] + vz_e[k] * vz_e[k];
            velocity = sqrt(v_sqr);
            energy   = 0.5 * E_MASS * v_sqr / EV_TO_J;
            energy_index = min( int(energy / DE_CS + 0.5), CS_RANGES-1);
            nu = sigma_tot_e[energy_index] * velocity;
            p_coll = 1 - exp(- nu * DT_E);                  // collision probability for electrons
            if (R01(MTgen) < p_coll) {                      // electron collision takes place
                collision_electron(x_e[k], &vx_e[k], &vy_e[k], &vz_e[k], energy_index);
                #pragma omp atomic
                N_e_coll++;
            }
        }
        
        if ((t % N_SUB) == 0) {                             // checking for occurrence of a collision for all ions in every N_SUB-th time steps (subcycling)
            for (k=0; k<N_i; k++){
                vx_a = RMB(MTgen);                          // pick velocity components of a random target gas atom
                vy_a = RMB(MTgen);
                vz_a = RMB(MTgen);
                gx   = vx_i[k] - vx_a;                       // compute the relative velocity of the collision partners
                gy   = vy_i[k] - vy_a;
                gz   = vz_i[k] - vz_a;
                g_sqr = gx * gx + gy * gy + gz * gz;
                g = sqrt(g_sqr);
                energy = 0.5 * MU_ARAR * g_sqr / EV_TO_J;
                energy_index = min( int(energy / DE_CS + 0.5), CS_RANGES-1);
                nu = sigma_tot_i[energy_index] * g;
                p_coll = 1 - exp(- nu * DT_I);              // collision probability for ions
                if (R01(MTgen)< p_coll) {                   // ion collision takes place
                    collision_ion (&vx_i[k], &vy_i[k], &vz_i[k], &vx_a, &vy_a, &vz_a, energy_index);
                    N_i_coll++;
                }
            }
        }
        
        if (measurement_mode) {
            
            // collect 'xt' data from the grid
            
            for (p=0; p<N_G; p++) {
                pot_xt   [p][t_index] += pot[p];
                efield_xt[p][t_index] += efield[p];
                ne_xt    [p][t_index] += e_density[p];
                ni_xt    [p][t_index] += i_density[p];
            }
        }
        
        if ((t % 1000) == 0) printf(" c = %8d  t = %8d  #e = %8d  #i = %8d\n", cycle,t,N_e,N_i);
    }
    fprintf(datafile,"%8d  %8d  %8d\n",cycle,N_e,N_i);
}

//-----------------------------------------------------------------//
// solve Poisson equation and calculate induced E-field for ICP    //
//-----------------------------------------------------------------//
void solve_Poisson_ICP(xvector rho1, double tt) {
    const double A =  1.0;
    const double B = -2.0;
    const double C =  1.0;
    const double S = 1.0 / (2.0 * DX);
    const double ALPHA = -DX * DX / EPSILON0;
    xvector g, w, f;
    xvector E_induced;
    int i;
    
    // Calculate induced E-field (E = -dB/dt)
    // For sinusoidal B-field: B = B_0 * sin(ωt)
    // Induced E-field: E_induced = -dB/dt = -ω * B_0 * cos(ωt)
    for (i = 0; i < N_G; i++) {
        double pos = i * DX;  // physical position in meters
        
        // Spatially varying induced E-field that mimics radial coil geometry
        // Using sin(π*pos/L) to create maximum field at center, zero at walls
        // The negative sign accounts for Faraday's law: E = -dB/dt
        E_induced[i] = -OMEGA * B_STRENGTH * cos(OMEGA * tt) * sin(PI * pos / L);
    }
    
    // Apply boundary conditions (both ends grounded in ICP)
    pot[0] = pot[N_G-1] = 0.0;
    
    // Solve Poisson equation for electrostatic potential
    for (i = 1; i <= N_G-2; i++) {
        f[i] = ALPHA * rho1[i];
    }
    f[1] -= pot[0];
    f[N_G-2] -= pot[N_G-1];
    
    // Thomas algorithm forward sweep
    w[1] = C/B;
    g[1] = f[1]/B;
    for (i = 2; i <= N_G-2; i++) {
        w[i] = C / (B - A * w[i-1]);
        g[i] = (f[i] - A * g[i-1]) / (B - A * w[i-1]);
    }
    
    // Thomas algorithm backward sweep
    pot[N_G-2] = g[N_G-2];
    for (i = N_G-3; i > 0; i--) {
        pot[i] = g[i] - w[i] * pot[i+1];
    }
    
    // Compute total electric field (electrostatic + induced)
    // Interior points: use central difference for electrostatic field
    for (i = 1; i <= N_G-2; i++) {
        efield[i] = (pot[i-1] - pot[i+1]) * S + E_induced[i];
    }
    
    // Boundary conditions for E-field (including surface charge effects)
    efield[0] = (pot[0] - pot[1]) * INV_DX - rho1[0] * DX / (2.0 * EPSILON0) + E_induced[0];
    efield[N_G-1] = (pot[N_G-2] - pot[N_G-1]) * INV_DX + rho1[N_G-1] * DX / (2.0 * EPSILON0) + E_induced[N_G-1];
}

//---------------------------------------------------------------------//
// simulation of one radiofrequency cycle (ICP)                        //
//---------------------------------------------------------------------//

void do_one_cycle_ICP (void){
    const double DV       = electrode_area * DX;
    const double FACTOR_W = WEIGHT / DV;
    const double FACTOR_E = DT_E / E_MASS * E_CHARGE;
    const double FACTOR_I = DT_I / AR_MASS * E_CHARGE;
    const double MIN_X    = 0.45 * L;                       // min. position for EEPF collection
    const double MAX_X    = 0.55 * L;                       // max. position for EEPF collection
    int      k, t, p, energy_index;
    double   g, g_sqr, gx, gy, gz, vx_a, vy_a, vz_a, e_x, energy, nu, p_coll, v_sqr, velocity;
    double   mean_v, c0, c1, c2, rate;
    bool     out;
    xvector  rho;
    int      t_index;
    
    for (t=0; t<N_T; t++){          // the RF period is divided into N_T equal time intervals (time step DT_E)
        Time += DT_E;               // update of the total simulated time
        t_index = t / N_BIN;        // index for XT distributions
        
        // step 1: compute densities at grid points
 // step 1: compute densities at grid points (electron density every time step)

 #pragma omp parallel
 {
     int tid = omp_get_thread_num();
     static double e_density_private[MAX_THREADS][N_G]; // Still risky, better to allocate outside!
     memset(e_density_private[tid], 0, sizeof(double) * N_G);
 
     #pragma omp for
     for (int k = 0; k < N_e; k++) {
         double c0 = x_e[k] * INV_DX;
         int p = (int)c0;
         if (p >= 0 && p < N_G - 1) {
             double w1 = (p + 1 - c0) * FACTOR_W;
             double w2 = (c0 - p) * FACTOR_W;
             e_density_private[tid][p]   += w1;
             e_density_private[tid][p+1] += w2;
         }
     }
 
     #pragma omp barrier
 
     #pragma omp single
     {
         for (int p = 0; p < N_G; p++) e_density[p] = 0;
 
         for (int t = 0; t < omp_get_num_threads(); t++) {
             for (int p = 0; p < N_G; p++) {
                 e_density[p] += e_density_private[t][p];
             }
         }
 
         e_density[0]     *= 2.0;
         e_density[N_G-1] *= 2.0;
 
         for (int p = 0; p < N_G; p++) cumul_e_density[p] += e_density[p];
     }
 }
 
 if ((t % N_SUB) == 0) {                                            // ion density - computed in every N_SUB-th time steps (subcycling) - CAN'T PARALLELISE!
    for(p=0; p<N_G; p++) i_density[p] = 0;
    for(k=0; k<N_i; k++){
        c0 = x_i[k] * INV_DX;
        p  = int(c0);
        i_density[p]   += (p + 1 - c0) * FACTOR_W;  
        i_density[p+1] += (c0 - p) * FACTOR_W;
    }
    i_density[0]     *= 2.0;
    i_density[N_G-1] *= 2.0;
}
for(p=0; p<N_G; p++) cumul_i_density[p] += i_density[p];


        // step 2: solve Poisson equation
        #pragma omp parallel for
        for(p=0; p<N_G; p++) rho[p] = E_CHARGE * (i_density[p] - e_density[p]);  // get charge density
        solve_Poisson_ICP(rho,Time);                                                 // compute potential and electric field
        
        // steps 3 & 4: move particles according to electric field interpolated to particle positions
        
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            int nthreads = omp_get_num_threads();
        
            // Private thread-local arrays
            double counter_e_xt_private[N_G] = {0.0};
            double ue_xt_private[N_G] = {0.0};
            double meanee_xt_private[N_G] = {0.0};
            double ioniz_rate_xt_private[N_G] = {0.0};
            double eepf_private[N_EEPF] = {0.0};
            double mean_energy_accu_center_private = 0.0;
            double mean_energy_counter_center_private = 0.0;
        
            #pragma omp for private(c0, p, c1, c2, e_x, mean_v, v_sqr, energy, energy_index, velocity, rate)
            for (int k = 0; k < N_e; k++) {  // move all electrons
                c0  = x_e[k] * INV_DX;
                p   = (int)c0;
                if (p < 0 || p >= N_G - 1) continue;  // guard against out-of-bounds
                c1  = p + 1.0 - c0;
                c2  = c0 - p;
                e_x = c1 * efield[p] + c2 * efield[p+1];
        
                if (measurement_mode) {
                    mean_v = vx_e[k] - 0.5 * e_x * FACTOR_E;
        
                    counter_e_xt_private[p]   += c1;
                    counter_e_xt_private[p+1] += c2;
        
                    ue_xt_private[p]   += c1 * mean_v;
                    ue_xt_private[p+1] += c2 * mean_v;
        
                    v_sqr  = mean_v * mean_v + vy_e[k] * vy_e[k] + vz_e[k] * vz_e[k];
                    energy = 0.5 * E_MASS * v_sqr / EV_TO_J;
        
                    meanee_xt_private[p]   += c1 * energy;
                    meanee_xt_private[p+1] += c2 * energy;
        
                    energy_index = (int)(energy / DE_CS + 0.5);
                    if (energy_index > CS_RANGES - 1) energy_index = CS_RANGES - 1;
                    velocity = sqrt(v_sqr);
                    rate = sigma[E_ION][energy_index] * velocity * DT_E * GAS_DENSITY;
        
                    ioniz_rate_xt_private[p]   += c1 * rate;
                    ioniz_rate_xt_private[p+1] += c2 * rate;
        
                    if ((MIN_X < x_e[k]) && (x_e[k] < MAX_X)) {
                        energy_index = (int)(energy / DE_EEPF);
                        if (energy_index < N_EEPF) {
                            eepf_private[energy_index] += 1.0;
                        }
                        mean_energy_accu_center_private += energy;
                        mean_energy_counter_center_private += 1.0;
                    }
                }
        
                // update velocity and position
                vx_e[k] -= e_x * FACTOR_E;
                x_e[k]  += vx_e[k] * DT_E;
            }
        
            // Reduction: sum private arrays into global arrays (only one thread at a time here)
            #pragma omp critical
            {
                for (int i = 0; i < N_G; i++) {
                    counter_e_xt[i][t_index] += counter_e_xt_private[i];
                    ue_xt[i][t_index] += ue_xt_private[i];
                    meanee_xt[i][t_index] += meanee_xt_private[i];
                    ioniz_rate_xt[i][t_index] += ioniz_rate_xt_private[i];
                }
                for (int i = 0; i < N_EEPF; i++) {
                    eepf[i] += eepf_private[i];
                }
                mean_energy_accu_center += mean_energy_accu_center_private;
                mean_energy_counter_center += mean_energy_counter_center_private;
            }
        }
        
        
        if ((t % N_SUB) == 0) {                       // move all ions in every N_SUB-th time steps (subcycling)
            #pragma omp parallel
            {
                int tid = omp_get_thread_num();
                int nthreads = omp_get_num_threads();
            
                // Private thread-local arrays
                double counter_i_xt_private[N_G] = {0.0};
                double ui_xt_private[N_G] = {0.0};
                double meanei_xt_private[N_G] = {0.0};
            
                #pragma omp for private(k, c0, p, c1, c2, e_x, mean_v, v_sqr, energy)
                for (k = 0; k < N_i; k++) {
                    c0  = x_i[k] * INV_DX;
                    p   = (int)c0;
                    if (p < 0 || p >= N_G - 1) continue;  // guard against out-of-bounds
                    c1  = p + 1 - c0;
                    c2  = c0 - p;
                    e_x = c1 * efield[p] + c2 * efield[p+1];
            
                    if (measurement_mode) {
                        mean_v = vx_i[k] + 0.5 * e_x * FACTOR_I;
            
                        counter_i_xt_private[p]   += c1;
                        counter_i_xt_private[p+1] += c2;
            
                        ui_xt_private[p]   += c1 * mean_v;
                        ui_xt_private[p+1] += c2 * mean_v;
            
                        v_sqr  = mean_v * mean_v + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                        energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
            
                        meanei_xt_private[p]   += c1 * energy;
                        meanei_xt_private[p+1] += c2 * energy;
                    }
            
                    // update velocity and position
                    vx_i[k] += e_x * FACTOR_I;
                    x_i[k]  += vx_i[k] * DT_I;
                }
            
                // Reduction: sum private arrays into global arrays
                #pragma omp critical
                {
                    for (int i = 0; i < N_G; i++) {
                        counter_i_xt[i][t_index] += counter_i_xt_private[i];
                        ui_xt[i][t_index] += ui_xt_private[i];
                        meanei_xt[i][t_index] += meanei_xt_private[i];
                    }
                }
            }
            
        }
        
        // step 5: check boundaries
        
        k = 0;
        while(k < N_e) {    // check boundaries for all electrons in every time step
            out = false;
            if (x_e[k] < 0) {N_e_abs_pow++; out = true;}    // the electron is out at the powered electrode
            if (x_e[k] > L) {N_e_abs_gnd++; out = true;}    // the electron is out at the grounded electrode
            if (out) {                                      // remove the electron, if out
                x_e [k] = x_e [N_e-1];
                vx_e[k] = vx_e[N_e-1];
                vy_e[k] = vy_e[N_e-1];
                vz_e[k] = vz_e[N_e-1];
                N_e--;
            } else k++;
        }
        
        if ((t % N_SUB) == 0) {   // check boundaries for all ions in every N_SUB-th time steps (subcycling)
            k = 0;
            while(k < N_i) {
                out = false;
                if (x_i[k] < 0) {       // the ion is out at the powered electrode
                    N_i_abs_pow++;
                    out    = true;
                    v_sqr  = vx_i[k] * vx_i[k] + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS *  v_sqr/ EV_TO_J;
                    energy_index = (int)(energy / DE_IFED);
                    if (energy_index < N_IFED) {ifed_pow[energy_index]++;}       // save IFED at the powered electrode
                }
                if (x_i[k] > L) {       // the ion is out at the grounded electrode
                    N_i_abs_gnd++;
                    out    = true;
                    v_sqr  = vx_i[k] * vx_i[k] + vy_i[k] * vy_i[k] + vz_i[k] * vz_i[k];
                    energy = 0.5 * AR_MASS * v_sqr / EV_TO_J;
                    energy_index = (int)(energy / DE_IFED);
                    if (energy_index < N_IFED) {ifed_gnd[energy_index]++;}        // save IFED at the grounded electrode
                }
                if (out) {  // delete the ion, if out
                    x_i [k] = x_i [N_i-1];
                    vx_i[k] = vx_i[N_i-1];
                    vy_i[k] = vy_i[N_i-1];
                    vz_i[k] = vz_i[N_i-1];
                    N_i--;
                } else k++;
            }
        }
        
        // step 6: collisions
        #pragma omp parallel for private(v_sqr, velocity, energy, energy_index, nu, p_coll)
        for (k=0; k<N_e; k++){                              // checking for occurrence of a collision for all electrons in every time step
            v_sqr = vx_e[k] * vx_e[k] + vy_e[k] * vy_e[k] + vz_e[k] * vz_e[k];
            velocity = sqrt(v_sqr);
            energy   = 0.5 * E_MASS * v_sqr / EV_TO_J;
            energy_index = min( int(energy / DE_CS + 0.5), CS_RANGES-1);
            nu = sigma_tot_e[energy_index] * velocity;
            p_coll = 1 - exp(- nu * DT_E);                  // collision probability for electrons
            if (R01(MTgen) < p_coll) {                      // electron collision takes place
                collision_electron(x_e[k], &vx_e[k], &vy_e[k], &vz_e[k], energy_index);
                #pragma omp atomic
                N_e_coll++;
            }
        }
        
        if ((t % N_SUB) == 0) {                             // checking for occurrence of a collision for all ions in every N_SUB-th time steps (subcycling)
            for (k=0; k<N_i; k++){
                vx_a = RMB(MTgen);                          // pick velocity components of a random target gas atom
                vy_a = RMB(MTgen);
                vz_a = RMB(MTgen);
                gx   = vx_i[k] - vx_a;                       // compute the relative velocity of the collision partners
                gy   = vy_i[k] - vy_a;
                gz   = vz_i[k] - vz_a;
                g_sqr = gx * gx + gy * gy + gz * gz;
                g = sqrt(g_sqr);
                energy = 0.5 * MU_ARAR * g_sqr / EV_TO_J;
                energy_index = min( int(energy / DE_CS + 0.5), CS_RANGES-1);
                nu = sigma_tot_i[energy_index] * g;
                p_coll = 1 - exp(- nu * DT_I);              // collision probability for ions
                if (R01(MTgen)< p_coll) {                   // ion collision takes place
                    collision_ion (&vx_i[k], &vy_i[k], &vz_i[k], &vx_a, &vy_a, &vz_a, energy_index);
                    N_i_coll++;
                }
            }
        }
        
        if (measurement_mode) {
            
            // collect 'xt' data from the grid
            
            for (p=0; p<N_G; p++) {
                pot_xt   [p][t_index] += pot[p];
                efield_xt[p][t_index] += efield[p];
                ne_xt    [p][t_index] += e_density[p];
                ni_xt    [p][t_index] += i_density[p];
            }
        }
        
        if ((t % 1000) == 0) printf(" c = %8d  t = %8d  #e = %8d  #i = %8d\n", cycle,t,N_e,N_i);
    }
    fprintf(datafile,"%8d  %8d  %8d\n",cycle,N_e,N_i);
}

//---------------------------------------------------------------------//
// save particle coordinates                                           //
//---------------------------------------------------------------------//

void save_particle_data(){
    double   d;
    FILE   * f;
    char fname[80];
    
    strcpy(fname,"picdata.bin");
    f = fopen(fname,"wb");
    fwrite(&Time,sizeof(double),1,f);
    d = (double)(cycles_done);
    fwrite(&d,sizeof(double),1,f);
    d = (double)(N_e);
    fwrite(&d,sizeof(double),1,f);
    fwrite(x_e, sizeof(double),N_e,f);
    fwrite(vx_e,sizeof(double),N_e,f);
    fwrite(vy_e,sizeof(double),N_e,f);
    fwrite(vz_e,sizeof(double),N_e,f);
    d = (double)(N_i);
    fwrite(&d,sizeof(double),1,f);
    fwrite(x_i, sizeof(double),N_i,f);
    fwrite(vx_i,sizeof(double),N_i,f);
    fwrite(vy_i,sizeof(double),N_i,f);
    fwrite(vz_i,sizeof(double),N_i,f);
    fclose(f);
    printf(">> eduPIC: data saved : %d electrons %d ions, %d cycles completed, time is %e [s]\n",N_e,N_i,cycles_done,Time);
}

//---------------------------------------------------------------------//
// load particle coordinates                                           //
//---------------------------------------------------------------------//

void load_particle_data(){
    double   d;
    FILE   * f;
    char fname[80];
    
    strcpy(fname,"picdata.bin");
    f = fopen(fname,"rb");
    if (f==NULL) {printf(">> eduPIC: ERROR: No particle data file found, try running initial cycle using argument '0'\n"); exit(0); }
    fread(&Time,sizeof(double),1,f);
    fread(&d,sizeof(double),1,f);
    cycles_done = int(d);
    fread(&d,sizeof(double),1,f);
    N_e = int(d);
    fread(x_e, sizeof(double),N_e,f);
    fread(vx_e,sizeof(double),N_e,f);
    fread(vy_e,sizeof(double),N_e,f);
    fread(vz_e,sizeof(double),N_e,f);
    fread(&d,sizeof(double),1,f);
    N_i = int(d);
    fread(x_i, sizeof(double),N_i,f);
    fread(vx_i,sizeof(double),N_i,f);
    fread(vy_i,sizeof(double),N_i,f);
    fread(vz_i,sizeof(double),N_i,f);
    fclose(f);
    printf(">> eduPIC: data loaded : %d electrons %d ions, %d cycles completed before, time is %e [s]\n",N_e,N_i,cycles_done,Time);
}

//---------------------------------------------------------------------//
// save density data                                                   //
//---------------------------------------------------------------------//

void save_density(void){
    FILE *f;
    double c;
    int m;
    
    f = fopen("density.dat","w");
    c = 1.0 / (double)(no_of_cycles) / (double)(N_T);
    for(m=0; m<N_G; m++){
        fprintf(f,"%8.5f  %12e  %12e\n",m * DX, cumul_e_density[m] * c, cumul_i_density[m] * c);
    }
    fclose(f);
}

//---------------------------------------------------------------------//
// save EEPF data                                                      //
//---------------------------------------------------------------------//

void save_eepf(void) {
    FILE   *f;
    int    i;
    double h,energy;
    
    h = 0.0;
    for (i=0; i<N_EEPF; i++) {h += eepf[i];}
    h *= DE_EEPF;
    f = fopen("eepf.dat","w");
    for (i=0; i<N_EEPF; i++) {
        energy = (i + 0.5) * DE_EEPF;
        fprintf(f,"%e  %e\n", energy, eepf[i] / h / sqrt(energy));
    }
    fclose(f);
}

//---------------------------------------------------------------------//
// save IFED data                                                      //
//---------------------------------------------------------------------//

void save_ifed(void) {
    FILE   *f;
    int    i;
    double h_pow,h_gnd,energy;
    
    h_pow = 0.0;
    h_gnd = 0.0;
    for (i=0; i<N_IFED; i++) {h_pow += ifed_pow[i]; h_gnd += ifed_gnd[i];}
    h_pow *= DE_IFED;
    h_gnd *= DE_IFED;
    mean_i_energy_pow = 0.0;
    mean_i_energy_gnd = 0.0;
    f = fopen("ifed.dat","w");
    for (i=0; i<N_IFED; i++) {
        energy = (i + 0.5) * DE_IFED;
        fprintf(f,"%6.2f %10.6f %10.6f\n", energy, (double)(ifed_pow[i])/h_pow, (double)(ifed_gnd[i])/h_gnd);
        mean_i_energy_pow += energy * (double)(ifed_pow[i]) / h_pow;
        mean_i_energy_gnd += energy * (double)(ifed_gnd[i]) / h_gnd;
    }
    fclose(f);
}

//--------------------------------------------------------------------//
// save XT data                                                       //
//--------------------------------------------------------------------//

void save_xt_1(xt_distr distr, char *fname) {
    FILE   *f;
    int    i, j;
    
    f = fopen(fname,"w");
    for (i=0; i<N_G; i++){
        for (j=0; j<N_XT; j++){
            fprintf(f,"%e  ", distr[i][j]);
        }
        fprintf(f,"\n");
    }
    fclose(f);
}

void norm_all_xt(void){
    double f1, f2;
    int    i, j;
    
    // normalize all XT data
    
    f1 = (double)(N_XT) / (double)(no_of_cycles * N_T);
    f2 = WEIGHT / ( electrode_area * DX) / (no_of_cycles * (PERIOD / (double)(N_XT)));
    
    for (i=0; i<N_G; i++){
        for (j=0; j<N_XT; j++){
            pot_xt[i][j]    *= f1;
            efield_xt[i][j] *= f1;
            ne_xt[i][j]     *= f1;
            ni_xt[i][j]     *= f1;
            if (counter_e_xt[i][j] > 0) {
                ue_xt[i][j]     =  ue_xt[i][j] / counter_e_xt[i][j];
                je_xt[i][j]     = -ue_xt[i][j] * ne_xt[i][j] * E_CHARGE;
                meanee_xt[i][j] =  meanee_xt[i][j] / counter_e_xt[i][j];
                ioniz_rate_xt[i][j] *= f2;
             } else {
                ue_xt[i][j]         = 0.0;
                je_xt[i][j]         = 0.0;
                meanee_xt[i][j]     = 0.0;
                ioniz_rate_xt[i][j] = 0.0;
            }
            if (counter_i_xt[i][j] > 0) {
                ui_xt[i][j]     = ui_xt[i][j] / counter_i_xt[i][j];
                ji_xt[i][j]     = ui_xt[i][j] * ni_xt[i][j] * E_CHARGE;
                meanei_xt[i][j] = meanei_xt[i][j] / counter_i_xt[i][j];
            } else {
                ui_xt[i][j]     = 0.0;
                ji_xt[i][j]     = 0.0;
                meanei_xt[i][j] = 0.0;
            }
            powere_xt[i][j] = je_xt[i][j] * efield_xt[i][j];
            poweri_xt[i][j] = ji_xt[i][j] * efield_xt[i][j];
        }
    }
}

void save_all_xt(void){
    char fname[80];
    
    strcpy(fname,"pot_xt.dat");     save_xt_1(pot_xt, fname);
    strcpy(fname,"efield_xt.dat");  save_xt_1(efield_xt, fname);
    strcpy(fname,"ne_xt.dat");      save_xt_1(ne_xt, fname);
    strcpy(fname,"ni_xt.dat");      save_xt_1(ni_xt, fname);
    strcpy(fname,"je_xt.dat");      save_xt_1(je_xt, fname);
    strcpy(fname,"ji_xt.dat");      save_xt_1(ji_xt, fname);
    strcpy(fname,"powere_xt.dat");  save_xt_1(powere_xt, fname);
    strcpy(fname,"poweri_xt.dat");  save_xt_1(poweri_xt, fname);
    strcpy(fname,"meanee_xt.dat");  save_xt_1(meanee_xt, fname);
    strcpy(fname,"meanei_xt.dat");  save_xt_1(meanei_xt, fname);
    strcpy(fname,"ioniz_xt.dat");   save_xt_1(ioniz_rate_xt, fname);
}

//---------------------------------------------------------------------//
// simulation report including stability and accuracy conditions       //
//---------------------------------------------------------------------//
Info check_and_save_info(double timedifsec) {
    FILE *f;
    double kT, e_max, v_max, c;
    int i, j;
    bool conditions_OK;
    Info info;

    
    info.density    = cumul_e_density[N_G / 2] / (double)(no_of_cycles) / (double)(N_T);  // e density @ center
    info.plas_freq  = E_CHARGE * sqrt(info.density / EPSILON0 / E_MASS);                       // e plasma frequency @ center
    info.meane      = mean_energy_accu_center / (double)(mean_energy_counter_center);     // e mean energy @ center
    kT         = 2.0 * info.meane * EV_TO_J / 3.0;                                        // k T_e @ center (approximate)
    info.sim_time   = (double)(no_of_cycles) / FREQUENCY;                                 // simulated time
    info.ecoll_freq = (double)(N_e_coll) / info.sim_time / (double)(N_e);                      // e collision frequency
    info.icoll_freq = (double)(N_i_coll) / info.sim_time / (double)(N_i);                      // ion collision frequency
    info.debye_length = sqrt(EPSILON0 * kT / info.density) / E_CHARGE;                         // e Debye length @ center    
    double skinDepth = c_light / info.plas_freq; //skin depth of plasma

    f = fopen("info.txt","w");
    fprintf(f,"########################## eduPIC simulation report ############################\n");
    fprintf(f,"Simulation parameters:\n");
    fprintf(f,"Gap distance                          = %12.3e [m]\n",  L);
    fprintf(f,"# of grid divisions                   = %12d\n",      N_G);
    fprintf(f,"Frequency                             = %12.3e [Hz]\n", FREQUENCY);
    fprintf(f,"# of time steps / period              = %12d\n",      N_T);
    fprintf(f,"# of electron / ion time steps        = %12d\n",      N_SUB);
    fprintf(f,"Voltage amplitude                     = %12.3e [V]\n",  VOLTAGE);
    fprintf(f,"Pressure (Ar)                         = %12.3e [Pa]\n", PRESSURE);
    fprintf(f,"Temperature                           = %12.3e [K]\n",  TEMPERATURE);
    fprintf(f,"Superparticle weight                  = %12.3e\n",      WEIGHT);
    fprintf(f,"# of simulation cycles in this run    = %12d\n",      no_of_cycles);
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fprintf(f,"Plasma characteristics:\n");
    fprintf(f,"Electron density @ center             = %12.3e [m^{-3}]\n", info.density);
    fprintf(f,"Plasma frequency @ center             = %12.3e [rad/s]\n",  info.plas_freq);
    fprintf(f,"Plasma Skin Depth         = %12.10e [m]\n", skinDepth);
    fprintf(f,"Debye length @ center                 = %12.3e [m]\n",      info.debye_length);
    fprintf(f,"Electron collision frequency          = %12.3e [1/s]\n",    info.ecoll_freq);
    fprintf(f,"Ion collision frequency               = %12.3e [1/s]\n",    info.icoll_freq);
    fprintf(f,"--------------------------------------------------------------------------------\n");
    fprintf(f,"Stability and accuracy conditions:\n");
    conditions_OK = true;
    c = info.plas_freq * DT_E;
    fprintf(f,"Plasma frequency @ center * DT_E      = %12.3f (OK if less than 0.20)\n", c);
    if (c > 0.2) {conditions_OK = false;}
    c = DX / info.debye_length;
    fprintf(f,"DX / Debye length @ center            = %12.3f (OK if less than 1.00)\n", c);
    if (c > 1.0) {conditions_OK = false;}
    c = max_electron_coll_freq() * DT_E;
    fprintf(f,"Max. electron coll. frequency * DT_E  = %12.3f (OK if less than 0.05)\n", c);
    if (c > 0.05) {conditions_OK = false;}
    c = max_ion_coll_freq() * DT_I;
    fprintf(f,"Max. ion coll. frequency * DT_I       = %12.3f (OK if less than 0.05)\n", c);
    if (c > 0.05) {conditions_OK = false;}
    if (conditions_OK == false){
        fprintf(f,"--------------------------------------------------------------------------------\n");
        fprintf(f,"** STABILITY AND ACCURACY CONDITION(S) VIOLATED - REFINE SIMULATION SETTINGS! **\n");
        fprintf(f,"--------------------------------------------------------------------------------\n");
        fclose(f);
        printf(">> eduPIC: ERROR: STABILITY AND ACCURACY CONDITION(S) VIOLATED!\n");
        printf(">> eduPIC: for details see 'info.txt' and refine simulation settings!\n");
    }
    else
    {
        // calculate maximum energy for which the Courant-Friedrichs-Levy condition holds:
        
        v_max = DX / DT_E;
        e_max = 0.5 * E_MASS * v_max * v_max / EV_TO_J;
        fprintf(f,"Max e- energy for CFL condition       = %12.3f [eV]\n", e_max);
        fprintf(f,"Check EEPF to ensure that CFL is fulfilled for the majority of the electrons!\n");
        fprintf(f,"--------------------------------------------------------------------------------\n");
        
        // saving of the following data is done here as some of the further lines need data
        // that are computed / normalized in these functions
        
        printf(">> eduPIC: saving diagnostics data\n");
        save_density();
        save_eepf();
        save_ifed();
        norm_all_xt();
        save_all_xt();
        fprintf(f,"Particle characteristics at the electrodes:\n");
        fprintf(f,"Ion flux at powered electrode         = %12.3e [m^{-2} s^{-1}]\n", N_i_abs_pow * WEIGHT /  electrode_area / (no_of_cycles * PERIOD));
        fprintf(f,"Ion flux at grounded electrode        = %12.3e [m^{-2} s^{-1}]\n", N_i_abs_gnd * WEIGHT /  electrode_area / (no_of_cycles * PERIOD));
        fprintf(f,"Mean ion energy at powered electrode  = %12.3e [eV]\n", mean_i_energy_pow);
        fprintf(f,"Mean ion energy at grounded electrode = %12.3e [eV]\n", mean_i_energy_gnd);
        fprintf(f,"Electron flux at powered electrode    = %12.3e [m^{-2} s^{-1}]\n", N_e_abs_pow * WEIGHT /  electrode_area / (no_of_cycles * PERIOD));
        fprintf(f,"Electron flux at grounded electrode   = %12.3e [m^{-2} s^{-1}]\n", N_e_abs_gnd * WEIGHT /  electrode_area / (no_of_cycles * PERIOD));
        fprintf(f,"--------------------------------------------------------------------------------\n");
        
        // calculate spatially and temporally averaged power absorption by the electrons and ions
        
        info.power_e = 0.0;
        info.power_i = 0.0;
        for (i=0; i<N_G; i++){
            for (j=0; j<N_XT; j++){
                info.power_e += powere_xt[i][j];
                info.power_i += poweri_xt[i][j];
            }
        }
        info.power_e /= (double)(N_XT * N_G);
        info.power_i /= (double)(N_XT * N_G);
        fprintf(f,"Absorbed power calculated as <j*E>:\n");
        fprintf(f,"Electron power density (average)      = %12.3e [W m^{-3}]\n", info.power_e);
        fprintf(f,"Ion power density (average)           = %12.3e [W m^{-3}]\n", info.power_i);
        fprintf(f,"Total power density(average)          = %12.3e [W m^{-3}]\n", info.power_e + info.power_i);
        fprintf(f,"--------------------------------------------------------------------------------\n");

        // Heating Efficiency calculation 
        double V_rms = VOLTAGE / sqrt(2.0);  // Convert peak to RMS voltage
        double ne = info.density;  // Electron density AT CENTRE(m^-3)
        double nu_m_e = info.ecoll_freq;  // Electron collision frequency (1/s)
        double nu_m_i = info.icoll_freq;  // Ion collision frequency (1/s)
        double power_density = info.power_e + info.power_i;  // Power density (W/m^3)
        
        // Compute plasma resistance
        info.R_plasma = (((E_MASS * info.ecoll_freq) + (AR_MASS * info.icoll_freq)) / (E_CHARGE * E_CHARGE * info.density)) * (L / electrode_area);
        fprintf(f,"Plasma Resistance        = %12.3e [Ohms]\n", info.R_plasma);
        fprintf(f,"--------------------------------------------------------------------------------\n");
        
        // Compute input power
        double P_in = (V_rms * V_rms) / info.R_plasma;  //50ohm for ideal matched load

        // Compute absorbed power
        double plasma_volume =  electrode_area * L;  // Plasma volume (m^3)
        double P_abs = power_density * plasma_volume;  // Absorbed power (W)

        // Compute efficiency
        info.heating_efficiency = (P_abs / P_in) * 100.0;  // Store the value in the global variable
        fprintf(f,"Total Input Power         = %12.3e [W m^{-3}]\n", P_in);
        fprintf(f,"--------------------------------------------------------------------------------\n");
        // Print to file
        fprintf(f,"Heating Efficiency        = %12.10e [%%]\n", info.heating_efficiency);
        fprintf(f,"--------------------------------------------------------------------------------\n");

        int hours   = timedifsec / 3600;
        int minutes = ((int)timedifsec % 3600) / 60;
        int seconds = ((int)timedifsec % 60);
        
        fprintf(f,"Total simulation run time             = %02d:%02d:%02d [hh:mm:ss]\n", hours, minutes, seconds);

        fprintf(f,"--------------------------------------------------------------------------------\n");

        fclose(f);    
        }
        return info;
    }

//----------------------------------------------------------------------------------------
//   Genetic Algorithm method of optimistion
//----------------------------------------------------------------------------------------

    double compute_heating_efficiency(Info info) {
        double V_rms = VOLTAGE / sqrt(2.0);
        double P_in = (V_rms * V_rms) / info.R_plasma;
        double plasma_volume = electrode_area * L;
        double P_abs = (info.power_e + info.power_i) * plasma_volume;

        return (P_abs / P_in) * 100.0;
    }

// Typedefs for modular optimisation
typedef double (*MutationFunc)(double);
typedef double (*InitFunc)();

// Store original GUI parameters during optimisation
double gui_power, gui_pressure, gui_gap, gui_temperature, gui_turns;

double crossover(double parent1, double parent2) {
    return (parent1 + parent2) / 2.0;
}
double random_temperature() {
    return MIN_TEMPERATURE + ((double)rand() / RAND_MAX) * (MAX_TEMPERATURE - MIN_TEMPERATURE);
}
double random_pressure() {
    return MIN_PRESSURE + ((double)rand() / RAND_MAX) * (MAX_PRESSURE - MIN_PRESSURE);
}
double random_voltage() {
    return MIN_VOLTAGE + ((double)rand() / RAND_MAX) * (MAX_VOLTAGE - MIN_VOLTAGE);
}
double random_gap() {
    return MIN_GAP + ((double)rand() / RAND_MAX) * (MAX_GAP - MIN_GAP);
}


double mutate_temperature(double temperature) {
    double base = gui_temperature;  // Use GUI value as baseline
    if (((double)rand() / RAND_MAX) < MUTATION_RATE) {
        double mutation = ((double)rand()/RAND_MAX)*100.0 - 50.0; // mutate by +/- 50K     //((double)rand() / RAND_MAX) * 0.2 - 0.1; // a ±10% variation
        temperature += mutation;     //base * (1.0 + mutation);
        temperature = fmax(MIN_TEMPERATURE, fmin(temperature, MAX_TEMPERATURE));
    }
    return temperature;
}
double mutate_pressure(double pressure) {
    double base = gui_pressure;  // Use GUI value as baseline
    if (((double)rand() / RAND_MAX) < MUTATION_RATE) {
        double mutation = ((double)rand()/RAND_MAX)* 6.0 - 3.0; // mutate by +/- 2Pa   //((double)rand() / RAND_MAX) * 0.2 - 0.1; // a ±10% variation
        pressure += mutation;     //base * (1.0 + mutation);
        pressure = fmax(MIN_PRESSURE, fmin(pressure, MAX_PRESSURE));
    }
    return pressure;
}
double mutate_voltage(double voltage) {
    double base = VOLTAGE;  // Use GUI value as baseline
    if (((double)rand() / RAND_MAX) < MUTATION_RATE) {
        double mutation = ((double)rand()/RAND_MAX)*40.0 - 20.0; // mutate by +/- 20V
        voltage += mutation;
        voltage = fmax(MIN_VOLTAGE, fmin(voltage, MAX_VOLTAGE));
    }
    return voltage;
}
double mutate_gap(double gap) {
    double base = gui_gap;  // Use GUI value as baseline
    if (((double)rand() / RAND_MAX) < MUTATION_RATE) {
        double mutation = ((double)rand()/RAND_MAX)*1.0 - 0.5; // mutate by +/- 0.5mm     //((double)rand() / RAND_MAX) * 0.2 - 0.1; // a ±10% variation
        gap += mutation;      // base * (1.0 + mutation);
        gap = fmax(MIN_GAP, fmin(gap, MAX_GAP));
    }
    return gap;
}

// Configuration structure for optimisation modes
typedef struct {
    const char* param_name;
    const char* filename;
    MutationFunc mutate;
    InitFunc init;
    double* param_ptr;
} ModeConfig;

// Modified logging function with parameter name
void log_heating_efficiency(int cycle, double efficiency, double param_value, const char* param_name) {
    FILE *log_file = fopen("heating_efficiency_log.txt", "a");
    if (log_file == NULL) {
        perror("Error opening log file");
        return;
    }
    fprintf(log_file, "Cycle: %d | Efficiency: %.10f %% | %s: %.6e\n", 
            cycle, efficiency, param_name, param_value);
    fclose(log_file);
}
void log_heating_efficiency_no_opt(int cycle, double efficiency) {
    FILE *log_file = fopen("heating_efficiency_log.txt", "a");
    if (log_file == NULL) {
        perror("Error opening log file");
        return;
    }
    fprintf(log_file, "Cycle: %d | Efficiency: %.10f %% ", 
            cycle, efficiency);
    fclose(log_file);
}

// Mode-specific configuration getter
ModeConfig get_mode_config(char mode) {
    ModeConfig config;
    switch(mode) {
        case 'T':
            config.param_name = "Temperature";
            config.filename = "optimised_temp.txt";
            config.mutate = mutate_temperature;
            config.init = random_temperature;
            config.param_ptr = &TEMPERATURE;
            break;
        case 'P':
            config.param_name = "Pressure";
            config.filename = "optimised_press.txt";
            config.mutate = mutate_pressure;
            config.init = random_pressure;
            config.param_ptr = &PRESSURE;
            break;
        case 'W':
            config.param_name = "Voltage";
            config.filename = "optimised_volt.txt";
            config.mutate = mutate_voltage;
            config.init = random_voltage;
            config.param_ptr = &VOLTAGE;
            break;
        case 'G':
            config.param_name = "Gap";
            config.filename = "optimised_gap.txt";
            config.mutate = mutate_gap;
            config.init = random_gap;
            config.param_ptr = &L;
            break;
        default:
            printf(">> Unknown mode '%c'\n", mode);
            exit(1);
    }
    return config;
}

// Modified genetic algorithm handler
double handle_genetic_algorithm(Info info, const ModeConfig* config) {
    
    // Preserve non-optimised parameters from GUI
    double fixed_params[3];
    int param_index = 0;
    
    for(int i=0; i<4; i++){
        if(config->param_ptr != &VOLTAGE + i) {
            fixed_params[param_index++] = *(&VOLTAGE + i);
        }
    }

    double population[POP_SIZE];
    double fitness[POP_SIZE];

    for (int i = 0; i < POP_SIZE; i++) {
        population[i] = config->init();
    }

    for (int gen = 0; gen < GENERATIONS; gen++) {
        for (int i = 0; i < POP_SIZE; i++) {
            double original = *(config->param_ptr);
            *(config->param_ptr) = population[i];
            double efficiency = compute_heating_efficiency(info);

            // Discard or penalize unphysical solutions
            if (efficiency > 100.0 || efficiency < 0.0) {
                fitness[i] = -1.0;  // Penalize bad efficiency (or use 0.0)
            } else {
                fitness[i] = efficiency;
            }
            *(config->param_ptr) = original;


            if (fitness[i] >= TARGET_EFFICIENCY) return population[i];
        }

        int best1 = 0, best2 = 1;
        if (fitness[best2] > fitness[best1]) { int t = best1; best1 = best2; best2 = t; }
        for (int i = 2; i < POP_SIZE; i++) {
            if (fitness[i] > fitness[best1]) {
                best2 = best1;
                best1 = i;
            } else if (fitness[i] > fitness[best2]) {
                best2 = i;
            }
        }

        for (int i = 0; i < POP_SIZE; i++) {
            double child = crossover(population[best1], population[best2]);
            population[i] = config->mutate(child);
        }
    }

    int best_idx = 0;
    for (int i = 1; i < POP_SIZE; i++) {
        if (fitness[i] > fitness[best_idx]) best_idx = i;
    }
    
    return population[best_idx];
}


// Modified optimization runner with simulation type selection
__declspec(dllexport)
void run_optimisation(const char* folder_path, char mode, char sim_type, int arg1, double power, double pressure, double gap, double temp, double turns, double min_power, double max_power, double min_pressure, double max_pressure, double min_temp, double max_temp, double min_l, double max_l ) {
    omp_set_num_threads(4);
    Info info;
    double current_param = 0.0;
    double current_efficiency = 0.0;
    
    if (mode == 'N') {

     // Set parameters from GUI inputs
    gui_power = POWER = power;
    gui_pressure = PRESSURE = pressure;
    gui_gap = L = gap;
    gui_temperature = TEMPERATURE = temp;
    gui_turns = N_TURNS = turns;
    
    // create folder if it doesn't exist
    if (access(folder_path, 0) == -1) {
        _mkdir(folder_path);
    }
    // move into the directory
    chdir(folder_path); // set working directory


    printf(">> [opt.c] Running simulation WITHOUT optimization.\n");
    measurement_mode = true;

    
    VOLTAGE = sqrt(2 * power * 50);
            printf(">> VOLTAGE = %.3f V.\n", VOLTAGE);


    time1 = clock() / CLOCKS_PER_SEC;

    set_electron_cross_sections_ar();
    set_ion_cross_sections_ar();
    calc_total_cross_sections();
    

    if (arg1 == 0) {
        // Initialization path (same as optimization mode)
        if (FILE *f = fopen("picdata.bin", "r")) {
            fclose(f);
            printf(">> Warning: Previous data exists.\n");
            exit(0);
        }

        datafile = fopen("conv.dat", "a");
        no_of_cycles = 1;
        cycle = 1;
        init(N_INIT);
        Time = 0;

        if (sim_type == 'C') {
            do_one_cycle_CCP();
        } else {
            do_one_cycle_ICP();
        }
        cycles_done = 1;

        fclose(datafile);
        save_particle_data();
        printf(">> [DLL] Initialization completed (non-optimized mode).\n");
    
    } else {
        no_of_cycles = arg1;
        load_particle_data();
        // Create run_1 folder to match GUI expectations
        if (access("run_1", F_OK) == -1) {
            mkdir("run_1");
        }
        chdir("run_1");

        // Full run without optimization loop
        datafile = fopen("conv.dat", "w");
        if (!datafile) {
            printf(">> Error: could not open conv.dat for writing\n");
            exit(1);
        }

        for (cycle = cycles_done + 1; cycle <= cycles_done + no_of_cycles; cycle++) {
            if (sim_type == 'C') {
                do_one_cycle_CCP();
            } else {
                do_one_cycle_ICP();
            }
        }

            fclose(datafile);
            save_particle_data();

            timedif1 = (clock() / CLOCKS_PER_SEC) - time1;
           
            info = check_and_save_info(timedif1);
            heating_efficiency = info.heating_efficiency;
            log_heating_efficiency_no_opt(cycle, heating_efficiency);
        
            timedifHR = ((clock() / CLOCKS_PER_SEC) - time1) / 3600.0;

            fclose(datafile);
            save_particle_data();
        printf(">> [DLL] Non-optimized simulation completed in %.3f hours.\n", timedifHR);
        chdir("..");
    }

    return;
} else {
    // create folder if it doesn't exist
    if (access(folder_path, 0) == -1) {
        _mkdir(folder_path);
    }

    // move into the directory
    chdir(folder_path); // set working directory



    ModeConfig config = get_mode_config(mode);
    
    MIN_POWER = min_power;
    MAX_POWER = max_power;

    MIN_PRESSURE = min_pressure;
    MAX_PRESSURE = max_pressure;

    MIN_L = min_l;
    MAX_L = max_l;

    MIN_TEMPERATURE = min_temp;
    MAX_TEMPERATURE = max_temp;
    printf(">> opt.c bounds are: min power = %.3f, max power = %.3f, min pressure = %.3f, max pressure = %.3f, min temp = %.1f, max temp = %.1f, min gap = %.4f, max gap = %.4f\n",
        MIN_POWER, MAX_POWER, MIN_PRESSURE, MAX_PRESSURE, MIN_TEMPERATURE, MAX_TEMPERATURE, MIN_L, MAX_L);
    
    MIN_VOLTAGE = sqrt(MIN_POWER*50);
    MAX_VOLTAGE = sqrt(MAX_POWER*50);
    MIN_GAP = MIN_L;
    MAX_GAP = MAX_L;

    // Set parameters from GUI inputs
    gui_power = POWER = power;
    gui_pressure = PRESSURE = pressure;
    gui_gap = L = gap;
    gui_temperature = TEMPERATURE = temp;
    gui_turns = N_TURNS = turns;
    
    double Vrms = sqrt(power * 50);
    VOLTAGE = sqrt(2) * Vrms;

    // Validate simulation type
    if(sim_type != 'C' && sim_type != 'I') {
        printf(">> Invalid simulation type! Use 'C' for CCP or 'I' for ICP\n");
        exit(1);
    }

    measurement_mode = (mode == 'm');
    measurement_mode = true;
    time1 = clock() / CLOCKS_PER_SEC;

    set_electron_cross_sections_ar();
    set_ion_cross_sections_ar();
    calc_total_cross_sections();


    if (arg1 == 0) {
        if (FILE *f = fopen("picdata.bin", "r")) {
            fclose(f);
            printf(">> Warning: Previous data exists.\n");
            exit(0);
        }
        datafile = fopen("conv.dat", "a");
        no_of_cycles = 1;
        cycle = 1;
        init(N_INIT);
        Time = 0;
        // Choose simulation type
        if (sim_type == 'C'){
            do_one_cycle_CCP();
        } else {
            do_one_cycle_ICP();
        } 
        cycles_done = 1;
        fclose(datafile);
        save_particle_data();
    } else {
        int folder_index = 1;
        char prev_folder[64] = ".";

        while (1) {
            no_of_cycles = arg1;
            load_particle_data();

            char folder_name[64];
            snprintf(folder_name, sizeof(folder_name), "run_%d", folder_index);
            if (access(folder_name, F_OK) == -1) {
                mkdir(folder_name);
            }
            chdir(folder_name);

            datafile = fopen("conv.dat", "w");
            if (!datafile) {
                printf(">> Error: could not open conv.dat in %s\n", folder_name);
                exit(1);
            }

            if (folder_index > 1) {
                char param_path[128];
                snprintf(param_path, sizeof(param_path), "../%s/%s", prev_folder, config.filename);
                
                FILE *f = fopen(param_path, "r");
                double val;
                if (f && fscanf(f, "%lf", &val) == 1) {
                    *(config.param_ptr) = val;
                    printf(">> Using optimised %s from %s: %.6f\n", 
                          config.param_name, param_path, val);
                } else {
                    printf(">> Failed to read parameter from %s\n", param_path);
                    exit(1);
                }
            }

            for (cycle = cycles_done + 1; cycle <= cycles_done + no_of_cycles; cycle++) {
                // Simulation type selection in main loop
                if (sim_type == 'C'){
                    do_one_cycle_CCP();
                } else {
                    do_one_cycle_ICP();
                } 
            }
            cycles_done += no_of_cycles;

            fclose(datafile);
            save_particle_data();

            timedif1 = (clock() / CLOCKS_PER_SEC) - time1;
           
            info = check_and_save_info(timedif1);
/*
            // Use actual simulated plasma resistance from this run
            R_plas = info.R_plasma;
            Vrms = sqrt(power * R_plas);
            VOLTAGE = sqrt(2) * Vrms;
*/
            heating_efficiency = info.heating_efficiency;
            current_param = *(config.param_ptr);
            log_heating_efficiency(cycle, heating_efficiency, current_param, config.param_name);

            if (heating_efficiency >= TARGET_EFFICIENCY) {
                printf(">> Target efficiency reached: %.2f%% in %s.\n", heating_efficiency, folder_name);
                break;
            } else {strncpy(info.param_name, config.param_name, sizeof(info.param_name));
            info.param_name[sizeof(info.param_name) - 1] = '\0';  // Ensure null-termination

            double new_value = handle_genetic_algorithm(info, &config);

            FILE *opt_file = fopen(config.filename, "w");
            if (opt_file) {
                fprintf(opt_file, "%.6f\n", new_value);
                fclose(opt_file);
                printf(">> Saved optimised %s: %.6f\n", config.param_name, new_value);
            } else {
                printf(">> Failed to write to %s\n", config.filename);
                exit(1);
            }

            strcpy(prev_folder, folder_name);
            folder_index++;
            chdir("..");
            cycles_done = 1;
            }    
        }
    }

/*
    timedif2 = (clock() / CLOCKS_PER_SEC) - time1;
    timedifHR = timedif2 / 3600.0;

    if (arg1 > 0) {
        info = check_and_save_info(timedif2);
        log_heating_efficiency(cycle, heating_efficiency, current_param, config.param_name);
    }
*/
    printf(">> [DLL] Simulation completed in %.3f hours.\n", timedifHR);
}
}