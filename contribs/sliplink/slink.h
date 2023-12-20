#pragma once

#include <map>

struct Bead
{
double XYZ[3];
double XYZold[3];
int index; // = global id of bead
int terminal; // true if bead is the right most
int index_in_chain; // 0 is left-most bead in chain, etc.
int id_ch;
int id_inf;
int id_sup;
double F_SP[3];
double F_SL[3];
};

extern std::map< int , int > g_bead_index;

extern struct Bead *bead;

struct SlipLink
{
double XYZsj[3];
double XYZaj[3];
double xj;
int truncxj;
double id;
int chainj;
int id_left_bead;
};

extern struct SlipLink *sl;

void sl_init(SlipLink *sl, Bead *bead, double b, double Ns, int zzz, int n_e, int n_ch, int n_m, long int *in_ran1, int *num_SL);
void chain_init(Bead *bead, double b, long int *in_ran1, int n_ch, int n_m, double box[3]);
void sl_destroy(int *destroy, int zzz);
void sl_renew(double b, Bead *bead, SlipLink *sl, double Ns, int zzz, int n_ch, int istep, int n_m, long int *in_ran1, int *destroy);
void move_beads(Bead *bead, SlipLink *sl, double h, double ksi, double kBT, long int *in_ran1, int n_ch, int n_m, int istep);
void move_sl(Bead *bead, SlipLink *sl, double b, double h, double kBT, double ksi, double Ns, int n_ch, int n_m, long int *in_ran1, int zzz, int *destroy);
void force_rouse(Bead *bead, double b, double kBT, int n_ch, int n_m, int istep);
void force_sliplink(Bead *bead, SlipLink *sl, double b, double kBT, int n_ch, int n_m, double Ns, int zzz, int istep);
void sl_analysis(Bead *bead, double ***XYZ0, double *MSD, double **COM0, double *P2, double *Rbmean, double ***Rb0, double *endtoend, double *Ree2, double *DRee, double **Ree0, int n_ch, int n_m, int istep);
void update_sl(Bead *bead, SlipLink *sl, int zzz, int n_m, int n_ch, int *num_SL, int istep);


double ran1(long *idum);
