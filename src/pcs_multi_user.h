
#include <gmp.h>
#include <omp.h>
#include <inttypes.h>
#include "pcs.h"

#define __NB_ENSEMBLES_MU__ 1
#define __NB_USERS__ 5 // right now also number of threads so 8 = max for this pc

//void combLin(point_t * R, mpz_t a, mpz_t b);

void pcs_mu_init(point_t P_init,
                   point_t Q_init[__NB_USERS__],
                   elliptic_curve_t E_init,
                   mpz_t n_init,
		   mpz_t *A_init,
                   uint8_t nb_bits_init,
                   uint8_t trailling_bits_init,
                   int type_struct,
                   int nb_threads,
                   uint8_t level);

long long int pcs_mu_run(mpz_t x_res,
                           int nb_threads,
                           int nb_collisions);

long long int pcs_mu_run_order(mpz_t x_res[__NB_USERS__],
                           int nb_threads);




void pcs_mu_clear();
