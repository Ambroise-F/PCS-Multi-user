
#include <gmp.h>
#include <omp.h>
#include <inttypes.h>
#include "pcs.h"

#define __NB_USERS__ 65535
//#define __NB_USERS__ 1


//#define SEED 0xE5CA1ADE
#define SEED (time(NULL))
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
			       int nb_threads,
			       unsigned long long int times[__NB_USERS__],
			       unsigned long int pts_per_users[__NB_USERS__],
			       mpz_t x_true[__NB_USERS__]);




void pcs_mu_clear();
