#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <sys/time.h>
#include <assert.h>
#include "pcs_elliptic_curve_operations.h"
#include "pcs_pollard_rho.h"
#include "pcs_storage.h"
#include "pcs.h"
#include "pcs_multi_user.h"

#define FF fflush(stdout)
#define verbose 0
//#define SEED (time(null))
//#define SEED 0xE5CA1ADE

elliptic_curve_t E;
point_t P;
point_t Q[__NB_USERS__]; // One Q per user
mpz_t n;
mpz_t *A;
mpz_t *B;
mpz_t X_res[__NB_USERS__]; // One x per user
point_t M[__NB_ENSEMBLES__]; // precomputed a*P values
uint8_t trailling_bits;
uint8_t nb_bits;

/** Determines whether a point is a distinguished one.
 *
 *  @param[in]	R				A point on an elliptic curve.
 *  @param[in]	trailling_bits	Number of trailling zero bits in a ditinguished point.
 *  @param[out]	q				The x-coordinate, without the trailling zeros.
 *  @return 	1 if the point is distinguished, 0 otherwise.
 */
int is_distinguished_mu(point_t R, int trailling_bits, mpz_t *q)
{
	int res;
	mpz_t r;
	mpz_inits(r, NULL);
	mpz_tdiv_qr_ui(*q, r, R.x, (unsigned long int)pow(2, trailling_bits));
	res=(mpz_sgn(r) == 0);
	mpz_clears(r, NULL);
	return (res);
}


/** Checks if the linear combination aP+bQ is equal to R or its inverse.
 *
 *  @param[in]	R	A point on an elliptic curve.
 *  @param[in]	a	a coefficient.
 *  @param[in]	b	b coefficient.
 *  @param[in]  user    user id
 *  @return 	1 if aP+bQ[user] = R, 0 if aP+bQ[user] = -R.
 */
int same_point_mu(point_t R, mpz_t a, mpz_t b, uint16_t user)
{
	int res;
	point_t S1, S2, S;
	mpz_inits(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
	double_and_add(&S1, P, a, E); // S1 = aP
	double_and_add(&S2, Q[user], b, E); // S2 = bQ[user]
	add(&S, S1, S2, E); // S2 = aP + bQ
	res=(mpz_cmp(R.y, S.y) == 0); //
	mpz_clears(S1.x, S1.y, S1.z, S2.x, S2.y, S2.z, S.x, S.y, S.z, NULL);
	return res;
}

/** Computes R = a * P  on E.
 *
 *  @param[out]	R	Resulting point.
 *  @param[in]	a	a coefficient.
 */
void lin_comb_mu(point_t * R, mpz_t a)
{
	double_and_add(R, P, a, E);
}

/** Checks if there is a collision.
 *
 */
//int is_collision(mpz_t x, mpz_t a1, mpz_t a2, int trailling_bits)
int is_collision_mu(mpz_t x, mpz_t b1, uint16_t userid1, mpz_t b2, uint16_t userid2, int trailling_bits, mpz_t x_true1, mpz_t x_true2)
{
	uint8_t r;
	mpz_t xDist_;
	int retval = 0;
	mpz_t a1, a2;
	point_t R;
	// mpz_t xR, xI; // debug
	
	
        
	point_init(&R);
	mpz_inits(a1, a2, xDist_, NULL);
	
	mpz_set_ui(a2, 0); // a1 = a2 = 0
	mpz_set_ui(a1, 0);

	if (userid1!=userid2){
	  if (mpz_cmp(x_true1,x_true2)==0)
	    {
	      printf("même clé : %d,%d\n",userid1,userid2);FF;
	    }
	}

	
	//recompute first a,b pair
	
	double_and_add(&R, Q[userid1], b1, E); // R = b1 * Qi
	
	while(!is_distinguished_mu(R, trailling_bits, &xDist_))
	{
		r = hash(R.y);
		compute_a(a1, A[r], n); // a1 = a1 + A[r] % n   <=> ajout de A[r] au coeff a1
		f(R, M[r], &R, E);      // R  = R  + M[r]       <=> calcul du point R = R + M[r] = R + A[r]*P
	}

	
	//recompute second a,b pair
	
	double_and_add(&R, Q[userid2], b2, E);
	
	while(!is_distinguished_mu(R, trailling_bits, &xDist_))
	{
		r = hash(R.y);
		compute_a(a2, A[r], n); 
		f(R, M[r], &R, E);
	}



	// debug
	/*
	if (userid1==2771)
	  {
	    gmp_printf("a1 = %-10Zd, b1 = %-10Zd, a2 = %-10Zd, b2 = %-10Zd, x1 = %-10Zd, x2 = %-10Zd\n",a1,b1,a2,b2,x_true1,x_true2);
	    gmp_printf("n = %Zd\n",n);

	    mpz_inits(xR, xI, NULL);
        
	    mpz_mul(xR,b2,x_true2);  // b2*x2
	    mpz_mmod(xR,xR,n);
	    gmp_printf("xR = %Zd\n",xR);
	    mpz_add(xR,xR,a2);  // +a2
	    mpz_mmod(xR,xR,n);
	    gmp_printf("xR = %Zd\n",xR);
	    mpz_sub(xR,xR,a1);  // -a1
	    mpz_mmod(xR,xR,n);
	    gmp_printf("xR = %Zd\n",xR);
	    mpz_invert(xI,b1,n);// b1^-1 mod n
	    gmp_printf("xI = %Zd\n",xI);
	    mpz_mul(xR,xR,xI);  // (b2*x2+a2-a1) * (b1^-1)
	    mpz_mmod(xR,xR,n);
	    gmp_printf("x = %Zd\n",xR);
	    
	
	    mpz_clears(xR,xI, NULL);
	  }
	
	// end debug
	*/
	if(userid1==userid2 && (mpz_cmp(b1, b2) != 0)) //two different pairs with the same Q, so collision
	{
          if(!same_point_mu(R, a1, b1,userid1)) //it's the inverse point // to be modified
		{	
			mpz_neg(a2, a2); 
			mpz_mmod(a2, a2, n);
			mpz_neg(b2, b2);
			mpz_mmod(b2, b2, n);
		}
		compute_x(x, a1, a2, b1, b2, n);
		retval = 1;
	}
        else if(userid1!=userid2) //two different Qs, so collision - userid2 has to be known
        {
          if(!same_point_mu(R, a1, b1,userid1)) //it's the inverse point // to be modified
		{	
			mpz_neg(a2, a2); 
			mpz_mmod(a2, a2, n);
			mpz_neg(b2, b2);
			mpz_mmod(b2, b2, n);
		}
          compute_x_2users(x, a1, a2, b1, b2, X_res[userid2], n);
          retval = 1;
        }
	point_clear(&R);
	mpz_clears(a1, a2, xDist_, NULL);
	return retval;
}

/** Initialize all variables needed to do a PCS algorithm.
 *
 */
void pcs_mu_init(point_t  P_init,
                 point_t Q_init[__NB_USERS__],
                 elliptic_curve_t E_init,
                 mpz_t n_init,
                 mpz_t *A_init,
                 uint8_t nb_bits_init,
                 uint8_t trailling_bits_init,
                 int type_struct,
                 int nb_threads,
                 uint8_t level)
{
  uint8_t i; //  __NB_ENSEMBLES__
  int j; // __NB_USERS__
  uint16_t user;

        
	point_init(&P);
	//point_init(&Q);
	curve_init(&E);
	mpz_init(n);
        
        
	mpz_set(P.x, P_init.x);
	mpz_set(P.y, P_init.y);
	mpz_set(P.z, P_init.z);
	
	//mpz_set(Q.x, Q_init.x);
	//mpz_set(Q.y, Q_init.y);
	//mpz_set(Q.z, Q_init.z);
	
	mpz_set(E.A, E_init.A);
	mpz_set(E.B, E_init.B);
	mpz_set(E.p, E_init.p);
	
	mpz_set(n, n_init);
	
	A = A_init;
        
        

        for(j=0; j<__NB_USERS__; j++) // Q init
          {
            mpz_set(Q[j].x,Q_init[j].x); 
            mpz_set(Q[j].y,Q_init[j].y);
            mpz_set(Q[j].z,Q_init[j].z);
            user = (uint16_t) j;
            
          }
	for(i=0; i<__NB_ENSEMBLES__; i++) // has to be after Q inits at index j since it is needed for lin_comb_mu
	  {
	    mpz_inits(M[i].x,M[i].y,M[i].z,NULL);
	    lin_comb_mu(&M[i],A[i]);
	  }
	
	trailling_bits = trailling_bits_init;
	nb_bits = nb_bits_init;
	
	struct_init_mu(type_struct, n, trailling_bits, nb_bits, nb_threads, level); // TODO : mu adaptation : done?
}





/** Run the PCS algorithm.
 *
 */
long long int pcs_mu_run_order(mpz_t x_res[__NB_USERS__], int nb_threads, unsigned long long int times[__NB_USERS__],unsigned long int pts_per_users[__NB_USERS__],mpz_t x_true[__NB_USERS__])
{
  point_t R;
  mpz_t b, b2;
  mpz_t x, xDist;
  uint8_t r;
  int trail_length;
  int col;
  int trail_length_max = pow(2, trailling_bits) * 20; // 20 * 1<<trailling_bits
  int collision_count = 0;
  unsigned long long int time1,time2;
  // int userid1;
  // int* userid2;
  uint16_t userid1,userid2;
  char xDist_str[50];
  struct timeval tv1,tv2;
  unsigned long int nb_pts;
  //unsigned long long int times[__NB_USERS__];
  //nb_threads = 1;
  //nb_threads = omp_get_max_threads();
  //nb_threads = __NB_USERS__; // for testing purposes : userid1 = thread number
  
  for (userid1=0; userid1<__NB_USERS__; userid1++)
    {
      pts_per_users[userid1] = 0;
      gettimeofday(&tv1,NULL);
#pragma omp parallel private(nb_pts,userid2, R, b, b2, x, r, xDist, xDist_str, trail_length,col) shared(collision_count, X_res, trail_length_max,pts_per_users) num_threads(nb_threads)
      {
	nb_pts = 0;
	col = 0;
	point_init(&R);
	mpz_inits(x, b, b2, xDist, NULL);
	
	
	//Initialize a starting point
	gmp_randstate_t r_state;
	gmp_randinit_default(r_state);
	gmp_randseed_ui(r_state, SEED * (userid1+1) * (omp_get_thread_num() + 1));
	mpz_urandomb(b, r_state, nb_bits); // random b
	double_and_add(&R, Q[userid1], b, E); // R = bQi
	trail_length = 0;
	collision_count = 0;
	while(collision_count < 1) 
	  {
	    if(is_distinguished_mu(R, trailling_bits, &xDist)) // xDist = R.x >> trailling_bits
	      {
		//printf("dist point!\n");fflush(stdout);
		
		userid2 = __NB_USERS__; // debug / useless
		nb_pts++;
		//printf(".");FF;
		
		if(struct_add_mu(b2, &userid2, b, userid1, xDist, xDist_str)) // ajout de b dans la mémoire, b2 = b d'un autre point avec collision 
		  {
		    //printf("added\n)");fflush(stdout);
		    if(is_collision_mu(x, b, userid1, b2, userid2, trailling_bits,x_true[userid1],x_true[userid2])) // si b et b2 forment une vraie collision
		      {
                      //printf("\nThread num %d :\n",omp_get_thread_num());
                      //printf("True collision %2hu - %2hu",userid1,userid2);
			if(verbose)
			  {
			    printf("coll(%hu-%hu);",userid1,userid2);FF;
			  }
			  //if(userid1!=userid2) printf(" ---- different origin");
			col = 1;
			//printf("\n");
			
                        #pragma omp critical
			{
			  collision_count++;
			  mpz_init_set(X_res[userid1],x);
			}
		      }
		  }
		//              else // pt distingué ajouté, pas de collision
		if(col==0)
		  { 
		    mpz_urandomb(b, r_state, nb_bits);
		    double_and_add(&R, Q[userid1], b, E); // new start, R = bQi
		    trail_length = 0;
		  }
	      }
	    else // R n'est pas un pt dist.
	      {

            
		r=hash(R.y); // y%20
		//printf("f...\n");FF;
		//printf("%d,%d\n",userid1,r);FF;
		//gmp_printf("%Zd\n",M[userid1][r].x);FF;
		f(R, M[r], &R, E); // 1 step (among 20) of the path
		//gmp_printf("%Zd.",M[userid1][r].x);
		//gmp_printf("%Zd-",M[userid1][r].y);FF;
		//printf("f - ok\n");FF;
              
		trail_length++;
		if(trail_length > trail_length_max)
		  {
		    mpz_urandomb(b, r_state, nb_bits); // new random start  
		    double_and_add(&R, Q[userid1], b, E);
		    trail_length = 0;
		  }

	      }
          
	  } // end while
	//printf("thread %d : %lu\n",omp_get_thread_num(),nb_pts);
	#pragma omp critical
	{
	  pts_per_users[userid1]+=nb_pts;
	}
	
	point_clear(&R);
	mpz_clears(b, b2, x, xDist, NULL);
	gmp_randclear(r_state);
      } // end omp parallel
      
      //printf("user %d : %lu pts\n",userid1,pts_per_users[userid1]); 
      gettimeofday(&tv2,NULL);
      // times.append(tv2-tv1)
      time1 = (tv1.tv_sec) * 1000000 + tv1.tv_usec;
      time2 = (tv2.tv_sec) * 1000000 + tv2.tv_usec;
      times[userid1] = time2 - time1;
      if (verbose && !userid1%((int)(__NB_USERS__/100)))
	{ 
	  printf("#");FF;
	}

      
    } // end for
  if (verbose)
    printf("\n");
  for (userid1=0; userid1<__NB_USERS__; userid1++)
    {
      mpz_init_set(x_res[userid1],X_res[userid1]);
    }
  return 0;
}









/** Free all variables used in the previous PCS run.
 *
 */
void pcs_mu_clear()
{
  
  uint32_t i;
  point_clear(&P);
  /*
    for (i=0;i<__NB_USERS__;i++)
    {
    point_clear(&Q[i]);
    //mpz_clears(Q[i].x, Q[i].y, Q[i].z, NULL);
    }
    free(Q);
    printf("cleared\n");
    printf("clearing M... ");fflush(stdout);
  */
  for(i = 0; i < __NB_ENSEMBLES__; i++)
    {
      point_clear(&M[i]);
      //mpz_clears(M[i].x, M[i].y, M[i].z, NULL);
    }
    
  curve_clear(&E);
  mpz_clear(n);
  struct_free_mu();
}
