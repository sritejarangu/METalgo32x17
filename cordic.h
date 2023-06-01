#ifndef _CORDIC_H_
#define _CORDIC_H_
#include "ap_fixed.h"

typedef unsigned int UINTYPE_12;
/*typedef ap_fixed<32,16>  THETA_TYPE;
typedef ap_fixed<32,16> COS_SIN_TYPE;*/
typedef ap_fixed<32,25,AP_RND_INF>  THETA_TYPE;
typedef ap_fixed<32,25,AP_RND_INF> COS_SIN_TYPE;
const ap_uint<5> NUM_ITERATIONS=8;
static  THETA_TYPE cordic_phase[10]={0.78539816339744828000,0.46364760900080609000,
 		                             0.24497866312686414000,0.12435499454676144000,0.06241880999595735000,
		                             0.03123983343026827700,0.01562372862047683100,0.00781234106010111110,
    	                             0.00390623013196697180,0.00195312251647881880
		                             };

void cordic (THETA_TYPE theta, COS_SIN_TYPE  &s, COS_SIN_TYPE  &c);
#endif
