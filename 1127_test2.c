#include <stdio.h>
#include <stdlib.h>
#include "ecc.h"

void test()
{
    ECC_PT P, R, Q, result_sum;
    
    ECC_pt_cpy(&P, &base_p256);
    ECC_pt_dbl(&R, &P);
    ECC_pt_cpy(&Q, &P);
    ECC_bn_sub_mod(&Q.y, &prime_p256, &P.y, &prime_p256);
    ECC_pt_add(&result_sum, &R, &Q);
  
    if (ECC_bn_cmp(&result_sum.x, &P.x) || ECC_bn_cmp(&result_sum.y, &P.y) ||
        result_sum.point_at_infinity != P.point_at_infinity) {
        printf("test failed: R + Q != P\n");
        printf("R = 2P, Q = -P\n");
        return;
    }
    printf("test: PASS\n");
}
void main()
{
    test();
}