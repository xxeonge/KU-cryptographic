#pragma once
#include <stdio.h>
#include <stdlib.h>

#define __GMP_ENABLE

#ifdef __GMP_ENABLE
#include "gmp.h"
#endif

#define ECC_PASS  1
#define ECC_FAIL -1

#define ECC_P256_WORD_NUM 20

typedef struct _ECC_BN
{
	unsigned int dat[ECC_P256_WORD_NUM];
	unsigned int len;
}ECC_BN;


// FFFFFFFF 00000001 000000000000000000000000 FFFFFFFF FFFFFFFF FFFFFFFF
static ECC_BN prime_p256 = {
  {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000,
   0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF}, 8
};


#ifdef __GMP_ENABLE
int ECC_ecc_bn_to_mpz(mpz_t c, ECC_BN* a);
int ECC_mpz_to_ecc_bn(ECC_BN* c, mpz_t a);
#endif

int ECC_bn_cpy(ECC_BN* c, ECC_BN* a);
int ECC_bn_cmp(ECC_BN* a, ECC_BN* b);

// 
int ECC_bn_add(ECC_BN* c, ECC_BN* a, ECC_BN* b);
int ECC_bn_sub(ECC_BN* c, ECC_BN* a, ECC_BN* b);
int ECC_bn_mul(ECC_BN* c, ECC_BN* a, ECC_BN* b);

//int ECC_bn_mod_p256(ECC_BN* c, ECC_BN* a, ECC_BN* p);

// 10/30
int ECC_bn_mod_p256(ECC_BN* c, ECC_BN* a);

int ECC_bn_add_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p);
int ECC_bn_sub_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p);