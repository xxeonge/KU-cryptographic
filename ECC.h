#include <stdio.h>
#include <stdlib.h>

//#define __GMP_ENABLE

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

typedef struct _ECC_PT
{
  ECC_BN x;
  ECC_BN y;
  int point_at_infinity; //1인 경우 무한원점으로 정의
}ECC_PT;



// FFFFFFFF 00000001 00000000 00000000 00000000 FFFFFFFF FFFFFFFF FFFFFFFF
static ECC_BN prime_p256 = {
  {0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0x00000000,
   0x00000000, 0x00000000, 0x00000001, 0xFFFFFFFF}, 8
};

// base x : 6b17d1f2 e12c4247 f8bce6e5 63a440f2 77037d81 2deb33a0 f4a13945 d898c296
// base y : 4fe342e2 fe1a7f9b 8ee7eb4a 7c0f9e16 2bce3357 6b315ece cbb64068 37bf51f5
static ECC_PT base_p256 = {
  {{ 0xd898c296, 0xf4a13945, 0x2deb33a0, 0x77037d81, 0x63a440f2, 0xf8bce6e5, 0xe12c4247, 0x6b17d1f2 }, 8 },
  {{ 0x37bf51f5, 0xcbb64068, 0x6b315ece, 0x2bce3357, 0x7c0f9e16, 0x8ee7eb4a, 0xfe1a7f9b, 0x4fe342e2 }, 8 },
  0
};

// order = FFFFFFFF 00000000 FFFFFFFF FFFFFFFF BCE6FAAD A7179E84 F3B9CAC2 FC632551

#ifdef __GMP_ENABLE
int ECC_ecc_bn_to_mpz(mpz_t c, ECC_BN* a);
int ECC_mpz_to_ecc_bn(ECC_BN* c, mpz_t a);
#endif


//int ECC_bn_add(ECC_BN* c, ECC_BN* a, ECC_BN *b);
//int ECC_bn_sub(ECC_BN* c, ECC_BN* a, ECC_BN* b);
//int ECC_bn_mul(ECC_BN* c, ECC_BN* a, ECC_BN* b);
//int ECC_bn_mod(ECC_BN* c, ECC_BN* a, ECC_BN* p);
//int ECC_bn_mod_p256(ECC_BN* c, ECC_BN* a);

int ECC_bn_cpy(ECC_BN* c, ECC_BN* a);
int ECC_bn_cmp(ECC_BN* a, ECC_BN* b);
int ECC_bn_1bit_rshift (ECC_BN* c, ECC_BN* a);
int ECC_bn_add_mod     (ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p);
int ECC_bn_sub_mod     (ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p);
int ECC_bn_mul_mod     (ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p);
int ECC_bn_mul_mod_p256(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p);
int ECC_bn_binary_inv  (ECC_BN* c, ECC_BN* a, ECC_BN* p);


int ECC_pt_init(ECC_PT* P);
int ECC_pt_cpy(ECC_PT* R, ECC_PT* P);

// R = P + Q
int ECC_pt_add(ECC_PT* R, ECC_PT* P, ECC_PT* Q);
// R = 2P
int ECC_pt_dbl(ECC_PT* R, ECC_PT* P);