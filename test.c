#include <stdio.h>
#include <stdlib.h>
#include “ECC.h”
void test_add_sub()
{
    int i;
    mpz_t ma, mb, mc, mp;
    ECC_BN a, b, c, tmp;
    gmp_randstate_t state;
    mpz_init(ma); mpz_init(mb); mpz_init(mc); mpz_init(mp);
    gmp_randinit_default(state);
    ECC_ecc_bn_to_mpz(mp, &prime_p256);
    for (i = 0; i < 100000; i++) {
        mpz_urandomm(ma, state, mp);
        mpz_urandomm(mb, state, ma);
        mpz_add(mc, ma, mb);
        ECC_mpz_to_ecc_bn(&a, ma);
        ECC_mpz_to_ecc_bn(&b, mb);
        ECC_mpz_to_ecc_bn(&tmp, mc);
        ECC_bn_add(&c, &a, &b);
        if (ECC_bn_cmp(&c, &tmp)) {
            printf(“ecc_bn_add: FAIL %d \n”, i);
            return;
        }
        mpz_sub(mc, ma, mb);
        ECC_mpz_to_ecc_bn(&tmp, mc);
        ECC_bn_sub(&c, &a, &b);
        if (ECC_bn_cmp(&c, &tmp)) {
            printf(“ecc_bn_sub: FAIL %d\n”, i);
            return;
        }
    }
    printf(“ecc_bn_add & sub: PASS \n”);
    mpz_clear(ma); mpz_clear(mb); mpz_clear(mc); mpz_clear(mp);
}
void test_add_sub_mod()
{
    int i;
    mpz_t ma, mb, mc, mp;
    ECC_BN a, b, c, tmp;
    gmp_randstate_t state;
    mpz_init(ma); mpz_init(mb); mpz_init(mc); mpz_init(mp);
    gmp_randinit_default(state);
    ECC_ecc_bn_to_mpz(mp, &prime_p256);
    for (i = 0; i < 100000; i++) {
        mpz_urandomm(ma, state, mp);
        mpz_urandomm(mb, state, mp);
        mpz_add(mc, ma, mb);
        mpz_mod(mc, mc, mp);
        ECC_mpz_to_ecc_bn(&a, ma);
        ECC_mpz_to_ecc_bn(&b, mb);
        ECC_mpz_to_ecc_bn(&tmp, mc);
        ECC_bn_add_mod(&c, &a, &b, &prime_p256);
        if (ECC_bn_cmp(&c, &tmp)) {
            printf(“ecc_bn_add_mod: FAIL %d\n”, i);
            return;
        }
        mpz_sub(mc, ma, mb);
        mpz_mod(mc, mc, mp);
        ECC_mpz_to_ecc_bn(&tmp, mc);
        ECC_bn_sub_mod(&c, &a, &b, &prime_p256);
        if (ECC_bn_cmp(&c, &tmp)) {
            printf(“ecc_bn_sub_mod: FAIL %d\n”, i);
            return;
        }
    }
    printf(“ecc_bn_add & sub_mod: PASS \n”);
    mpz_clear(ma); mpz_clear(mb); mpz_clear(mc); mpz_clear(mp);
}
// a,b 길이는 가변--> 랜덤값 나누기 9  a.len =  a.dat[i]%9 10만번 반복
void test_mul()
{
    int i;
    mpz_t ma, mb, mc, mp;
    ECC_BN a, b, c, tmp;
    gmp_randstate_t state;
    mpz_init(ma); mpz_init(mb); mpz_init(mc); mpz_init(mp);
    gmp_randinit_default(state);
    ECC_ecc_bn_to_mpz(mp, &prime_p256);

    for (i = 0; i < 1000000; i++) {
        mpz_urandomm(ma, state, mp);
        mpz_urandomm(mb, state, mp);
        mpz_mul(mc, ma, mb);
        ECC_mpz_to_ecc_bn(&a, ma);
        ECC_mpz_to_ecc_bn(&b, mb);
        ECC_mpz_to_ecc_bn(&tmp, mc);
        ECC_bn_mul(&c, &a, &b);
        if (ECC_bn_cmp(&c, &tmp)) {
            printf(“test_mul: FAIL at iteration %d\n”, i);
            return;
        }
    }
    mpz_clear(ma); mpz_clear(mb); mpz_clear(mc); mpz_clear(mp);
    gmp_randclear(state);
}

void test_mod() {

    int i;
    mpz_t ma, mb, mc, md, mp;

    ECC_BN a, b, c, d, tmp;

    gmp_randstate_t state;

    mpz_init(ma); mpz_init(mb); mpz_init(mc); mpz_init(md); mpz_init(mp);

    gmp_randinit_default(state);

    ECC_ecc_bn_to_mpz(mp, &prime_p256);

    for (i = 0; i < 100000; i++) {
        mpz_urandomm(ma, state, mp);
        mpz_urandomm(mb, state, mp);

        mpz_mul(mc, ma, mb);
        mpz_mod(md, mc, mp);

        ECC_mpz_to_ecc_bn(&a, ma);
        ECC_mpz_to_ecc_bn(&b, mb);

        ECC_bn_mul(&c, &a, &b);
        ECC_bn_mod_p256(&d, &c);
        ECC_mpz_to_ecc_bn(&tmp, md);

        if (ECC_bn_cmp(&d, &tmp)) {
            printf(“test_mod: FAIL  %d\n”, i);
            gmp_printf(“mpz: %Zx\n”, md);
            return;
        }
    }
    printf(“test_mod: PASS\n”);

    mpz_clear(ma); mpz_clear(mc); mpz_clear(md); mpz_clear(mb); mpz_clear(mp);

    gmp_randclear(state);
}

void test_binary_inv() {
    int i;
    mpz_t ma, mc, mp, one, product;
    ECC_BN a, c, tmp;
    gmp_randstate_t state;

    mpz_init(ma);
    mpz_init(mc);
    mpz_init(mp);
    mpz_init_set_ui(one, 1); 
    mpz_init(product);
    gmp_randinit_default(state);

    ECC_ecc_bn_to_mpz(mp, &prime_p256);

    for (i = 0; i < 100000; i++) {
        mpz_urandomm(ma, state, mp);
        

        ECC_mpz_to_ecc_bn(&a, ma);

        if (ECC_bn_binary_inv(&c, &a) != ECC_PASS) {
            printf("ECC_bn_binary_inv: FAIL - Could not find inverse %d\n", i);
            return;
        }

        ECC_ecc_bn_to_mpz(mc, &c);

        mpz_mul(product, ma, mc);
        mpz_mod(product, product, mp);

        if (mpz_cmp(product, one) != 0) {
            printf("ECC_bn_binary_inv: FAIL  %d\n", i);
            gmp_printf("ma: %Zx\n", ma);
            gmp_printf("mc (inv): %Zx\n", mc);
            gmp_printf("product: %Zx\n", product);
            return;
        }
    }

    printf("ECC_bn_binary_inv: PASS\n");

    mpz_clear(ma);
    mpz_clear(mc);
    mpz_clear(mp);
    mpz_clear(one);
    mpz_clear(product);
    gmp_randclear(state);
}

void test_ECC_bn_binary_inv() {
    int i;
    mpz_t ma, mc, mp, one, product;
    ECC_BN a, c, tmp;
    gmp_randstate_t state;

    mpz_init(ma);
    mpz_init(mc);
    mpz_init(mp);
    mpz_init_set_ui(one, 1);
    mpz_init(product);
    gmp_randinit_default(state);

    ECC_ecc_bn_to_mpz(mp, &prime_p256);

    for (i = 0; i < 100000; i++) {
        mpz_urandomm(ma, state, mp);
        
        ECC_mpz_to_ecc_bn(&a, ma);

        if (ECC_bn_binary_inv(&c, &a) != ECC_PASS) {
            printf("ECC_bn_binary_inv: FAIL - Could not find inverse %d\n", i);
            return;
        }

        ECC_ecc_bn_to_mpz(mc, &c);

        mpz_mul(product, ma, mc);
        mpz_mod(product, product, mp);

        if (mpz_cmp(product, one) != 0) {
            printf("ECC_bn_binary_inv: FAIL %d\n", i);
            gmp_printf("ma: %Zx\n", ma);
            gmp_printf("mc (inv): %Zx\n", mc);
            gmp_printf("product: %Zx\n", product);
            return;
        }
    }

    printf("ECC_bn_binary_inv: PASS\n");

    mpz_clear(ma);
    mpz_clear(mc);
    mpz_clear(mp);
    mpz_clear(one);
    mpz_clear(product);
    gmp_randclear(state);
}


void main()
{
    // test_add_sub();
    // test_add_sub_mod();
    // test_mul();
    // test_mod();
    test_ECC_bn_binary_inv();
}