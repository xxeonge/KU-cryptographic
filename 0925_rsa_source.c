#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"


void main(void) {
    
    mpz_t p, q, n, e, d;
    gmp_randstate_t state; 

    // 변수 생성
    mpz_init(p);
    mpz_init(q); 
    mpz_init(n); 
    mpz_init(e); 
    mpz_init(d);  
    gmp_randinit_default(state); 
        
    // 키생성
    // 1. e=65337, 2^16+1=0x10001
    e-> _mp_d[0] = 0x1001;
    e-> _mp_size = 1;

    mpz_set_ui(e, 0x1001);

    // 2. p, q는 1024 비트 랜덤 소수
    while (1) {
        mpz_urandomb(p, state, 1024);
        if (p-> _mp_size == 32) {
            p-> _mp_d[31] = _mp_d[31] | 0x80000000;
            p-> _mp_d[ 0] = _mp_d[ 0] | 0x00000001;

            // p-1 sub_ui 사용 가능
            p-> _mp_d[ 0] = _mp_d[ 0] | 0xfffffffe;
            mpx_gcd(n, p, e);

            // 값 원상복귀
            p-> _mp_d[ 0] = _mp_d[ 0] | 0x00000001;

            //if (!mpz_cmp_ui(n, 1))
            if( n-> _mp_size == 1) && ( n->_mp_d[0]==1) {
                if (mpz_probab_prime_p(p, 56)) {
                    break;
                }
            }
        }
    }

    // 3. n = pq 계산
    mpz_mul(n, p, q);

    // 4. d = e^(-1) mod phi(n) 계산, {phi(n)=(p-1)(q-1)}
    mpz_t phi_n;
    mpz_init(phi_n);

    mpz_sub_ui(p, p, 1);  // p - 1
    mpz_sub_ui(q, q, 1);  // q - 1
    mpz_mul(phi_n, p, q); 

    if (mpz_invert(d, e, phi_n) == 0) {
        printf("역원 계산 실패\n");
        return;
    }

    // 메세지 선택 : m
    mpz_t m;
    mpz_init(m); 
    mpz_set_ui(m, 42);

    // 암호화: m^e mod n = c
    mpz_t c;
    mpz_init(c);
    mpz_powm(c, m, e, n);

    // 복호화: c^d mod n = out
    mpz_t out;
    mpz_init(out);
    mpz_powm(out, c, d, n);

    // 테스트: m?=output 출력
    if (mpz_cmp(m, out) == 0) {
        printf("성공: 원본 메시지와 복호화된 메시지가 일치합니다.\n");
    } else {
        printf("실패: 원본 메시지와 복호화된 메시지가 일치하지 않습니다.\n");
    }

    // 변수 제거
    mpz_clear(p);
    mpz_clear(q);
    mpz_clear(n);
    mpz_clear(e);
    mpz_clear(d);
    gmp_randclear(state);
}