// #include <stdio.h>
// #include <stdlib.h>
// #include "gmp.h"

// void main(void) {

//     mpz_t a, b, c; //아무것도 할당안됨 alloc, size, d

//     mpz_init(a); // 메모리 할당 - init 지나가면 alloc 1, size 0, d는 쓰레기 값으로 차있음
//     mpz_init(b);
//     mpz_init(c);

//     // c = a * b
//     a-> _mp_d[0] = 0xffffffff;
//     a-> _mpz_size = 1;

//     b-> _mp_d[0] = 0xfffffffe;
//     b-> _mp_size = 1;

//     mpz_mul(c, a, b);

//     mpz_clear(a); // 동적할당을 하면 free를 항상 해줘야 함 (메모리 누수 방지)
//     mpz_clear(b);
//     mpz_clear(c);

//     return;
// }

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void main(void) {
    mpz_t n, d, q, r, cmp_n;
    gmp_randstate_t state;
    
    // mpz 자료형 초기화
    mpz_init(n);
    mpz_init(d);
    mpz_init(q);
    mpz_init(r);
    mpz_init(cmp_n);
    
    // 난수 생성을 위한 시드 생성
    gmp_randinit_default(state);
    gmp_randseed_ui(state, time(NULL));
    
    // 2048비트 난수 n 생성
    mpz_urandomb(n, state, 2048);
    
    // 1024비트 난수 d 생성
    mpz_urandomb(d, state, 1024);
    
    // n을 d로 나누어 몫 q와 나머지 r을 구함
    mpz_fdiv_qr(q, r, n, d);
    
    // cmp_n = q*d + r
    mpz_mul(cmp_n, q, d);
    mpz_add(cmp_n, cmp_n, r);
    
    // 계산에 사용된 모든 데이터 출력
    printf("n = ");
    mpz_out_str(stdout, 10, n);
    
    printf("\nd = ");
    mpz_out_str(stdout, 10, d);
    
    printf("\nq = ");
    mpz_out_str(stdout, 10, q);
    
    printf("\nr = ");
    mpz_out_str(stdout, 10, r);

    printf("\ncmp_n = ");
    mpz_out_str(stdout, 10, cmp_n);
    printf("\n\n");
    
    // n과 cmp_n 비교
    if (mpz_cmp(n, cmp_n) == 0) {
        printf("n과 cmp_n은 같습니다.\n");
    } else {
        printf("n과 cmp_n은 다릅니다.\n");
    }
    
    // 자료형 해제
    mpz_clear(n);
    mpz_clear(d);
    mpz_clear(q);
    mpz_clear(r);
    mpz_clear(cmp_n);
    gmp_randclear(state);
    
    return;
}
