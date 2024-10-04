#include <stdio.h>
#include <stdlib.h>
#include "gmp.h"


//ctrl+f5  : 실행
//ctrl+f10 : 디버깅모드 커서까지 실행 
//shift+f5 : 디버깅 종료


void main(void)
{
	mpz_t           p, q, n, e, d, phi_n, m, c, out;
	gmp_randstate_t state;

	// 변수 생성
	mpz_init(p); mpz_init(q); mpz_init(d); mpz_init(n); mpz_init(e);
	mpz_init(phi_n), mpz_init(m), mpz_init(c), mpz_init(out);
	gmp_randinit_default(state);

	// 키생성

	//1. e=65537 = 2^16+1=0x10001
	e->_mp_d[0] = 0x10001;
	e->_mp_size = 1;
	mpz_set_ui(e, 0x10001);

	//2. p,q 는 1024 비트 랜덤 소수
	while (1) {
		mpz_urandomb(p, state, 1024);
		if (p->_mp_size == 32) {
			p->_mp_d[31] = p->_mp_d[31] | 0x80000000;
			p->_mp_d[0] = p->_mp_d[0] | 0x00000001;

			//p-1 
			p->_mp_d[0] = p->_mp_d[0] & 0xfffffffe;
			mpz_gcd(n, p, e);
			p->_mp_d[0] = p->_mp_d[0] | 0x00000001;

			if ((n->_mp_size == 1) && (n->_mp_d[0] == 1)) {
				//if( !mpz_cmp_ui(n,1) ){

				if (mpz_probab_prime_p(p, 56)) {
					break;
				}
			}
		}
	}

	while (1) {
		mpz_urandomb(q, state, 1024);
		if (q->_mp_size == 32) {
			q->_mp_d[31] = q->_mp_d[31] | 0x80000000;
			q->_mp_d[0] = q->_mp_d[0] | 0x00000001;

			//q-1 
			q->_mp_d[0] = q->_mp_d[0] & 0xfffffffe;
			mpz_gcd(n, q, e);
			q->_mp_d[0] = q->_mp_d[0] | 0x00000001;

			if ((n->_mp_size == 1) && (n->_mp_d[0] == 1)) {
				//if( !mpz_cmp_ui(n,1) ){

				if (mpz_probab_prime_p(q, 56)) {
					break;
				}
			}
		}
	}

	//3. n = p*q 계산 , mpz_mul()
	mpz_mul(n, p, q);

	//4. d = e^(-1) mod phi(n) 계산 {phi(n)=(p-1)(q-1)}  mpz_invert(), mpz_sub_ui...
	mpz_sub_ui(p, p, 1);
	mpz_sub_ui(q, q, 1);
	mpz_mul(phi_n, p, q);

	// +1 add
	mpz_add_ui(p, p, 1);
	mpz_add_ui(q, q, 1);
	
	mpz_invert(d, e, phi_n);

	//메시지 선택 : m  , mpz_urandomb, mpz_urandomm
	mpz_urandomm(m, state, n);

	//암호화 : m^e mod n =c,  mpz_powm
	mpz_powm(c, m, e, n);

	//복호화 : c^d mod n = out,  mpz_powm
	mpz_powm(out, c, d, n);

	// 테스트 : m?=out 출력, mpz_cmp 
	if (mpz_cmp(m, out) == 0) {
		printf("성공: m==0\n");
	}
	else {
		printf("실패: m!=0\n");
	}

	//gmp_print("p = %Zx", p);
	//gmp_print("q = %Zxn", q);
	//gmp_print("qn = %Zx", n);

	// 변수 제거
	mpz_clear(p); mpz_clear(q); mpz_clear(d); mpz_clear(n); mpz_clear(e);
	mpz_clear(phi_n), mpz_clear(m), mpz_clear(c), mpz_clear(out), gmp_randclear(state);
	gmp_randclear(state);
}