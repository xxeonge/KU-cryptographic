#include "ECC.h"


#ifdef __GMP_ENABLE
int ECC_ecc_bn_to_mpz(mpz_t c, ECC_BN* a)
{
	int i;

	if (c->_mp_alloc < a->len) {
		mpz_realloc2(c, a->len << 5);
	}

	for (i = 0; i < a->len; i++) {
		c->_mp_d[i] = a->dat[i];
	}

	c->_mp_size = a->len;

	return ECC_PASS;
}

int ECC_mpz_to_ecc_bn(ECC_BN* c, mpz_t a)
{
	int i;

	if (a->_mp_size < 0)
		return ECC_FAIL;

	if (a->_mp_size > ECC_P256_WORD_NUM)
		return ECC_FAIL;

	for (i = 0; i < a->_mp_size; i++) {
		c->dat[i] = a->_mp_d[i];
	}

	c->len = a->_mp_size;

	return ECC_PASS;
}
#endif

int ECC_bn_cpy(ECC_BN* c, ECC_BN* a)
{
	int i;

	for (i = 0; i < a->len; i++) {
		c->dat[i] = a->dat[i];
	}

	c->len = a->len;

	return ECC_PASS;
}

// a>b => 1 , a<b => -1 , a=b => 0
int ECC_bn_cmp(ECC_BN* a, ECC_BN* b)
{
	int i;

	if (a->len > b->len)
		return 1;

	if (a->len < b->len)
		return -1;

	for (i = a->len - 1; i >= 0; i--) {
		if (a->dat[i] > b->dat[i]) {
			return 1;
		}
		else {
			if (a->dat[i] < b->dat[i]) {
				return -1;
			}
		}
	}

	return 0;
}


// len(a)>=len(b), a,b>=0
int ECC_bn_add(ECC_BN* c, ECC_BN* a, ECC_BN* b)
{
	int i, carry = 0;
	ECC_BN out;

	for (i = 0; i < b->len; i++) {
		out.dat[i] = a->dat[i] + b->dat[i] + carry;
		if (carry)
			carry = a->dat[i] >= (~b->dat[i]);
		else
			carry = a->dat[i] > (~b->dat[i]);
	}

	for (; i < a->len; i++) {
		out.dat[i] = a->dat[i] + carry;
		carry = out.dat[i] < carry;
	}
	out.dat[i] = carry;

	if (carry)
		out.len = a->len + 1;
	else
		out.len = a->len;

	ECC_bn_cpy(c, &out);

	return ECC_PASS;
}


// a>=b, a,b>=0
int ECC_bn_sub(ECC_BN* c, ECC_BN* a, ECC_BN* b)
{
	int i, borrow = 0;
	ECC_BN out;

	for (i = 0; i < b->len; i++) {
		out.dat[i] = a->dat[i] - b->dat[i] - borrow;
		borrow = (borrow) ? (a->dat[i] <= b->dat[i]) : (a->dat[i] < b->dat[i]);
	}

	for (; i < a->len; i++) {
		out.dat[i] = a->dat[i] - borrow;
		borrow = a->dat[i] < borrow;
	}

	out.len = a->len;

	while ((out.len > 0) && (out.dat[out.len - 1] == 0))
	{
		out.len--;
	}

	ECC_bn_cpy(c, &out);

	return ECC_PASS;
}

// c = a * b (a, b in GF(p))
int ECC_bn_mul(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p) 
{
    int                     i, j;
    unsigned long long int carry=0;
    ECC_BN                   out;

    // a, b != 0을 의미
    if ((a->len==0) || (b->len == 0)) {
        c->len = 0;
        return ECC_PASS;
    }

    for (i=0; i<ECC_P256_WORD_NUM; i++)
        out.dat[i]=0;

    for(i=0; i< b->len; i++)
    {
        carry=0;
        for (j=0; j< a->len; j++){
            carry = (unsigned long long int)a->dat[j] * b->dat[i] + (unsigned long long int)out.dat[i+j]
            + (unsigned long long int)carry;

            out.dat[i+j] = (unsigned long int)carry;
            carry = (unsigned long long int)carry >> 32;
        }
        out.dat[i+j] = (unsigned long int)carry;
    }

    // 길이 결정
	out.len = i + j;
	while ((out.len > 0) && (out.dat[out.len - 1] == 0)) {
		out.len--;
	}
	if (out.len == 0) out.len = 1;

    // 결과 cpy
    ECC_bn_cpy(c, &out);

    return ECC_PASS;
}


// a,b in GF(p)
int ECC_bn_add_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p)
{
	ECC_BN out;

	if ((a->len < 0) || (b->len < 0))
		return ECC_FAIL;

	if (a->len >= b->len) ECC_bn_add(&out, a, b);
	else ECC_bn_add(&out, b, a);

	if (ECC_bn_cmp(&out, p) >= 0)
	{
		ECC_bn_sub(c, &out, p);
	}
	else ECC_bn_cpy(c, &out);

	return ECC_PASS;

}

// 10/30 

int ECC_bn_mod_p256(ECC_BN* c, ECC_BN* a) {
    int i;
    ECC_BN in, out, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, result;

    for (i = 0; i < a->len; i++)
        in.dat[i] = a->dat[i];

    for (; i < ECC_P256_WORD_NUM; i++)
        in.dat[i] = 0;

	// 0으로 초기화 (memset or for문)

	tmp1 = tmp2 = tmp3 = tmp4 = tmp5 = tmp6 = tmp7 = tmp8 = tmp9 = {0};

    // tmp1 = s2 = (c15, c14, c13, c12, c11, 0, 0, 0)
    // tmp2 = s3 = (0, c15, c14, c13, c12, 0, 0, 0)
    tmp1.dat[0] = tmp1.dat[1] = tmp1.dat[2] = 
    tmp2.dat[0] = tmp2.dat[1] = tmp2.dat[2] = tmp2.dat[7] = 0;
    
    tmp1.dat[3] = in.dat[11];
    tmp1.dat[4] = tmp2.dat[3] = in.dat[12];
    tmp1.dat[5] = tmp2.dat[4] = in.dat[13];
    tmp1.dat[6] = tmp2.dat[5] = in.dat[14];
    tmp1.dat[7] = tmp2.dat[6] = in.dat[15];

    tmp1.len = 8;
    tmp2.len = 8;

    while ((tmp1.len > 0) && (tmp1.dat[tmp1.len - 1] == 0)) 
        tmp1.len--;
    while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) 
        tmp2.len--;

    // out = tmp1 + tmp2
    ECC_bn_add_mod(&out, &tmp1, &tmp2, &prime_p256);


	//비트 처리
	ECC_bn_add_mod(&tmp1, &tmp1, &tmp1, &prime_p256);


    // for (i = out.len - 1; i > 0; i--)
    //     out.dat[i] = (out.dat[i] << 1) | (out.dat[i - 1] >> 31);
    // out.dat[0] = out.dat[0] << 1;

    // ----------------------------------------------------
    // s1 + 2s2 + 2s3 + s4 + s5 - s6 - s7 - s8 - s9 mod p256

    // s4 정의: (c15, c14, 0, 0, 0, c10, c9, c8)
    tmp3.dat[0] = in.dat[8];
    tmp3.dat[1] = in.dat[9];
    tmp3.dat[2] = in.dat[10];
    tmp3.dat[3] = in.dat[15];
    tmp3.dat[4] = in.dat[14];
    tmp3.dat[5] = tmp3.dat[6] = tmp3.dat[7] = 0;
    tmp3.len = 8;

    // s5 정의: (c8, c13, c15, c14, c13, c11, c10, c9)
    tmp4.dat[0] = in.dat[9];
    tmp4.dat[1] = in.dat[10];
    tmp4.dat[2] = in.dat[11];
    tmp4.dat[3] = in.dat[13];
    tmp4.dat[4] = in.dat[15];
    tmp4.dat[5] = in.dat[8];
    tmp4.dat[6] = in.dat[12];
    tmp4.dat[7] = in.dat[14];
    tmp4.len = 8;

    // s6 정의: (c10, c8, 0, 0, 0, c13, c12, c11)
    tmp5.dat[0] = in.dat[10];
    tmp5.dat[1] = in.dat[8];
    tmp5.dat[5] = in.dat[13];
    tmp5.dat[6] = in.dat[12];
    tmp5.dat[7] = in.dat[11];
    tmp5.len = 8;

    // s7 정의: (c11, c9, 0, 0, c15, c14, c13, c12)
    tmp6.dat[0] = in.dat[11];
    tmp6.dat[1] = in.dat[9];
    tmp6.dat[4] = in.dat[15];
    tmp6.dat[5] = in.dat[14];
    tmp6.dat[6] = in.dat[13];
    tmp6.dat[7] = in.dat[12];
    tmp6.len = 8;

    // s8 정의: (c12, 0, c10, c9, c8, c15, c14, c13)
    tmp7.dat[1] = in.dat[10];
    tmp7.dat[2] = in.dat[9];
    tmp7.dat[3] = in.dat[8];
    tmp7.dat[5] = in.dat[15];
    tmp7.dat[6] = in.dat[14];
    tmp7.dat[7] = in.dat[13];
    tmp7.len = 8;

    // s9 정의: (c13, 0, c11, c10, c9, 0, c15, c14)
    tmp8.dat[1] = in.dat[11];
    tmp8.dat[2] = in.dat[10];
    tmp8.dat[3] = in.dat[9];
    tmp8.dat[6] = in.dat[15];
    tmp8.dat[7] = in.dat[14];
    tmp8.len = 8;

    tmp9 = in;


	//길이보정 tmp별 길이
	for (i = result.len - 1; i > 0; i--) {
			if (result.dat[i] != 0) break;
			result.len--;
		}


    // result = s1 + 2s2 + 2s3 + s4 + s5 - s6 - s7 - s8 - s9 mod p256
    ECC_bn_add_mod(&result, &out, &tmp3, &prime_p256);  // + s4
    ECC_bn_add_mod(&result, &result, &tmp4, &prime_p256); // + s5
    ECC_bn_sub_mod(&result, &result, &tmp5, &prime_p256); // - s6
    ECC_bn_sub_mod(&result, &result, &tmp6, &prime_p256); // - s7
    ECC_bn_sub_mod(&result, &result, &tmp7, &prime_p256); // - s8
    ECC_bn_sub_mod(&result, &result, &tmp8, &prime_p256); // - s9
    ECC_bn_add_mod(&result, &result, &tmp9, &prime_p256); // + s1

    for (i = result.len - 1; i > 0; i--) {
        result.dat[i] = (result.dat[i] << 1) | (result.dat[i - 1] >> 31);
    }
    result.dat[0] = result.dat[0] << 1;

   ECC_bn_copy(c, &result);

    return ECC_PASS;
}









// a,b in GF(p)
int ECC_bn_sub_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p)
{
	ECC_BN out;

	if ((a->len < 0) || (b->len < 0))
		return ECC_FAIL;

	if (ECC_bn_cmp(a, b) >= 0) {
		ECC_bn_sub(c, a, b);
	}
	else {
		ECC_bn_sub(&out, b, a);
		ECC_bn_sub(c, p, &out);
		
	}

	return ECC_PASS;
}