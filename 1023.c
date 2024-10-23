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