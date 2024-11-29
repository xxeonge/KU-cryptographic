#include "ecc.h"


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

int ECC_bn_1bit_rshift(ECC_BN* c, ECC_BN* a)
{
  int i;
  
  if (a->len == 0) {
	c->len = 0;
	return ECC_PASS;
  }

  for (i = 0; i < (a->len - 1); i++) {
	c->dat[i] = (a->dat[i+1]<<31) | (a->dat[i]>>1);
  }
  c->dat[i] = a->dat[i] >> 1;

  if (c->dat[i])
	c->len = a->len;
  else
	c->len = a->len - 1;

  return ECC_PASS;
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
	borrow = (borrow) ? (a->dat[i]<=b->dat[i]) :(a->dat[i]<b->dat[i]);
  }

  for (; i < a->len; i++) {
	out.dat[i] = a->dat[i] - borrow;
	borrow = a->dat[i] < borrow;
  }

  out.len = a->len; 

  while ((out.len>0) && (out.dat[out.len-1]==0)) out.len--;

  ECC_bn_cpy(c, &out);

  return ECC_PASS;
}

// c = a * b (a,b in GF(p))
int ECC_bn_mul(ECC_BN* c, ECC_BN* a, ECC_BN* b)
{
  int                    i, j;
  unsigned long long int carry=0;
  ECC_BN                 out;

  if ((a->len == 0) || (b->len == 0)) {
	c->len = 0;
	return ECC_PASS;
  }

  for (i = 0; i < ECC_P256_WORD_NUM; i++)
	out.dat[i] = 0;

  for (i = 0; i < b->len; i++)
  {
	carry = 0;
	for (j = 0; j < a->len; j++) {
	  carry = (unsigned long long int)a->dat[j] * b->dat[i] + (unsigned long long int)out.dat[i + j]
		    + (unsigned long long int)carry;

	  out.dat[i + j] = (unsigned long int)carry;
	  carry = (unsigned long long int)carry >> 32;
	}
	out.dat[i + j] = (unsigned long int)carry;
  }
 
  // 길이 결정
  out.len = a->len + b->len;
  if (out.dat[out.len - 1] == 0)
	out.len--;

  // 결과 cpy
  ECC_bn_cpy(c, &out);

  return ECC_PASS;
}


int ECC_bn_mod_p256(ECC_BN* c, ECC_BN* a)
{
  int    i;
  ECC_BN in, tmp1, tmp2;

  for (i = 0; i < a->len; i++)
	in.dat[i] = a->dat[i];

  for (     ; i < ECC_P256_WORD_NUM; i++)
	in.dat[i] = 0;

  // 2(s2+s3)
  //                7    6    5    4   3   2  1  0
  // tmp1 = s2 = (c15, c14, c13, c12, c11, 0, 0, 0)
  // tmp2 = s3 = (0,   c15, c14, c13, c12, 0, 0, 0)
  tmp1.dat[0] = tmp1.dat[1] = tmp1.dat[2] = tmp2.dat[0] = tmp2.dat[1] = tmp2.dat[2] = tmp2.dat[7] = 0;
  tmp1.dat[3] = in.dat[11];
  tmp1.dat[4] = tmp2.dat[3] = in.dat[12];
  tmp1.dat[5] = tmp2.dat[4] = in.dat[13];
  tmp1.dat[6] = tmp2.dat[5] = in.dat[14];
  tmp1.dat[7] = tmp2.dat[6] = in.dat[15];
  tmp1.len = 8;
  tmp2.len = 7;

  while ((tmp1.len > 0) && (tmp1.dat[tmp1.len - 1] == 0)) tmp1.len--;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_add_mod(&tmp1, &tmp1, &tmp2, &prime_p256);
  ECC_bn_add_mod(&tmp1, &tmp1, &tmp1, &prime_p256);

  // 2(s2+s3) + s1
  //             7   6   5   4   3   2   1   0
  // in = s1 = (c7, c6, c5, c4, c3, c2, c1, c0)  
  if(in.len > 8)
	in.len = 8;
    
  ECC_bn_add_mod(&tmp1, &tmp1, &in, &prime_p256);

  // 2(s2+s3) + s1 + s4
  //                7    6  5  4  3   2    1   0
  // tmp2 = s4 = (c15, c14, 0, 0, 0, c10, c9, c8)  
  tmp2.dat[3] = tmp2.dat[4] = tmp2.dat[5] = 0;
  tmp2.dat[0] = in.dat[8];
  tmp2.dat[1] = in.dat[9];
  tmp2.dat[2] = in.dat[10];
  tmp2.dat[6] = in.dat[14];
  tmp2.dat[7] = in.dat[15];
  tmp2.len = 8;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_add_mod(&tmp1, &tmp1, &tmp2, &prime_p256);

  // 2(s2+s3) + s1 + s4 + s5
  //               7    6    5    4    3    2    1   0
  // tmp2 = s5 = (c8, c13, c15, c14, c13, c11, c10, c9)  
  tmp2.dat[0] = in.dat[ 9];
  tmp2.dat[1] = in.dat[10];
  tmp2.dat[2] = in.dat[11];
  tmp2.dat[3] = in.dat[13];
  tmp2.dat[4] = in.dat[14];
  tmp2.dat[5] = in.dat[15];
  tmp2.dat[6] = in.dat[13];
  tmp2.dat[7] = in.dat[ 8];
  tmp2.len = 8;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_add_mod(&tmp1, &tmp1, &tmp2, &prime_p256);
  
  // 2(s2+s3) + s1 + s4 + s5 - s6
  //                7   6  5  4  3    2    1    0
  // tmp2 = s6 = (c10, c8, 0, 0, 0, c13, c12, c11)  
  tmp2.dat[3] = tmp2.dat[4] = tmp2.dat[5] = 0;
  tmp2.dat[0] = in.dat[11];
  tmp2.dat[1] = in.dat[12];
  tmp2.dat[2] = in.dat[13];
  tmp2.dat[6] = in.dat[ 8];
  tmp2.dat[7] = in.dat[10];
  tmp2.len = 8;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_sub_mod(&tmp1, &tmp1, &tmp2, &prime_p256);

  // 2(s2+s3) + s1 + s4 + s5 - s6 - s7
  //                7   6  5  4    3    2    1    0
  // tmp2 = s7 = (c11, c9, 0, 0, c15, c14, c13, c12)  
  tmp2.dat[4] = tmp2.dat[5] = 0;
  tmp2.dat[0] = in.dat[12];
  tmp2.dat[1] = in.dat[13];
  tmp2.dat[2] = in.dat[14];
  tmp2.dat[3] = in.dat[15];
  tmp2.dat[6] = in.dat[ 9];
  tmp2.dat[7] = in.dat[11];
  tmp2.len = 8;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_sub_mod(&tmp1, &tmp1, &tmp2, &prime_p256);

  // 2(s2+s3) + s1 + s4 + s5 - s6 - s7 - s8
  //                7   6   5   4   3    2    1    0
  // tmp2 = s8 = (c12, 0, c10, c9, c8, c15, c14, c13)  
  tmp2.dat[0] = in.dat[13];
  tmp2.dat[1] = in.dat[14];
  tmp2.dat[2] = in.dat[15];
  tmp2.dat[3] = in.dat[ 8];
  tmp2.dat[4] = in.dat[ 9];
  tmp2.dat[5] = in.dat[10];
  tmp2.dat[6] = 0;
  tmp2.dat[7] = in.dat[12];
  tmp2.len = 8;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_sub_mod(&tmp1, &tmp1, &tmp2, &prime_p256);

  // 2(s2+s3) + s1 + s4 + s5 - s6 - s7 - s8 - s9
  //                7   6    5   4   3   2    1    0
  // tmp2 = s9 = (c13,  0, c11, c10, c9, 0, c15, c14)  
  tmp2.dat[0] = in.dat[14];
  tmp2.dat[1] = in.dat[15];
  tmp2.dat[2] = tmp2.dat[6] = 0;
  tmp2.dat[3] = in.dat[ 9];
  tmp2.dat[4] = in.dat[10];
  tmp2.dat[5] = in.dat[11];
  tmp2.dat[7] = in.dat[13];
  tmp2.len = 8;
  while ((tmp2.len > 0) && (tmp2.dat[tmp2.len - 1] == 0)) tmp2.len--;

  ECC_bn_sub_mod(c, &tmp1, &tmp2, &prime_p256);

  return ECC_PASS;
}




int ECC_bn_mul_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p) 
{
  ECC_bn_mul(c, a, b);

  // 모듈러감산 추가요

  return ECC_PASS;
}

int ECC_bn_mul_mod_p256(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p)
{
  ECC_bn_mul(c, a, b);
  ECC_bn_mod_p256(c, c);

  return ECC_PASS;
}

// c= a^(-1) mod p
int ECC_bn_binary_inv(ECC_BN* c, ECC_BN* a, ECC_BN* p)
{
  ECC_BN u, v, x1, x2;

  if (a->len == 0) {
	return ECC_FAIL;
  }

  ECC_bn_cpy(&u, a);
  ECC_bn_cpy(&v, p);
  x1.dat[0] = 1;
  x1.len    = 1;
  x2.len    = 0;

  while (!((u.len==1)&&(u.dat[0]==1)) && !((v.len == 1) && (v.dat[0] == 1))) {
	while ( !(u.dat[0]&1) ) {
	  ECC_bn_1bit_rshift(&u,&u);
	  if (x1.dat[0] & 1) {
		ECC_bn_add(&x1, p, &x1);
	  }
	  ECC_bn_1bit_rshift(&x1, &x1);
	}

	while (!(v.dat[0] & 1)) {
	  ECC_bn_1bit_rshift(&v, &v);
	  if (x2.dat[0] & 1) {
		ECC_bn_add(&x2, p, &x2);
	  }
	  ECC_bn_1bit_rshift(&x2, &x2);
	}

	if (ECC_bn_cmp(&u, &v) >= 0) {
	  ECC_bn_sub(&u, &u, &v);
	  ECC_bn_sub_mod(&x1, &x1, &x2, p);
	}
	else {
	  ECC_bn_sub(&v, &v, &u);
	  ECC_bn_sub_mod(&x2, &x2, &x1, p);
	}	
  }

  if ((u.len == 1) && u.dat[0] == 1) {
	ECC_bn_cpy(c, &x1);
  }
  else {
	ECC_bn_cpy(c, &x2);
  }

  return ECC_PASS;
}

// a,b in GF(p)
int ECC_bn_add_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p)
{
  ECC_BN out;

  if ((a->len < 0) || (b->len < 0))
	return ECC_FAIL;
  
  if (a->len >= b->len) ECC_bn_add(&out, a, b);
  else                  ECC_bn_add(&out, b, a);

  if(ECC_bn_cmp(&out,p)>=0) 
	ECC_bn_sub(c, &out, p);
  else                      
	ECC_bn_cpy(c, &out);

  return ECC_PASS;
}

// a,b in GF(p)
int ECC_bn_sub_mod(ECC_BN* c, ECC_BN* a, ECC_BN* b, ECC_BN* p)
{
  ECC_BN out;

  if ((a->len < 0) || (b->len < 0))
	return ECC_FAIL;

  if (ECC_bn_cmp(a, b)>=0) {
	ECC_bn_sub(c, a, b);
  }
  else {
	ECC_bn_sub(&out, p, b);

	if(out.len>=a->len) ECC_bn_add(c, &out, a);
	else                ECC_bn_add(c, a, &out);
  }

  return ECC_PASS;
}




int ECC_pt_init(ECC_PT* P)
{
  P->point_at_infinity = 1;

  return ECC_PASS;
}

int ECC_pt_cpy(ECC_PT* R, ECC_PT* P)
{
  ECC_bn_cpy(&R->x, &P->x);
  ECC_bn_cpy(&R->y, &P->y);

  R->point_at_infinity = P->point_at_infinity;

  return ECC_PASS;
}

// R = P + Q
int ECC_pt_add(ECC_PT* R, ECC_PT* P, ECC_PT* Q)
{
  if (P->point_at_infinity == 1) {
	ECC_pt_cpy(R, Q);
	return ECC_PASS;
  }

  if (Q->point_at_infinity == 1) {
	ECC_pt_cpy(R, P);
	return ECC_PASS;
  }

  if (ECC_bn_cmp(&P->x, &Q->x) == 0) {
	if (ECC_bn_cmp(&P->y, &Q->y) == 0) {
	  ECC_pt_dbl(R, P);
	}
	else {
	  R->point_at_infinity = 1;
	}
  }
  else {
	// add 계산
    ECC_BN m, temp, temp_t;

		// m = (y2 - y1)/(x2 - x1)
		ECC_bn_sub_mod(&temp, &Q->y, &P->y, &prime_p256);
		ECC_bn_sub_mod(&m, &Q->x, &P->x, &prime_p256);
		ECC_bn_binary_inv(&m, &m, &prime_p256);

		ECC_bn_mul_mod_p256(&m, &temp, &m, &prime_p256);
		// x3 = m^2 - x1 - x2
		ECC_bn_mul_mod_p256(&temp, &m, &m, &prime_p256);
		ECC_bn_sub_mod(&temp, &temp, &P->x, &prime_p256);
		ECC_bn_sub_mod(&temp_t, &temp, &Q->x, &prime_p256);

    
		// y3 = m(x1-x3) - y1
		ECC_bn_sub_mod(&temp, &P->x, &temp_t, &prime_p256);
		ECC_bn_mul_mod_p256(&temp, &m, &temp, &prime_p256);
		ECC_bn_sub_mod(&R->y, &temp, &P->y, &prime_p256);
    ECC_bn_cpy(&R->x, &temp_t);
    
		R->point_at_infinity = 0;
  }

  return ECC_PASS;
}

// R = 2P
int ECC_pt_dbl(ECC_PT* R, ECC_PT* P)
{
  if (P->point_at_infinity == 1) {
	R->point_at_infinity = 1;
	return ECC_PASS;
  }

  // dbl 계산
  ECC_BN m, temp, out_x;
	ECC_BN x_squared, t_x_squared;

	ECC_bn_mul_mod_p256(&x_squared, &P->x, &P->x, &prime_p256);    
	ECC_bn_add_mod(&t_x_squared, &x_squared, &x_squared, &prime_p256); 
	ECC_bn_add_mod(&t_x_squared, &t_x_squared, &x_squared, &prime_p256); 

  ECC_BN three = {{3}, 1};
	ECC_bn_sub_mod(&temp, &t_x_squared, &three, &prime_p256);
	ECC_bn_add_mod(&m, &P->y, &P->y, &prime_p256); 
	ECC_bn_binary_inv(&m, &m, &prime_p256);    
	ECC_bn_mul_mod_p256(&m, &temp, &m, &prime_p256); 

	ECC_bn_mul_mod_p256(&temp, &m, &m, &prime_p256);
	ECC_bn_add_mod(&out_x, &P->x, &P->x, &prime_p256);        
	ECC_bn_sub_mod(&out_x, &temp, &out_x, &prime_p256);       

	ECC_bn_sub_mod(&temp, &P->x, &out_x, &prime_p256);
	ECC_bn_mul_mod_p256(&temp, &m, &temp, &prime_p256);
	ECC_bn_sub_mod(&R->y, &temp, &P->y, &prime_p256);

  ECC_bn_cpy(&R->x, &out_x);
	R->point_at_infinity = 0;

  return ECC_PASS;
}