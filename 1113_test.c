int ECC_bn_binary_inv(ECC_BN* c, ECC_BN* a)
{
    ECC_BN u, v, x1, x2;

    // a가 0인 경우 처리
    if (a->len == 0) {
        return ECC_PASS;
    }

    // u와 v 초기화 (u = a, v = p)
    ECC_bn_cpy(&u, a);
    ECC_bn_cpy(&v, &prime_p256);

    // 초기값 설정 (x1 = 1, x2 = 0)
    x1.dat[0] = 1;
    x1.len = 1;
    x2.len = 0;

    // 반복문: u = 1 and v != 1일 때까지 반복
    while (!(u.len == 1 && u.dat[0] == 1) && !(v.len == 1 && v.dat[0] == 1)) {

        // u가 짝수일 때 처리
        while (u.len > 0 && (u.dat[0] & 1) == 0) {
            ECC_bn_1bit_rshift(&u, &u);  // u = u / 2
            if ((x1.dat[0] & 1) == 0) {
                ECC_bn_1bit_rshift(&x1, &x1);  // x1 = x1 / 2
            } else {
                ECC_bn_add_mod(&x1, &x1, &prime_p256);  // x1 = (x1 + p) / 2
                ECC_bn_1bit_rshift(&x1, &x1);
            }
        }

        // v가 짝수일 때 처리
        while (v.len > 0 && (v.dat[0] & 1) == 0) {
            ECC_bn_1bit_rshift(&v, &v);  // v = v / 2
            if ((x2.dat[0] & 1) == 0) {
                ECC_bn_1bit_rshift(&x2, &x2);  // x2 = x2 / 2
            } else {
                ECC_bn_add_mod(&x2, &x2, &prime_p256);  // x2 = (x2 + p) / 2
                ECC_bn_1bit_rshift(&x2, &x2);
            }
        }

        // u >= v이면 u - v, x1 - x2 (mod p) 수행
        if (ECC_bn_cmp(&u, &v) >= 0) {
            ECC_bn_sub(&u, &u, &v);  // u = u - v
            ECC_bn_sub_mod(&x1, &x1, &x2, &prime_p256);  // x1 = x1 - x2 mod p
        } else {
            // v > u이면 v - u, x2 - x1 (mod p) 수행
            ECC_bn_sub(&v, &v, &u);  // v = v - u
            ECC_bn_sub_mod(&x2, &x2, &x1, &prime_p256);  // x2 = x2 - x1 mod p
        }
    }

        // If u = 1 then return(x1 mod p); else return(x2 mod p)
    if (u.len == 1 && u.dat[0] == 1) {
        ECC_bn_cpy(c, &x1);
    } else {
        ECC_bn_cpy(c, &x2);
    }
    return ECC_PASS;
}

