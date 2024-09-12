#include <stdio.h>
#include <stdlib.h>
#include "KISA_SHA256.h"

int main() {
    // 메시지 선언 및 초기화 (0x0102030405)
    uint8_t message[5] = {0x01, 0x02, 0x03, 0x04, 0x05};
    
    SHA256_INFO sha_info;

    uint8_t digest[SHA256_DIGEST_VALUELEN];

    SHA256_Init(&sha_info);

    SHA256_Process(&sha_info, message, sizeof(message));

    SHA256_Close(&sha_info, digest);

    // 결과 출력
    printf("SHA256(message) = ");
    for (int i = 0; i < SHA256_DIGEST_VALUELEN; i++) {
        printf("%02x", digest[i]);
    }
    printf("\n");

    return 0;
}


unsinged char msg[5]  = {0};
usinged char hash_out[SHA256_DIGEST_VALUELEN] = { 0 };
unsinged int msg_len = 0;
SHA256_INFO info;

msg[0] = 0x01;

SHA256_Encrpyt(msg,msg_len, hash_out);

printf("hash_out0=");
for (int i =0; i <SHA256)


SHA256_Process(&info, &msg[0], 1);

SHA256_Close(&info, hash_out);

printf("hash_out1");
for (int 0; i<SHA256)