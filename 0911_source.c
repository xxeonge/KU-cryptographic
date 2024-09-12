#include <stdio.h>
#include <stdlib.h>
#include "KISA_SHA256.h"
void main(void)
{
  unsigned char msg[5] = { 0 };
  unsigned char hash_out[SHA256_DIGEST_VALUELEN] = { 0 };
  unsigned int  msg_len = 0;
  SHA256_INFO   info;
  //msg = 0x0102030405
  msg[0] = 0x01;
  msg[1] = 0x02;
  msg[2] = 0x03;
  msg[3] = 0x04;
  msg[4] = 0x05;
  msg_len = 5;
  SHA256_Encrpyt(msg, msg_len, hash_out);
  printf("hash_out0=");
  for (int i = 0; i < SHA256_DIGEST_VALUELEN; i++) {
    printf("%02x", hash_out[i]);
  }
  printf("\n");

  SHA256_Init(&info);
  //SHA256_Process(&info, msg,msg_len);
  SHA256_Process(&info, &msg[0], 1);
  SHA256_Process(&info, &msg[1], 2);
  SHA256_Process(&info, &msg[3], 2);
  SHA256_Close(&info, hash_out);
  printf("hash_out1=");
  for (int i = 0; i < SHA256_DIGEST_VALUELEN; i++) {
    printf("%02x", hash_out[i]);
  }
  printf("\n");


}