#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

int main(int argc, char **argv) {

  FILE * pFile;
  unsigned char* buffer;
  clock_t t0, t1;

  t0 = clock();

  // read buffer	
  const int BUFF_SIZE=1024;
  buffer = (unsigned char*)malloc(BUFF_SIZE);
  
  // array (histogram) with occurances of different byte patterns
  // each item index represents the value of a byte (0...255)
  // each item value represents the # of occurances of that byte
  long byteValues[256];
  for (int i = 0; i < 256; i++) byteValues[i] = 0;

  // total bytes
  int total = 0;

  // reading the file
  if (argc==1) pFile=stdin;
  else pFile = fopen(argv[1], "rb");
  int n;
  do { // binary read of file in chunks of BUFF_SIZE
    n = fread(buffer, 1, BUFF_SIZE, pFile);
    for (int i = 0; i < n; i++) {
      // update occurance of each byte read from buffer
      byteValues[(int) buffer[i]]++;
      total++;
    }
  } while (n != 0);
  fclose(pFile);

  //  entropy calculation
  double entropy = 0.0;
  for (int i = 0; i < 256; i++) {
    if (byteValues[i] != 0) {
      double r = (double) byteValues[i] / (double) total;
      entropy += -r * log2(r);
    }
  }
  t1 = clock();
  entropy/=8;  // divide by number of bits in each word
  printf("Entropy = %2.5f\n", entropy);
  printf("Exec time: %g milliseconds.\n", (double) (t1 - t0) * 1000/ (double) CLOCKS_PER_SEC);

  return 0;
}
