// perform an XOR of testfile with keyfile
// the key file is repeated if it is shorter than the testfile

#include <cstdio>

int main(int argc, char *argv[]){ 

  FILE *inp, *key=0;

  if (argc<=2) {
    printf("Usage: xor testfile keyfile\n"); 
    return 1;
  }
  inp=fopen(argv[1],"r");
  key=fopen(argv[2],"r");
  
  unsigned char c;

  while(1){
    c=fgetc(inp);          // read in 1 byte at a time
    if (feof(inp)) break;
    char x=fgetc(key);
    if (feof(key)) {
      rewind(key);
      x=fgetc(key);
    }
    c=c^x;
    fputc(c,stdout);
  }
  return 0;
}
