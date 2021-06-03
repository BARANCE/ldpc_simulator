void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);
unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

/*
genrand_int32() generates unsigned 32-bit integers.
genrand_int31() generates unsigned 31-bit integers.
genrand_real1() generates uniform real in [0,1] (32-bit resolution).
genrand_real2() generates uniform real in [0,1) (32-bit resolution).
genrand_real3() generates uniform real in (0,1) (32-bit resolution).
genrand_res53() generates uniform real in [0,1) with 53-bit resolution.
 */
