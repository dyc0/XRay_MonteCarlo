#include <iostream>
#include <immintrin.h>

int main()
{
    double x = 0;
    double y = 1;
    double z = 2;
    double r = 3;

    __m256d firsts = _mm256_set_pd(x, y, z, r);
    std::cout << std::endl;
}