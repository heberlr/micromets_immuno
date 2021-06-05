#include <iostream>
#include <cstdio>
#include <omp.h>

int main(int argc, char *argv[])
{
FILE * pFile = fopen (argv[2],"w");
#pragma omp parallel for
    for (int i = 0; i < omp_get_num_threads(); ++i)
    {
        fprintf(pFile,"C OpenMP. Process-thread: %s-%d\n", argv[1], i);
    }
    return 0;
}
