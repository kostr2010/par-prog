#include <cassert>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    assert(argc == 2 /* argument STEPS needed */);

    int steps = atoi(argv[1]);
    long double euler = 0, fac = 1;

    for (int i = steps; i >= 1; i--) {
        fac *= i;
        euler += fac;
    }

    euler /= fac;

    printf("euler:      %.20llf\n", euler);
    printf("true euler: 2.71828182845904523536028747135266249775724709369995...\n");

    return 0;
}