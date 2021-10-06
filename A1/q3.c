#include <stdio.h>
#include <math.h>

void run_single_precision() {
    float eapprox, n;
    for (int i = 0; i <= 12; i++) {
        n = pow(10, i);
        eapprox = pow((1 + 1 / n), n);
        printf("%13.0f %13.10f %13.10f %11.3e %14.11f\n",
            n,
            eapprox,
            exp(1),
            (exp(1) - eapprox) / exp(1),
            (1 + 1 / n)
        );
    }
}

void run_double_precision() {
    double eapprox, n;
    for (int i = 0; i <= 12; i++) {
        n = pow(10, i);
        eapprox = pow((1 + 1 / n), n);
        printf("%13.0f %13.10f %13.10f %11.3e %14.11f\n",
            n,
            eapprox,
            exp(1),
            (exp(1) - eapprox) / exp(1),
            (1 + 1 / n)
        );
    }
}

int main() {
    run_single_precision();
    printf("\n");
    run_double_precision();
    return 0;
}