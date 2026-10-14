#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Function to compute the location kernel
// Forward declaration of lkern from aws.c
extern double lkern(int kern, double xsq);

// Function to compute the sum of location weights
static double sofw(double bw, int kern, double wght[2]) {
    double h2 = bw * bw;
    int ih3 = (int)floor(bw * wght[1]);
    int ih2 = (int)floor(bw * wght[0]);
    int ih1 = (int)floor(bw);
    int dlw1 = 2 * ih1 + 1;
    int dlw2 = 2 * ih2 + 1;
    int dlw3 = 2 * ih3 + 1;
    int clw1 = (dlw1 + 1) / 2;
    int clw2 = (dlw2 + 1) / 2;
    int clw3 = (dlw3 + 1) / 2;
    double sw = 0.0;
    double sw2 = 0.0;
    
    for (int j1 = 1; j1 <= dlw1; j1++) {
        double z1 = (clw1 - j1);
        z1 = z1 * z1;
        if (wght[0] > 0.0) {
            ih2 = (int)floor(sqrt(h2 - z1) * wght[0]);
            for (int j2 = clw2 - ih2; j2 <= clw2 + ih2; j2++) {
                double z2 = (clw2 - j2) / wght[0];
                z2 = z1 + z2 * z2;
                if (wght[1] > 0.0) {
                    ih3 = (int)floor(sqrt(h2 - z2) * wght[1]);
                    for (int j3 = clw3 - ih3; j3 <= clw3 + ih3; j3++) {
                        double z3 = (clw3 - j3) / wght[1];
                        double z = lkern(kern, (z3 * z3 + z2) / h2);
                        sw += z;
                        sw2 += z * z;
                    }
                } else {
                    double z = lkern(kern, z2 / h2);
                    sw += z;
                    sw2 += z * z;
                }
            }
        } else {
            double z = lkern(kern, z1 / h2);
            sw += z;
            sw2 += z * z;
        }
    }
    return sw * sw / sw2;
}

// Main function to compute bandwidth
void cgethani(double x, double y, int kern, double value, double wght[2], double eps, double *bw) {
    if (x >= y) {
        *bw = 1.0;
        return;
    }
    double fw1 = sofw(x, kern, wght);
    double fw2 = sofw(y, kern, wght);
    
    while (fw1 > value) {
        x = x * x / y;
        fw1 = sofw(x, kern, wght);
    }
    
    while (fw2 <= value) {
        y = y * y / x;
        fw2 = sofw(y, kern, wght);
    }
    
    while (fmin(fw2 / value, value / fw1) > 1.0 + eps) {
        double z = x + (value - fw1) / (fw2 - fw1) * (y - x);
        double fw3 = sofw(z, kern, wght);
        if (fw3 <= value) {
            x = z;
            fw1 = fw3;
        }
        if (fw3 >= value) {
            y = z;
            fw2 = fw3;
        }
    }
    
    if (fw2 / value > value / fw1) {
        *bw = x + (value - fw1) / (fw2 - fw1) * (y - x);
    } else {
        *bw = y - (fw2 - value) / (fw2 - fw1) * (y - x);
    }
}