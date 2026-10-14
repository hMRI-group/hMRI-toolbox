#include <stdio.h>
#include <math.h>
#include <stdbool.h>

/**
 * Computes the kernel weight for a given kernel type and squared distance.
 * 
 * This function implements various kernel functions commonly used in
 * non-parametric density estimation and smoothing. The kernel type is
 * specified by the `kern` parameter, and the squared distance is provided
 * as `xsq`.
 * 
 * @param kern The kernel type identifier:
 *             - 1: Triangular kernel
 *             - 2: Epanechnikov kernel (default)
 *             - 3: Quadratic kernel
 *             - 4: Cubic kernel
 *             - 5: Gaussian kernel
 * @param xsq The squared distance from the center of the kernel.
 * 
 * @return The kernel weight for the given squared distance.
 *         Returns 0.0 if `xsq >= 1.0`.
 */
double lkern(int kern, double xsq) {
    double z;
    
    if (xsq >= 1.0) {
        return 0.0;
    } else if (kern == 1) {
        if (xsq <= 0.5) {
            return 1.0;
        } else {
            return 2.0 * (1.0 - xsq);
        }
    } else if (kern == 2) {
        return 1.0 - xsq;
    } else if (kern == 3) {
        z = 1.0 - xsq;
        return z * z;
    } else if (kern == 4) {
        z = 1.0 - xsq;
        return z * z * z;
    } else if (kern == 5) {
        return exp(-xsq * 8.0);
    } else {
        // Default: Epanechnikov kernel
        return 1.0 - xsq;
    }
}

/**
 * Computes the Kullback-Leibler (KL) divergence between two multivariate
 * Gaussian distributions.
 * 
 * This function calculates the squared Mahalanobis distance between two
 * vectors `thi` and `thj`, weighted by the inverse covariance matrix `si2`.
 * The KL divergence is a measure of the difference between two probability
 * distributions and is used in the adaptive weights smoothing (AWS) to
 * determine the similarity between voxels.
 * 
 * @param thi Pointer to the first vector (mean of the first Gaussian).
 * @param thj Pointer to the second vector (mean of the second Gaussian).
 * @param si2 Pointer to the inverse covariance matrix (stored in lower
 *            triangular form as a vector).
 * @param nv The number of dimensions (variables) in the vectors.
 * 
 * @return The squared Mahalanobis distance (KL divergence) between the
 *         two vectors.
 */
double KLdistsi(double *thi, double *thj, double *si2, int nv) {
    double z = 0.0;
    double zdk;
    int k, l, m = 0;

    for (k = 0; k < nv; k++) {
        zdk = thi[k] - thj[k];
        if (k > 0) {
            for (l = 0; l < k; l++) {
                z += 2.0 * (thi[l] - thj[l]) * zdk * si2[m];
                m++;
            }
        }
        z += zdk * zdk * si2[m];
        m++;
    }

    return z;
}

/**
 * Performs Patchwise Vector Adaptive Weights Smoothing (PVAWS) on voxel data.
 * 
 * This function implements the PVAWS algorithm, which adaptively smooths voxel
 * data while preserving edges and boundaries between different tissue types.
 * It uses a kernel-based smoothing approach with adaptive weights to reduce
 * noise while maintaining structural integrity. The adaptive weighting scheme
 * (AWS) adjusts the smoothing based on the similarity between neighboring voxels,
 * measured using the Kullback-Leibler (KL) divergence.
 * 
 * @param y Pointer to the input voxel data (e.g., MRI signal intensities).
 * @param pos Pointer to the voxel positions or indices.
 * @param nv Number of variables (e.g., dimensions in the data).
 * @param nvd Number of variables in the inverse covariance matrix.
 * @param n1 Number of voxels in the first dimension (x-axis).
 * @param n2 Number of voxels in the second dimension (y-axis).
 * @param n3 Number of voxels in the third dimension (z-axis).
 * @param hakt Bandwidth parameter for the kernel smoothing.
 * @param lambda Regularization parameter for the adaptive weighting.
 * @param theta Pointer to the initial parameter estimates for each voxel.
 * @param bi Pointer to the initial local sum of weights for each voxel.
 * @param bin Pointer to the output local sum of weights after smoothing.
 * @param thnew Pointer to the output parameter estimates after smoothing.
 * @param invcov Pointer to the inverse covariance matrices for each voxel.
 * @param spmin Plateau of the kernel function starts here.
 * @param lwght Pointer to the location weights for the kernel.
 * @param wght Pointer to the voxel size ratios (w1, w2, w3).
 * @param swjy Pointer to the intermediate weighted sums for the parameters.
 * @param np1 Size of the patches in the first dimension.
 * @param np2 Size of the patches in the second dimension.
 * @param np3 Size of the patches in the third dimension.
 */
void pvaws(
    double *y,       // 1
    int *pos,        // 2
    int nv,          // 3
    int nvd,         // 4
    int n1,          // 5
    int n2,          // 6
    int n3,          // 7
    double hakt,     // 8
    double lambda,   // 9
    double *theta,   // 10
    double *bi,      // 11
    double *bin,     // 12
    double *thnew,   // 13
    double *invcov,  // 14
    double spmin,    // 16
    double *lwght,   // 17
    double *wght,    // 18
    double *swjy,    // 19   
    int np1,         // 20
    int np2,         // 21
    int np3          // 22
) {
    int ih1, ih2, ih3, i1, i2, i3, j1, j2, j3, jw1, jw2, jw3, jwind3, jwind2;
    int iind, jind, jind3, jind2, clw1, clw2, clw3, dlw1, dlw2, dlw3, dlw12, n12, k; 
    int ip1, ip2, ip3, nph1, nph2, nph3, ipind, jp1, jp2, jp3, jpind;
    int iindp, jindp, ipindp, jpindp;
    double sij, swj, z1, z2, z3, wj, w1, w2, sijp, hakt2, spf;
    bool aws;

    hakt2 = hakt * hakt;
    spf = 1.0 / (1.0 - spmin);
    aws = lambda < 1e35;

    w1 = wght[0];      // ratio of y-extension to the x-extension of the voxel
    w2 = wght[1];      // ratio of z-extension to the x-extension of the voxel
    ih1 = (int)floor(hakt);       // rounding down the actual bandwidth of this step to the nearest integer
    ih2 = (int)floor(hakt / w1);  // rounding down the actual bandwidth of this step in y-direction to the nearest integer
    ih3 = (int)floor(hakt / w2);  // rounding down the actual bandwidth of this step in z-direction to the nearest integer
    if (n2 == 1) ih2 = 0;  // can this actually happen in practice? If n2 is 1, then the data is 1D and there is no bandwidth in y-direction, so set ih2 to 0
    if (n3 == 1) ih3 = 0;  // can this actually happen in practice? If n3 is 1, then the data is 2D and there is no bandwidth in z-direction, so set ih3 to 0

    // define some variables in connection with the kernel indexing and the bandwidth
    clw1 = ih1;
    clw2 = ih2;
    clw3 = ih3;
    dlw1 = ih1 + clw1 + 1; // the number of points in the kernel in x-direction (from -ih1 to ih1, plus the center point, plus the points in the other direction)
                           // we could also just set dlw1 to 2*ih1+1
    dlw2 = ih2 + clw2 + 1; // the number of points in the kernel in y-direction (from -ih2 to ih2, plus the center point, plus the points in the other direction)
                           // we could also just set dlw2 to 2*ih2+1
    dlw3 = ih3 + clw3 + 1; // the number of points in the kernel in z-direction (from -ih3 to ih3, plus the center point, plus the points in the other direction)
                           // we could also just set dlw3 to 2*ih3+1
    dlw12 = dlw1 * dlw2;   // the number of points in the kernel in x and y direction
    nph1 = (np1 - 1) / 2;  // this is actually the same as clw1 or ih1
    nph2 = (np2 - 1) / 2;  // this is actually the same as clw2 or ih2
    nph3 = (np3 - 1) / 2;  // this is actually the same as clw3 or ih3

    n12 = n1 * n2;         // the number of voxels in one z-slice

    z2 = 0.0;  // this is the temporary variable for the squared distance in y-direction                                             
    z3 = 0.0;  // this is the temporary variable for the squared distance in z-direction

    // First calculate location weights for all possible locations within the bandwidth of the kernel. 
    // This is done in a nested loop over the possible locations in z, y and x direction (j3, j2, j1). 
    // The location weights are stored in the array lwght, which is indexed by the location 
    // in the kernel (j1, j2, j3) and the bandwidth (clw1, clw2, clw3). 
    // The location weights are calculated using the function lkern, which takes as input 
    // the squared distance from the center of the kernel (z1^2 + z2^2 + z3^2) divided 
    // by the squared bandwidth (hakt^2). The location weights are used later to weight 
    // the contributions of the neighboring voxels to the estimation of theta and bi for each voxel.
    for (j3 = -clw3; j3 <= clw3; j3++) {
        if (n3 > 1) {
            z3 = j3 * w2;
            z3 = z3 * z3;
            ih2 = (int)floor(sqrt(hakt2 - z3) / w1);
            jind3 = (j3 + clw3) * dlw12;
        } else {
            jind3 = 0;
        }
        for (j2 = -ih2; j2 <= ih2; j2++) {
            if (n2 > 1) {
                z2 = j2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = (int)floor(sqrt(hakt2 - z2));
                jind2 = jind3 + (j2 + clw2) * dlw1;
            } else {
                jind2 = 0;
            }
            for (j1 = -ih1; j1 <= ih1; j1++) {
                jind = j1 + clw1 + jind2; //minimal value is jind=0, maximal value is jind=dlw1*dlw2*dlw3-1
                z1 = j1;
                lwght[jind] = lkern(2, (z1 * z1 + z2) / hakt2);
            }
        }
    }

    // main AWS loop over all voxels in the data (iind) to calculate the new estimates of theta and bi for each voxel.
    for (iind = 0; iind < n1 * n2 * n3; iind++) {
        iindp = pos[iind];
        if (iindp == 0) continue;

        // This is the conversion from the linear index iind
        // to the 3D indices i1, i2, i3 for the voxel location.
        // The linear index iind runs from 0 to n1*n2*n3-1.
        // The 3D indices i1, i2, i3 run from 1 to n1, 1 to n2, and 1 to n3, respectively.
        i1 = (iind + 1) % n1;
        if (i1 == 0) i1 = n1;
        i2 = ((iind + 1 - i1) / n1 + 1) % n2;
        if (i2 == 0) i2 = n2;
        i3 = (iind + 1 - i1 - (i2 - 1) * n1) / n12 + 1;

        swj = 0.0;
        for (k = 0; k < nv; k++) {
            swjy[k] = 0.0;
        }

        for (jw3 = -clw3; jw3 <= clw3; jw3++) {
            j3 = jw3 + i3;
            if (j3 < 1 || j3 > n3) continue;
            jwind3 = (jw3 + clw3) * dlw12;
            jind3 = (j3 - 1) * n12;
            z3 = jw3 * w2;
            z3 = z3 * z3;
            if (n2 > 1) ih2 = (int)floor(sqrt(hakt2 - z3) / w1);
            for (jw2 = -ih2; jw2 <= ih2; jw2++) {
                j2 = jw2 + i2;
                if (j2 < 1 || j2 > n2) continue;
                jwind2 = jwind3 + (jw2 + clw2) * dlw1;
                jind2 = (j2 - 1) * n1 + jind3;
                z2 = jw2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = (int)floor(sqrt(hakt2 - z2));
                for (jw1 = -ih1; jw1 <= ih1; jw1++) {
                    j1 = jw1 + i1;
                    if (j1 < 1 || j1 > n1) continue;
                    jind = j1 + jind2;
                    jindp = pos[jind - 1];
                    if (jindp == 0) continue;
                    wj = lwght[jw1 + clw1 + jwind2];
                    if (aws) {
                        sij = 0.0;
                        for (ip1 = i1 - nph1; ip1 <= i1 + nph1; ip1++) {
                            if (ip1 <= 0 || ip1 > n1) continue;
                            jp1 = ip1 + jw1;
                            if (jp1 <= 0 || jp1 > n1) continue;
                            for (ip2 = i2 - nph2; ip2 <= i2 + nph2; ip2++) {
                                if (ip2 <= 0 || ip2 > n2) continue;
                                jp2 = ip2 + jw2;
                                if (jp2 <= 0 || jp2 > n2) continue;
                                for (ip3 = i3 - nph3; ip3 <= i3 + nph3; ip3++) {
                                    if (sij > 1.0) continue;
                                    if (ip3 <= 0 || ip3 > n3) continue;
                                    ipind = ip1 + (ip2 - 1) * n1 + (ip3 - 1) * n12;
                                    ipindp = pos[ipind - 1];
                                    if (ipindp == 0) continue;
                                    jp3 = ip3 + jw3;
                                    if (jp3 <= 0 || jp3 > n3) continue;
                                    jpind = jp1 + (jp2 - 1) * n1 + (jp3 - 1) * n12;
                                    jpindp = pos[jpind - 1];
                                    if (jpindp == 0) continue;
                                    sijp = KLdistsi(&theta[(jpindp - 1) * nv],
                                                    &theta[(ipindp - 1) * nv],
                                                    &invcov[(ipindp - 1) * nvd],
                                                    nv);
                                    sij = fmax(sij, bi[ipindp - 1] / lambda * sijp);
                                }
                            }
                        }
                        if (sij >= 1.0) continue;
                        if (sij > spmin) wj = wj * (1.0 - spf * (sij - spmin));
                    }
                    swj += wj;
                    for (k = 0; k < nv; k++) {
                        swjy[k] += wj * y[(jindp - 1) * nv + k];
                    }
                }
            }
        }
        for (k = 0; k < nv; k++) {
            thnew[(iindp - 1) * nv + k] = swjy[k] / swj;
        }
        bin[iindp - 1] = swj;
    }
}


void pvawslast(
    double *y,       // 1
    double *yd,      // 2
    int *pos,        // 3
    int nv,          // 4
    int nvd,         // 5
    int nd,          // 6
    int n1,          // 7
    int n2,          // 8
    int n3,          // 9
    double hakt,     // 10
    double lambda,   // 11
    double *theta,   // 12
    double *bi,      // 13
    double *bin,     // 14
    double *thnew,   // 15
    double *ydnew,   // 16
    double *invcov,  // 17
    int ncores,      // 18 OMP not needed for pvawsme, but we keep it as an argument for consistency with pvaws
    double spmin,    // 19
    double *lwght,   // 20 
    double *wght,    // 21
    double *swjy,    // 22
    double *swjd,    // 23
    int np1,         // 24
    int np2,         // 25
    int np3          // 26
) {

    int ih1, ih2, ih3, i1, i2, i3, j1, j2, j3, jw1, jw2, jw3, jwind3, jwind2;
    int iind, jind, jind3, jind2, clw1, clw2, clw3, dlw1, dlw2, dlw3, dlw12, n12, k; 
    int iindp, jindp, ipindp, jpindp;
    int ip1, ip2, ip3, nph1, nph2, nph3, ipind, jp1, jp2, jp3, jpind;
    double sij, swj, z1, z2, z3, wj, hakt2, w1, w2, sijp, spf;
    bool aws;

    hakt2 = hakt * hakt;
    spf = 1.0 / (1.0 - spmin);
    aws = lambda < 1e35;

    w1 = wght[0];
    w2 = wght[1];
    ih1 = (int)floor(hakt);
    ih2 = (int)floor(hakt / w1);
    ih3 = (int)floor(hakt / w2);
    if (n2 == 1) ih2 = 0;
    if (n3 == 1) ih3 = 0;

    clw1 = ih1;
    clw2 = ih2;
    clw3 = ih3;
    dlw1 = ih1 + clw1 + 1;
    dlw2 = ih2 + clw2 + 1;
    dlw3 = ih3 + clw3 + 1;
    dlw12 = dlw1 * dlw2;
    nph1 = (np1 - 1) / 2;
    nph2 = (np2 - 1) / 2;
    nph3 = (np3 - 1) / 2;

    n12 = n1 * n2;

    z2 = 0.0;
    z3 = 0.0;

    for (j3 = -clw3; j3 <= clw3; j3++) {
        if (n3 > 1) {
            z3 = j3 * w2;
            z3 = z3 * z3;
            ih2 = (int)floor(sqrt(hakt2 - z3) / w1);
            jind3 = (j3 + clw3) * dlw12;
        } else {
            jind3 = 0;
        }
        for (j2 = -ih2; j2 <= ih2; j2++) {
            if (n2 > 1) {
                z2 = j2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = (int)floor(sqrt(hakt2 - z2));
                jind2 = jind3 + (j2 + clw2) * dlw1;
            } else {
                jind2 = 0;
            }
            for (j1 = -ih1; j1 <= ih1; j1++) {
                jind = j1 + clw1 + jind2;
                z1 = j1;
                lwght[jind] = lkern(2, (z1 * z1 + z2) / hakt2);
            }
        }
    }

    // Call to rchkusr() - assuming it's a placeholder for some user-defined function
    // rchkusr();

    for (iind = 0; iind < n1 * n2 * n3; iind++) {
        iindp = pos[iind];
        if (iindp == 0) continue;

        // This is the conversion from the linear index iind
        // to the 3D indices i1, i2, i3 for the voxel location.
        // The linear index iind runs from 0 to n1*n2*n3-1.
        // The 3D indices i1, i2, i3 run from 1 to n1, 1 to n2, and 1 to n3, respectively.
        i1 = (iind + 1) % n1;
        if (i1 == 0) i1 = n1;
        i2 = ((iind + 1 - i1) / n1 + 1) % n2;
        if (i2 == 0) i2 = n2;
        i3 = (iind + 1 - i1 - (i2 - 1) * n1) / n12 + 1;

        swj = 0.0;
        for (k = 0; k < nv; k++) {
            swjy[k] = 0.0;
        }
        for (k = 0; k < nd; k++) {
            swjd[k] = 0.0;
        }

        for (jw3 = -clw3; jw3 <= clw3; jw3++) {
            j3 = jw3 + i3;
            if (j3 < 1 || j3 > n3) continue;
            jwind3 = (jw3 + clw3) * dlw12;
            jind3 = (j3 - 1) * n12;
            z3 = jw3 * w2;
            z3 = z3 * z3;
            if (n2 > 1) ih2 = (int)floor(sqrt(hakt2 - z3) / w1);
            for (jw2 = -ih2; jw2 <= ih2; jw2++) {
                j2 = jw2 + i2;
                if (j2 < 1 || j2 > n2) continue;
                jwind2 = jwind3 + (jw2 + clw2) * dlw1;
                jind2 = (j2 - 1) * n1 + jind3;
                z2 = jw2 * w1;
                z2 = z3 + z2 * z2;
                ih1 = (int)floor(sqrt(hakt2 - z2));
                for (jw1 = -ih1; jw1 <= ih1; jw1++) {
                    j1 = jw1 + i1;
                    if (j1 < 1 || j1 > n1) continue;
                    jind = j1 + jind2;
                    jindp = pos[jind - 1];
                    if (jindp == 0) continue;
                    wj = lwght[jw1 + clw1 + jwind2];
                    if (aws) {
                        sij = 0.0;
                        for (ip1 = i1 - nph1; ip1 <= i1 + nph1; ip1++) {
                            if (ip1 <= 0 || ip1 > n1) continue;
                            jp1 = ip1 + jw1;
                            if (jp1 <= 0 || jp1 > n1) continue;
                            for (ip2 = i2 - nph2; ip2 <= i2 + nph2; ip2++) {
                                if (ip2 <= 0 || ip2 > n2) continue;
                                jp2 = ip2 + jw2;
                                if (jp2 <= 0 || jp2 > n2) continue;
                                for (ip3 = i3 - nph3; ip3 <= i3 + nph3; ip3++) {
                                    if (sij > 1.0) continue;
                                    if (ip3 <= 0 || ip3 > n3) continue;
                                    ipind = ip1 + (ip2 - 1) * n1 + (ip3 - 1) * n12;
                                    ipindp = pos[ipind - 1];
                                    if (ipindp == 0) continue;
                                    jp3 = ip3 + jw3;
                                    if (jp3 <= 0 || jp3 > n3) continue;
                                    jpind = jp1 + (jp2 - 1) * n1 + (jp3 - 1) * n12;
                                    jpindp = pos[jpind - 1];
                                    if (jpindp == 0) continue;
                                    sijp = KLdistsi(&theta[(jpindp - 1) * nv],
                                                    &theta[(ipindp - 1) * nv],
                                                    &invcov[(ipindp - 1) * nvd],
                                                    nv);
                                    sij = fmax(sij, bi[ipindp - 1] / lambda * sijp);
                                }
                            }
                        }
                        if (sij >= 1.0) continue;
                        if (sij > spmin) wj = wj * (1.0 - spf * (sij - spmin));
                    }
                    swj += wj;
                    for (k = 0; k < nv; k++) {
                        swjy[k] += wj * y[(jindp - 1) * nv + k];
                    }
                    for (k = 0; k < nd; k++) {
                        swjd[k] += wj * yd[(jindp - 1) * nd + k];
                    }
                }
            }
        }
        for (k = 0; k < nv; k++) {
            thnew[(iindp - 1) * nv + k] = swjy[k] / swj;
        }
        for (k = 0; k < nd; k++) {
            ydnew[(iindp - 1) * nd + k] = swjd[k] / swj;
        }
        bin[iindp - 1] = swj;
    }
}


