#ifndef VIENNA_RNA_PACKAGE_CONFIG_H
#define VIENNA_RNA_PACKAGE_CONFIG_H

/* version number */
#define VRNA_VERSION  2.4.1

/*
 * The following pre-processor definitions specify whether
 * or not certain features were activated upon build-time
 */

/*
 * Build with deactivated C11 Features
 *
 * It this feature is missing, the next line defines
 * 'VRNA_DISABLE_C11_FEATURES'
 */


/*
 * Build with OpenMP support
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_OPENMP'
 */
#define VRNA_WITH_OPENMP

/*
 * Build with single precision partition function
 *
 * If this feature is present, the next line defines
 * 'USE_FLOAT_PF'
 */


/*
 * Build with Boustrophedon speedup in stochastic backtracking
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_BOUSTROPHEDON'
 */
#define VRNA_WITH_BOUSTROPHEDON

/*
 * Build with JSON input/output support
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_JSON_SUPPORT'
 */
#define VRNA_WITH_JSON_SUPPORT

/*
 * Build with Support Vector Machine (SVM) Z-score feature in RNALfold
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_SVM'
 */
#define VRNA_WITH_SVM

/*
 * Build with GSL minimizers
 *
 * If this feature is present, the next line defines
 * 'VRNA_WITH_GSL'
 */


/*
 * Build with colored TTY output
 *
 * If this feature is missing, the next line defines
 * 'VRNA_WITHOUT_TTY_COLORS'
 */


/*
 * Build with Link Time Optimization support
 *
 * If this feature is enabled, the next line defines
 * 'VRNA_WITH_LTO'
 */
#define VRNA_WITH_LTO

#endif
