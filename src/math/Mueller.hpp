//==============================================================================
// base routines for treating Mueller and Johns matrices < 11.10.2002 >
//==============================================================================
#ifndef __Mueller_hpp_
#define __Mueller_hpp_

#ifndef MATRIX_HPP
#include "matrix.hpp"
#endif

#include "JonesMatrix.h"
#include "Matrix4x4.h"

/** @addtogroup Tracer Beam splitting algorithm
 * @{
 */

/** @addtogroup AxFunc Axillary functions
 * @{
 */

	/**
	@brief The function returns Mueller matrix calculated from Jones matrix.

	The dimension of Mueller matrix is 4x4.
	The dimension of Jones matrix must be 2x2.
	@param in: complex-value matrix (matrixC)
	@return real-value matrix
	*/
matrix Mueller(const matrixC & in);

matrix Mueller(const Matrix2x2c& in);

/// Add the Mueller matrix corresponding to a 2x2 Jones matrix directly to a
/// contiguous row-major 4x4 cell. This avoids constructing the legacy
/// dynamically allocated matrixC and matrix objects in direction-grid loops.
inline void AddMuellerFromJones(double *out,
                                const complex &j00, const complex &j01,
                                const complex &j10, const complex &j11,
                                double weight)
{
    const double j00r = real(j00), j00i = imag(j00);
    const double j01r = real(j01), j01i = imag(j01);
    const double j10r = real(j10), j10i = imag(j10);
    const double j11r = real(j11), j11i = imag(j11);

    const double a11 = j00r*j00r + j00i*j00i;
    const double a12 = j01r*j01r + j01i*j01i;
    const double a21 = j10r*j10r + j10i*j10i;
    const double a22 = j11r*j11r + j11i*j11i;

    double A1 = a11 + a21;
    double A2 = a12 + a22;
    out[0] += ((A1 + A2) * 0.5) * weight;
    out[1] += ((A1 - A2) * 0.5) * weight;
    A1 = a11 - a21;
    A2 = a12 - a22;
    out[4] += ((A1 + A2) * 0.5) * weight;
    out[5] += ((A1 - A2) * 0.5) * weight;

    double c1r = j00r*j01r + j00i*j01i;
    double c1i = j00i*j01r - j00r*j01i;
    double c2r = j11r*j10r + j11i*j10i;
    double c2i = j11i*j10r - j11r*j10i;
    out[2] += (-c1r - c2r) * weight;
    out[3] += ( c2i - c1i) * weight;
    out[6] += ( c2r - c1r) * weight;
    out[7] += (-c1i - c2i) * weight;

    c1r = j00r*j10r + j00i*j10i;
    c1i = j00i*j10r - j00r*j10i;
    c2r = j11r*j01r + j11i*j01i;
    c2i = j11i*j01r - j11r*j01i;
    out[8] += (-c1r - c2r) * weight;
    out[9] += ( c2r - c1r) * weight;
    out[12] += ( c1i - c2i) * weight;
    out[13] += ( c2i + c1i) * weight;

    c1r = j00r*j11r + j00i*j11i;
    c1i = j00i*j11r - j00r*j11i;
    c2r = j01r*j10r + j01i*j10i;
    c2i = j01i*j10r - j01r*j10i;
    out[10] += ( c1r + c2r) * weight;
    out[11] += ( c1i - c2i) * weight;
    out[14] += (-c1i - c2i) * weight;
    out[15] += ( c1r - c2r) * weight;
}

/**
 @brief Right multiplication of matrix \b m by <i> rotation matrix </i> with cos(f)=cs, sin(f)=sn.

 Rotation matrix
  @code
	| 1	0	0	0 |
	| 0	cs	sn	0 |
	| 0	-sn	cs	0 |
	| 0	0	0	1 |
 @endcode
*/
void RightRotateMueller(matrix&, double, double);

/**
 @brief Left multiplication of matrix \b m by <i> rotation matrix </i> with cos(f)=cs, sin(f)=sn.

 Rotation matrix
  @code
	| 1	0	0	0 |
	| 0	cs	sn	0 |
	| 0	-sn	cs	0 |
	| 0	0	0	1 |
 @endcode
*/
void LeftRotateMueller(matrix&, double, double);

/**
 @brief Multiplication of matrix \b m by left and right <i> rotation matrixes </i>

 Right rotation matrix
  @code
	| 1	0	0	0 |
	| 0	cs	sn	0 |
	| 0	-sn	cs	0 |
	| 0	0	0	1 |
 @endcode

  Left rotation matrix
  @code
	| 1	 0		0	0 |
	| 0	 _cs	_sn	0 |
	| 0	-(_sn)	_cs	0 |
	| 0	 0		0	1 |
 @endcode
*/
void RotateMueller(matrix& m, double _cs, double _sn, double cs, double sn);

/**
 @brief Calculate Mueller matrix in forward direction for randomly oriented particle. <b>See equation  below.</b>

 Mueller matrix in forward direction for randomly oriented particle
  @code
	| M11	0			0			0 	|
	| 0		(M22+M33)/2	0			0 	|
	| 0		0			(M22+M33)/2	0 	|
	| 0		0			0			M44 |
 @endcode

*/
void ForwardScattering(matrix &);

/**
 @brief Calculate Mueller matrix in backward direction for randomly oriented particle. <b>See equation  below. </b>

 Mueller matrix in backward direction for randomly oriented particle
  @code
	| M11	0			0			0 	|
	| 0		(M22-M33)/2	0			0 	|
	| 0		0			(M33-M22)/2	0 	|
	| 0		0			0			M44 |
 @endcode

*/
void BackwardScattering(matrix &);

#endif

/** @}
*/

/** @}
*/
