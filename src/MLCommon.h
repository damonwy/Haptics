#pragma once
#ifndef MUSCLEMASS_SRC_MLCOMMON_H_
#define MUSCLEMASS_SRC_MLCOMMON_H_
#include <json.hpp>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen\Dense>
//#include <Eigen/Sparse>
#include <Eigen/Eigenvalues> 

#define GLM_FORCE_RADIANS
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
#include <complex>
#include <unsupported/Eigen/CXX11/Tensor>

#define MAX(a, b) ((a)>(b)?(a):(b))
#define ToRadian(x) (double)(((x) * M_PI / 180.0))

#include <omp.h>

//typedef Eigen::Triplet<double> T;



typedef Eigen::Matrix<int, 6, 1> Vector6i;
typedef Eigen::Matrix<float, 3, 1> Vector3f;
typedef Eigen::Matrix<int, 3, 1> Vector3i;
typedef Eigen::Matrix<double, 2, 1> Vector2d;

typedef Eigen::Matrix<double, 3, 1> Vector3d;
typedef Eigen::Matrix<double, 4, 1> Vector4d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<int, 9, 1> Vector9i;
typedef Eigen::Matrix<double, 9, 1> Vector9d;
typedef Eigen::Matrix<double, 12, 1> Vector12d;
typedef Eigen::Matrix<double, 81, 1> Vector81d;
typedef Eigen::Matrix<double, 108, 1> Vector108d;

typedef Eigen::Matrix<double, 3, 3> Matrix3d;
typedef Eigen::Matrix<double, 4, 4> Matrix4d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
typedef Eigen::Matrix<double, 9, 9> Matrix9d;
typedef Eigen::Matrix<double, 12, 12> Matrix12d;

typedef Eigen::Matrix<double, 3, 4> Matrix3x4d;
typedef Eigen::Matrix<double, 3, 6> Matrix3x6d;
typedef Eigen::Matrix<double, 3, 12> Matrix3x12d;

typedef Eigen::Matrix<double, 4, 2> Matrix4x2d;
typedef Eigen::Matrix<double, 4, 3> Matrix4x3d;

typedef Eigen::Matrix<double, 5, 6> Matrix5x6d;

typedef Eigen::Matrix<double, 6, 2> Matrix6x2d;
typedef Eigen::Matrix<double, 6, 3> Matrix6x3d;

typedef Eigen::Matrix<double, 9, 12> Matrix9x12d;

typedef Eigen::TensorFixedSize<double, Eigen::Sizes<4, 4, 6>> Tensor4x4x6d;
typedef Eigen::TensorFixedSize<double, Eigen::Sizes<6, 2, 2>> Tensor6x2x2d;


enum Integrator { REDMAX_EULER, REDUCED_ODE45, REDMAX_ODE45 };
enum Axis {X_AXIS, Y_AXIS, Z_AXIS};

struct Energy {
	double K;
	double V;
};
enum SparseSolver {CG, CG_ILUT, QR, BICG,BICG_ILUT, SLDLT, LU, PARDISO_LU, PARDISO_LDLT, MINRES_SOLVER, GMRES_SOLVER, SUPER_LU
};

template<typename T>
using  MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename Scalar, int rank, typename sizeType>
inline auto Tensor_to_Matrix(const Eigen::Tensor<Scalar, rank> &tensor, const sizeType rows, const sizeType cols)
{
	return Eigen::Map<const MatrixType<Scalar>>(tensor.data(), rows, cols);
}

template<typename Scalar, typename... Dims>
inline auto Matrix_to_Tensor(const MatrixType<Scalar> &matrix, Dims... dims)
{
	constexpr int rank = sizeof... (Dims);
	return Eigen::TensorMap<Eigen::Tensor<const Scalar, rank>>(matrix.data(), { dims... });
}



inline Eigen::Vector3d findOrthonormalVector(Eigen::Vector3d input) {
	// Find a unit vector that is orthogonal to an input vector v

	// find smallest abs component of v
	int smallestIndex = 0;
	for (int dim = 1; dim<3; dim++)
		if (fabs(input[dim]) < fabs(input[smallestIndex]))
			smallestIndex = dim;

	Eigen::Vector3d axis;
	axis.setZero();
	axis(smallestIndex) = 1.0;

	// this cross-product will be non-zero (as long as v is not zero)
	Eigen::Vector3d result = input.cross(axis).normalized();
	return result;

}

inline Eigen::Matrix3d gs3(Eigen::Matrix3d A) {
	Eigen::Matrix3d R;
	Eigen::Vector3d a0 = A.col(0);
	Eigen::Vector3d a1 = A.col(1);
	Eigen::Vector3d r0 = a0 / a0.norm();
	Eigen::Vector3d r1 = a1 - r0.dot(a1) * r0;
	r1 = r1 / r1.norm();
	R.col(0) = r0;
	R.col(1) = r1;
	R.col(2) = r0.cross(r1);
	return R;
}

inline Matrix6d hooke(double young, double poisson) {
	double v = poisson;
	double v1 = 1 - v;
	double v2 = 1 - 2 * v;
	double s = young / ((1 + v)*(1 - 2 * v));

	Matrix6d E;
	E <<
		v1, v, v, 0, 0, 0,
		v, v1, v, 0, 0, 0,
		v, v, v1, 0, 0, 0,
		0, 0, 0, v2, 0, 0,
		0, 0, 0, 0, v2, 0,
		0, 0, 0, 0, 0, v2;

	E = s * E;
	return E;
}

inline static double hypot2(double x, double y){
	return sqrt(x*x + y*y);
}

// Eigen types to/from GLM types
inline glm::mat3 eigen_to_glm(const Eigen::Matrix3d &m) {
	return glm::make_mat3x3((const double *)m.data());
}

inline glm::mat4 eigen_to_glm(const Eigen::Matrix4d &m) {
	return glm::make_mat4x4((const double *)m.data());
}

inline glm::vec3 eigen_to_glm(const Eigen::Vector3d &v) {
	return glm::vec3(v[0], v[1], v[2]);
}

inline Eigen::Vector3d glm_to_eigen(const glm::vec3 &v) {
	return Eigen::Vector3d(v[0], v[1], v[2]);
}

inline Eigen::Matrix3d glm_to_eigen(const glm::mat3 &m) {
	const float *_m = glm::value_ptr(m);
	float m_cp[3 * 3];
	memcpy(m_cp, _m, sizeof(m_cp));
	Eigen::Map<Eigen::Matrix3f> result(m_cp);
	return result.cast<double>();
}

inline Eigen::Matrix4d glm_to_eigen(const glm::mat4 &m) {
	const float *_m = glm::value_ptr(m);
	float m_cp[4 * 4];
	memcpy(m_cp, _m, sizeof(m_cp));
	Eigen::Map<Eigen::Matrix4f> result(m_cp);
	return result.cast<double>();
}


inline bool rayTriangleIntersects(Eigen::Vector3d v1, Eigen::Vector3d v2, Eigen::Vector3d v3, Eigen::Vector3d dir, Eigen::Vector3d pos, double &t, double &u, double &v) {

	Eigen::Vector3d e1 = v2 - v1;
	Eigen::Vector3d e2 = v3 - v1;

	// Calculate planes normal vector
	//cross product
	Eigen::Vector3d pvec = dir.cross(e2);

	//dot product
	double det = e1.dot(pvec);

	// Ray is parallel to plane
	if (det <1e-8 && det > -1e-8) {
		return false;
	}

	double inv_det = 1 / det;

	// Distance from v1 to ray pos
	Eigen::Vector3d tvec = pos - v1;
	u = (tvec.dot(pvec))*inv_det;
	if (u < 0 || u > 1) {
		return false;
	}

	Eigen::Vector3d qvec = tvec.cross(e1);
	v = dir.dot(qvec) * inv_det;
	if (v<0 || u + v>1) {
		return false;
	}

	t = e2.dot(qvec) * inv_det;
	if (t > 1e-8) { return true; }
	return false;
}

inline static void tql2(double V[3][3], double d[3], double e[3]) {

	//  This is derived from the Algol procedures tql2, by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	for (int i = 1; i < 3; i++) {
		e[i - 1] = e[i];
	}
	e[3 - 1] = 0.0;

	double f = 0.0;
	double tst1 = 0.0;
	double eps = pow(2.0, -52.0);
	for (int l = 0; l < 3; l++) {

		// Find small subdiagonal element

		tst1 = MAX(tst1, fabs(d[l]) + fabs(e[l]));
		int m = l;
		while (m <3) {
			if (fabs(e[m]) <= eps*tst1) {
				break;
			}
			m++;
		}

		// If m == l, d[l] is an eigenvalue,
		// otherwise, iterate.

		if (m > l) {
			int iter = 0;
			do {
				iter = iter + 1;  // (Could check iteration count here.)

								  // Compute implicit shift

				double g = d[l];
				double p = (d[l + 1] - g) / (2.0 * e[l]);
				double r = hypot2(p, 1.0);
				if (p < 0) {
					r = -r;
				}
				d[l] = e[l] / (p + r);
				d[l + 1] = e[l] * (p + r);
				double dl1 = d[l + 1];
				double h = g - d[l];
				for (int i = l + 2; i < 3; i++) {
					d[i] -= h;
				}
				f = f + h;

				// Implicit QL transformation.

				p = d[m];
				double c = 1.0;
				double c2 = c;
				double c3 = c;
				double el1 = e[l + 1];
				double s = 0.0;
				double s2 = 0.0;
				for (int i = m - 1; i >= l; i--) {
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = hypot2(p, e[i]);
					e[i + 1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i + 1] = h + s * (c * g + s * d[i]);

					// Accumulate transformation.

					for (int k = 0; k < 3; k++) {
						h = V[k][i + 1];
						V[k][i + 1] = s * V[k][i] + c * h;
						V[k][i] = c * V[k][i] - s * h;
					}
				}
				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;

				// Check for convergence.

			} while (fabs(e[l]) > eps*tst1);
		}
		d[l] = d[l] + f;
		e[l] = 0.0;
	}

	// Sort eigenvalues and corresponding vectors.

	for (int i = 0; i < 3 - 1; i++) {
		int k = i;
		double p = d[i];
		for (int j = i + 1; j < 3; j++) {
			if (d[j] < p) {
				k = j;
				p = d[j];
			}
		}
		if (k != i) {
			d[k] = d[i];
			d[i] = p;
			for (int j = 0; j <3; j++) {
				p = V[j][i];
				V[j][i] = V[j][k];
				V[j][k] = p;
			}
		}
	}
}

inline static void tred2(double V[3][3], double d[3], double e[3]) {

	//  This is derived from the Algol procedures tred2 by
	//  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
	//  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
	//  Fortran subroutine in EISPACK.

	for (int j = 0; j < 3; j++) {
		d[j] = V[3 - 1][j];
	}

	// Householder reduction to tridiagonal form.

	for (int i = 3 - 1; i > 0; i--) {

		// Scale to avoid under/overflow.

		double scale = 0.0;
		double h = 0.0;
		for (int k = 0; k < i; k++) {
			scale = scale + fabs(d[k]);
		}
		if (scale == 0.0) {
			e[i] = d[i - 1];
			for (int j = 0; j < i; j++) {
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
				V[j][i] = 0.0;
			}
		}
		else {

			// Generate Householder vector.

			for (int k = 0; k < i; k++) {
				d[k] /= scale;
				h += d[k] * d[k];
			}
			double f = d[i - 1];
			double g = sqrt(h);
			if (f > 0) {
				g = -g;
			}
			e[i] = scale * g;
			h = h - f * g;
			d[i - 1] = f - g;
			for (int j = 0; j < i; j++) {
				e[j] = 0.0;
			}

			// Apply similarity transformation to remaining columns.

			for (int j = 0; j < i; j++) {
				f = d[j];
				V[j][i] = f;
				g = e[j] + V[j][j] * f;
				for (int k = j + 1; k <= i - 1; k++) {
					g += V[k][j] * d[k];
					e[k] += V[k][j] * f;
				}
				e[j] = g;
			}
			f = 0.0;
			for (int j = 0; j < i; j++) {
				e[j] /= h;
				f += e[j] * d[j];
			}
			double hh = f / (h + h);
			for (int j = 0; j < i; j++) {
				e[j] -= hh * d[j];
			}
			for (int j = 0; j < i; j++) {
				f = d[j];
				g = e[j];
				for (int k = j; k <= i - 1; k++) {
					V[k][j] -= (f * e[k] + g * d[k]);
				}
				d[j] = V[i - 1][j];
				V[i][j] = 0.0;
			}
		}
		d[i] = h;
	}

	// Accumulate transformations.

	for (int i = 0; i < 3 - 1; i++) {
		V[3 - 1][i] = V[i][i];
		V[i][i] = 1.0;
		double h = d[i + 1];
		if (h != 0.0) {
			for (int k = 0; k <= i; k++) {
				d[k] = V[k][i + 1] / h;
			}
			for (int j = 0; j <= i; j++) {
				double g = 0.0;
				for (int k = 0; k <= i; k++) {
					g += V[k][i + 1] * V[k][j];
				}
				for (int k = 0; k <= i; k++) {
					V[k][j] -= g * d[k];
				}
			}
		}
		for (int k = 0; k <= i; k++) {
			V[k][i + 1] = 0.0;
		}
	}
	for (int j = 0; j < 3; j++) {
		d[j] = V[3 - 1][j];
		V[3 - 1][j] = 0.0;
	}
	V[3 - 1][3 - 1] = 1.0;
	e[0] = 0.0;
}

inline void eigen_decomposition(double A[3][3], double V[3][3], double d[3]) {
	double e[3];
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			V[i][j] = A[i][j];
		}
	}
	tred2(V, d, e);
	tql2(V, d, e);
}

inline void eigen_sym(Eigen::Matrix3d &a, Eigen::Vector3d &eig_val, Eigen::Matrix3d &eig_vec) {
	double A[3][3] = { { a(0,0), a(0,1), a(0,2) },
	{ a(1,0), a(1,1), a(1,2) },
	{ a(2,0), a(2,1), a(2,2) } };
	double V[3][3];
	double d[3];
	eigen_decomposition(A, V, d);

	eig_val = Eigen::Vector3d(d[2], d[1], d[0]);

	eig_vec.col(0) = Eigen::Vector3d(V[0][2], V[1][2], V[2][2]);
	eig_vec.col(1) = Eigen::Vector3d(V[0][1], V[1][1], V[2][1]);
	eig_vec.col(2) = Eigen::Vector3d(V[0][0], V[1][0], V[2][0]);
}

inline int SVD(Eigen::Matrix3d &F,
	Eigen::Matrix3d &U,
	Eigen::Vector3d &Sigma,
	Eigen::Matrix3d &V,
	double sv_eps,
	int modifiedSVD) {

	// Adapted from Jernej Barbic's code
	// https://github.com/starseeker/VegaFEM/blob/master/libraries/minivector/mat3d.cpp

	// The code handles the following special situations:

	//---------------------------------------------------------
	// 1. det(V) == -1
	//    - multiply the first column of V by -1
	//---------------------------------------------------------
	// 2. An entry of Sigma is near zero
	//---------------------------------------------------------
	// (if modifiedSVD == 1) :
	// 3. negative determinant (Tet is inverted in solid mechanics).
	//    - check if det(U) == -1
	//    - If yes, then negate the minimal element of Sigma
	//      and the corresponding column of U
	//---------------------------------------------------------

	// form F^T F and do eigendecomposition

	Eigen::Matrix3d normalEq;
	normalEq.noalias() = F.transpose() * F;
	Eigen::Vector3d eigenValues;
	Eigen::Matrix3d eigenVectors;

	eigen_sym(normalEq, eigenValues, eigenVectors);
	V = eigenVectors;

	// Handle situation:
	// 1. det(V) == -1
	//    - multiply the first column of V by -1
	if (V.determinant() < -0.0000001) {
		// convert V into a rotation (multiply column 1 by -1)
		V.col(0) *= -1.0;
	}

	Sigma(0) = (eigenValues(0) > 0.0) ? sqrt(eigenValues(0)) : 0.0;
	Sigma(1) = (eigenValues(1) > 0.0) ? sqrt(eigenValues(1)) : 0.0;
	Sigma(2) = (eigenValues(2) > 0.0) ? sqrt(eigenValues(2)) : 0.0;

	// compute inverse of singular values
	// also check if singular values are close to zero
	Eigen::Vector3d SigmaInverse;
	SigmaInverse(0) = (Sigma(0) > sv_eps) ? (1.0 / Sigma(0)) : 0.0;
	SigmaInverse(1) = (Sigma(1) > sv_eps) ? (1.0 / Sigma(1)) : 0.0;
	SigmaInverse(2) = (Sigma(2) > sv_eps) ? (1.0 / Sigma(2)) : 0.0;

	// compute U using the formula:
	// U = F * V * diag(SigmaInverse)
	U.noalias() = F * V;
	Eigen::Matrix3d Utemp;
	Utemp.noalias() = U * SigmaInverse.asDiagonal();
	U = Utemp;

	// In theory, U is now orthonormal, U^T U = U U^T = I .. it may be a rotation or a reflection, depending on F.
	// But in practice, if singular values are small or zero, it may not be orthonormal, so we need to fix it.
	// Handle situation:
	// 2. An entry of Sigma is near zero
	// ---------------------------------------------------------
	if ((Sigma(0) < sv_eps) && (Sigma(1) < sv_eps) && (Sigma(2) < sv_eps))
	{
		// extreme case, all singular values are small, material has collapsed almost to a point
		// see [Irving 04], p. 4
		U.setIdentity();
	}
	else
	{
		// handle the case where two singular values are small, but the third one is not
		// handle it by computing two (arbitrary) vectors orthogonal to the eigenvector for the large singular value
		int done = 0;
		for (int dim = 0; dim<3; dim++)
		{
			int dimA = dim;
			int dimB = (dim + 1) % 3;
			int dimC = (dim + 2) % 3;
			if ((Sigma(dimB) < sv_eps) && (Sigma(dimC) < sv_eps))
			{
				// only the column dimA can be trusted, columns dimB and dimC correspond to tiny singular values
				Eigen::Vector3d tmpVec1 = U.col(dimA); // column dimA
				Eigen::Vector3d tmpVec2;
				tmpVec2 = findOrthonormalVector(tmpVec1);
				Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2).normalized();
				U.col(dimB) = tmpVec2;
				U.col(dimC) = tmpVec3;

				if (U.determinant() < -0.0000001)
				{
					U.col(dimB) *= -1.0;
				}
				done = 1;
				break; // out of for
			}
		}


		// handle the case where one singular value is small, but the other two are not
		// handle it by computing the cross product of the two eigenvectors for the two large singular values
		if (!done)
		{
			for (int dim = 0; dim<3; dim++)
			{
				int dimA = dim;
				int dimB = (dim + 1) % 3;
				int dimC = (dim + 2) % 3;

				if (Sigma(dimA) < sv_eps)
				{
					// columns dimB and dimC are both good, but column dimA corresponds to a tiny singular value
					Eigen::Vector3d tmpVec1 = U.col(dimB); // column dimB
					Eigen::Vector3d tmpVec2 = U.col(dimC); // column dimC
					Eigen::Vector3d tmpVec3 = tmpVec1.cross(tmpVec2).normalized();
					U.col(dimA) = tmpVec3;

					if (U.determinant() < -0.0000001)
					{
						U.col(dimA) *= -1.0;
					}

					done = 1;
					break; // out of for
				}
			}
		}

		if ((!done) && (modifiedSVD == 1))
		{
			// Handle situation:
			// 3. negative determinant (Tet is inverted in solid mechanics)
			//    - check if det(U) == -1
			//    - If yes, then negate the minimal element of Sigma
			//      and the corresponding column of U

			double detU = U.determinant();
			if (detU < 0.0)
			{
				// negative determinant
				// find the smallest singular value (they are all non-negative)
				int smallestSingularValueIndex = 0;
				for (int dim = 1; dim<3; dim++)
					if (Sigma(dim) < Sigma(smallestSingularValueIndex))
						smallestSingularValueIndex = dim;

				// negate the smallest singular value
				Sigma(smallestSingularValueIndex) *= -1.0;
				U.col(smallestSingularValueIndex) *= -1.0;

			}
		}
	}

	return 0;
}

//void eigen_sym(Eigen::Matrix3d &a, Eigen::Vector3d &eig_val, Eigen::Matrix3d &eig_vec) {
//	Eigen::EigenSolver<Eigen::Matrix3d> es(a);
//
//	for (int i = 0; i < 3; i++) {
//		std::complex<double> ev = es.eigenvalues()[i];
//		eig_val(i) = ev.real();
//		Eigen::Vector3cd v = es.eigenvectors().col(i);
//		eig_vec(0, i) = v(0).real();
//		eig_vec(1, i) = v(1).real();
//		eig_vec(2, i) = v(2).real();
//	}
//
//
//}

const int MIN_ITERATOR_NUM = 4;
inline int getThreadsNumber(int n, int min_n) {
	int ncore = omp_get_num_procs();
	int max_tn = n / min_n;
	int tn = max_tn > 2 * ncore ? 2 * ncore : max_tn;
	if (tn < 1) {
		tn = 1;
	}
	return tn;
}

#endif // MUSCLEMASS_SRC_MLCOMMON_H_