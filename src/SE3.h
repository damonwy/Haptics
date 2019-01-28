#pragma once

#ifndef MUSCLEMASS_SRC_SE3_H_
#define MUSCLEMASS_SRC_SE3_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

#define THRESH 1.0e-9

class SE3 {

public:
	static Eigen::Matrix4d inverse(const Eigen::Matrix4d &E);
	static Matrix3x6d gamma(const Eigen::Vector3d &r);
	static Matrix6d adjoint(const Eigen::Matrix4d &E);
	static Eigen::Matrix3d bracket3(const Eigen::Vector3d &a);
	static Eigen::Matrix4d bracket6(const Vector6d &a);
	static Eigen::Vector3d unbracket3(const Eigen::Matrix3d &A);
	static Vector6d unbracket6(const Eigen::Matrix4d &A);
	static Eigen::Matrix4d integrate(const Eigen::Matrix4d &E0, const Eigen::VectorXd &phi, double h);
	static Matrix6d dAddt(const Eigen::Matrix4d &E, const Eigen::VectorXd &phi);
	static Vector6d inertiaCuboid(Eigen::Vector3d whd, double density);
	static Eigen::Matrix3d aaToMat(Eigen::Vector3d axis, double angle);
	static Eigen::Matrix4d RpToE(Eigen::Matrix3d R, Eigen::Vector3d p);
	static void EToRp(const Eigen::Matrix4d &E, Eigen::Matrix3d &R, Eigen::Vector3d &p);
	static Matrix6d ad(Vector6d phi);
	static Eigen::Matrix3d exp(const Vector3d &phi);
	static Eigen::Matrix4d exp(const Vector6d &phi);
	static Eigen::Matrix4d exp(const Eigen::Matrix4d &phi);
	static Vector6d log(const Eigen::Matrix4d &A);
	static bool reparam(Eigen::VectorXd &w);
};



#endif // MUSCLEMASS_SRC_SE3_H_