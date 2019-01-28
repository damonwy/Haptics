#include "rmpch.h"
#include "SE3.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Matrix4d SE3::inverse(const Matrix4d &E)
{
	Matrix4d Einv = Matrix4d::Identity();
	Matrix3d R = E.block<3, 3>(0, 0);
	Vector3d p = E.block<3, 1>(0, 3);
	Matrix3d Rt = R.transpose();
	Einv.block<3, 3>(0, 0) = Rt;
	Einv.block<3, 1>(0, 3) = -Rt * p;
	return Einv;
}

Matrix3x6d SE3::gamma(const Eigen::Vector3d &r)
{
	// Gets the 3x6 Gamma matrix, for computing the point velocity
	Matrix3x6d G = Matrix3x6d::Zero();
	G.block<3, 3>(0, 0) = bracket3(r).transpose();
	G.block<3, 3>(0, 3) = Matrix3d::Identity();
	return G;
}

Matrix6d SE3::adjoint(const Matrix4d &E)
{
	// Gets the adjoint transform
	Matrix6d Ad = Matrix6d::Zero();
	Matrix3d R = E.block<3, 3>(0, 0);
	Vector3d p = E.block<3, 1>(0, 3);
	Ad.block(0, 0, 3, 3) = R;
	Ad.block(3, 0, 3, 3) = bracket3(p) * R;
	Ad.block(3, 3, 3, 3) = R;
	return Ad;
}

Matrix6d SE3::dAddt(const Matrix4d &E, const VectorXd &phi) 
{
	// Gets the time derivative of the adjoint
	Matrix6d dA;
	dA.setZero();
	Matrix3d R = E.block<3, 3>(0, 0);
	Vector3d p = E.block<3, 1>(0, 3);
	Vector3d w = phi.segment<3>(0);
	Vector3d v = phi.segment<3>(3);
	Matrix3d wbrac = bracket3(w);
	Matrix3d vbrac = bracket3(v);
	Matrix3d pbrac = bracket3(p);
	Matrix3d Rwbrac = R * wbrac;
	dA.block<3, 3>(0, 0) = Rwbrac;
	dA.block<3, 3>(3, 3) = Rwbrac;
	dA.block<3, 3>(3, 0) = R * vbrac + pbrac * Rwbrac;

	return dA;
}

Matrix3d SE3::bracket3(const Vector3d &a)
{
	// Gets S = [x], the skew symmetric matrix
	Matrix3d A = Matrix3d::Zero();
	A(0, 1) = -a(2);
	A(0, 2) = a(1);
	A(1, 0) = a(2);
	A(1, 2) = -a(0);
	A(2, 0) = -a(1);
	A(2, 1) = a(0);
	return A;
}

Matrix4d SE3::bracket6(const Vector6d &a)
{
	Matrix4d A = Matrix4d::Zero();
	A.block<3, 3>(0, 0) = bracket3(a.segment<3>(0));
	A.block<3, 1>(0, 3) = a.segment<3>(3);
	return A;
}

Vector3d SE3::unbracket3(const Matrix3d &A)
{
	// Gets ]x[ = S, the vector corresponding to a skew symmetric matrix
	Vector3d a;
	a(0) = A(2, 1);
	a(1) = A(0, 2);
	a(2) = A(1, 0);
	return a;
}

Vector6d SE3::unbracket6(const Matrix4d &A)
{
	Vector6d a;
	a.segment<3>(0) = unbracket3(A.block<3, 3>(0, 0));
	a(3) = A(0, 3);
	a(4) = A(1, 3);
	a(5) = A(2, 3);
	return a;
}

Matrix4d SE3::integrate(const Matrix4d &E0, const VectorXd &phi, double h)
{
	Matrix3d I = Matrix3d::Identity();
	Vector3d w = phi.segment<3>(0);
	Vector3d v = phi.segment<3>(3);
	Matrix4d phib = Matrix4d::Identity();
	phib.block<3, 1>(0, 3) = h*v;
	double wlen = w.norm();
	if (wlen > 1e-10) {
		w /= wlen;
		v /= wlen;
		// Rodrigues formula %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		double wX = w(0);
		double wY = w(1);
		double wZ = w(2);
		double c = cos(wlen * h);
		double s = sin(wlen * h);
		double c1 = 1.0 - c;
		Matrix3d R;
		R << c + wX * wX * c1, -wZ * s + wX * wY * c1, wY * s + wX * wZ * c1,
			wZ * s + wX * wY * c1, c + wY * wY * c1, -wX * s + wY * wZ * c1,
			-wY * s + wX * wZ * c1, wX * s + wY * wZ * c1, c + wZ * wZ * c1;
		//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		Matrix3d A = I - R;
		Vector3d cc = w.cross(v);
		Vector3d d = A * cc;
		double wv = w.dot(v);
		Vector3d p = (wv * wlen * h) * w + d;
		phib.block<3, 3>(0, 0) = R;
		phib.block<3, 1>(0, 3) = p;
		//cout << phib << endl;
	}
	return E0 * phib;
}

Vector6d SE3::inertiaCuboid(Eigen::Vector3d whd, double density) {
	// Gets the diagonal inertia of a cuboid with (width, height, depth)
	Vector6d m;
	m.setZero();
	double volume = whd(0) * whd(1) * whd(2);
	double mass = density * volume;
	m(0) = (1.0 / 12.0) * mass * (whd(1) * whd(1) + whd(2) * whd(2));
	m(1) = (1.0 / 12.0) * mass * (whd(2) * whd(2) + whd(0) * whd(0));
	m(2) = (1.0 / 12.0) * mass * (whd(0) * whd(0) + whd(1) * whd(1));
	m(3) = mass;
	m(4) = mass;
	m(5) = mass;
	return m;
}


Matrix3d SE3::aaToMat(Vector3d axis, double angle) {
	// Creates a rotation matrix from an (axis, angle) pair
	Matrix3d R;
	R.setIdentity();

	double mag = axis.norm();
	if (mag > THRESH) {
		mag = 1.0 / mag;
		axis = axis * mag;
		if (abs(axis(0)) < THRESH && abs(axis(1)) < THRESH) {
			// Rotation about Z
			if (axis(2) < 0.0) {
				angle = -angle;
			}
			double sinTheta = sin(angle);
			double cosTheta = cos(angle);

			R(0, 0) = cosTheta;
			R(0, 1) = -sinTheta;
			R(1, 0) = sinTheta;
			R(1, 1) = cosTheta;

		}
		else if (abs(axis(1)) < THRESH && abs(axis(2)) < THRESH) {
			// Rotation about X
			if (axis(0) < 0.0) {
				angle = -angle;
			}
			double sinTheta = sin(angle);
			double cosTheta = cos(angle);
			R(1, 1) = cosTheta;
			R(1, 2) = -sinTheta;
			R(2, 1) = sinTheta;
			R(2, 2) = cosTheta;
		}
		else if (abs(axis(2)) < THRESH && abs(axis(0)) < THRESH) {
			// Rotation about Y
			if (axis(1) < 0.0) {
				angle = -angle;
			}
			double sinTheta = sin(angle);
			double cosTheta = cos(angle);
			R(0, 0) = cosTheta;
			R(0, 2) = sinTheta;
			R(2, 0) = -sinTheta;
			R(2, 2) = cosTheta;
		}
		else {
			// General rotation
			double sinTheta = sin(angle);
			double cosTheta = cos(angle);
			double t = 1.0 - cosTheta;
			double xz = axis(0) * axis(2);
			double xy = axis(0) * axis(1);
			double yz = axis(1) * axis(2);
			R(0, 0) = t * axis(0) * axis(0) + cosTheta;
			R(0, 1) = t * xy - sinTheta * axis(2);
			R(0, 2) = t * xz + sinTheta * axis(1);
			R(1, 0) = t * xy + sinTheta * axis(2);
			R(1, 1) = t * axis(1) * axis(1) + cosTheta;
			R(1, 2) = t * yz - sinTheta * axis(0);
			R(2, 0) = t * xz - sinTheta * axis(1);
			R(2, 1) = t * yz + sinTheta * axis(0);
			R(2, 2) = t * axis(2) * axis(2) + cosTheta;
		}
	}
	return R;
}

Matrix4d SE3::RpToE(Matrix3d R, Vector3d p) {
	Matrix4d E;
	E.setIdentity();
	E.block<3, 3>(0, 0) = R;
	E.block<3, 1>(0, 3) = p;
	return E;
}

void SE3::EToRp(const Eigen::Matrix4d &E, Eigen::Matrix3d &R, Eigen::Vector3d &p) {
	R = E.block<3, 3>(0, 0);
	p = E.block<3, 1>(0, 3);
}

Matrix6d SE3::ad(Vector6d phi) {
	Matrix6d a;
	a.setZero();

	Vector3d w = phi.segment<3>(0);
	Vector3d v = phi.segment<3>(3);

	Matrix3d W = bracket3(w);
	a.block<3, 3>(0, 0) = W;
	a.block<3, 3>(3, 3) = W;
	a.block<3, 3>(3, 0) = bracket3(v);
	return a;
}

Matrix4d SE3::exp(const Vector6d & phi)
{
	Matrix3d I = Matrix3d::Identity();
	Vector3d w = phi.segment<3>(0);
	double wlen = w.norm();
	Matrix3d R = Matrix3d::Identity();;
	if (wlen > 1e-9) {
		w = w / wlen;
		// Rodrigues forumula ----------------------
		double wX = w(0);
		double wY = w(1);
		double wZ = w(2);
		double c = cos(wlen);
		double s = sin(wlen);
		double c1 = 1 - c;
		Matrix3d fillR;
		fillR <<
			c + wX*wX*c1, -wZ*s + wX*wY*c1, wY*s + wX*wZ*c1,
			wZ*s + wX*wY*c1, c + wY*wY*c1, -wX*s + wY*wZ*c1,
			-wY*s + wX*wZ*c1, wX*s + wY*wZ*c1, c + wZ*wZ*c1;
		R = fillR;
		//------------------------------------------
	}
	Matrix4d E = Matrix4d::Identity();
	E.block<3, 3>(0, 0) = R;
	// Translational part
	Eigen::Vector3d v = phi.segment<3>(3);
	if (wlen > 1e-9) {
		v = v / wlen;
		Matrix3d A = I - R;
		Vector3d cc = w.cross(v);
		Vector3d d = A * cc;
		double wv = w.transpose() * v;
		Vector3d p = (wv * wlen) * w + d;
		E.block<3, 1>(0, 3) = p;
	}
	else {
		E.block<3, 1>(0, 3) = v;
	}

	return E;
}

Matrix3d SE3::exp(const Vector3d &phi) {
	Matrix3d I = Matrix3d::Identity();
	Vector3d w = phi.segment<3>(0);
	double wlen = w.norm();
	Matrix3d R = Matrix3d::Identity();;
	if (wlen > 1e-9) {
		w = w / wlen;
		// Rodrigues forumula ----------------------
		double wX = w(0);
		double wY = w(1);
		double wZ = w(2);
		double c = cos(wlen);
		double s = sin(wlen);
		double c1 = 1 - c;
		Matrix3d fillR;
		fillR <<
			c + wX*wX*c1, -wZ*s + wX*wY*c1, wY*s + wX*wZ*c1,
			wZ*s + wX*wY*c1, c + wY*wY*c1, -wX*s + wY*wZ*c1,
			-wY*s + wX*wZ*c1, wX*s + wY*wZ*c1, c + wZ*wZ*c1;
		R = fillR;
		//------------------------------------------
	}
	return R;
}


Matrix4d SE3::exp(const Matrix4d &phi) {
	// Convert from skew symmetric matrix to vector
	Vector6d phi_ = unbracket6(phi);
	Matrix4d E = exp(phi_);
	return E;
}

Vector6d SE3::log(const Matrix4d & A)
{
	Matrix3d R = A.block<3, 3>(0, 0);
	Vector3d p = A.block<3, 1>(0, 3);
	Vector6d phi;

	double cosTheta = 0.5*(R.trace() - 1);
	double theta = acos(cosTheta);

	if (cosTheta > 0.999999999999)
		theta = 0;
	else if (cosTheta < -0.999999999999)
		theta = M_PI;

	if (abs(theta) < 1e-8)
	{
		phi.segment<3>(0) << 0.0, 0.0, 0.0;
		phi.segment<3>(3) = p;
	}
	else
	{
		double sinTheta = sin(theta);
		Matrix3d wBracket = theta / (2 * sinTheta)*(R - R.transpose());
		phi.segment<3>(0) = unbracket3(wBracket);
		Matrix3d V = Eigen::Matrix3d::Identity() + ((1 - cosTheta) / (theta * theta) * wBracket) + ((theta - sinTheta) / (theta * theta * theta) * wBracket*wBracket);
		phi.segment<3>(3) = V.colPivHouseholderQr().solve(p);
	}

	return phi;
}

bool SE3::reparam(VectorXd &w) {
	bool flag = false;
	double wnorm = w.norm();
	while (wnorm > 1.5 * M_PI) {
		flag = true;
		double a = (1.0 - 2.0 * M_PI / wnorm);
		w *= a;
		wnorm = abs(a * wnorm);
	}
	return flag;
}