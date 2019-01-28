#pragma once
// Body A rigid body connected through a joint to a parent

#ifndef REDUCEDCOORD_SRC_BODY_H_
#define REDUCEDCOORD_SRC_BODY_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class Shape;
class Joint;
class Program;
class MatrixStack;
class ConstraintPrescBody;
class ConstraintPrescBodyAttachPoint;

typedef Eigen::Triplet<double> T;

class Body 
{
public:
	Body();
	Body(double density);
	virtual ~Body() {}

	void setDamping(double damping) { m_damping = damping; }
	void setTransform(Matrix4d E);	
	void setSides(Vector3d sides) { m_sides = sides; }
	void setJoint(std::shared_ptr<Joint> joint) { m_joint = joint; };
	void setAttachedColor(Vector3f color) { m_attached_color = color; }
	void setDrawingOption(bool drawing) { m_isDrawing = drawing; }
	std::string getName() const { return m_name; };
	std::shared_ptr<Joint> getJoint() const { return m_joint; };
	Matrix4d getEndPoint() { return (E_wi * E_ie); }
	Matrix4d getBodyByEndPoint(Matrix4d E_we) { return (E_we * E_ei); }
	Vector3d getBodyVelocityByEndPointVelocity(Vector3d v_e);
	void computeInertia();
	void countDofs(int &nm);
	int countM(int &nm, int data);
	void computeMass(Eigen::MatrixXd &M);
	void computeMassSparse(std::vector<T> &M_);
	void computeGrav(Vector3d grav, Eigen::VectorXd &f);
	void computeMassGrav(Vector3d grav, Eigen::MatrixXd &M, Eigen::VectorXd &f);
	void computeForceDamping(Eigen::VectorXd &f, Eigen::MatrixXd &D);
	void computeForceDampingSparse(Eigen::VectorXd &f, std::vector<T> &D_);
	void computeEnergies(Vector3d grav, Energy &energies);
	void setColor(Vector3f color) { m_body_color = color; };
	void load(const std::string &RESOURCE_DIR, std::string box_shape);
	void init(int &nm);
	void update();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;

	double m_density;					// Mass/volume
	Vector6d I_i;						// Inertia at body center
	Matrix6d M_i;						// Inertia Matrix at body center
	Matrix4d E_ji;						// Where the body is wrt joint
	Matrix4d E_ij;						// Where the joint is wrt body
	Matrix4d E_ie;						// Where the end point is wrt body
	Matrix4d E_ei;						// Where the body is wrt the end point
	Matrix4d E_wi;						// Where the body is wrt world
	Matrix4d E_iw;						// Where the world is wrt body
	Matrix4d E_ip;						// Where the parent is wrt body
	Matrix6d Ad_ji;						// Adjoint of E_ji
	Matrix6d Ad_ij;						// Adjoint of E_ij
	Matrix6d Ad_iw;						// Adjoint of E_iw
	Matrix6d Ad_wi;						// Adjoint of E_wi
	Matrix6d Ad_ip;						// Adjoint of E_ip
	Matrix6d Addot_wi;					// Adjoint dot of E_wi
	Vector6d phi;						// Twist at body center
	Vector6d phidot;					// Acceleration at body center
	Vector6d wext_i;					// External wrench in body space (not used by redmax)
	Matrix6d Kmdiag;					// Maximal stiffness block diagonal term
	Matrix6d Dmdiag;					// Maximal damping block diagonal term

	Vector6d fcor;
	Vector6d fgrav;
	Matrix3d R_wi;
	Matrix3d R_iw;
	bool   m_isDrawing;
	double m_damping;					// Viscous damping
	std::shared_ptr<Joint> m_joint;		// Joint to parent
	int idxM;							// Maximal indices
	std::shared_ptr<Body> next;			// Next body in traversal order

	std::shared_ptr<Body> m_parent;
	std::shared_ptr<ConstraintPrescBody> presc;						// Presribed motion constraint
	std::vector<std::shared_ptr<ConstraintPrescBodyAttachPoint> > m_presc_attach_points;

	Vector3f m_attached_color;
	Vector3f m_sliding_color;

	Vector3f m_body_color;
	Vector3d m_sides;
	void toggleDrawing(bool isDrawing) { m_isDrawing = isDrawing; }

protected:
	std::shared_ptr<Shape> bodyShape;
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;
	virtual void computeInertia_() {}
	std::string m_name;
	int m_UID;
};

#endif // REDUCEDCOORD_SRC_BODY_H_
