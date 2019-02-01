#pragma once

#ifndef MUSCLEMASS_SRC_WORLD_H_
#define MUSCLEMASS_SRC_WORLD_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"

class Joint;
class JointRevolute;
class JointRevoluteHyperReduced;
class JointUniversalXY;

class Body;
class MatrixStack;
class Program;
class Constraint;
class ConstraintJointLimit;
class ConstraintPrescJoint;
class ConstraintPrescBody;
class ConstraintPrescBodyAttachPoint;
class ConstraintNull;
class Spring;
class SpringNull;
class SpringDamper;
class Deformable;
class DeformableSpring;
class DeformableNull;
class JointNull;
class JointFixed;
class Comp;
class CompNull;
class CompSphere;
class CompCylinder;
class CompDoubleCylinder;
class WrapObst;
class WrapSphere;
class WrapCylinder;
class WrapDoubleCylinder;
class Muscle;
class MuscleSpring;
class MuscleNull;
class Line;
class Node;
class JointFree;

enum WorldType {
	SERIAL_CHAIN,
	DIFF_REVOLUTE_AXES,
	BRANCHING,
	SHPERICAL_JOINT,
	LOOP,
	JOINT_TORQUE,
	JOINT_LIMITS,
	EQUALITY_CONSTRAINED_ANGLES,
	EQUALITY_AND_LOOP, 
	HYBRID_DYNAMICS,
	EXTERNAL_WORLD_FORCE,
	JOINT_STIFFNESS,
	SPRINGS,
	COMPONENT,
	WRAP_SPHERE,
	WRAP_CYLINDER,
	WRAP_DOUBLECYLINDER,
	SPRING_DAMPER,
	FREEJOINT,
	TEST_MAXIMAL_HYBRID_DYNAMICS,
	TEST_REDUCED_HYBRID_DYNAMICS,
	TEST_HYPER_REDUCED_COORDS,
	TEST_JOINT_UNIVERSAL,
	TEST_CONSTRAINT_PRESC_BODY_ATTACH_POINT,
	MUSCLE_INERTIA
};

struct Floor {
	float y;
	Eigen::Vector2f xrange;
	Eigen::Vector2f zrange;
	Floor(float y_, Eigen::Vector2f x, Eigen::Vector2f z):
		y(y_), xrange(x), zrange(z)
	{

	}
};

class World
{
public:
	World();
	World(WorldType type);
	virtual ~World() {}

	std::shared_ptr<Body> addBody(
		double density, 
		Vector3d sides, 
		Vector3d p, 
		Matrix3d R, 
		const std::string &RESOURCE_DIR, 
		std::string file_name);

	std::shared_ptr<JointRevolute> addJointRevolute(
		std::shared_ptr<Body> body, 
		Vector3d axis, 
		Vector3d p, 
		Matrix3d R, 
		double q, 
		const std::string &RESOURCE_DIR,
		std::shared_ptr<Joint> parent=nullptr);

	std::shared_ptr<JointRevoluteHyperReduced> addJointRevoluteHyperReduced(
		std::shared_ptr<Body> body,
		Vector3d axis,
		std::shared_ptr<Joint> friend_joint,
		double scalar,
		Vector3d p,
		Matrix3d R,
		double q,
		const std::string &RESOURCE_DIR,
		std::shared_ptr<Joint> parent = nullptr);

	std::shared_ptr<JointUniversalXY> addJointUniversalXY(
		std::shared_ptr<Body> body,
		Vector3d p,
		Matrix3d R,
		const std::string &RESOURCE_DIR,
		std::shared_ptr<Joint> parent = nullptr);

	std::shared_ptr<JointFixed> addJointFixed(
		std::shared_ptr<Body> body, 
		Vector3d p, 
		Matrix3d R, 
		double q, 
		std::shared_ptr<Joint> parent = nullptr);
	
	std::shared_ptr<JointFree> addJointFree(
		std::shared_ptr<Body> body,
		Vector3d p,
		Matrix3d R,
		Vector6d q0,
		Vector6d qdot0,
		const std::string &RESOURCE_DIR,
		std::shared_ptr<Joint> parent = nullptr
	);

	std::shared_ptr<ConstraintJointLimit> addConstraintJointLimit(
		std::shared_ptr<Joint> joint, 
		double ql, 
		double qu);
	
	std::shared_ptr<DeformableSpring> addDeformableSpring(
		double mass, 
		int n_points, 
		std::shared_ptr<Body> body0, 
		Vector3d r0, 
		std::shared_ptr<Body> body1, 
		Vector3d r1);

	std::shared_ptr<ConstraintNull> addConstraintNull();
	std::shared_ptr<DeformableNull> addDeformableNull();
	std::shared_ptr<JointNull> addJointNull();
	std::shared_ptr<CompNull> addCompNull();
	std::shared_ptr<WrapObst> addWrapNull();
	std::shared_ptr<SpringNull> addSpringNull();

	std::shared_ptr<CompSphere> addCompSphere(
		double r, 
		std::shared_ptr<Body> parent, 
		Matrix4d E, 
		const std::string &RESOURCE_DIR);
	
	std::shared_ptr<CompCylinder> addCompCylinder(
		double r, 
		std::shared_ptr<Body> parent, 
		Matrix4d E, 
		Vector3d z, 
		Vector3d o, 
		const std::string &RESOURCE_DIR, 
		std::string shape);
	
	std::shared_ptr<CompDoubleCylinder> addCompDoubleCylinder(
		double rA, 
		std::shared_ptr<Body> parentA, 
		Matrix4d EA, 
		Vector3d z_a, 
		Vector3d o_a,
		double rB, 
		std::shared_ptr<Body> parentB, 
		Matrix4d EB,
		Vector3d z_b, 
		Vector3d o_b,
		const std::string &RESOURCE_DIR, 
		std::string shapeA, 
		std::string shapeB);
	
	std::shared_ptr<WrapSphere> addWrapSphere(
		std::shared_ptr<Body> b0, 
		Vector3d r0, 
		std::shared_ptr<Body> b1, 
		Vector3d r1, 
		std::shared_ptr<CompSphere> compSphere, 
		int num_points, 
		const std::string &RESOURCE_DIR);

	std::shared_ptr<WrapCylinder> addWrapCylinder(
		std::shared_ptr<Body> b0, 
		Vector3d r0, 
		std::shared_ptr<Body> b1, 
		Vector3d r1, 
		std::shared_ptr<CompCylinder> compCylinder, 
		int num_points, 
		const std::string &RESOURCE_DIR);

	std::shared_ptr<WrapDoubleCylinder> addWrapDoubleCylinder(
		std::shared_ptr<Body> b0,
		Vector3d r0,
		std::shared_ptr<Body> b1,
		Vector3d r1,
		Vector3d u,
		Vector3d v,
		Vector3d z_u,
		Vector3d z_v,
		std::shared_ptr<CompDoubleCylinder> compDoubleCylinder,
		int num_points,
		const std::string &RESOURCE_DIR);
	std::shared_ptr<ConstraintPrescJoint> addConstraintPrescJoint(
		std::shared_ptr<Joint> j);

	std::shared_ptr<ConstraintPrescBody> addConstraintPrescBody(
		std::shared_ptr<Body> b,
		Eigen::VectorXi dof
	);

	std::shared_ptr<ConstraintPrescBodyAttachPoint> addConstraintPrescBodyAttachPoint(
		std::shared_ptr<Body> b,
		Eigen::Vector3d r,
		Eigen::VectorXi dof
	);

	void sceneTestReducedHD(double t);
	void sceneTestMaximalHD(double t);
	void sceneTestHyperReduced(double t);
	void sceneAttachPoint(double t);

	void setMaximalPrescStates(std::shared_ptr<Body> b, Vector3d vt_w, Vector3d vtdot_w, Vector3d wt_i, Vector3d wtdot_i);
	void setMaximalPrescStates(int index_body, Vector3d vt_w, Vector3d wt_i);
	void setMaximalPrescAttachPointStates(int index_body, int index_point, Vector3d vt_w);
	
	void setReducedPrescStates(std::shared_ptr<Joint> j, Eigen::VectorXd q, Eigen::VectorXd dq);
	void setReducedPrescStates(std::shared_ptr<Joint> j, double q, double dq);

	void setListOfReducedPrescStates(Eigen::VectorXi rcon, Eigen::VectorXd q, Eigen::VectorXd dq);
	void setListOfMaximalPrescStates(Eigen::VectorXi mcon, Vector3d vt_w, Vector3d vtdot_w, Vector3d wt_i, Vector3d wtdot_i);
	void computeTargetQ(double t0, double t1, double t, double angle, double q0, double &q, double &dq);
	
	void deactivateListOfPrescConstraints(Eigen::VectorXi mcon, Eigen::VectorXi rcon);
	Energy computeEnergy();

	void load(const std::string &RESOURCE_DIR);
	void init();
	void update();
	
	void draw(
		std::shared_ptr<MatrixStack> MV, 
		const std::shared_ptr<Program> prog, 
		const std::shared_ptr<Program> progSimple, 
		const std::shared_ptr<Program> progSoft, 
		std::shared_ptr<MatrixStack> P);

	void drawFloor(
		Floor f,
		std::shared_ptr<MatrixStack> MV, 
		const std::shared_ptr<Program> progSimple, 
		std::shared_ptr<MatrixStack> P);

	void setTime(double t) { m_t = t; }
	double getTime() const { return m_t; }
	double getH() const { return m_h; }
	void incrementTime() { m_t += m_h; }

	void setGrav(Eigen::Vector3d grav) { m_grav = grav; }
	Vector3d getGrav() const { return m_grav; }

	std::shared_ptr<Body> getBody(int uid);
	std::shared_ptr<Body> getBody(const std::string &name);
	std::shared_ptr<Constraint> getConstraint(int uid);
	std::shared_ptr<Constraint> getConstraint(const std::string &name);

	std::shared_ptr<Body> getBody0() const { return m_bodies[0]; }
	std::shared_ptr<Joint> getJoint0() const { return m_joints[0]; }
	std::shared_ptr<Deformable> getDeformable0() const { return m_deformables[0]; }
	std::shared_ptr<Constraint> getConstraint0() const { return m_constraints[0]; }
	std::shared_ptr<Spring> getSpring0() const { return m_springs[0]; }

	Vector2d getTspan() const { return m_tspan; }
	int getNsteps();

	int nm;
	int nr;
	int nR;
	int nem;
	int ner;
	int ne;
	int nim;
	int nir;
	int m_countS;
	int m_countCM;

	int m_nbodies;
	int m_njoints;
	int m_ndeformables;
	int m_nsprings;
	int m_nconstraints;
	int m_ncomps;
	int m_nwraps;
	int m_nmuscles;
	
	int m_dense_nm;
	int m_dense_nr;
	int m_ntets;

	int nsegments;
	WorldType m_type;

private:
	Energy m_energy;		// the energy in current state
	Energy m_energy0;		// the energy in initial state

	Eigen::Vector3d m_grav;
	double m_stiffness;
	double m_damping;
	double m_t;
	double m_h;
	Eigen::Vector2d m_tspan;	

	double m_Hexpected;		// used to check correctness
	
	// These are the actual objects that are created
	std::vector<std::shared_ptr<Body>> m_bodies;
	std::vector<std::shared_ptr<Comp>> m_comps;
	std::vector<std::shared_ptr<WrapObst>> m_wraps;
	std::vector<std::shared_ptr<Muscle>> m_muscles;
	std::vector<std::shared_ptr<Joint>> m_joints;
	std::vector<std::shared_ptr<Deformable>> m_deformables;
	std::vector<std::shared_ptr<Spring>> m_springs;
	std::vector<std::shared_ptr<Constraint>> m_constraints;
	std::vector<std::shared_ptr<Line> > m_lines;
	std::vector<Floor> m_floors;

	// These are extra data structures that point to the objects
	typedef std::map<std::string, std::shared_ptr<Body>> MapBodyName;
	typedef std::map<std::string, std::shared_ptr<Constraint>> MapConstraintName;
	typedef std::map<int, std::shared_ptr<Body>> MapBodyUID;
	typedef std::map<int, std::shared_ptr<Constraint>> MapConstraintUID;

	MapBodyName m_bodyName;
	MapBodyUID m_bodyUID;
	MapConstraintName m_constraintName;
	MapConstraintUID m_constraintUID;
};

#endif // MUSCLEMASS_SRC_WORLD_H_