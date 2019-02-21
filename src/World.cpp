#include "rmpch.h"
#include "World.h"

#include "Joint.h"
#include "JointNull.h"
#include "JointFixed.h"
#include "JointRevolute.h"
#include "JointRevoluteHyperReduced.h"
#include "JointSphericalExp.h"
#include "JointUniversal.h"
#include "JointUniversalXY.h"

#include "JointFree.h"
#include "JointTranslational.h"

#include "Node.h"
#include "Body.h"
#include "BodyCuboid.h"

#include "JsonEigen.h"

#include "ConstraintJointLimit.h"
#include "ConstraintNull.h"
#include "ConstraintLoop.h"
#include "ConstraintAttachSpring.h"
#include "ConstraintPrescBody.h"
#include "ConstraintPrescJoint.h"
#include "ConstraintPrescBodyAttachPoint.h"

#include "Deformable.h"
#include "DeformableSpring.h"
#include "DeformableNull.h"

#include "Spring.h"
#include "SpringNull.h"
#include "SpringDamper.h"

#include "Comp.h"
#include "CompNull.h"
#include "CompSphere.h"
#include "CompCylinder.h"
#include "CompDoubleCylinder.h"

#include "WrapNull.h"
#include "WrapSphere.h"
#include "WrapCylinder.h"
#include "WrapDoubleCylinder.h"
#include "Vector.h"

#include "Muscle.h"
#include "MuscleNull.h"
#include "MuscleSpring.h"

#include "Line.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

World::World() :
nm(0), nr(0), nR(0), nem(0), ner(0), ne(0), nim(0), nir(0), m_countS(0),m_countCM(0), m_nbodies(0), m_njoints(0), m_ndeformables(0), m_nsprings(0), m_ncomps(0), m_nwraps(0), m_nconstraints(0), m_nmuscles(0)
{
	m_energy.K = 0.0;
	m_energy.V = 0.0;
}

World::World(WorldType type) :
nm(0), nr(0), nR(0), nem(0), ner(0), ne(0), nim(0), nir(0), m_countS(0),m_countCM(0),  m_nbodies(0), m_njoints(0), m_ndeformables(0), m_nsprings(0), m_ncomps(0), m_nwraps(0), m_type(type), m_nconstraints(0), m_nmuscles(0)

{
	m_energy.K = 0.0;
	m_energy.V = 0.0;
}

void World::load(const std::string &RESOURCE_DIR) {

	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	json js;
	i >> js;
	i.close();

	double density;
	Eigen::Vector3d sides;
	Matrix4d E;
	Vector3d p;

	switch (m_type)
	{
	case SERIAL_CHAIN:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 3; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		//auto con0 = make_shared<ConstraintPrescJoint>(m_joints[0], REDMAX_EULER);
		//m_nconstraints++;
		//m_constraints.push_back(con0);
		break;
	}
	case DIFF_REVOLUTE_AXES:
		break;
	case BRANCHING:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		Vector3d sides_0;
		sides_0 << 1.0, 10.0, 1.0;
		Vector3d sides_1;
		sides_1 << 20.0, 1.0, 1.0;

		auto b0 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b1 = addBody(density, sides_1, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box20_1_1.obj");
		auto b2 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b3 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b4 = addBody(density, sides, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
		auto b5 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b6 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b7 = addBody(density, sides, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
		auto b8 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b9 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");

		auto j0 = addJointRevolute(b0, Vector3d::UnitX(), Vector3d(0.0, 15.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
		auto j1 = addJointRevolute(b1, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j0);
		auto j2 = addJointRevolute(b2, Vector3d::UnitX(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j1);
		auto j3 = addJointRevolute(b3, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j1);
		auto j4 = addJointRevolute(b4, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j2);
		auto j5 = addJointRevolute(b5, Vector3d::UnitX(), Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j4);
		auto j6 = addJointRevolute(b6, Vector3d::UnitY(), Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j4);
		auto j7 = addJointRevolute(b7, Vector3d::UnitY(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j3);
		auto j8 = addJointRevolute(b8, Vector3d::UnitX(), Vector3d(-5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j7);
		auto j9 = addJointRevolute(b9, Vector3d::UnitY(), Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), M_PI / 4.0, RESOURCE_DIR, j7);
	}
	break;
	case SHPERICAL_JOINT:
		break;
	case LOOP:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		Vector3d sides_0;
		sides_0 << 1.0, 10.0, 1.0;
		Vector3d sides_1;
		sides_1 << 20.0, 1.0, 1.0;

		auto b0 = addBody(density, sides_1, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box20_1_1.obj");
		auto b1 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b2 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");
		auto b3 = addBody(density, sides_1, Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box20_1_1.obj");
		auto b4 = addBody(density, sides_0, Vector3d(0.0, -5.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box1_10_1.obj");

		auto j0 = addJointRevolute(b0, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
		auto j1 = addJointRevolute(b1, Vector3d::UnitZ(), Vector3d(-10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j0);
		auto j2 = addJointRevolute(b2, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j0);
		auto j3 = addJointRevolute(b3, Vector3d::UnitZ(), Vector3d(0.0, -10.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j1);
		auto j4 = addJointRevolute(b4, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, j3);
		j4->m_qdot(0) = 5.0;

		auto constraint = make_shared<ConstraintLoop>(b2, b3);
		m_constraints.push_back(constraint);
		constraint->setPositions(Vector3d(0.0, -5.0, 0.0), Vector3d(10.0, 0.0, 0.0));
		m_nconstraints++;

	}
	break;
	case JOINT_TORQUE:
		break;
	case JOINT_LIMITS:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 10; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}

			// Init constraints
			if (i > 0) {
				addConstraintJointLimit(m_joints[i], -M_PI / 4, M_PI / 4);
			}
		}
	}

	break;
	case EQUALITY_CONSTRAINED_ANGLES:
		break;
	case EQUALITY_AND_LOOP:
		break;
	case HYBRID_DYNAMICS:
		break;
	case EXTERNAL_WORLD_FORCE:
		break;
	case JOINT_STIFFNESS:
	{	
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, 0.0, 0.0;
		Eigen::from_json(js["sides"], sides);
		
		m_stiffness = 1.0e4;
		m_damping = 1.0e3;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 3; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
			m_joints[i]->setStiffness(m_stiffness);
			m_joints[i]->setDamping(m_damping);

		}

		m_joints[0]->m_qdot(0) = 1.0;

		break;

	}

		break;
	case SPRINGS:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		m_stiffness = 5.0e3;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		// Init springs
		auto deformable0 = addDeformableSpring(sides.x()*sides.y()*sides.z() * density, m_stiffness, 3, nullptr, Vector3d(10.0 * m_nbodies + 10.0, 10.0, 0.0), m_bodies[m_nbodies - 1], Vector3d(5.0, 0.0, 0.0), RESOURCE_DIR);
		deformable0->setStiffness(m_stiffness);
		auto deformable1 = addDeformableSpring(sides.x()*sides.y()*sides.z() * density, m_stiffness, 2, m_bodies[0], Vector3d(0.0, 0.0, 0.0), m_bodies[m_nbodies - 1], Vector3d(0.0, 0.0, 0.0), RESOURCE_DIR);
		deformable1->setStiffness(m_stiffness);
		for (int i = 0; i < (int)m_deformables.size(); i++) {
			m_deformables[i]->load(RESOURCE_DIR);
		}
	}
	break;

	case COMPONENT:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);

				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				joint->m_qdot(0) = -5.0;
			}
		}
		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(3.0, 1.0, 0.0));

		auto compSphere = addCompSphere(1.0, m_bodies[0], E, RESOURCE_DIR);
		auto compCylinder = addCompCylinder(1.0, m_bodies[1], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0), RESOURCE_DIR, "obstacle.obj");
		auto compDoubleCylinder = addCompDoubleCylinder(
			0.5, m_bodies[0], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0),
			0.5, m_bodies[2], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0),
			RESOURCE_DIR, "obstacle.obj", "obstacle.obj");
	}
	break;

	case WRAP_SPHERE:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				joint->m_qdot(0) = -0.0;
			}
		}

		auto compSphere = addCompSphere(1.0, m_bodies[1], Matrix4d::Identity(), RESOURCE_DIR);
		auto wrapSphere = addWrapSphere(m_bodies[0], Vector3d(1.0, 0.0, 0.0), m_bodies[2], Vector3d(1.0, 0.0, 0.0), compSphere, 20, RESOURCE_DIR);

	}
	break;

	case WRAP_CYLINDER:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				//joint->m_qdot(0) = -5.0;
			}
		}

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(0.0, 1.0, 0.0));
		auto compCylinder = addCompCylinder(0.5, m_bodies[1], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0), RESOURCE_DIR, "obstacle.obj");
		auto wrapCylinder = addWrapCylinder(m_bodies[0], Vector3d(1.0, 0.0, 0.0), m_bodies[2], Vector3d(1.0, 0.0, 0.0), compCylinder, 20, RESOURCE_DIR);

	}
	break;

	case WRAP_DOUBLECYLINDER:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				joint->m_qdot(0) = -0.0;
			}
		}

		Matrix4d E = SE3::RpToE(SE3::aaToMat(Vector3d(1.0, 0.0, 0.0), 0.0), Vector3d(3.0, 1.0, 0.0));
		auto compDoubleCylinder = addCompDoubleCylinder(
			0.5, m_bodies[0], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0),
			0.5, m_bodies[2], E, Vector3d(0.0, 0.0, 1.0), Vector3d(0.0, 0.0, 0.0), 
			RESOURCE_DIR, "obstacle.obj", "obstacle.obj");

		auto wrapDoubleCylinder = addWrapDoubleCylinder(
			m_bodies[0], Vector3d(-5.0, 0.5, 0.0), 
			m_bodies[2], Vector3d(5.0, 0.5, 0.0), 
			Vector3d(0.0, 0.0, 0.0), Vector3d(0.0, 0.0, 0.0), 
			Vector3d(0.0, 0.0, -1.0), Vector3d(0.0, 0.0, 1.0), 
			compDoubleCylinder, 20, RESOURCE_DIR);

	}
	break;

	case SPRING_DAMPER:
	{
		m_h = 1.0e-1;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			// Inits joints
			if (i == 0) {
				addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		m_joints[1]->m_q(0) = - M_PI / 2.0;
		m_joints[2]->m_q(0) = M_PI / 2.0;

		auto spring = make_shared<SpringDamper>(m_bodies[0], Vector3d(-2.0, -0.5, 0.0), m_bodies[1], Vector3d(2.0, -0.5, 0.0));
		m_springs.push_back(spring);
		m_nsprings++;
		spring->setStiffness(1.0e5);
		spring->setDamping(1.0e3);
		spring->load(RESOURCE_DIR);

	}
	break;

	case FREEJOINT: 
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -1.0, 0.0;
		Eigen::from_json(js["cube_sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 1; i++) {

			auto body = addBody(density, sides, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box_1_1_1.obj");
			// Inits joints
			if (i == 0) {
				Vector6d qdot0;
				qdot0 << 0.2, 0.4, 0.6, 0.0, 0.0, 3.0;
				addJointFree(body, Vector3d::Zero(), Matrix3d::Identity(), Vector6d::Zero(),qdot0, RESOURCE_DIR);
				
			}
			else {
				
			}
		}
		break;
	}

	case TEST_REDUCED_HYBRID_DYNAMICS:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 3; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}

		auto con0 = make_shared<ConstraintPrescJoint>(m_joints[0], REDMAX_EULER);
		m_nconstraints++;
		m_constraints.push_back(con0);


		break;
	}
	case TEST_MAXIMAL_HYBRID_DYNAMICS:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -980, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 4; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -ToRadian(90)), 0.0, RESOURCE_DIR);
			}
			else if (i == 1) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), ToRadian(90)), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}
			else if (i == 2) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), ToRadian(90)), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				
			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -ToRadian(90)), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}
		Vector3i dof;
		dof << 2, 3, 4;
		auto con0 = make_shared<ConstraintPrescBody>(m_bodies[3], dof, REDMAX_EULER);
		m_nconstraints++;
		m_constraints.push_back(con0);
		break;
	}

	case TEST_HYPER_REDUCED_COORDS:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -9.8, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 2; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");

			//// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				//addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
				addJointRevoluteHyperReduced(body, Vector3d::UnitZ(), m_joints[i-1], 0.8, Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}
		VectorXi dof(1);
		dof << 2;
		addConstraintPrescJoint(m_joints[0]);		
		addConstraintPrescBody(m_bodies[1], dof);
	}
	break;
	case TEST_CONSTRAINT_PRESC_BODY_ATTACH_POINT:
	{
		m_h = 1.0e-2;
		density = 1.0;
		m_grav << 0.0, -98, 0.0;
		Eigen::from_json(js["sides"], sides);
		//m_nbodies = 5;
		//m_njoints = 5;
		m_Hexpected = 10000; // todo
		m_tspan << 0.0, 5.0;
		m_t = 0.0;
		// Inits rigid bodies
		for (int i = 0; i < 4; i++) {

			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "box10_1_1.obj");
			// Inits joints
			if (i == 0) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(0.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -M_PI / 2.0), 0.0, RESOURCE_DIR);
			}
			else if (i == 1) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), M_PI / 2.0), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}
			else if (i == 2) {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), M_PI / 2.0), 0.0, RESOURCE_DIR, m_joints[i - 1]);

			}
			else {
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), SE3::aaToMat(Vector3d(0.0, 0.0, 1.0), -M_PI / 2.0), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}
		VectorXi dof(2);
		dof << 0, 1;
		auto con0 = make_shared<ConstraintPrescBodyAttachPoint>(m_bodies[3], Vector3d(5.0, 0.0, 0.0), dof, REDMAX_EULER);
		m_nconstraints++;
		m_constraints.push_back(con0);
		break;
	}

	case MUSCLE_INERTIA:
	{
		m_h = 1.0e-2;
		m_tspan << 0.0, 50.0;
		m_t = 0.0;
		density = 1.0;
		m_grav << 0.0, -9.81, 0.0;
		Eigen::from_json(js["sides"], sides);

		for (int i = 0; i < 2; i++) {
			auto body = addBody(density, sides, Vector3d(5.0, 0.0, 0.0), Matrix3d::Identity(), RESOURCE_DIR, "cylinder_9.obj");

			// Inits joints
			if (i == 0) {
				//addJointFixed(body, Vector3d(0.0, 0.0, 0.0), Matrix3d::Identity(), 0.0);
				addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR);
			}
			else {
				auto joint = addJointRevolute(body, Vector3d::UnitZ(), Vector3d(10.0, 0.0, 0.0), Matrix3d::Identity(), 0.0, RESOURCE_DIR, m_joints[i - 1]);
			}
		}
		//m_joints[0]->m_q(0) = -M_PI / 2.0;
		//m_joints[1]->m_q(0) = -M_PI / 2.0;
		//addDeformableSpring(10.0, 0.0, 2, m_bodies[0], Vector3d(-5.0, 0.0, 0.0), m_bodies[1], Vector3d(0.0, 0.0, 0.0), RESOURCE_DIR);
		/*auto spring = make_shared<SpringDamper>(m_bodies[0], Vector3d(-5.0, 0.0, 0.0), m_bodies[1], Vector3d(0.0, 0.0, 0.0));
		m_springs.push_back(spring);
		m_nsprings++;
		spring->setStiffness(0.0);
		spring->setDamping(0.0);
		spring->load(RESOURCE_DIR);*/
		std::vector<std::shared_ptr<Body>> related_bodies;
		related_bodies.push_back(m_bodies[0]);
		related_bodies.push_back(m_bodies[1]);

		auto muscle = make_shared<MuscleSpring>(related_bodies, 3);
		muscle->setMass(10.0);
		m_muscles.push_back(muscle);
		m_nmuscles++;
		muscle->setAttachments(m_bodies[0], Vector3d(-5.0, 0.0, 0.0), m_bodies[1], Vector3d(0.0, 0.0, 0.0));
		muscle->load(RESOURCE_DIR);
	}
	break;

default:
		break;
	}
}

shared_ptr<ConstraintPrescJoint> World::addConstraintPrescJoint(shared_ptr<Joint> j) {
	auto con = make_shared<ConstraintPrescJoint>(j, REDMAX_EULER);
	m_nconstraints++;
	m_constraints.push_back(con);
	return con;
}

shared_ptr<ConstraintPrescBody> World::addConstraintPrescBody(shared_ptr<Body> b, VectorXi dof) {
	auto con = make_shared<ConstraintPrescBody>(b, dof, REDMAX_EULER);
	m_nconstraints++;
	m_constraints.push_back(con);
	return con;
}

shared_ptr<ConstraintPrescBodyAttachPoint> World::addConstraintPrescBodyAttachPoint(shared_ptr<Body> b, Vector3d r, VectorXi dof) {
	auto con = make_shared<ConstraintPrescBodyAttachPoint>(b, r, dof, REDMAX_EULER);
	m_nconstraints++;
	m_constraints.push_back(con);
	return con;
}

shared_ptr<Body> World::addBody(double density, Vector3d sides, Vector3d p, Matrix3d R, const string &RESOURCE_DIR, string file_name) {
	auto body = make_shared<BodyCuboid>(density, sides);
	Matrix4d E = SE3::RpToE(R, p);
	body->setTransform(E);
	body->load(RESOURCE_DIR, file_name);
	m_bodies.push_back(body);
	m_nbodies++;
	return body;
}

shared_ptr<JointRevolute> World::addJointRevolute(shared_ptr<Body> body, 
	Vector3d axis, 
	Vector3d p, 
	Matrix3d R, 
	double q, 
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointRevolute>(body, axis, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q(0) = q;
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointUniversalXY> World::addJointUniversalXY(shared_ptr<Body> body,
	Vector3d p,
	Matrix3d R,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointUniversalXY>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointRevoluteHyperReduced> World::addJointRevoluteHyperReduced(shared_ptr<Body> body,
	Vector3d axis,
	shared_ptr<Joint> friend_joint,
	double scalar,
	Vector3d p,
	Matrix3d R,
	double q,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) {
	auto joint = make_shared<JointRevoluteHyperReduced>(body, axis, friend_joint, scalar, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q(0) = q;
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointFree> World::addJointFree(
	shared_ptr<Body> body,
	Vector3d p,
	Matrix3d R,
	Vector6d q0,
	Vector6d qdot0,
	const string &RESOURCE_DIR,
	shared_ptr<Joint> parent) 
{
	auto joint = make_shared<JointFree>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	joint->m_q = q0;
	joint->m_qdot = qdot0;
	joint->load(RESOURCE_DIR, "sphere2.obj");
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<JointFixed> World::addJointFixed(shared_ptr<Body> body, Vector3d p, Matrix3d R, double q, shared_ptr<Joint> parent) {
	auto joint = make_shared<JointFixed>(body, parent);
	Matrix4d E = SE3::RpToE(R, p);
	joint->setJointTransform(E);
	
	m_joints.push_back(joint);
	m_njoints++;
	return joint;
}

shared_ptr<ConstraintJointLimit> World::addConstraintJointLimit(shared_ptr<Joint> joint, double ql, double qu) {
	auto constraint = make_shared<ConstraintJointLimit>(joint);
	m_constraints.push_back(constraint);
	constraint->setLimits(ql, qu);
	m_nconstraints++;
	return constraint;
}

shared_ptr<DeformableSpring> World::addDeformableSpring(double mass, double k, int n_points, shared_ptr<Body> body0, Vector3d r0, shared_ptr<Body> body1, Vector3d r1, const string &RESOURCE_DIR) {

	auto deformable = make_shared<DeformableSpring>(n_points, m_countS, m_countCM);
	m_deformables.push_back(deformable);
	deformable->setStiffness(k);
	deformable->setMass(mass);
	deformable->setAttachments(body0, r0, body1, r1);
	deformable->load(RESOURCE_DIR);
	m_ndeformables++;
	return deformable;
}

shared_ptr<CompSphere> World::addCompSphere(double r, shared_ptr<Body> parent, Matrix4d E, const string &RESOURCE_DIR) {
	auto comp = make_shared<CompSphere>(parent, r);
	m_comps.push_back(comp);
	comp->setTransform(E);
	comp->load(RESOURCE_DIR);
	m_ncomps++;
	return comp;
}

shared_ptr<CompCylinder> World::addCompCylinder(double r, shared_ptr<Body> parent, Matrix4d E, Vector3d z, Vector3d o, const string &RESOURCE_DIR, string shape) {
	auto comp = make_shared<CompCylinder>(parent, r);
	m_comps.push_back(comp);
	comp->setTransform(E);
	auto z_axis = make_shared<Vector>();
	z_axis->dir0 = z;
	comp->setZAxis(z_axis);
	auto origin = make_shared<Node>();
	origin->x0 = o;
	comp->setOrigin(origin);
	comp->load(RESOURCE_DIR, shape);

	m_ncomps++;
	return comp;
}

shared_ptr<CompDoubleCylinder> World::addCompDoubleCylinder(
	double rA, shared_ptr<Body> parentA, Matrix4d EA, Vector3d z_a, Vector3d o_a,
	double rB, shared_ptr<Body> parentB, Matrix4d EB, Vector3d z_b, Vector3d o_b,
	const string &RESOURCE_DIR, string shapeA, string shapeB) {
	auto comp = make_shared<CompDoubleCylinder>(parentA, rA, parentB, rB);
	m_comps.push_back(comp);
	comp->setTransformA(EA);
	comp->setTransformB(EB);
	auto  za = make_shared<Vector>();
	za->dir0 = z_a;
	comp->setZAxisA(za);
	auto zb = make_shared<Vector>();
	zb->dir0 = z_b;
	comp->setZAxisB(zb);
	auto originA = make_shared<Node>();
	originA->x0 = o_a;
	auto originB = make_shared<Node>();
	originB->x0 = o_b;
	comp->setOriginA(originA);
	comp->setOriginB(originB);
	comp->load(RESOURCE_DIR, shapeA, shapeB);
	m_ncomps++;
	return comp;
}

shared_ptr<ConstraintNull> World::addConstraintNull() {
	auto constraint = make_shared<ConstraintNull>();
	m_nconstraints++;
	m_constraints.push_back(constraint);
	return constraint;
}

shared_ptr<JointNull> World::addJointNull() {
	auto joint = make_shared<JointNull>();
	m_njoints++;
	m_joints.push_back(joint);
	return joint;
}

shared_ptr<DeformableNull> World::addDeformableNull() {

	auto deformable = make_shared<DeformableNull>();
	m_ndeformables++;
	m_deformables.push_back(deformable);
	return deformable;
}

shared_ptr<CompNull> World::addCompNull() {
	auto comp = make_shared<CompNull>();
	m_ncomps++;
	m_comps.push_back(comp);
	return comp;
}

shared_ptr<WrapObst> World::addWrapNull() {
	auto wrap = make_shared<WrapNull>();
	m_nwraps++;
	m_wraps.push_back(wrap);
	return wrap;

}

shared_ptr<SpringNull> World::addSpringNull() {
	auto spring = make_shared<SpringNull>();
	m_nsprings++;
	m_springs.push_back(spring);
	return spring;
}

std::shared_ptr<MuscleNull> World::addMuscleNull()
{
	auto muscle = make_shared<MuscleNull>();
	m_nmuscles++;
	m_muscles.push_back(muscle);
	return muscle;
}

shared_ptr<WrapSphere> World::addWrapSphere(shared_ptr<Body> b0, Vector3d r0, shared_ptr<Body> b1, Vector3d r1, shared_ptr<CompSphere> compSphere, int num_points, const string &RESOURCE_DIR) {
	auto P = make_shared<Node>();
	P->x0 = r0;
	P->setParent(b0);
	auto S = make_shared<Node>();
	S->x0 = r1;
	S->setParent(b1);

	auto wrapSphere = make_shared<WrapSphere>(P, S, compSphere, num_points);
	m_nwraps++;
	wrapSphere->load(RESOURCE_DIR);
	m_wraps.push_back(wrapSphere);
	return wrapSphere;
}

shared_ptr<WrapCylinder> World::addWrapCylinder(shared_ptr<Body> b0, Vector3d r0, shared_ptr<Body> b1, Vector3d r1, shared_ptr<CompCylinder> compCylinder, int num_points, const string &RESOURCE_DIR) {
	auto P = make_shared<Node>();
	P->x0 = r0;
	P->setParent(b0);

	auto S = make_shared<Node>();
	S->x0 = r1;
	S->setParent(b1);

	auto wrapCylinder = make_shared<WrapCylinder>(P, S, compCylinder, num_points);
	m_nwraps++;
	wrapCylinder->load(RESOURCE_DIR);
	m_wraps.push_back(wrapCylinder);

	return wrapCylinder;
}

shared_ptr<WrapDoubleCylinder> World::addWrapDoubleCylinder(shared_ptr<Body> b0, Vector3d r0, shared_ptr<Body> b1, Vector3d r1, Vector3d u, Vector3d v, Vector3d z_u, Vector3d z_v, shared_ptr<CompDoubleCylinder> compDoubleCylinder, int num_points, const string &RESOURCE_DIR) {
	auto z_axis_U = make_shared<Vector>();
	z_axis_U->dir0 = z_u;
	auto z_axis_V = make_shared<Vector>();
	z_axis_V->dir0 = z_v;
	compDoubleCylinder->setZAxisA(z_axis_U);
	compDoubleCylinder->setZAxisB(z_axis_V);

	auto origin_U = make_shared<Node>();
	origin_U->x0 = u;
	compDoubleCylinder->setOriginA(origin_U);
	auto origin_V = make_shared<Node>();
	origin_V->x0 = v;
	compDoubleCylinder->setOriginB(origin_V);

	auto P = make_shared<Node>();
	P->x0 = r0;
	P->setParent(b0);
	auto S = make_shared<Node>();
	S->x0 = r1;
	S->setParent(b1);
	auto wrapDoubleCylinder = make_shared<WrapDoubleCylinder>(P, S, compDoubleCylinder, num_points);
	wrapDoubleCylinder->load(RESOURCE_DIR);
	m_wraps.push_back(wrapDoubleCylinder);
	m_nwraps++;
	return wrapDoubleCylinder;
}


void World::init() {
	for (int i = 0; i < m_nbodies; i++) {
		m_bodies[i]->init(nm);
		if (i < m_nbodies - 1) {
			m_bodies[i]->next = m_bodies[i + 1];
		}
	}
	//cout << m_bodies[m_nbodies - 1]->idxM << endl;

	nm = 0;
	/*for (int i = 0; i < m_njoints; i++) {
	m_joints[i]->init(nm, nr);
	}*/

	//joint ordering
	// todo
	for (int i = m_njoints - 1; i > -1; i--) {
		m_joints[i]->countHRDofs(nR);
	}

	for (int i = m_njoints - 1; i > -1; i--) {		
		m_joints[i]->init(nm, nr);
	}

	m_dense_nm = nm;
	m_dense_nr = nr;
	// Until now, save nm, nr for later Sparse Jacobian Computation, 
	// only dense_nm, dense_nr 


	for (int i = 0; i < m_njoints; i++) {
		if (i < m_njoints - 1) {
			m_joints[i]->next = m_joints[i + 1];
		}
		if (i > 0) {
			m_joints[i]->prev = m_joints[i - 1];
		}
	}

	for (int i = 0; i < m_ncomps; ++i) {
		m_comps[i]->init();
		if (i < m_ncomps - 1) {
			m_comps[i]->next = m_comps[i + 1];
		}
	}

	for (int i = 0; i < m_nsprings; ++i) {
		m_springs[i]->init();
		if (i < m_nsprings - 1) {
			m_springs[i]->next = m_springs[i + 1];
		}
	}

	for (int i = 0; i < m_nwraps; ++i) {
		m_wraps[i]->init();
		if (i < m_nwraps - 1) {
			m_wraps[i]->next = m_wraps[i + 1];
		}
	}

	if (m_njoints == 0) {
		addJointNull();
	}

	if (m_ncomps == 0) {
		addCompNull();
	}

	if (m_nwraps == 0) {
		addWrapNull();
	}

	if (m_nsprings == 0) {
		addSpringNull();
	}

	m_joints[0]->update();
	m_comps[0]->update();
	m_wraps[0]->update();
	m_springs[0]->update();
	
	for (int i = 0; i < m_ndeformables; i++) {
		m_deformables[i]->countDofs(nm, nr);

		m_deformables[i]->init();
		// Create attachment constraints
		auto constraint = make_shared<ConstraintAttachSpring>(m_deformables[i]);
		m_constraints.push_back(constraint);
		m_nconstraints++;
		if (i < m_ndeformables - 1) {
			m_deformables[i]->next = m_deformables[i + 1];
		}
	}

	
	if (m_ndeformables == 0) {
		addDeformableNull();
	}

	for (int i = 0; i < m_nmuscles; i++) {

		m_muscles[i]->init();
		if (i < m_nmuscles - 1) {
			m_muscles[i]->next = m_muscles[i + 1];
		}
	}

	if (m_nmuscles == 0) {
		addMuscleNull();
	}

	int tet = nm;


	tet = nm - tet;
	tet = tet / 3;
	m_ntets = tet;


	// init constraints

	for (int i = 0; i < m_nconstraints; i++) {
		m_constraints[i]->countDofs(nem, ner, nim, nir);
		m_constraints[i]->init();
		if (i < m_nconstraints - 1) {
			m_constraints[i]->next = m_constraints[i + 1];
		}
	}

	if (m_nconstraints == 0) {
		addConstraintNull();
	}
}

void World::update() {
	m_comps[0]->update();
	m_wraps[0]->update();
	m_springs[0]->update();
	m_muscles[0]->update();
}

int World::getNsteps() {
	// Computes the number of results
	int nsteps = int((m_tspan(1) - m_tspan(0)) / m_h);
	return nsteps;
}

void World::drawFloor(Floor f, shared_ptr<MatrixStack> MV, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) {
	// Draw grid
	progSimple->bind();
	glUniformMatrix4fv(progSimple->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(progSimple->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	glLineWidth(2.0f);
	float x0 = f.xrange(0);
	float x1 = f.xrange(1);
	float z0 = f.zrange(0);
	float z1 = f.zrange(1);
	int gridSize = 10;
	glLineWidth(1.0f);
	glBegin(GL_LINES);
	for (int i = 1; i < gridSize; ++i) {
		if (i == gridSize / 2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		}
		else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float x = x0 + i / (float)gridSize * (x1 - x0);
		glVertex3f(x, f.y, z0);
		glVertex3f(x, f.y, z1);
	}
	for (int i = 1; i < gridSize; ++i) {
		if (i == gridSize / 2) {
			glColor3f(0.1f, 0.1f, 0.1f);
		}
		else {
			glColor3f(0.8f, 0.8f, 0.8f);
		}
		float z = z0 + i / (float)gridSize * (z1 - z0);
		glVertex3f(x0, f.y, z);
		glVertex3f(x1, f.y, z);
	}
	glEnd();
	glColor3f(0.4f, 0.4f, 0.4f);
	glBegin(GL_LINE_LOOP);
	glVertex3f(x0, f.y, z0);
	glVertex3f(x1, f.y, z0);
	glVertex3f(x1, f.y, z1);
	glVertex3f(x0, f.y, z1);
	glEnd();
	progSimple->unbind();
}

void World::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, const shared_ptr<Program> progSoft, shared_ptr<MatrixStack> P) {
	m_bodies[0]->draw(MV, prog, P);
	m_joints[0]->draw(MV, prog, progSimple, P);
	m_deformables[0]->draw(MV, prog, progSimple, P);
	m_comps[0]->draw(MV, prog, P);
	m_wraps[0]->draw(MV, prog, progSimple, P);
	m_springs[0]->draw(MV, prog, progSimple, P);
	m_muscles[0]->draw(MV, prog, progSimple, P);
	for (int i = 0; i < (int)m_floors.size(); ++i) {
		drawFloor(m_floors[i], MV, progSimple, P);
	}
}

Energy World::computeEnergy() {
	m_energy.K = 0.0;
	m_energy.V = 0.0;

	m_joints[0]->computeEnergies(m_grav, m_energy);
	m_deformables[0]->computeEnergies(m_grav, m_energy);
	m_springs[0]->computeEnergies(m_grav, m_energy);
	m_muscles[0]->computeEnergies(m_grav, m_energy);

	if (m_t == 0.0) {
		m_energy0 = m_energy;
	}

	m_energy.V -= m_energy0.V;

	return m_energy;
}

shared_ptr<Body> World::getBody(int uid) {
	MapBodyUID::const_iterator it = m_bodyUID.find(uid);
	return (it == m_bodyUID.end() ? NULL : it->second);
}

shared_ptr<Body> World::getBody(const string &name) {
	MapBodyName::const_iterator it = m_bodyName.find(name);
	return (it == m_bodyName.end() ? NULL : it->second);
}

shared_ptr<Constraint> World::getConstraint(int uid) {
	MapConstraintUID::const_iterator it = m_constraintUID.find(uid);
	return (it == m_constraintUID.end() ? NULL : it->second);
}

shared_ptr<Constraint> World::getConstraint(const string &name) {
	MapConstraintName::const_iterator it = m_constraintName.find(name);
	return (it == m_constraintName.end() ? NULL : it->second);
}

void World::sceneTestReducedHD(double t) {
	m_joints[0]->presc->m_q[0] = 0.0;
	m_joints[0]->presc->m_qdot[0] = 0.0;
	m_joints[0]->presc->m_qddot[0] = 0.0;

}

void World::sceneTestMaximalHD(double t) {
	Matrix4d E = m_bodies[3]->E_wi;
	Matrix3d R = E.topLeftCorner(3, 3);
	Vector6d phi = m_bodies[3]->phi;
	Vector3d vt_w, wt_i, vtdot_w, wtdot_i;

	if (t < 2.0) {
		vt_w.setZero();
		wt_i << 0.0, 0.0, t;
		vtdot_w.setZero();
		wtdot_i << 0.0, 0.0, 1.0;

	}
	else if (t < 4.0) {
		double t_ = t - 4.0;
		vt_w.setZero();
		wt_i << 0.0, 0.0, -t_;
		vtdot_w.setZero();
		wtdot_i << 0.0, 0.0, -1.0;

	}
	else if (t < 6.0) {
		double t_ = t - 4.0;
		vt_w << -2 * t_, 0.0, 0.0;
		wt_i << 0.0, 0.0, -t_;
		vtdot_w << -2.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, -1.0;
	}
	else if (t < 8.0) {
		double t_ = t - 8.0;
		vt_w << 2 * t_, 0.0, 0.0;
		wt_i << 0.0, 0.0, t_;
		vtdot_w << 2.0, 0.0, 0.0;
		wtdot_i << 0.0, 0.0, 1.0;
	}
	else {
		vt_w.setZero();
		wt_i.setZero();
		vtdot_w.setZero();
		wtdot_i.setZero();
	}

	Vector3d vt_i = R.transpose() * vt_w;
	m_bodies[3]->presc->m_qdot.segment<3>(0) = wt_i;
	m_bodies[3]->presc->m_qdot.segment<3>(3) = vt_i;
	m_bodies[3]->presc->m_qddot.segment<3>(0) = wtdot_i;
	m_bodies[3]->presc->m_qddot.segment<3>(3) = R.transpose() * vtdot_w - phi.segment<3>(0).cross(vt_i);
}


void World::setMaximalPrescStates(shared_ptr<Body> b, Vector3d vt_w, Vector3d vtdot_w, Vector3d wt_i, Vector3d wtdot_i) {
	b->presc->setActive();
	Matrix4d E = b->E_wi;
	Matrix3d R = E.topLeftCorner(3, 3);
	Vector6d phi = b->phi;

	Vector3d vt_i = R.transpose() * vt_w;
	b->presc->m_qdot.segment<3>(0) = wt_i;
	b->presc->m_qdot.segment<3>(3) = vt_i;
	b->presc->m_qddot.segment<3>(0) = wtdot_i;
	b->presc->m_qddot.segment<3>(3) = R.transpose() * vtdot_w - phi.segment<3>(0).cross(vt_i);
}

void World::setMaximalPrescStates(int index_body, Vector3d vt_w, Vector3d wt_i) {
	auto b = m_bodies[index_body];
	b->presc->setActive();
	Matrix4d E = b->E_wi;
	Matrix3d R = E.topLeftCorner(3, 3);
	Vector6d phi = b->phi;

	Vector3d vt_i = R.transpose() * vt_w;
	b->presc->m_qdot.segment<3>(0) = wt_i;
	b->presc->m_qdot.segment<3>(3) = vt_i;
}

void World::setMaximalPrescAttachPointStates(int index_body, int index_point, Vector3d vt_w) {
	auto b = m_bodies[index_body];
	auto con = b->m_presc_attach_points[index_point];
	con->setActive();

	con->m_qdot = vt_w;

}


void World::setReducedPrescStates(shared_ptr<Joint> j, VectorXd q, VectorXd dq) {
	j->presc->setActive();
	j->presc->m_q = q;
	j->presc->m_qdot = dq;

}

void World::setReducedPrescStates(shared_ptr<Joint> j, double q, double dq) {
	j->presc->setActive();
	j->presc->m_q[0] = q;
	j->presc->m_qdot[0] = dq;
}

void World::computeTargetQ(double t0, double t1, double t, double angle, double q0, double &q, double &dq) {

	double a = 7.0;
	double b = angle;
	double s = 2 * ((t - t0) / (t1 - t0) - 0.5);
	q = q0 + b / (1 + exp(-a * s));
	double T = t - t0;
	double TT = t0 - t1;
	double w = 2 * T;
	double P = w / TT + 1;
	double e = exp(a * P);
	double Q = T / TT + 1;
	double f = exp(a * Q);
	dq = -(2 * a*b*e) / (TT * (f + 1) *(f + 1));
}

void World::deactivateListOfPrescConstraints(Eigen::VectorXi mcon, Eigen::VectorXi rcon) {

	for (int i = 0; i < (int)mcon.size(); ++i) {
		if (mcon(i) > -1) {
			m_bodies[mcon(i)]->presc->setInactive();
		}	
	}

	for (int i = 0; i < (int)rcon.size(); ++i) {
		if (rcon(i) > -1) {
			m_joints[rcon(i)]->presc->setInactive();
		}	
	}
}

void World::setListOfReducedPrescStates(Eigen::VectorXi rcon, VectorXd q, VectorXd dq) {
	for (int i = 0; i < (int)rcon.size(); ++i) {
		setReducedPrescStates(m_joints[rcon(i)], q, dq);
	}
}

void World::setListOfMaximalPrescStates(Eigen::VectorXi mcon, Vector3d vt_w, Vector3d vtdot_w, Vector3d wt_i, Vector3d wtdot_i) {
	for (int i = 0; i < (int)mcon.size(); ++i) {
		setMaximalPrescStates(m_bodies[mcon(i)], vt_w, vtdot_w, wt_i, wtdot_i);
	}
}


void World::sceneTestHyperReduced(double t) {
	if (t < 2.0) {
		
		//computeTargetQ(0.0, 2.0, t, M_PI / 8.0, 0.0, q, dq);
		//setReducedPrescStates(m_joints[0], q, dq);
		Vector3d vt_w, wt_i, vtdot_w, wtdot_i;
		vt_w.setZero();
		setMaximalPrescStates(m_bodies[1], vt_w, vt_w, vt_w, vt_w);
	}
	if (t < 4.0) {

		//computeTargetQ(2.0, 4.0, t, -M_PI / 8.0, M_PI / 8.0, q, dq);
		//setReducedPrescStates(m_joints[0], q, dq);

	}
}

void World::sceneAttachPoint(double t) {

	Vector3d vt_w, wt_i;
	vt_w.setZero();
	if (t < 105.0) {
		vt_w << 0, 0.0, 0.0;
		setMaximalPrescAttachPointStates(3, 0, vt_w);
	}
}
