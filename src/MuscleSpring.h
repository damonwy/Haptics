#pragma once

#include "Muscle.h"

class Body;
class Node;
class World;

class MuscleSpring : public Muscle
{
public:
	MuscleSpring() {}
	MuscleSpring(std::vector<std::shared_ptr<Body>>, int n_nodes);
	virtual ~MuscleSpring() {}

	void setStiffness(double K) { m_K = K; }
	void setMass(double mass);
	void setAttachments(std::shared_ptr<Body> body0, Vector3d r0, std::shared_ptr<Body> body1, Vector3d r1);

protected:
	void init_();
	void update_();
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;

	void computeMass_(Eigen::MatrixXd &M);
	void computeMassSparse_(std::vector<T> &M_);
	void computeForce_(Vector3d grav, Eigen::VectorXd &f);
	void computeForceDamping_(Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd &D);
	void computeForceDampingSparse_(Vector3d grav, Eigen::VectorXd &f, std::vector<T> &D_);

	void computeEnergies_(Vector3d grav, Energy &ener);
	void computeJMJ_(Eigen::MatrixXd &J, std::shared_ptr<World> world);
	void computeJMJSparse_(std::vector<T> &J_);
	void computeJMJdotqdot_(Eigen::VectorXd & f, const Eigen::VectorXd &qdot, std::shared_ptr<World> world);

	std::shared_ptr<Body> m_body0;
	std::shared_ptr<Body> m_body1;
};