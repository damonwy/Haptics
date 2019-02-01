#pragma once

#include "Muscle.h"

class Body;
class Node;

class MuscleSpring : public Muscle
{
public:
	MuscleSpring() {}
	MuscleSpring(std::vector<std::shared_ptr<Body>>, int n_nodes);
	virtual ~MuscleSpring() {}

	void setStiffness(double K) { m_K = K; }
	void setMass(double mass) { m_mass = mass; }
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
	void computeJacobian_(Eigen::MatrixXd &J);
	void computeJacobianSparse_(std::vector<T> &J_);

	std::shared_ptr<Body> m_body0;
	std::shared_ptr<Body> m_body1;
};