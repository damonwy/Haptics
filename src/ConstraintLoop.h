#pragma once
// ConstraintLoop Spherical constraint between two bodies

#include "Constraint.h"

class Body;

class ConstraintLoop : public Constraint
{
public:
	ConstraintLoop();
	ConstraintLoop(std::shared_ptr<Body> bodyA, std::shared_ptr<Body> bodyB);
	void setPositions(Eigen::Vector3d xA, Eigen::Vector3d xB) { m_xA = xA; m_xB = xB; }
	
	Eigen::Vector3d m_xA;			// Local position wrt A
	Eigen::Vector3d m_xB;			// Local position wrt B

	std::shared_ptr<Body> m_bodyA;
	std::shared_ptr<Body> m_bodyB;

protected:
	void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);


};