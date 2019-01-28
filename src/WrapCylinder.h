#pragma once
#ifndef MUSCLEMASS_SRC_WRAPCYLINDER_H_
#define MUSCLEMASS_SRC_WRAPCYLINDER_H_

/*
* WrapCylinder.hpp
*
* Obstacle Set Algorithm Simulation for Cylinder Obstacles
*
*/

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Vector;
class CompCylinder;
class Node;

#include "WrapObst.h"

class WrapCylinder : public WrapObst
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
	std::shared_ptr<Vector> m_vec_z;      // Cylinder Positive z axis
	std::shared_ptr<CompCylinder> m_compCylinder;

public:
	WrapCylinder();
	WrapCylinder(const std::shared_ptr<Node> &P, const std::shared_ptr<Node> &S, const std::shared_ptr<CompCylinder> compCylinder, const int num_points);

	void compute();	
	Eigen::MatrixXd getPoints(int num_points) const;
	void init();
	void load(const std::string &RESOURCE_DIR);
	
protected:
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	void update_();
};

#endif // MUSCLEMASS_SRC_WRAPCYLINDER_H_