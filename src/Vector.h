#pragma once

#ifndef MUSCLEMASS_SRC_VECTOR_H_
#define MUSCLEMASS_SRC_VECTOR_H_

#include <vector>
#include <memory>
#include "Node.h"

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

class Node;
class Program;
class MatrixStack;
class Body;

class Vector
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

	Vector();
	Vector(std::shared_ptr<Node> p, std::shared_ptr<Body> body, Eigen::Vector3d dir);

	virtual ~Vector() {}
	void reset();
	void setP(std::shared_ptr<Node> p) { this->m_p = p; }
	void draw(std::shared_ptr<MatrixStack> MV, std::shared_ptr<MatrixStack> P, const std::shared_ptr<Program> p) const;
	void update(Eigen::Matrix4d E);
	void update();

	Eigen::Vector3d dir0; // initial direction
	Eigen::Vector3d dir;  // current direction

	bool fixed;

private:
	std::shared_ptr<Node> m_p;	// starting point
	std::shared_ptr<Body> m_body; // attached body
};

#endif // MUSCLEMASS_SRC_VECTOR_H_