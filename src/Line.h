#pragma once

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include "MLCommon.h"
class Node;
class Body;

class Line {

public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	Line() {}
	Line(std::shared_ptr<Node> n0, std::shared_ptr<Node> n1);

	Line(Vector3d x0, Vector3d x1);
	void draw();
	void setBody(std::shared_ptr<Body> b) { m_body = b; }
	std::shared_ptr<Body> getBody() { return m_body; };
	bool isInLine(Vector3d x);
	bool isInLine(std::shared_ptr<Node> n);
	void addSampleNodes(int n, std::vector<std::shared_ptr<Node>> &nodes);

protected:
	std::shared_ptr<Node> m_n0;
	std::shared_ptr<Node> m_n1;
	Vector3d m_x0;
	Vector3d m_x1;
	std::shared_ptr<Body> m_body;
};