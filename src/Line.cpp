#include "rmpch.h"
#include "Line.h"
#include "Node.h"

using namespace std;
using namespace Eigen;

Line::Line(std::shared_ptr<Node> n0, std::shared_ptr<Node> n1):
m_n0(n0), m_n1(n1)
{
	m_x0 = n0->x;
	m_x1 = n1->x;
}

Line::Line(Vector3d x0, Vector3d x1):
m_x0(x0), m_x1(x1)
{

}

void Line::draw() {
}

void Line::addSampleNodes(int n, std::vector<std::shared_ptr<Node>> &nodes) {
	int idx = (int)nodes.size();
	for (int s = 0; s < n; s++) {
		double f = double(s + 1) / double(n + 1);
		auto node = make_shared<Node>();
		node->i = idx + s;
		node->x = f * m_x1 + (1.0 - f) * m_x0;
		nodes.push_back(node);
	}
}


bool Line::isInLine(Vector3d x) {
	double diff = (x - m_x0).norm() + (x - m_x1).norm() - (m_x0 - m_x1).norm();
	if (diff < 0.0001) {
		return true;
	}
	else {
		return false;
	}
}

bool Line::isInLine(std::shared_ptr<Node> n) {
	double diff = (n->x - m_x0).norm() + (n->x - m_x1).norm() - (m_x0 - m_x1).norm();

	if (diff < 0.0001) {
		return true;
	}
	else {
		return false;
	}
}