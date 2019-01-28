#pragma once
#include "Comp.h"

class Node;

class CompSphere : public Comp
{
public:
	CompSphere();
	CompSphere(std::shared_ptr<Body> parent, double r);
	virtual ~CompSphere() {}

	void load(const std::string &RESOURCE_DIR);
	void init();
	void setTransform(Eigen::Matrix4d E);
	double getRadius() { return m_r; }
	std::shared_ptr<Node> getOrigin() { return m_O; }

protected:
	double m_r;
	std::shared_ptr<Node> m_O;
	Eigen::Matrix4d E_wi;	// Where the component is wrt world
	Eigen::Matrix4d E_ji;	// Where the component is wrt body

	std::shared_ptr<Shape> m_shape;
	std::shared_ptr<Body> m_parent;	
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;
	void update_();

};