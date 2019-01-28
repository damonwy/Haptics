#pragma once
#include "Comp.h"

class Vector;
class Node;

class CompCylinder : public Comp
{
public:
	CompCylinder();
	CompCylinder(std::shared_ptr<Body> parent, double r);
	virtual ~CompCylinder() {}

	void load(const std::string &RESOURCE_DIR, std::string shape);
	void init();
	void setTransform(Matrix4d E);
	double getRadius() { return m_r; }
	std::shared_ptr<Vector> getZAxis() { return m_Z; }
	std::shared_ptr<Node> getOrigin() { return m_O; }
	void setZAxis(std::shared_ptr<Vector> Z) { m_Z = Z; }
	void setOrigin(std::shared_ptr<Node> O) { m_O = O; }

protected:
	double m_r;
	double m_h;
	std::shared_ptr<Vector> m_Z;	//Z axis
	std::shared_ptr<Node> m_O;      // Origin
	Matrix4d E_wi;	// Where the component is wrt world
	Matrix4d E_ji;	// Where the component is wrt body
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, std::shared_ptr<MatrixStack> P)const;
	void update_();

	std::shared_ptr<Shape> m_shape;
	std::shared_ptr<Body> m_parent;
};