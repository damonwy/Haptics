#pragma once
#ifndef MUSCLEMASS_SRC_WRAPDOUBLECYLINDER_H_
#define MUSCLEMASS_SRC_WRAPDOUBLECYLINDER_H_

/*
* WrapDoubleCylinder.hpp
*
* Obstacle Set Algorithm Simulation for Double Cylinder Obstacles
*
*/

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include "WrapObst.h"
#include "Node.h"
#include "Vector.h"

class Node;
class Vector;
class CompDoubleCylinder;

class WrapDoubleCylinder : public WrapObst
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:

	Eigen::MatrixXd
		M_U,          // Obstacle Coord Transformation Matrix for U
		M_V;          // Obstacle Coord Transformation Matrix for V

	double
		m_radius_U,     // U Cylinder Radius
		m_radius_V;     // V Cylinder Radius

	Status
		status_U,     // U Wrapping Status
		status_V;     // V Wrapping Status

	std::shared_ptr<Node> m_point_U;	// U Cylinder Origin

	std::shared_ptr<Vector> m_vec_z_U;		// Z axis of Cylinder U(direction matters)

	std::shared_ptr<Node> m_point_V;		// V Cylinder Origin

	std::shared_ptr<Vector> m_vec_z_V;		// Z axis of Cylinder V(direction matters)
	std::shared_ptr<Node> m_point_g;
	std::shared_ptr<Node> m_point_h;
	std::shared_ptr<CompDoubleCylinder> m_compDoubleCylinder;

public:

	WrapDoubleCylinder();
    virtual ~WrapDoubleCylinder(){}
	WrapDoubleCylinder(const std::shared_ptr<Node> &P,
		const std::shared_ptr<Node> &S,
		const std::shared_ptr<CompDoubleCylinder> compDoubleCylinder,
		const int num_points);

	void compute();
	Eigen::MatrixXd getPoints(int num_points) const;

	void init();
	void load(const std::string &RESOURCE_DIR);
	
protected:
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const;
	void update_();

};

#endif // MUSCLEMASS_SRC_WRAPDOUBLECYLINDER_H_
