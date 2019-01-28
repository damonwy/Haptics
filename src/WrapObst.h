#pragma once

#ifndef MUSCLEMASS_SRC_WRAPOBST_H_
#define MUSCLEMASS_SRC_WRAPOBST_H_

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <complex>
#include <vector>
#include <memory>
#include <string>
#include "Node.h"
#include "Body.h"
#include "MatrixStack.h"
#include "Program.h"

enum Status { wrap, inside_radius, no_wrap, empty_status};
enum Type { none, sphere, cylinder, double_cylinder };

class WrapObst
{
protected:
	std::shared_ptr<Node> m_point_P;	// Bounding-Fixed Via Point 1
	std::shared_ptr<Node> m_point_S;	// Bounding-Fixed Via Point 2
	std::shared_ptr<Node> m_point_O;	// Obstacle Center Point

	std::shared_ptr<Node> m_point_q;	// Obstacle Via Point 1 in Obstacle Frame
	std::shared_ptr<Node> m_point_t;	// Obstacle Via Point 2 in Obstacle Frame

	Eigen::MatrixXd M;	// Obstacle Coord Transformation Matrix
	Status m_status;		// Wrapping Status
	Type m_type;			// Obstacle Type
	double m_path_length;	// Wrapping Path Length
	double m_radius;			// Obstacle sphere radius
	int m_num_points;					// Number of points
	Eigen::MatrixXd m_arc_points;		// Each col stores the pos of a point 
	virtual void update_() {}
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P)const {}

public:
	std::shared_ptr<WrapObst> next;
	// Set muscle origin point
	void setOrigin(const std::shared_ptr<Node> &P) {
		m_point_P = P;
	}

	// Set muscle insertion point
	void setInsertion(const std::shared_ptr<Node> &S) {
		m_point_S = S;
	}

	WrapObst() {

		m_point_P = std::make_shared<Node>();
		m_point_P->x.setZero(); 
		m_point_S = std::make_shared<Node>();
		m_point_S->x.setZero();
		m_point_O = std::make_shared<Node>();
		m_point_O->x.setZero();
		m_point_q = std::make_shared<Node>();
		m_point_q->x.setZero();
		m_point_t = std::make_shared<Node>();
		m_point_t->x.setZero();

		M = Eigen::MatrixXd(3, 3);
		m_status = empty_status;
		m_path_length = 0.0;
		m_radius = 0.0;
		m_type = none;
	}

	// Constructor
	WrapObst(const std::shared_ptr<Node> &P,
		const std::shared_ptr<Node> &S, int num_points) :
		m_point_P(P), m_point_S(S), m_num_points(num_points)
	{
		m_point_O = std::make_shared<Node>();
		m_point_O->x.setZero();
		m_point_q = std::make_shared<Node>();
		m_point_q->x.setZero();
		m_point_t = std::make_shared<Node>();
		m_point_t->x.setZero();

		M.resize(3, 3); 
		M.setZero();
		m_status = empty_status;
		m_path_length = 0.0;
		m_radius = 0.0;
		m_type = none;
	}

	// wrap calculation
	virtual void compute() {}
	virtual void init() {}
	virtual void load(const std::string &RESOURCE_DIR) {}

	void update() {
		update_();
		if (next != nullptr) {
			next->update();
		}
	}

	virtual double getLength() { return this->m_path_length; }
	virtual Status getStatus() { return this->m_status; }
	double getRadius() { return this->m_radius; }
	virtual Eigen::MatrixXd getPoints() { return m_arc_points; }

	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P)const {
		draw_(MV, prog, progSimple, P);
		if (next != nullptr) {
			next->draw(MV, prog, progSimple, P);
		}
	}


};


#endif // MUSCLEMASS_SRC_WRAPOBST_H_
