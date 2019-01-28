#pragma once

#ifndef MUSCLEMASS_SRC_WRAPSHPERE_H_
#define MUSCLEMASS_SRC_WRAPSHPERE_H_

/*
* WrapSphere.hpp
*
* Obstacle Set Algorithm Simulation for Sphere Obstacles
*
*/

#include "WrapObst.h"

class CompSphere;

class WrapSphere : public WrapObst
{
public:
	WrapSphere();
	WrapSphere(const std::shared_ptr<Node> &P, const std::shared_ptr<Node> &S, const std::shared_ptr<CompSphere> &compSphere, int num_points);
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void compute();
	Eigen::MatrixXd getPoints(int num_points);

private:
	std::shared_ptr<CompSphere> m_compSphere;

protected:
	void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P)const;
	void update_();
};

#endif //MUSCLEMASS_SRC_WRAPSHPERE_H_
