#pragma once
#ifndef MUSCLEMASS_SRC_SCENE_H_
#define MUSCLEMASS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <json.hpp>
#include "MLCommon.h"

class Node;
class MatrixStack;
class Program;
class Shape;
class Solver;
class Spring;
class Joint;
class Vector;
class Stepper;
class World;
struct Solution;

class Scene
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW	
	Scene();
	virtual ~Scene();
	
	void load(const std::string &RESOURCE_DIR);
	void init();
	void reset();
	void step();
	void solve();
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, const std::shared_ptr<Program> prog3, std::shared_ptr<MatrixStack> P) const;
	double getTime() const { return t; }
	Eigen::VectorXd y;
private:
	int count;
	double t;
	double h;
	int time_step;
	int search_idx;
	int drawHz;
	double drawH;
	double tk;

	Eigen::Vector3d grav;

	nlohmann::json js;
	std::shared_ptr<World> m_world;
	std::shared_ptr<Solver> m_solver;
	std::shared_ptr<Solution> m_solution;

};

#endif // MUSCLEMASS_SRC_SCENE_H_