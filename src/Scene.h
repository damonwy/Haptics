#pragma once
#ifndef MUSCLEMASS_SRC_SCENE_H_
#define MUSCLEMASS_SRC_SCENE_H_

#include <vector>
#include <memory>
#include <string>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <nlohmann/json.hpp>
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
struct TrainingData;

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
	void saveEnergyData(int num_steps);
	void saveTrainingData(int num_steps, int gap);


	Eigen::VectorXd y;
	std::vector<Energy> m_energy_vector;
	std::vector<double> m_time_vector;
	std::vector<TrainingData> m_training_data_vector;

private:
	int count;
	double t;
	double h;
	int search_idx;
	int drawHz;
	double drawH;

	Eigen::Vector3d grav;

	nlohmann::json js;
	std::shared_ptr<World> m_world;
	std::shared_ptr<Solver> m_solver;
	std::shared_ptr<Solution> m_solution;
};

class Observer
{
public:
	Observer(std::vector<Eigen::VectorXd>& states, std::vector<double>& times): 
		m_states(states), m_times(times)
	{
	}

	void operator()(const Eigen::VectorXd& x, double t)
	{
		m_states.push_back(x);
		m_times.push_back(t);
	}

private:
	std::vector<Eigen::VectorXd>& m_states;
	std::vector<double>& m_times;
};

#endif // MUSCLEMASS_SRC_SCENE_H_
