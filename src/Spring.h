#pragma once
#ifndef MUSCLEMASS_SRC_SPRING_H_
#define MUSCLEMASS_SRC_SPRING_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class MatrixStack;
class Program;
class Node;
class Body;

typedef Eigen::Triplet<double> T;


class Spring {

public:
	Spring();
	virtual ~Spring() {}

	void computeEnergies(Vector3d grav, Energy &ener);
	void computeForceStiffnessDamping(Eigen::VectorXd &f, Eigen::MatrixXd &K, Eigen::MatrixXd &D);
	void computeForceStiffnessDampingSparse(Eigen::VectorXd &f, std::vector<T> &K_, std::vector<T> &D_);

	void computeStiffnessProd(Eigen::VectorXd x, Eigen::VectorXd &y);
	void computeDampingProd(Eigen::VectorXd x, Eigen::VectorXd &y);
	void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
	void init();
	void update();	

	virtual void load(const std::string &RESOURCE_DIR) {}
	virtual void init_() {}
	virtual void update_() {}

	std::shared_ptr<Spring> next;
	std::vector<std::shared_ptr<Node>> m_nodes;

protected:
	virtual void computeForceStiffnessDampingSparse_(Eigen::VectorXd &f, std::vector<T> &K_, std::vector<T> &D_) {}
	virtual void computeForceStiffnessDamping_(Eigen::VectorXd &f, Eigen::MatrixXd &K, Eigen::MatrixXd &D) {}
	virtual void computeStiffnessProd_(Eigen::VectorXd x, Eigen::VectorXd &y) {}
	virtual void computeDampingProd_(Eigen::VectorXd x, Eigen::VectorXd &y) {}
	virtual void computeEnergies_(Vector3d grav, Energy &ener) {}
	virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {}

};

#endif // MUSCLEMASS_SRC_SPRING_H_