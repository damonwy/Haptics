#pragma once

#include "Solver.h"

class SolverDense : public Solver {
public:
	SolverDense() {}
	SolverDense(std::shared_ptr<World> world, Integrator integrator);
	virtual ~SolverDense() {}
	std::shared_ptr<Solution> solve();
	Eigen::VectorXd dynamics(Eigen::VectorXd y);
	void initMatrix(int nm, int nr, int nem, int ner, int nim, int nir);

private:
	Eigen::MatrixXd Mm;

	Eigen::MatrixXd MDKr_;
	Eigen::MatrixXd K;
	Eigen::MatrixXd Km;

	Eigen::VectorXd fm;
	Eigen::MatrixXd J;

	Eigen::MatrixXd Jdot;
	Eigen::VectorXd q0;
	Eigen::VectorXd q1;
	Eigen::VectorXd qdot0;
	Eigen::VectorXd qdot1;
	Eigen::VectorXd qddot;

	Eigen::MatrixXd Mr;

	Eigen::MatrixXd Dm; // nm x nm
	Eigen::VectorXd tmp; // nm x 1
	Eigen::MatrixXd Dr;

	Eigen::MatrixXd Kr;
	Eigen::VectorXd fr;
	Eigen::VectorXd fr_;

	Eigen::MatrixXd Gm;
	Eigen::MatrixXd Gmdot;

	Eigen::VectorXd gm;
	Eigen::VectorXd gmdot;
	Eigen::VectorXd gmddot;

	Eigen::MatrixXd Gr;
	Eigen::MatrixXd Grdot;

	Eigen::VectorXd gr;
	Eigen::VectorXd grdot;
	Eigen::VectorXd grddot;

	Eigen::MatrixXd G;

	Eigen::VectorXd g;
	Eigen::VectorXd gdot;
	Eigen::VectorXd rhsG;

	Eigen::MatrixXd Cm;
	Eigen::MatrixXd Cmdot;
	Eigen::VectorXd cm;
	Eigen::VectorXd cmdot;
	Eigen::VectorXd cmddot;

	Eigen::MatrixXd Cr;
	Eigen::MatrixXd Crdot;
	Eigen::VectorXd cr;
	Eigen::VectorXd crdot;
	Eigen::VectorXd crddot;

	Eigen::VectorXd rhsC;

	Eigen::MatrixXd C;
	Eigen::VectorXd c;

	std::vector<int> rowsM;
	std::vector<int> rowsR;
};