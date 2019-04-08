#pragma once

#include "Solver.h"
#include "KKTSolver.h"


class SolverSparse : public Solver, public std::enable_shared_from_this<SolverSparse> {
public:
	SolverSparse() {}
	SolverSparse(std::shared_ptr<World> world, Integrator integrator, SparseSolver solver) : Solver(world, integrator), m_sparse_solver(solver) {}
	Eigen::VectorXd dynamics(Eigen::VectorXd y);
	Eigen::VectorXd dynamics_matlab(Eigen::VectorXd y);
	void test(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double);
	void initMatrix(int nm, int nr, int nem, int ner, int nim, int nir);
	Eigen::MatrixXd J_dense;	// dense_nm x dense_nr
	Eigen::MatrixXd Jdot_dense;
	Eigen::MatrixXd JMJ_mi;

	Eigen::SparseMatrix<double> Mm_sp;
	Eigen::SparseMatrix<double> Mr_sp;
	Eigen::VectorXd q0;
	Eigen::VectorXd q1;
	Eigen::VectorXd qdot0;
	Eigen::VectorXd qdot1;
	Eigen::VectorXd qddot;
	void exportTrainingData();

private:
	bool isCollided;
	SparseSolver m_sparse_solver;
	std::vector<T> Mm_;

	Eigen::SparseMatrix<double> MDKr_sp;
	Eigen::MatrixXd I_matlab;
	Eigen::MatrixXd Im_matlab;
	Eigen::MatrixXd Ir_matlab;

	Eigen::SparseMatrix<double> MDKr_sp_tp;

	Eigen::SparseMatrix<double> K_sp;
	std::vector<T> K_;

	Eigen::SparseMatrix<double> Km_sp;
	std::vector<T> Km_;

	Eigen::VectorXd fm;
	

	Eigen::SparseMatrix<double> J_sp;
	Eigen::SparseMatrix<double> J_t_sp;
	std::vector<T> J_;
	std::vector<T> J_pre;
	Eigen::SparseMatrix<double> Jdot_sp;
	std::vector<T> Jdot_;


	Eigen::VectorXd rhs;
	Eigen::VectorXd guess;

	int J_vec_idx;
	Eigen::SparseMatrix<double> Mr_sp_temp;
	Eigen::SparseMatrix<double> Dm_sp;
	std::vector<T> Dm_;
	Eigen::VectorXd tmp; // nm x 1
	Eigen::SparseMatrix<double> Dr_sp;
	std::vector<T> Dr_;

	Eigen::SparseMatrix<double> Kr_sp;
	std::vector<T> Kr_;
	Eigen::VectorXd fr;
	Eigen::VectorXd fr_;

	Eigen::SparseMatrix<double> Gm_sp;
	std::vector<T> Gm_;

	Eigen::SparseMatrix<double> Gmdot_sp;
	std::vector<T> Gmdot_;

	Eigen::VectorXd gm;
	Eigen::VectorXd gmdot;
	Eigen::VectorXd gmddot;

	Eigen::SparseMatrix<double> Gr_sp;
	std::vector<T> Gr_;

	Eigen::SparseMatrix<double> Grdot_sp;
	std::vector<T> Grdot_;

	Eigen::VectorXd gr;
	Eigen::VectorXd grdot;
	Eigen::VectorXd grddot;

	Eigen::VectorXd g;
	Eigen::VectorXd gdot;
	Eigen::VectorXd gddot;
	Eigen::VectorXd rhsG;

	Eigen::SparseMatrix<double> Cm_sp;
	std::vector<T> Cm_;
	Eigen::SparseMatrix<double> Cmdot_sp;
	std::vector<T> Cmdot_;

	Eigen::VectorXd cm;
	Eigen::VectorXd cmdot;
	Eigen::VectorXd cmddot;

	Eigen::SparseMatrix<double> Cr_sp;
	std::vector<T> Cr_;
	Eigen::SparseMatrix<double> Crdot_sp;
	std::vector<T> Crdot_;

	Eigen::VectorXd cr;
	Eigen::VectorXd crdot;
	Eigen::VectorXd crddot;

	Eigen::VectorXd rhsC;

	Eigen::MatrixXd C;
	Eigen::SparseMatrix<double> C_sp;
	std::vector<T> C_;
	Eigen::VectorXd c;

	std::vector<int> rowsM;
	std::vector<int> rowsR;
	std::vector<int> rowsEM;
	std::vector<int> rowsER;
	Eigen::MatrixXd G;

	Eigen::SparseMatrix<double, Eigen::RowMajor> G_sp;
	Eigen::SparseMatrix<double> G_sp_tp;
	Eigen::SparseMatrix<double> lhs_left_tp;
	Eigen::SparseMatrix<double> lhs_left;
	Eigen::SparseMatrix<double> lhs_right_tp;
	Eigen::SparseMatrix<double> lhs_right;
	Eigen::SparseMatrix<double> zero;
	Eigen::SparseMatrix<double> LHS_sp;

	int m_dense_nm;
	int m_dense_nr;

	Eigen::MatrixXd Cm;
	Eigen::MatrixXd Cmdot;

	Eigen::MatrixXd Cr;
	Eigen::MatrixXd Crdot;
	//
	Eigen::SparseLU<Eigen::SparseMatrix<double> > solver;
	Eigen::SparseMatrix<double> D_sp;

	// Muscle Inertia Matrix
	Eigen::VectorXd Jf_mi;
	Eigen::VectorXd fvm;
	Eigen::MatrixXd JMJdot_mi;

	Eigen::VectorXd m_fk;
	Eigen::VectorXd m_fk_matlab;



};