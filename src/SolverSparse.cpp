#include "rmpch.h"
//#include <Eigen/PardisoSupport>
//#include <Eigen/SuperLUSupport>
//#include <unsupported/Eigen/src/IterativeSolvers/MINRES.h>
//#include <unsupported\Eigen\src\IterativeSolvers\Scaling.h>
//#include <unsupported\Eigen\src\IterativeSolvers\GMRES.h>

#include "SolverSparse.h"
#include "World.h"
#include "Body.h"
#include "Joint.h"
#include "Spring.h"
#include "Muscle.h"
#include "MuscleSpring.h"
#include "SpringDamper.h"
#include "Deformable.h"
#include "DeformableSpring.h"
#include "ConstraintJointLimit.h"
#include "ConstraintLoop.h"
#include "ConstraintAttachSpring.h"
#include "QuadProgMosek.h"

using namespace std;
using namespace Eigen;

//#define CHECK_ENERGY

void SolverSparse::initMatrix(int nm, int nr, int nem, int ner, int nim, int nir) {
	ni = nim + nir;
	int nre = nr + ne;

	fm.setZero();
	fr.setZero();
	fr_.setZero();
	tmp.setZero();
	
	Mr_sp.resize(nr, nr);
	//Mr_sp.data().squeeze();
	Mr_sp_temp.resize(nr, nr);

	//Mr_sp_temp.data().squeeze();

	MDKr_sp.resize(nr, nr);
	//MDKr_sp.data().squeeze();
	MDKr_sp_tp.resize(nr, nr);
	//MDKr_sp_tp.data().squeeze();
	lhs_left_tp.resize(nr, nre);
	//lhs_left_tp.data().squeeze();
	lhs_right_tp.resize(ne, nre);
	//lhs_right_tp.data().squeeze();
	lhs_left.resize(nre, nr);
	//lhs_left.data().squeeze();
	lhs_right.resize(nre, ne);
	//lhs_right.data().squeeze();
	LHS_sp.resize(nre, nre);
	//LHS_sp.data().squeeze();

	Kr_sp.resize(nr, nr);
	//Kr_sp.data().squeeze();
	Kr_.clear();

	Dr_sp.resize(nr, nr);
	//Dr_sp.data().squeeze();
	Dr_.clear();
	
	Dm_sp.resize(nm, nm);
	//Dm_sp.data().squeeze();
	Dm_.clear();

	K_sp.resize(nm, nm);
	//K_sp.data().squeeze();
	K_.clear();

	Km_sp.resize(nm, nm);
	//Km_sp.data().squeeze();
	Km_.clear();


	J_dense.setZero();
	Jdot_dense.setZero();

	J_sp.resize(nm, nr);
	J_t_sp.resize(nr, nm);
	//J_sp.data().squeeze();
	//J_.clear();

	Jdot_sp.resize(nm, nr);
	//Jdot_sp.data().squeeze();
	Jdot_.clear();
	
	//Gm_sp.resize(nem, nm);
	//Gm_sp.data().squeeze();
	Gm_.clear();

	//Gmdot_sp.resize(nem, nm);
	//Gmdot_sp.data().squeeze();
	Gmdot_.clear();

	gm.setZero();
	gmdot.setZero();
	gmddot.setZero();

	//Gr_sp.resize(ner, nr);
	//Gr_sp.data().squeeze();
	Gr_.clear();

	//Grdot_sp.resize(ner, nr);
	//Grdot_sp.data().squeeze();
	Grdot_.clear();

	gr.setZero();
	grdot.setZero();
	grddot.setZero();

	g.setZero();
	gdot.setZero();
	gddot.setZero();

	rhsG.setZero();
	rhs.setZero();
	guess.setZero();

	Cm_sp.resize(nim, nm);
	//Cm_sp.data().squeeze();
	Cm_.clear();

	Cmdot_sp.resize(nim, nm);
	//Cmdot_sp.data().squeeze();
	Cmdot_.clear();

	cm.resize(nim);
	cm.setZero();
	cmdot.resize(nim);
	cmdot.setZero();
	cmddot.resize(nim);
	cmddot.setZero();

	Cr_sp.resize(nir, nr);
	//Cr_sp.data().squeeze();
	Cr_.clear();

	Crdot_sp.resize(nir, nr);
	//Crdot_sp.data().squeeze();
	Crdot_.clear();

	cr.resize(nir);
	cr.setZero();
	crdot.resize(nir);
	crdot.setZero();
	crddot.resize(nir);
	crddot.setZero();

	//G_sp.resize(ne, nr);
	//G_sp.data().squeeze();
	//G_sp_tp.resize(nr, ne);
	//G_sp_tp.data().squeeze();

	Cm.resize(nim, nm);
	Cm.setZero();
	Cmdot.resize(nim, nm);
	Cmdot.setZero();
	Cr.resize(nir, nr);
	Cr.setZero();
	Crdot.resize(nir, nr);
	Crdot.setZero();

	JMJ_mi.setZero();
	Jf_mi.setZero();
	fvm.setZero();
}

VectorXd SolverSparse::dynamics(VectorXd y)
{
	cout << "Input:" << endl << y << endl << endl;

	//SparseMatrix<double, RowMajor> G_sp;
		if (step == 0) {
			// constant during simulation
			isCollided = false;
			nr = m_world->nr;
			nm = m_world->nm;

			nem = m_world->nem;
			ner = m_world->ner;
			ne = ner + nem;

			nim = m_world->nim;
			nir = m_world->nir;
			ni = nim + nir;

			m_ntets = m_world->m_ntets;

			Mm_sp.resize(nm, nm);
			Mm_sp.data().squeeze();
			Mm_.clear();
			Mm_.reserve(nm);  
			
			m_dense_nm = m_world->m_dense_nm;
			m_dense_nr = m_world->m_dense_nr;

			yk.resize(2 * nr);
			ydotk.resize(2 * nr);

			fm.resize(nm);
			fr.resize(nr);
			fr_.resize(nr);

			tmp.resize(nm);

			J_dense.resize(m_dense_nm, m_dense_nr);
			Jdot_dense.resize(m_dense_nm, m_dense_nr);
			

			gr.resize(ner);
			grdot.resize(ner);
			grddot.resize(ner);
			g.resize(ne);
			rhsG.resize(ne);
			rhs.resize(nr + ne);
			guess.resize(nr + ne);
			gdot.resize(ne);
			gddot.resize(ne);
			gm.resize(nem);
			gmdot.resize(nem);
			gmddot.resize(nem);
			Grdot_sp.resize(ner, nr);
			Gm_sp.resize(nem, nm);
			Gmdot_sp.resize(nem, nm);
			Gr_sp.resize(ner, nr);

			body0 = m_world->getBody0();
			joint0 = m_world->getJoint0();
			deformable0 = m_world->getDeformable0();
			constraint0 = m_world->getConstraint0();
			spring0 = m_world->getSpring0();
			muscle0 = m_world->getMuscle0();

			t = m_world->getTspan()(0);
			h = m_world->getH();
			hsquare = h * h;
			this->grav = m_world->getGrav();
			zero.resize(ne, ne);

			deformable0->computeJacobianSparse(J_);		
			J_vec_idx = static_cast<int>(J_.size());

			JMJ_mi.resize(nr, nr);
			Jf_mi.resize(nr);
			fvm.resize(nr);
			// JMJdot_mi.resize
		}

#ifdef CHECK_ENERGY
		Energy ener = m_world->computeEnergy();
		cout << "V" << ener.V << endl;
		cout << "K" << ener.K << endl;
		cout << " sum " << ener.V + ener.K << endl;

#endif // CHECK_ENERGY

		nim = m_world->nim;
		nir = m_world->nir;
		
		ni = nim + nir;	

		q0 = y.segment(0, nr);
		qdot0 = y.segment(nr, nr);

		initMatrix(nm, nr, nem, ner, nim, nir);

		if (step == 0) {
			body0->computeMassSparse(Mm_);
			deformable0->computeMassSparse(Mm_);
			Mm_sp.setFromTriplets(Mm_.begin(), Mm_.end());
		}

		double t_i = m_world->getTime();
		switch (m_world->m_type) {
		case TEST_MAXIMAL_HYBRID_DYNAMICS:
			m_world->sceneTestMaximalHD(t_i);
			break;
		case TEST_HYPER_REDUCED_COORDS:
			m_world->sceneTestHyperReduced(t_i);
			break;
		case TEST_CONSTRAINT_PRESC_BODY_ATTACH_POINT:
			m_world->sceneAttachPoint(t_i);
		default:
			break;
		}


		body0->computeGrav(grav, fm);
		body0->computeForceDampingSparse(tmp, Dm_);
		deformable0->computeForce(grav, fm);
		deformable0->computeForceDampingSparse(grav, tmp, Dm_);
		joint0->computeForceStiffnessSparse(fr, Kr_);
		joint0->computeForceDampingSparse(tmp, Dr_);

		//// First get dense jacobian (only a small part of the matrix)
		joint0->computeJacobian(J_dense, Jdot_dense);
	
		//// Push back the dense part
		if (step == 0) {
			for (int i = 0; i < J_dense.rows(); ++i) {
				for (int j = 0; j < J_dense.cols(); ++j) {
					J_.push_back(T(i, j, J_dense(i, j)));
					Jdot_.push_back(T(i, j, Jdot_dense(i, j)));
				}
			}
		}
		else {
			int idx = J_vec_idx;
			for (int i = 0; i < J_dense.rows(); ++i) {
				for (int j = 0; j < J_dense.cols(); ++j) {
					J_[idx] = (T(i, j, J_dense(i, j)));
					Jdot_.push_back(T(i, j, Jdot_dense(i, j)));
					idx++;
				}
			}
		}

		spring0->computeForceStiffnessDampingSparse(fm, Km_, Dm_);
		//cout <<"JMJ_B:" << JMJ_mi << endl;
		muscle0->computeJMJ(JMJ_mi, m_world);
		//cout << "JMJ_A:" << JMJ_mi << endl;
		muscle0->computeForce(grav, Jf_mi);
		//cout << "Jf_mi" << Jf_mi << endl;
		muscle0->computeJMJdotqdot(fvm, qdot0, m_world);
		cout << "fk" << fvm << endl;

		Km_sp.setFromTriplets(Km_.begin(), Km_.end());
		Dm_sp.setFromTriplets(Dm_.begin(), Dm_.end());
		Dr_sp.setFromTriplets(Dr_.begin(), Dr_.end());
		K_sp.setFromTriplets(K_.begin(), K_.end()); // check
		//cout << "K" << K_.size() << endl;

		Kr_sp.setFromTriplets(Kr_.begin(), Kr_.end());
		J_sp.setFromTriplets(J_.begin(), J_.end()); // check

		Jdot_sp.setFromTriplets(Jdot_.begin(), Jdot_.end());

		J_t_sp = J_sp.transpose();

		Mr_sp = J_t_sp * (Mm_sp - hsquare * K_sp) * J_sp + JMJ_mi;
		// 
		//cout << "I " << endl << MatrixXd(Mr_sp) << endl;
		//cout << "JMJ_mi "<< JMJ_mi + MatrixXd(Mr_sp) << endl << endl;
		//Mr_sp_temp = Mr_sp.transpose();
		//Mr_sp += Mr_sp_temp;
		//Mr_sp *= 0.5;
	
		fr_ = /*Mr_sp * qdot0 + h* */ (J_t_sp * (fm - Mm_sp * Jdot_sp * qdot0) + fr + Jf_mi + fvm); 
		//
		
		cout << "fr " << endl << fr_ << endl;

		MDKr_sp = Mr_sp + J_t_sp * (h * Dm_sp - hsquare * Km_sp) * J_sp + h * Dr_sp - hsquare * Kr_sp;
		//cout << MatrixXd(MDKr_sp) << endl << endl;
		//cout << "fm" << fm << endl << endl;
		//cout << "fr" << fr << endl << endl;
		//cout << "Jdot_sp" << MatrixXd(Jdot_sp) << endl << endl;
		//cout << "Mr_sp"<< endl << MatrixXd(Mr_sp) << endl << endl;
		//cout << "J_sp" << endl << MatrixXd(J_sp) << endl << endl;
		//cout << "Dm_sp" << endl << MatrixXd(Dm_sp) << endl << endl;
		//cout << "Km_sp" << endl << MatrixXd(Km_sp) << endl << endl;
		//cout << "fr_"<< (fr_) << endl << endl;
		//cout <<"fm"<< fm << endl << endl;
	
		Mr_sp = J_t_sp * (Mm_sp - hsquare * K_sp) * J_sp;

		//Mr_sp_temp = Mr_sp.transpose();
		//Mr_sp += Mr_sp_temp;
		//Mr_sp *= 0.5;

		/*	MatrixXd Mrnew = MatrixXd(JmR.transpose() * (Mm_sp - hsquare * K_sp) * JmR);
		VectorXd fnew = fr.segment<1>(0);
		VectorXd qdot0new = qdot0.segment<1>(0);

		VectorXd f_new = Mrnew * qdot0new + h * (JMR.transpose() * (fm - Mm_sp * Jdot_sp * qdot0) + fnew);
		MatrixXd MDKr_new = Mrnew + JMR.transpose() * (h * Dm_sp - hsquare * Km_sp) * JMR;
		qdot1 = JrR * (MDKr_new.ldlt().solve(f_new));
		cout << qdot1 << endl;*/

		if (ne > 0) {
			constraint0->computeJacEqMSparse(Gm_, Gmdot_, gm, gmdot, gmddot);		
			constraint0->computeJacEqRSparse(Gr_, Grdot_, gr, grdot, grddot);
			
			rowsEM.clear();
			rowsER.clear();
			constraint0->getEqActiveList(rowsEM, rowsER);
			nem = static_cast<int>(rowsEM.size());
			ner = static_cast<int>(rowsER.size());
			ne = nem + ner;

			if (ne > 0) {
				Eigen::VectorXi m_rowsEM = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsEM.data(), rowsEM.size());
				Eigen::VectorXi m_rowsER = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsER.data(), rowsER.size());
				Gm_sp.setFromTriplets(Gm_.begin(), Gm_.end());
				Gmdot_sp.setFromTriplets(Gmdot_.begin(), Gmdot_.end());
				Gr_sp.setFromTriplets(Gr_.begin(), Gr_.end());
				Grdot_sp.setFromTriplets(Grdot_.begin(), Grdot_.end());

				MatrixXd m_Gm = MatrixXd(Gm_sp)(m_rowsEM, Eigen::placeholders::all);
				//cout << m_Gm << endl << endl;
				

				MatrixXd m_Gr = MatrixXd(Gr_sp)(m_rowsER, Eigen::placeholders::all);
				VectorXd m_gm = gm(m_rowsEM);
				//cout << m_gm << endl << endl;

				VectorXd m_gr = gr(m_rowsER);
				VectorXd m_gmdot = gmdot(m_rowsEM);
				VectorXd m_grdot = grdot(m_rowsER);
				VectorXd m_gmddot = gmddot(m_rowsEM);
				VectorXd m_grddot = grddot(m_rowsER);
				MatrixXd GmJ = m_Gm * MatrixXd(J_sp);
				//cout << "J" << endl << MatrixXd(J_sp) << endl;

				G.resize(GmJ.rows() + m_Gr.rows(), m_Gr.cols());
				///G_sp.resize(GmJ.rows() + m_Gr.rows(), m_Gr.cols());
				G << GmJ, m_Gr;
				//G_sp = G.sparseView();
				rhsG.resize(G.rows());
				VectorXd g(G.rows());
				g << m_gm, m_gr;
				VectorXd gdot(G.rows());
				gdot << m_gmdot, m_grdot;
				rhsG = -  gdot - 5.0 * g;
			}

			//Gm_sp.setFromTriplets(Gm_.begin(), Gm_.end());
			//Gmdot_sp.setFromTriplets(Gmdot_.begin(), Gmdot_.end());
			//Gr_sp.setFromTriplets(Gr_.begin(), Gr_.end());
			//Grdot_sp.setFromTriplets(Grdot_.begin(), Grdot_.end());

			//sparse_to_file_as_dense(Gm_sp * J_sp, "Gm_sp * J_sp");
			//G_sp.topRows(nem) = Gm_sp * J_sp;
			//G_sp.bottomRows(ner) = Gr_sp;
			
			//sparse_to_file_as_dense(G_sp, "G_sp");
			/*g.segment(0, nem) = gm;
			g.segment(nem, ner) = gr;
			gdot.segment(0, nem) = gmdot;
			gdot.segment(nem, ner) = grdot;
			gddot.segment(0, nem) = gmddot;
			gddot.segment(nem, ner) = grddot;
			rhsG = - gdot - 5.0 * g ;*/
		}

		if (ni > 0) {
			// Check for active inequality constraint
			//constraint0->computeJacIneqMSparse(Cm, Cmdot, cm, cmdot, cmddot);
			//constraint0->computeJacIneqRSparse(Cr, Crdot, cr, crdot, crddot);
			constraint0->computeJacIneqM(Cm, Cmdot, cm, cmdot, cmddot);
			constraint0->computeJacIneqR(Cr, Crdot, cr, crdot, crddot);

			rowsR.clear();
			rowsM.clear();

			constraint0->getActiveList(rowsM, rowsR);
			nim = static_cast<int>(rowsM.size());
			nir = static_cast<int>(rowsR.size());
			ni = nim + nir;

			if (ni > 0) {
				Eigen::VectorXi m_rowsM = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsM.data(), rowsM.size());
				Eigen::VectorXi m_rowsR = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(rowsR.data(), rowsR.size());

				MatrixXd m_Cm = Cm(m_rowsM, Eigen::placeholders::all);
				MatrixXd m_Cr = Cr(m_rowsR, Eigen::placeholders::all);
				VectorXd m_cm = cm(m_rowsM);
				VectorXd m_cr = cr(m_rowsR);
				VectorXd m_cmdot = cmdot(m_rowsM);
				VectorXd m_crdot = crdot(m_rowsR);

				MatrixXd CmJ = m_Cm * MatrixXd(J_sp);
				C.resize(CmJ.rows() + m_Cr.rows(), m_Cr.cols());
				C << CmJ, m_Cr;
				rhsC.resize(C.rows());
				VectorXd c(C.rows());
				c << m_cm, m_cr;
				VectorXd cdot(C.rows());
				cdot << m_cmdot, m_crdot;
				rhsC = -cdot - 5.0 * c;
			}
		}

		if (ne == 0 && ni == 0) {	// No constraints
			ConjugateGradient< SparseMatrix<double> > cg;
			cg.setMaxIterations(100000);
			cg.setTolerance(1e-10);
			cg.compute(MDKr_sp);
			if (cg.info() != Success) {
				// solve() failed

				cout << "solve failed" << endl << endl;
				exit(1);
			}

			//qdot1 = cg.solveWithGuess(fr_, qdot0);
			qddot = cg.solve(fr_);

			/*cout << MatrixXd(MDKr_sp) << endl << endl;
			cout << fr_ << endl << endl;
			cout << qdot1 << endl;*/
		}
		else if (ne > 0 && ni == 0) {  // Just equality
			//int rows = nr + ne;
			//int cols = nr + ne;

			/*
			//Slow

			MatrixXd LHS(rows, cols);		
			LHS.setZero();
			
			MatrixXd MDKr_ = MatrixXd(MDKr_sp);
			MatrixXd G = MatrixXd(G_sp);

			LHS.block(0, 0, nr, nr) = MDKr_;
			LHS.block(0, nr, nr, ne) = G.transpose();
			LHS.block(nr, 0, ne, nr) = G;
			SparseMatrix<double> LHS_sp(rows, cols);
			LHS_sp = LHS.sparseView(1e-8);
			
			*/
			int nre = nr + ne;
			lhs_left_tp.resize(nr, nre);
			lhs_right_tp.resize(ne, nre);
			lhs_left.resize(nre, nr);
			lhs_right.resize(nre, ne);

			LHS_sp.resize(nre, nre);
			guess.segment(0, nr) = qdot0;
			SparseMatrix<double> Gtemp = G.sparseView();
			SparseMatrix<double> Gtemp_tp = Gtemp.transpose();
			// Assemble sparse matrices
			MDKr_sp_tp = MDKr_sp.transpose();
			//G_sp_tp = G_sp.transpose();

			// Combine MDKr' and G' by column
			lhs_left_tp.leftCols(nr) = MDKr_sp_tp;
			//lhs_left_tp.rightCols(ne) = G_sp_tp;
			lhs_left_tp.rightCols(ne) = Gtemp_tp;

			// Combine G and Z by column
			//lhs_right_tp.leftCols(nr) = G_sp;
			//lhs_right_tp.rightCols(ne) = zero;
			lhs_right_tp.leftCols(nr) = Gtemp;
			zero.resize(ne, ne);
			lhs_right_tp.rightCols(ne) = zero;


			lhs_left = lhs_left_tp.transpose();  // rows x nr
			lhs_right = lhs_right_tp.transpose(); // rows x ne

			LHS_sp.leftCols(nr) = lhs_left;
			LHS_sp.rightCols(ne) = lhs_right;
			
			rhs.resize(nre);
			rhs.segment(0, nr) = fr_;
			rhs.segment(nr, ne) = rhsG;


			//cout << MatrixXd(LHS_sp) << endl << endl;
			//cout << rhs << endl << endl;

			//VectorXd sol = LHS.ldlt().solve(rhs);
			//qdot1 = sol.segment(0, nr);
			//VectorXd l = sol.segment(nr, sol.rows() - nr);
			switch (m_sparse_solver)
			{
			case CG: 
				{
					ConjugateGradient< SparseMatrix<double>, Lower | Upper> cg;
					//cg.setMaxIterations(2000);
					cg.setTolerance(1e-3);
					cg.compute(LHS_sp);
					qdot1 = cg.solveWithGuess(rhs, guess).segment(0, nr);
					
					//std::cout << "#iterations:     " << cg.iterations() << std::endl;
					//std::cout << "estimated error: " << cg.error() << std::endl;
					break;
				}
			case SLDLT:
				{
					SimplicialLDLT<SparseMatrix<double>, Lower, NaturalOrdering<int> > sldlt;

					sldlt.compute(LHS_sp);
					qdot1 = sldlt.solve(rhs).segment(0, nr);
					break;
				}			
			case LU: 
				{
					LHS_sp.makeCompressed();
					if (step == 0) {
						solver.analyzePattern(LHS_sp);
					}
					solver.compute(LHS_sp);
					if (solver.info() != Success) {
						// decomposition failed

						cout << "decomposition failed" << endl << endl;
						exit(1);
					}

					solver.factorize(LHS_sp);
					qdot1 = solver.solve(rhs).segment(0, nr);
					cout << MatrixXd(LHS_sp) << endl << endl;
					cout << rhs << endl << endl;
					cout << qdot1 << endl;
					break;
				}		
			default:
				{
					shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
					program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
					program_->setParamInt(MSK_IPAR_LOG, 10);
					program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
					program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
					program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
					program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
					program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
					program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
					program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);
					program_->setNumberOfVariables(nr);
					program_->setObjectiveMatrix(MDKr_sp);

					program_->setObjectiveVector(-fr_);

					program_->setNumberOfEqualities(ne);
					program_->setEqualityMatrix(G_sp);

					program_->setEqualityVector(rhsG);

					bool success = program_->solve();
                    assert(success == true);
					VectorXd sol = program_->getPrimalSolution();
					qdot1 = sol.segment(0, nr);	
					VectorXd l = program_->getDualEquality();
					constraint0->scatterForceEqM(MatrixXd(Gm_sp.transpose()), l.segment(0, nem) / h);
					constraint0->scatterForceEqR(MatrixXd(Gr_sp.transpose()), l.segment(nem, l.rows() - nem) / h);
				}
				break;
			}

		}
		else if (ne == 0 && ni > 0) {  // Just inequality
			shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
			program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
			program_->setParamInt(MSK_IPAR_LOG, 10);
			program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);

			program_->setNumberOfVariables(nr);
			program_->setObjectiveMatrix(MDKr_sp);
			program_->setObjectiveVector(-fr_);
			program_->setNumberOfInequalities(ni);
			program_->setInequalityMatrix(C.sparseView());

			VectorXd cvec(ni);
			cvec.setZero();
			program_->setInequalityVector(cvec);

			bool success = program_->solve();
            assert(success == true);
			VectorXd sol = program_->getPrimalSolution();
			qdot1 = sol.segment(0, nr);
		}
		else {  // Both equality and inequality
			shared_ptr<QuadProgMosek> program_ = make_shared <QuadProgMosek>();
			program_->setParamInt(MSK_IPAR_OPTIMIZER, MSK_OPTIMIZER_INTPNT);
			program_->setParamInt(MSK_IPAR_LOG, 10);
			program_->setParamInt(MSK_IPAR_LOG_FILE, 1);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_DFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_INFEAS, 1e-10);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_MU_RED, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_NEAR_REL, 1e3);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_PFEAS, 1e-8);
			program_->setParamDouble(MSK_DPAR_INTPNT_QO_TOL_REL_GAP, 1e-8);
			program_->setNumberOfVariables(nr);

			program_->setObjectiveMatrix(MDKr_sp);
			program_->setObjectiveVector(-fr_);
			program_->setNumberOfInequalities(ni);
			program_->setInequalityMatrix(C.sparseView());
			program_->setNumberOfEqualities(ne);
			VectorXd cvec(ni);
			cvec.setZero();

			program_->setInequalityVector(cvec);
			program_->setEqualityMatrix(G_sp);

			VectorXd gvec(ne);
			gvec.setZero();
			program_->setEqualityVector(rhsG);

			bool success = program_->solve();
            assert(success == true);
			VectorXd sol = program_->getPrimalSolution();
			qdot1 = sol.segment(0, nr);
		}

		VectorXd ydot_matlab = dynamics_matlab(y);
		VectorXd qddot_matlab = ydot_matlab.segment(nr, nr);
		VectorXd qdot1_matlab = qddot_matlab * h + qdot0;
		VectorXd q1_matlab = q0 + h * qdot1_matlab;
		
		//qddot = (qdot1 - qdot0) / h;
		qdot1 = qddot * h + qdot0;
		q1 = q0 + h * qdot1;
		yk.segment(0, nr) = q1;
		yk.segment(nr, nr) = qdot1;

		ydotk.segment(0, nr) = qdot1;
		ydotk.segment(nr, nr) = qddot;

		cout <<"qddot "<< endl << qddot << endl;
		cout << "qddot_matlab"<<endl << qddot_matlab << endl;
		cout << "diff qddot" << (qddot - qddot_matlab).norm() << endl;
		
		// Use the result of matlab for comparison
		q1_matlab = q0 + h * qdot0;
		//qdot1_matlab = qdot0;

		yk.segment(0, nr) = q1_matlab;
		yk.segment(nr, nr) = qdot1_matlab;
		ydotk.segment(0, nr) = qdot0;
		ydotk.segment(nr, nr) = qddot_matlab;

		//cout << "yk" << endl << yk << endl;
		//cout << "ydotk" << endl << ydotk << endl;

		joint0->scatterDofs(yk, nr);
		joint0->scatterDDofs(ydotk, nr);
		//joint0->reparam();
		//joint0->gatherDofs(yk, nr);

		deformable0->scatterDofs(yk, nr);
		deformable0->scatterDDofs(ydotk, nr);

		step++;
		return yk;
}

Eigen::VectorXd SolverSparse::dynamics_matlab(Eigen::VectorXd y) {
	VectorXd y_ = y;

	// Parameters
	double M = 10.0;
	double L = 10.0;
	double l = 1.0;
	double r = 0.5;
	double g = 9.81;
	double mu0 = 1.0;
	double mu1 = 1.0;
	double mum = 1.0;

	// Unpack data
	VectorXd q0 = y.segment(0, nr);
	q0(0) = y_(1);
	q0(1) = y_(0);

	q0(0) += M_PI / 2.0;

	VectorXd qdot0 = y.segment(nr, nr);
	qdot0(0) = y_(3);
	qdot0(1) = y_(2);

	double c0 = cos(q0(0));
	double s0 = sin(q0(0));
	double c1 = cos(q0(1));
	double s1 = sin(q0(1));
	double c01 = cos(q0(0) + q0(1));
	double s01 = sin(q0(0) + q0(1));

	// Compute inertia
	Matrix2d I0, I1, Im, I, dIdtheta1;
	I0 << l * l, 0.0, 0.0, 0.0;
	I0 *= mu0 / 3.0;
	//cout << "I0 matlab" << endl << I0 << endl;

	I1 << 1 + 3.0 * l *(l + c1), 1 + 1.5*l*c1, 1 + 1.5*l*c1, 1;
	I1 *= mu1 / 3.0;
	//cout << "I1 matlab" << endl << I1 << endl;

	Im << l * l + r * r + 2 * l*r*c1, r*r + l * r*c1, r*r + l * r*c1, r*r;
	Im *= mum / 3.0;
	I = M * L * L * (I0 + I1 + Im);
	//cout << "I matlab:" << endl << I << endl;

	Matrix2d temp0, temp1;
	temp0 << 1.0, 0.5, 0.5, 0.0;
	temp1 << 2.0, 1.0, 1.0, 0.0;

	dIdtheta1 = -s1 * (mu1*l*temp0 + mum / 3.0*l*r*temp1);
	cout << "dIdtheta matlab" << endl << dIdtheta1 << endl;

	Vector2d fk = -qdot0(1)*dIdtheta1*qdot0;
	cout << "test1" << endl << fk << endl;
	double var = 0.5*qdot0.transpose()*dIdtheta1*qdot0;
	fk(1) += var;
	cout << "test2" << endl << var << endl;

	Vector2d fv0, fv1;
	fv0 << l * s0, 0.0;
	fv0 *= -0.5 * mu0 * g;
	fv1 << 2 * l*s0 + s01, s01;
	fv1 *= -0.5*mu1*g;

	Vector2d fvm;
	fvm << l * s0 + r * s01, r * s01;
	fvm *= -0.5 * mum * g;

	Vector2d f;
	f = M * L * L * fk + M * L * (fv0 + fv1 + fvm);
	
	cout << "fk matlab:" << endl << M * L*L*fk << endl;
	cout << "f matlab:" << endl << f << endl;
	Vector2d qddot = I.ldlt().solve(f);
	VectorXd ydot(4);
	ydot << qdot0, qddot;

	ydot(0) = qdot0(1);
	ydot(1) = qdot0(0);
	ydot(2) = qddot(1);
	ydot(3) = qddot(0);

	return ydot;
}


void SolverSparse::test(const Eigen::VectorXd &x, Eigen::VectorXd &dxdt, const double) {
	//dynamics(x);
	//dxdt = ydotk;
}


