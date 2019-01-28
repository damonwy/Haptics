#include "rmpch.h"
#include "Constraint.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Constraint::Constraint() {

}

Constraint::Constraint(int _nconEM, int _nconER, int _nconIM, int _nconIR) :
nconEM(_nconEM),
nconER(_nconER),
nconIM(_nconIM),
nconIR(_nconIR),
activeM(false),
activeR(false),
activeEM(true),
activeER(true)
{

}

void Constraint::init() {
	init_();
	if (next != nullptr) {
		next->init();
	}
}

void Constraint::countDofs(int &nem, int &ner, int &nim, int &nir) {
	// Counts DOFs
	idxEM = nem;
	idxER = ner;
	idxIM = nim;
	idxIR = nir;

	nem += nconEM;
	ner += nconER;
	nim += nconIM;
	nir += nconIR;
}

void Constraint::getActiveList(std::vector<int> &listM, std::vector<int> &listR) {
	// Gets list of active inequality indices
	if (activeM) {
		listM.push_back(idxIM);
	}
	if (activeR) {
		listR.push_back(idxIR);
	}
	if (next != nullptr) {
		next->getActiveList(listM, listR);
	}
}

void Constraint::getEqActiveList(vector<int> &listEqM, vector<int> &listEqR) {
	// Gets list of active equality indices
	if (activeEM) {
		for (int i = 0; i < nconEM; i++) {
			listEqM.push_back(idxEM + i);
		}
	}
	if (activeER) {
		for (int i = 0; i < nconER; i++) {
			listEqR.push_back(idxER + i);
		}
	}
	if (next != nullptr) {
		next->getEqActiveList(listEqM, listEqR);
	}

}

void Constraint::computeJacEqM(MatrixXd &Gm, MatrixXd &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {

	computeJacEqM_(Gm, Gmdot, gm, gmdot, gmddot);
	if (next != nullptr) {
		next->computeJacEqM(Gm, Gmdot, gm, gmdot, gmddot);
	}
}

void Constraint::computeJacEqR(MatrixXd &Gr, MatrixXd &Grdot, VectorXd &gr, VectorXd &grdot, VectorXd &grddot) {

	computeJacEqR_(Gr, Grdot, gr, grdot, grddot);
	if (next != nullptr) {
		next->computeJacEqR(Gr, Grdot, gr, grdot, grddot);
	}
}

void Constraint::computeJacIneqM(MatrixXd &Cm, MatrixXd &Cmdot, VectorXd &cm, VectorXd &cmdot, VectorXd &cmddot) {
	computeJacIneqM_(Cm, Cmdot, cm, cmdot, cmddot);
	if (next != nullptr) {
		next->computeJacIneqM(Cm, Cmdot, cm, cmdot, cmddot);
	}
}

void Constraint::computeJacIneqR(MatrixXd &Cr, MatrixXd &Crdot, VectorXd &cr, VectorXd &crdot, VectorXd &crddot) {
	computeJacIneqR_(Cr, Crdot, cr, crdot, crddot);
	if (next != nullptr) {
		next->computeJacIneqR(Cr, Crdot, cr, crdot, crddot);
	}
}

void Constraint::computeJacEqMSparse(vector<T> &Gm, vector<T> &Gmdot, VectorXd &gm, VectorXd &gmdot, VectorXd &gmddot) {
	computeJacEqMSparse_(Gm, Gmdot, gm, gmdot, gmddot);
	if (next != nullptr) {
		next->computeJacEqMSparse(Gm, Gmdot, gm, gmdot, gmddot);
	}
}
void Constraint::computeJacEqRSparse(vector<T> &Gr, vector<T> &Grdot, VectorXd &gr, VectorXd &grdot, VectorXd &grddot) {
	computeJacEqRSparse_(Gr, Grdot, gr, grdot, grddot);
	if (next != nullptr) {
		next->computeJacEqRSparse(Gr, Grdot, gr, grdot, grddot);
	}
}

void Constraint::computeJacIneqMSparse(vector<T> &Cm, vector<T> &Cmdot, VectorXd &cm, VectorXd &cmdot, VectorXd &cmddot) {
	computeJacIneqMSparse_(Cm, Cmdot, cm, cmdot, cmddot);
	if (next != nullptr) {
		next->computeJacIneqMSparse(Cm, Cmdot, cm, cmdot, cmddot);
	}
}
void Constraint::computeJacIneqRSparse(vector<T> &Cr, vector<T> &Crdot, VectorXd &cr, VectorXd &crdot, VectorXd &crddot) {
	computeJacIneqRSparse_(Cr, Crdot, cr, crdot, crddot);
	if (next != nullptr) {
		next->computeJacIneqRSparse(Cr, Crdot, cr, crdot, crddot);
	}
}


void Constraint::scatterForceEqM(Eigen::MatrixXd Gmt, Eigen::VectorXd lm) {
	if (nconEM > 0) {
		int rows = idxQ.cols() * idxQ.rows();
		MatrixXd temp(rows, nconEM);
		temp.setZero();
		for (int i = 0; i < idxQ.cols(); i++) {
			temp.block(6 * i, 0, 6, nconEM) = Gmt.block(idxQ(0, i), idxEM, 6, nconEM);
		}
		fcon.resize(idxQ.cols() * idxQ.rows());
		fcon = -temp * lm.segment(idxEM, nconEM);
		//fcon = -Gmt.block(idxQ, idxEM, nQ, nconEM) * lm.segment(idxEM, nconEM);
	}
	else {
		fcon.resize(idxQ.cols() * idxQ.rows());
		fcon.setZero();
	}
	scatterForceEqM_();
	if (next != nullptr) {
		next->scatterForceEqM(Gmt, lm);
	}
}

void Constraint::scatterForceEqR(Eigen::MatrixXd Grt, Eigen::VectorXd lr) {
	if (nconER > 0) {
		int rows = idxQ.cols() * idxQ.rows();
		MatrixXd temp(rows, nconER);
		temp.setZero();
		for (int i = 0; i < idxQ.cols(); i++) {
			temp.block(6 * i, 0, 6, nconER) = Grt.block(idxQ(0, i), idxER, 6, nconER);
		}
		fcon = -temp * lr.segment(idxER, nconER);
		//fcon = -Grt.block(idxQ, idxER, nQ, nconER) * lr.segment(idxER, nconER);
	}
	else {
		fcon.resize(idxQ.cols() * idxQ.rows());
		fcon.setZero();
	}
	scatterForceEqR_();
	if (next != nullptr) {
		next->scatterForceEqR(Grt, lr);
	}
}

void Constraint::scatterForceIneqR(Eigen::MatrixXd Crt, Eigen::VectorXd lr) {
	if (nconIR > 0) {
		fcon = -Crt.block(idxQ(0), idxIR, idxQ.rows(), nconIR) * lr.segment(idxIR, nconIR);
	}
	else {
		fcon.resize(idxQ.rows());
		fcon.setZero();
	}
	scatterForceIneqR_();
	if (next != nullptr) {
		next->scatterForceIneqR(Crt, lr);
	}
}

void Constraint::scatterForceIneqM(Eigen::MatrixXd Cmt, Eigen::VectorXd lm) {
	if (nconIM > 0) {
		fcon = -Cmt.block(idxQ(0), idxIM, idxQ.rows(), nconIM) * lm.segment(idxEM, nconIM);
	}
	else {
		fcon.resize(idxQ.rows());
		fcon.setZero();
	}
	scatterForceIneqM_();
	if (next != nullptr) {
		next->scatterForceIneqM(Cmt, lm);
	}
}

void Constraint::ineqEventFcn(vector<double> &value, vector<int> &isterminal, vector<int> &direction) {
	ineqEventFcn_(value, isterminal, direction);
	if (next != nullptr) {
		next->ineqEventFcn(value, isterminal, direction);
	}
}

void Constraint::ineqProjPos() {
	ineqProjPos_();

	if (next != nullptr) {
		next->ineqProjPos();
	}
}