// Constraint Generic constraint
// A constraint can be applied in reduced or maximal coordinates.
// Reduced: keeping a quaternion to be of unit length.
// Maximal: holding two bodies together in world space.

#pragma once
#ifndef REDUCEDCOORD_SRC_CONSTRAINT_H_
#define REDUCEDCOORD_SRC_CONSTRAINT_H_

#include <vector>
#include <memory>

#define EIGEN_DONT_ALIGN_STATICALLY
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "MLCommon.h"

class Body;

typedef Eigen::Triplet<double> T;


class Constraint
{

public:
	Constraint();
	Constraint(int _nconEM, int _nconER, int _nconIM, int _nconIR);
	virtual ~Constraint() {}

	void computeJacEqM(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);	
	void computeJacEqR(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacIneqM(Eigen::MatrixXd &Cm, Eigen::MatrixXd &Cmdot, Eigen::VectorXd &cm, Eigen::VectorXd &cmdot, Eigen::VectorXd &cmddot);
	void computeJacIneqR(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr, Eigen::VectorXd &crdot, Eigen::VectorXd &crddot);

	void computeJacEqMSparse(std::vector<T> &Gm, std::vector<T> &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot);
	void computeJacEqRSparse(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot);
	void computeJacIneqMSparse(std::vector<T> &Cm, std::vector<T> &Cmdot, Eigen::VectorXd &cm, Eigen::VectorXd &cmdot, Eigen::VectorXd &cmddot);
	void computeJacIneqRSparse(std::vector<T> &Cr, std::vector<T> &Crdot, Eigen::VectorXd &cr, Eigen::VectorXd &crdot, Eigen::VectorXd &crddot);

	void init();
	void countDofs(int &nem, int &ner, int &nim, int &nir);
	virtual void setActive() {}
	virtual void setInactive() {}

	void getActiveList(std::vector<int> &listM, std::vector<int> &listR);

	void getEqActiveList(std::vector<int> &listEqM, std::vector<int> &listEqR);
	void scatterForceEqM(Eigen::MatrixXd Gmt, Eigen::VectorXd lm);
	void scatterForceEqR(Eigen::MatrixXd Grt, Eigen::VectorXd lr);
	void scatterForceIneqR(Eigen::MatrixXd Crt, Eigen::VectorXd lr);
	void scatterForceIneqM(Eigen::MatrixXd Cmt, Eigen::VectorXd lm);
	void ineqEventFcn(std::vector<double> &value, std::vector<int> &isterminal, std::vector<int> &direction);
	void ineqProjPos();

	int nconEM;								// Number of maximal equality constraints
	int nconER;								// Number of reduced equality constraints
	int nconIM;								// Number of maximal inequality constraints
	int nconIR;								// Number of reduced inequality constraints
	int nQ;								
	int idxEM;								// Maximal equality constraint indices
	int idxER;								// Reduced equality constraint indices
	int idxIM;								// Maximal inequality constraint indices
	int idxIR;								// Reduced inequality constraint indices
	Eigen::MatrixXi idxQ;					// Associated DOF indices
	int nidxQ;

	bool activeM;							// Whether the maximal inequality constraint is active
	bool activeR;							// Whether the reduced inequality constraint is active
	bool activeEM;							// Whether the maximal equality constraint is active
	bool activeER;							// Whether the reduced equality constraint is active
	Eigen::VectorXd fcon;					// Computed constraint force
	std::shared_ptr<Constraint> next;		// Next constraint in traversal order
	
protected:
	std::string m_name;
	void scatterForceEqM_() {}
	void scatterForceEqR_() {}
	void scatterForceIneqR_() {}
	void scatterForceIneqM_() {}

	virtual void computeJacEqM_(Eigen::MatrixXd &Gm, Eigen::MatrixXd &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot) {}
	virtual void computeJacEqR_(Eigen::MatrixXd &Gr, Eigen::MatrixXd &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot) {}
	virtual	void computeJacIneqM_(Eigen::MatrixXd &Cm, Eigen::MatrixXd &Cmdot, Eigen::VectorXd &cm, Eigen::VectorXd &cmdot, Eigen::VectorXd &cmddot) {}
	virtual void computeJacIneqR_(Eigen::MatrixXd &Cr, Eigen::MatrixXd &Crdot, Eigen::VectorXd &cr, Eigen::VectorXd &crdot, Eigen::VectorXd &crddot) {}

	virtual void computeJacEqMSparse_(std::vector<T> &Gm, std::vector<T> &Gmdot, Eigen::VectorXd &gm, Eigen::VectorXd &gmdot, Eigen::VectorXd &gmddot) {}
	virtual void computeJacEqRSparse_(std::vector<T> &Gr, std::vector<T> &Grdot, Eigen::VectorXd &gr, Eigen::VectorXd &grdot, Eigen::VectorXd &grddot) {}
	virtual	void computeJacIneqMSparse_(std::vector<T> &Cm, std::vector<T> &Cmdot, Eigen::VectorXd &cm, Eigen::VectorXd &cmdot, Eigen::VectorXd &cmddot) {}
	virtual void computeJacIneqRSparse_(std::vector<T> &Cr, std::vector<T> &Crdot, Eigen::VectorXd &cr, Eigen::VectorXd &crdot, Eigen::VectorXd &crddot) {}

	virtual void ineqEventFcn_(std::vector<double> &value, std::vector<int> &isterminal, std::vector<int> &direction) {}
	virtual void ineqProjPos_() {}
	virtual void init_() {}
};



#endif // REDUCEDCOORD_SRC_CONSTRAINT_H_