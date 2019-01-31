#include "rmpch.h"
#include "Deformable.h"
#include "Node.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Deformable::Deformable():
	m_K(0.0), m_damping(0.0), m_mass(1.0){
}

Deformable::Deformable(int &countS, int &countCM):
	m_K(0.0), m_damping(0.0), m_mass(1.0) {
	countS++;
	countCM++;
}

void Deformable::init() {
	init_();
	if (next != nullptr) {
		next->init();
	}
}

void Deformable::countDofs(int &nm, int &nr) {
	countDofs_(nm, nr);
	if (next != nullptr) {
		next->countDofs(nm, nr);
	}
}

void Deformable::gatherDofs(VectorXd &y, int nr) {
	gatherDofs_(y, nr);
	if (next != nullptr) {
		next->gatherDofs(y, nr);
	}
}

void Deformable::gatherDDofs(VectorXd &ydot, int nr) {
	gatherDDofs_(ydot, nr);
	if (next != nullptr) {
		next->gatherDofs(ydot, nr);
	}
}

void Deformable::scatterDofs(VectorXd &y, int nr) {
	scatterDofs_(y, nr);
	if (next != nullptr) {
		next->scatterDofs(y, nr);
	}
}

void Deformable::scatterDDofs(VectorXd &ydot, int nr) {
	scatterDDofs_(ydot, nr);
	if (next != nullptr) {
		next->scatterDDofs(ydot, nr);
	}	
}

void Deformable::computeJacobian(MatrixXd &J) {
	computeJacobian_(J);
	if (next != nullptr) {
		next->computeJacobian(J);
	}
}

void Deformable::computeJacobianSparse(std::vector<T> &J_){
	computeJacobianSparse_(J_);
	if (next != nullptr) {
		next->computeJacobianSparse(J_);
	}
}

void Deformable::computeMass(MatrixXd &M) {
	computeMass_(M);
	if (next != nullptr) {
		next->computeMass(M);
	}
}

void Deformable::computeMassSparse(vector<T> &M_) {
	computeMassSparse_(M_);
	if (next != nullptr) {
		next->computeMassSparse(M_);
	}
}

void Deformable::computeForce(Vector3d grav, VectorXd &f) {
	computeForce_(grav, f);
	if (next != nullptr) {
		next->computeForce(grav, f);
	}
}

void Deformable::computeForceDamping(Vector3d grav, VectorXd &f, MatrixXd &D) {
	computeForceDamping_(grav, f, D);
	if (next != nullptr) {
		next->computeForceDamping(grav, f, D);
	}
}

void Deformable::computeForceDampingSparse(Vector3d grav, VectorXd &f, vector<T> &D_) {
	computeForceDampingSparse_(grav, f, D_);
	if (next != nullptr) {
		next->computeForceDampingSparse(grav, f, D_);
	}
}

void Deformable::computeEnergies(Vector3d grav, Energy &ener) {
	computeEnergies_(grav, ener);
	if (next != nullptr) {
		next->computeEnergies(grav, ener);
	}
}

void Deformable::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	draw_(MV, prog, progSimple, P);
	if (next != nullptr)
	{
		next->draw(MV, prog, progSimple, P);
	}
}
