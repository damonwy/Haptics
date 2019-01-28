#include "rmpch.h"
#include "Spring.h"
#include "Deformable.h"
#include "Node.h"
#include "Body.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Spring::Spring(){

}

void Spring::init() {
	init_();
	if (next != nullptr) {
		next->init();
	}
}

void Spring::update() {
	update_();
	if (next != nullptr) {
		next->update();
	}
}

void Spring::computeForceStiffnessDamping(VectorXd &f, MatrixXd &K, MatrixXd &D) {
	computeForceStiffnessDamping_(f, K, D);

	if (next != nullptr) {
		next->computeForceStiffnessDamping(f, K, D);
	}

}

void Spring::computeForceStiffnessDampingSparse(Eigen::VectorXd &f, std::vector<T> &K_, std::vector<T> &D_) {
	computeForceStiffnessDampingSparse_(f, K_, D_);
	if (next != nullptr) {
		next->computeForceStiffnessDampingSparse(f, K_, D_);
	}
}

void Spring::computeStiffnessProd(VectorXd x, VectorXd &y) {
	// Computes y=K*x
	computeStiffnessProd_(x, y);
	if (next != nullptr) {
		next->computeStiffnessProd(x, y);
	}
}

void Spring::computeDampingProd(VectorXd x, VectorXd &y) {
	// Computes y=D*x
	computeDampingProd_(x, y);
	if (next != nullptr) {
		next->computeDampingProd(x, y);
	}
}

void Spring::computeEnergies(Vector3d grav, Energy &ener) {
	computeEnergies_(grav, ener);
	if (next != nullptr) {
		next->computeEnergies(grav, ener);
	}
}

void Spring::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	draw_(MV, prog, progSimple, P);
	if (next != nullptr)
	{
		next->draw(MV, prog, progSimple, P);
	}
}