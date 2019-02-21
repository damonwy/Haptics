#include "rmpch.h"
#include "Deformable.h"
#include "Node.h"
#include "Body.h"
#include "Muscle.h"
#include "World.h"
#include "Solver.h"
#include "SolverSparse.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Muscle::Muscle() {

}

Muscle::Muscle(std::vector<std::shared_ptr<Body>> bodies, int n_nodes):
m_bodies(bodies), m_n_nodes(n_nodes)
{
	m_n_bodies = static_cast<int>(bodies.size());
	for (int i = 0; i < n_nodes; i++) {
		auto node = make_shared<Node>();
		node->r = 0.2;
		node->x = Vector3d::Zero();
		node->v = Vector3d::Zero();
		node->a = Vector3d::Zero();
		node->m_J.resize(3, m_n_bodies);
		node->m_Jdot.resize(3, m_n_bodies);
		node->m_Jdot_lhs.resize(3, m_n_bodies);
		Matrix3x2d temp;
		temp.setZero();
		node->m_dJdq.reserve(m_n_bodies);

		// Now it's 2
		for (int i = 0; i < m_n_bodies; ++i) {
			node->m_dJdq.push_back(temp);
		}
			
		m_nodes.push_back(node);
	}
}

void Muscle::load(const std::string & RESOURCE_DIR)
{
	for (int i = 0; i < m_n_nodes; i++) {
		m_nodes[i]->load(RESOURCE_DIR);
	}
}

void Muscle::init() {
    init_();
    if (next != nullptr) {
        next->init();
    }
}

void Muscle::init_() {
	for (int i = 0; i < m_n_nodes; i++) {
		m_nodes[i]->init();
	}

}

void Muscle::update() {
	update_();
	if (next != nullptr) {
		next->update();
	}
}

void Muscle::gatherDofs(Eigen::VectorXd &y, int nr) {
    gatherDofs_(y, nr);
    if (next != nullptr) {
        next->gatherDofs(y, nr);
    }
}

void Muscle::gatherDDofs(Eigen::VectorXd &ydot, int nr) {
    gatherDDofs_(ydot, nr);
    if (next != nullptr) {
        next->gatherDofs(ydot, nr);
    }
}

void Muscle::scatterDofs(Eigen::VectorXd &y, int nr) {
    scatterDofs_(y, nr);
    if (next != nullptr) {
        next->scatterDofs(y, nr);
    }
}

void Muscle::scatterDDofs(Eigen::VectorXd &ydot, int nr) {
    scatterDDofs_(ydot, nr);
    if (next != nullptr) {
        next->scatterDDofs(ydot, nr);
    }
}

void Muscle::computeJacobian(Eigen::MatrixXd &J) { 
	computeJacobian_(J);
	if (next != nullptr) {
		next->computeJacobian(J);
	}
}

void Muscle::computeJacobianSparse(std::vector<T> &J_) { 
	computeJacobianSparse_(J_);
	if (next != nullptr) {
		next->computeJacobianSparse(J_);
	}
}

void Muscle::computeJMJ(Eigen::MatrixXd &JMJ, std::shared_ptr<World> world)
{
	computeJMJ_(JMJ, world);
	if (next != nullptr) {
		next->computeJMJ(JMJ, world);
	}
}

void Muscle::computeJMJSparse(std::vector<T>& JMJ_, std::shared_ptr<World> world)
{
	computeJMJSparse_(JMJ_, world);
	if (next != nullptr) {
		next->computeJMJSparse(JMJ_, world);
	}
}

void Muscle::computeMass(Eigen::MatrixXd &M) { 
	computeMass_(M);
	if (next != nullptr) {
		next->computeMass(M);
	}
}

void Muscle::computeForce(Vector3d grav, Eigen::VectorXd &f) { 
	computeForce_(grav, f);
	if (next != nullptr) {
		next->computeForce(grav, f);
	}
}

void Muscle::computeJMJdotqdot(Eigen::VectorXd & f, const Eigen::VectorXd &qdot, std::shared_ptr<World> world, std::shared_ptr<SolverSparse> solver)
{
	computeJMJdotqdot_(f, qdot, world, solver);
	if (next != nullptr) {
		next->computeJMJdotqdot(f, qdot, world, solver);
	}
}

void Muscle::computeMassSparse(std::vector<T> &M_) { 
    computeMassSparse_(M_);
    if (next != nullptr) {
        next->computeMassSparse(M_);
    }
}

void Muscle::computeForceDamping(Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd &D) { 
    computeForceDamping_(grav, f, D);
    
    if (next != nullptr) {
        next->computeForceDamping(grav, f, D);
    }
}

void Muscle::computeForceDampingSparse(Vector3d grav, Eigen::VectorXd &f, std::vector<T> &D_) { 
    computeForceDampingSparse_(grav, f, D_);
    if (next != nullptr) {
        next->computeForceDampingSparse(grav, f, D_);
    }
}

void Muscle::computeEnergies(Eigen::Vector3d grav, Energy &ener) { 
    computeEnergies_(grav, ener);
    if (next != nullptr) {
        next->computeEnergies(grav, ener);
    }
}

void Muscle::draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const { 
    draw_(MV, prog, progSimple, P);
    if (next != nullptr)
    {
        next->draw(MV, prog, progSimple, P);
    }
}
