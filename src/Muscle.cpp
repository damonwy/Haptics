#include "rmpch.h"
#include "Deformable.h"
#include "Node.h"
#include "Body.h"
#include "Muscle.h"

Muscle::Muscle() {

}

Muscle::Muscle(std::vector<std::shared_ptr<Body>> bodies):
m_bodies(bodies)
{
    
}

void Muscle::init() {
    init_();
    if (next != nullptr) {
        next->init();
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
    
}

void Muscle::computeJacobianSparse(std::vector<T> &J_) { 
    
}

void Muscle::computeMass(Eigen::MatrixXd &M) { 
   
}

void Muscle::computeForce(Vector3d grav, Eigen::VectorXd &f) { 
   
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
