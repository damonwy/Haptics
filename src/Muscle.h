#ifndef HAPTICS_SRC_MUSCLE_H_
#define HAPTICS_SRC_MUSCLE_H_

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

class Muscle {
    
public:
    Muscle();
    Muscle(std::vector<std::shared_ptr<Body>>);
    virtual ~Muscle() {}
    
    void setDamping(double damping) { m_damping = damping; }
    
    virtual void load(const std::string &RESOURCE_DIR) {}
    void init();
    void gatherDofs(Eigen::VectorXd &y, int nr);
    void gatherDDofs(Eigen::VectorXd &ydot, int nr);
    void scatterDofs(Eigen::VectorXd &y, int nr);
    void scatterDDofs(Eigen::VectorXd &ydot, int nr);
    
    void computeJacobian(Eigen::MatrixXd &J);
    void computeJacobianSparse(std::vector<T> &J_);
    void computeMass(Eigen::MatrixXd &M);
    void computeForce(Vector3d grav, Eigen::VectorXd &f);
    void computeMassSparse(std::vector<T> &M_);
    
    void computeForceDamping(Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd &D);
    void computeForceDampingSparse(Vector3d grav, Eigen::VectorXd &f, std::vector<T> &D_);
    
    void computeEnergies(Eigen::Vector3d grav, Energy &ener);
    
    void draw(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const;
    
    virtual void init_() {}
    virtual void draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> progSimple, std::shared_ptr<MatrixStack> P) const {}
    
    virtual void gatherDofs_(Eigen::VectorXd &y, int nr) {}
    virtual void gatherDDofs_(Eigen::VectorXd &ydot, int nr) {}
    virtual void scatterDofs_(Eigen::VectorXd &y, int nr) {}
    virtual void scatterDDofs_(Eigen::VectorXd &ydot, int nr) {}
    
    virtual void computeMass_(Eigen::MatrixXd &M) {}
    virtual void computeMassSparse_(std::vector<T> &M_) {}
    virtual void computeForce_(Vector3d grav, Eigen::VectorXd &f) {}
    virtual void computeForceDamping_(Vector3d grav, Eigen::VectorXd &f, Eigen::MatrixXd &D) {}
    virtual void computeForceDampingSparse_(Vector3d grav, Eigen::VectorXd &f, std::vector<T> &D_) {}
    virtual void computeEnergies_(Vector3d grav, Energy &ener) {}
    virtual void computeJacobian_(Eigen::MatrixXd &J) {}
    virtual void computeJacobianSparse_(std::vector<T> &J_) {}
    
    std::shared_ptr<Muscle> next;
    std::string m_name;
    int m_uid;
    
    std::vector<std::shared_ptr<Node>> m_nodes;
    double m_K;
    double m_damping;
    std::vector<std::shared_ptr<Body>> m_bodies;
    Eigen::Vector3d m_r0;
    Eigen::Vector3d m_r1;
    
    double m_mass;
};

#endif // HAPTICS_SRC_MUSCLE_H_
