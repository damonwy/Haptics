#include "rmpch.h"
#include "Body.h"
#include "Joint.h"
#include "ConstraintPrescBody.h"
#include "ConstraintPrescBodyAttachPoint.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

Body::Body() {
	
}

Body::Body(double density):
m_density(density)
{
	m_isDrawing = true;
	m_damping = 0.0;
	I_i.setIdentity();
	E_ji.setIdentity();
	E_ij.setIdentity();
	E_ie.setIdentity();
	E_wi.setIdentity();
	R_wi.setIdentity();
	E_iw.setIdentity();
	R_iw.setIdentity();
	E_ip.setIdentity();
	Ad_ji.setIdentity();
	Ad_ij.setIdentity();
	Ad_iw.setIdentity();
	Ad_wi.setIdentity();
	Ad_ip.setIdentity();
	Addot_wi.setIdentity();
	phi.setZero();
	phidot.setZero();
	wext_i.setZero();
	fgrav.setZero();
	fcor.setZero();
	m_body_color << 0.8f, 0.7f, 0.7f;
	m_attached_color << static_cast<float>((rand() % 255)/255.0), static_cast<float>((rand() % 255)/255.0), static_cast<float>((rand() % 255)/255.0);
	m_sliding_color << static_cast<float>((rand() % 255) / 255.0), static_cast<float>((rand() % 255) / 255.0), static_cast<float>((rand() % 255) / 255.0);
	presc = nullptr;
}

void Body::load(const string &RESOURCE_DIR, string box_shape) {
	//read a JSON file
	//ifstream i(RESOURCE_DIR + "input.json");
	//json js;
	//i >> js;
	//i.close();
	//string box_shape = js[box_shape];

	// Inits shape
	bodyShape = make_shared<Shape>();
	bodyShape->loadMesh(RESOURCE_DIR + box_shape);
}

void Body::init(int &nm) {
	bodyShape->init();
	countDofs(nm);
}

void Body::setTransform(Eigen::Matrix4d E) {
	// Sets the transform of this body wrt parent joint
	E_ji = E;
	Vector3d p = E.block<3, 1>(0, 3);
	E_ie = E;
	E_ie.block<3, 1>(0, 3) = -p;
	E_ij = SE3::inverse(E_ji);
	Ad_ji = SE3::adjoint(E_ji);
	Ad_ij = SE3::adjoint(E_ij);
}

Vector3d Body::getBodyVelocityByEndPointVelocity(Vector3d v_we) {
	Vector3d v_ew = - v_we;
	Matrix6d Ad_ie = SE3::adjoint(E_ie);
	Matrix3d R_ie = Ad_ie.block<3, 3>(0, 0);
	Vector3d v_iw = R_ie * v_ew;
	return(-v_iw);
}

void Body::update() {
	computeInertia();
	// Updates this body's transforms and velocities
	E_wi = m_joint->E_wj * E_ji;
	E_iw = SE3::inverse(E_wi);
	Ad_wi = SE3::adjoint(E_wi);
	Ad_iw = SE3::adjoint(E_iw);
	E_ip = Matrix4d::Identity();
	
	if (m_joint->getParent() != nullptr) {
		m_parent = m_joint->getParent()->getBody();
		E_ip = E_iw * m_parent->E_wi;
	}
	Ad_ip = SE3::adjoint(E_ip);

	// Body velocity
	phi = Ad_ij * m_joint->V;
	Addot_wi = SE3::dAddt(E_wi, phi);
	phidot = Ad_ij * m_joint->Vdot;
}

void Body::computeEnergies(Eigen::Vector3d grav, Energy &energies) {
	// Compute kinetic and potential energies
	energies.K = energies.K + 0.5 * this->phi.transpose() * I_i.asDiagonal() * this->phi;
	energies.V = energies.V - I_i(5) * grav.transpose() * E_wi.block<3, 1>(0, 3);
}

void Body::computeInertia() {
	// Computes inertia at body and joint
	computeInertia_();
	m_joint->computeInertia();
}

void Body::countDofs(int &nm) {
	// Counts maximal DOFs
	idxM = countM(nm, 6);
}

int Body::countM(int &nm, int data) {
	nm = nm + data;
	return (nm - data);
}

void Body::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) const
{
	draw_(MV, prog, P);
	if (next != nullptr) {
		next->draw(MV, prog, P);
	}
}

void Body::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, shared_ptr<MatrixStack> P) const {
	prog->bind();
	if (bodyShape && m_isDrawing) {
		glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
		glUniform3f(prog->getUniform("lightPos1"), 66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_1"), 0.6f);
		glUniform3f(prog->getUniform("lightPos2"), -66.0f, 25.0f, 25.0f);
		glUniform1f(prog->getUniform("intensity_2"), 0.2f);
		glUniform1f(prog->getUniform("s"), 300.0f);
		glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
		glUniform3fv(prog->getUniform("kd"), 1, this->m_body_color.data());
		glUniform3f(prog->getUniform("ks"), 1.0f, 0.9f, 0.8f);
		MV->pushMatrix();
		MV->multMatrix(eigen_to_glm(E_wi));

		glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
		bodyShape->draw(prog);
		MV->popMatrix();
	}
	prog->unbind();
}

void Body::computeMass(MatrixXd &M) {
	// Computes maximal mass matrix and force vector	
	M.block<6, 6>(idxM, idxM) = M_i;

	if (next != nullptr) {
		next->computeMass(M);
	}
}

void Body::computeMassGrav(Vector3d grav, MatrixXd &M, VectorXd &f) {
	// Computes maximal mass matrix and force vector	
	M.block<6, 6>(idxM, idxM) = M_i;

	fcor = SE3::ad(phi).transpose() * M_i * phi;
	R_wi = E_wi.block<3, 3>(0, 0);
	R_iw = R_wi.transpose();

	fgrav.setZero();
	fgrav.segment<3>(3) = M_i(3, 3) * R_iw * grav; // wrench in body space
	
	f.segment<6>(idxM) = fcor + fgrav;
	// External wrench: used only by recurse (not redmax) to accumulate the wrenches
	// to be applied to the joint in rhdPas2(). For redmax, the array of wrenches is
	// used to add wrenches to the bodies directly.

	this->wext_i.setZero();
	this->Kmdiag.setZero();
	this->Dmdiag.setZero();
	
	// Joint torque
	// This is how we would apply a joint torque using maximal coordinates. 
	// It's much easier in reduced, so we'll do that instead. See Joint::computeForce().

	//if (!m_joint->presc) {
	//	
	//	Vector6d tau = m_joint->m_S * (m_joint->m_tau - m_joint->m_K * m_joint->m_q);
	//	f.segment<6>(idxM) += Ad_ji.transpose() * tau;
	//	// Also apply to parent
	//	if (m_joint->getParent() != nullptr) {
	//		m_parent = m_joint->getParent()->getBody();
	//		int idxM_P = m_parent->idxM;
	//		Matrix4d E_jp = E_ji * E_iw * m_parent->E_wi; // this joint -> parent body
	//		f.segment<6>(idxM_P) -= SE3::adjoint(E_jp).transpose() * tau;
	//	}
	//}

	if (next != nullptr) {
		next->computeMassGrav(grav, M, f);
	}
}

void Body::computeMassSparse(vector<T> &M_) {
	// Computes maximal mass matrix (just once)
	for (int i = 0; i < 6; ++i) {
		M_.push_back(T(idxM + i, idxM + i, I_i(i)));
	}

	if (next != nullptr) {
		next->computeMassSparse(M_);
	}
}

void Body::computeGrav(Vector3d grav, Eigen::VectorXd &f) {
	fcor = SE3::ad(phi).transpose() * M_i * phi;
	//cout << "fcor: " << fcor << endl;

	R_wi = E_wi.block<3, 3>(0, 0);
	R_iw = R_wi.transpose();
	fgrav.setZero();
	fgrav.segment<3>(3) = M_i(3, 3) * R_iw * grav; // wrench in body space
	f.segment<6>(idxM) = fcor + fgrav;
	this->wext_i.setZero();
	this->Kmdiag.setZero();
	this->Dmdiag.setZero();
	if (next != nullptr) {
		next->computeGrav(grav, f);
	}
}


void Body::computeForceDamping(Eigen::VectorXd &f, Eigen::MatrixXd &D) {
	// Computes maximal damping force vector and matrix
	if (m_damping > 0.0) {
		Vector6d fi = -m_damping * phi;
		Matrix6d Di = m_damping * Matrix6d::Identity();
		f.segment<6>(this->idxM) += fi;
		D.block<6, 6>(this->idxM, this->idxM) += Di;
		// Used by recursive algorithm
		this->wext_i += fi;
		this->Dmdiag += Di;
	}

	if (next != nullptr) {
		next->computeForceDamping(f, D);
	}
}

void Body::computeForceDampingSparse(Eigen::VectorXd &f, std::vector<T> &D_) {
	// Computes maximal damping force vector and matrix
	if (m_damping > 0.0) {
		Vector6d fi = -m_damping * phi;
		Matrix6d Di = m_damping * Matrix6d::Identity();
		f.segment<6>(this->idxM) += fi;

		for (int i = 0; i < 6; ++i) {
			D_.push_back(T(idxM + i, idxM + i, m_damping));
		}

		// Used by recursive algorithm
		this->wext_i += fi;
		this->Dmdiag += Di;
	}

	if (next != nullptr) {
		next->computeForceDampingSparse(f, D_);
	}
}