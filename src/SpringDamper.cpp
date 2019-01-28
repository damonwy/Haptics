#include "rmpch.h"
#include "SpringDamper.h"

#include "Body.h"
#include "Node.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;

SpringDamper::SpringDamper() {

}


SpringDamper::SpringDamper(shared_ptr<Body> body0, Vector3d r0, shared_ptr<Body> body1, Vector3d r1):
Spring(),
m_body0(body0), m_body1(body1), 
m_r0(r0), m_r1(r1),
m_K(1.0), m_L(0.0), m_damping(1.0)
{
	for (int i = 0; i < 2; i++) {
		auto node = make_shared<Node>();
		node->r = 0.2;
		node->x = Vector3d::Zero();
		node->v = Vector3d::Zero();
		node->a = Vector3d::Zero();
		
		m_nodes.push_back(node);
	}

	m_nodes[0]->setParent(body0);
	m_nodes[1]->setParent(body1);
}

void SpringDamper::init_() {

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->init();
	}

	Matrix4d E0, E1;
	if (m_body0 == nullptr) {
		E0 = Matrix4d::Identity();
	}
	else {
		E0 = m_body0->E_wi;
	}

	if (m_body1 == nullptr) {
		E1 = Matrix4d::Identity();
	}
	else {
		E1 = m_body1->E_wi;
	}
	Vector4d r0, r1;
	r0.segment<3>(0) = m_r0;
	r1.segment<3>(0) = m_r1;
	r0(3) = 1.0;
	r1(3) = 1.0;
	Vector4d x0 = E0 * r0;
	Vector4d x1 = E1 * r1;

	// Set the nodal positions
	m_nodes[0]->x0 = x0.segment<3>(0);
	m_nodes[1]->x0 = x1.segment<3>(0);
	
}

void SpringDamper::load(const string &RESOURCE_DIR) {

	for (int i = 0; i < (int)m_nodes.size(); i++) {
		m_nodes[i]->load(RESOURCE_DIR);
	}

}

void SpringDamper::update_() {
	
	m_nodes[0]->update();
	m_nodes[1]->update();
	
}

void SpringDamper::computeEnergies_(Vector3d grav, Energy &ener) {
	Matrix4d E0, E1;
	E0.setIdentity();
	E1.setIdentity();

	if (m_body0 != nullptr) {
		E0 = m_body0->E_wi;
	}

	if (m_body1 != nullptr) {
		E1 = m_body1->E_wi;
	}

	Vector4d temp0, temp1;
	temp0 << m_r0, 1.0;
	temp1 << m_r1, 1.0;
	Vector3d x0_w = E0.block<3, 4>(0, 0)*temp0;
	Vector3d x1_w = E1.block<3, 4>(0, 0)*temp1;

	m_l = (x1_w - x0_w).norm();
	if (m_L == 0.0) {
		m_L = m_l;
	}

	double e = (m_l - m_L) / m_L;
	ener.V = ener.V + 0.5 * m_K * e * e;
	
}

void SpringDamper::computeStiffnessProd_(VectorXd x, VectorXd &y) {
	Vector12d f_;
	Matrix12d K_, D_;
	computeFKD(f_, K_, D_);
	int idxM0, idxM1;
	if (m_body0 != nullptr) {
		idxM0 = m_body0->idxM;
		Matrix6d K00 = K_.block<6, 6>(0, 0);
		y.segment<6>(idxM0) += K00 * x.segment<6>(idxM0);
	}

	if (m_body1 != nullptr) {
		idxM1 = m_body1->idxM;
		Matrix6d K11 = K_.block<6, 6>(6, 6);
		y.segment<6>(idxM1) += K11 * x.segment<6>(idxM1);
	}

	if (m_body0 != nullptr && m_body1 != nullptr) {
		Matrix6d K01 = K_.block<6, 6>(0, 6);
		Matrix6d K10 = K_.block<6, 6>(6, 0);

		y.segment<6>(idxM0) += K01 * x.segment<6>(idxM1);
		y.segment<6>(idxM1) += K10 * x.segment<6>(idxM0);
		
	}
}

void SpringDamper::computeDampingProd_(VectorXd x, VectorXd &y) {
	Vector12d f_;
	Matrix12d K_, D_;
	computeFKD(f_, K_, D_);
	int idxM0, idxM1;
	if (m_body0 != nullptr) {
		idxM0 = m_body0->idxM;
		Matrix6d D00 = D_.block<6, 6>(0, 0);
		y.segment<6>(idxM0) += D00 * x.segment<6>(idxM0);
	}

	if (m_body1 != nullptr) {
		idxM1 = m_body1->idxM;
		Matrix6d D11 = D_.block<6, 6>(6, 6);
		y.segment<6>(idxM1) += D11 * x.segment<6>(idxM1);
	}

	if (m_body0 != nullptr && m_body1 != nullptr) {
		Matrix6d D01 = D_.block<6, 6>(0, 6);
		Matrix6d D10 = D_.block<6, 6>(6, 0);

		y.segment<6>(idxM0) += D01 * x.segment<6>(idxM1);
		y.segment<6>(idxM1) += D10 * x.segment<6>(idxM0);

	}
}

void SpringDamper::computeForceStiffnessDamping_(VectorXd &f, MatrixXd &K, MatrixXd &D) {
	Vector12d f_;
	Matrix12d K_, D_;
	computeFKD(f_, K_, D_);

	int idxM0, idxM1;

	if (m_body0 != nullptr) {
		idxM0 = m_body0->idxM;
		Vector6d f0 = f_.segment<6>(0);
		Matrix6d K00 = K_.block<6, 6>(0, 0);
		Matrix6d D00 = D_.block<6, 6>(0, 0);
		f.segment<6>(idxM0) += f0;
		K.block<6, 6>(idxM0, idxM0) += K00;
		D.block<6, 6>(idxM0, idxM0) += D00;
		// TODO RECURSIVE DYM
	}

	if (m_body1 != nullptr) {
		idxM1 = m_body1->idxM;
		Vector6d f1 = f_.segment<6>(6);
		Matrix6d K11 = K_.block<6, 6>(6, 6);
		Matrix6d D11 = D_.block<6, 6>(6, 6);
		f.segment<6>(idxM1) += f1;
		K.block<6, 6>(idxM1, idxM1) += K11;
		D.block<6, 6>(idxM1, idxM1) += D11;
	}

	if (m_body0 != nullptr && m_body1 != nullptr) {
		Matrix6d K01 = K_.block<6, 6>(0, 6);
		Matrix6d K10 = K_.block<6, 6>(6, 0);

		Matrix6d D01 = D_.block<6, 6>(0, 6);
		Matrix6d D10 = D_.block<6, 6>(6, 0);

		K.block<6, 6>(idxM0, idxM1) += K01;
		K.block<6, 6>(idxM1, idxM0) += K10;

		D.block<6, 6>(idxM0, idxM1) += D01;
		D.block<6, 6>(idxM1, idxM0) += D10;
	}

}

void SpringDamper::computeForceStiffnessDampingSparse_(VectorXd &f, std::vector<T> &K, std::vector<T> &D) {
	Vector12d f_;
	Matrix12d K_, D_;
	computeFKD(f_, K_, D_);

	int idxM0, idxM1;

	if (m_body0 != nullptr) {
		idxM0 = m_body0->idxM;
		Vector6d f0 = f_.segment<6>(0);
		Matrix6d K00 = K_.block<6, 6>(0, 0);
		Matrix6d D00 = D_.block<6, 6>(0, 0);
		f.segment<6>(idxM0) += f0;

		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				K.push_back(T(idxM0 + i, idxM0 + j, K00(i, j)));
				D.push_back(T(idxM0 + i, idxM0 + j, D00(i, j)));
			}
		}
	}

	if (m_body1 != nullptr) {
		idxM1 = m_body1->idxM;
		Vector6d f1 = f_.segment<6>(6);
		Matrix6d K11 = K_.block<6, 6>(6, 6);
		Matrix6d D11 = D_.block<6, 6>(6, 6);
		f.segment<6>(idxM1) += f1;

		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				K.push_back(T(idxM1 + i, idxM1 + j, K11(i, j)));
				D.push_back(T(idxM1 + i, idxM1 + j, D11(i, j)));
			}
		}
	}

	if (m_body0 != nullptr && m_body1 != nullptr) {
		Matrix6d K01 = K_.block<6, 6>(0, 6);
		Matrix6d K10 = K_.block<6, 6>(6, 0);

		Matrix6d D01 = D_.block<6, 6>(0, 6);
		Matrix6d D10 = D_.block<6, 6>(6, 0);

		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				K.push_back(T(idxM0 + i, idxM1 + j, K01(i, j)));
				D.push_back(T(idxM0 + i, idxM1 + j, D01(i, j)));
				K.push_back(T(idxM1 + i, idxM0 + j, K10(i, j)));
				D.push_back(T(idxM1 + i, idxM0 + j, D10(i, j)));
			}
		}
	}
}

void SpringDamper::computeFKD(Vector12d &f, Matrix12d &K, Matrix12d &D) {
	Matrix4d E0, E1;
	E0.setIdentity();
	E1.setIdentity();

	Vector6d phi0, phi1;
	phi0.setZero();
	phi1.setZero();

	if (m_body0 != nullptr) {
		E0 = m_body0->E_wi;
		phi0 = m_body0->phi;
	}

	if (m_body1 != nullptr) {
		E1 = m_body1->E_wi;
		phi1 = m_body1->phi;
	}

	Vector4d temp0, temp1;
	temp0 << m_r0, 1.0;
	temp1 << m_r1, 1.0;

	Vector3d x0_w, x1_w, dx_w;
	x0_w = E0.block<3, 4>(0, 0)*temp0;
	x1_w = E1.block<3, 4>(0, 0)*temp1;
	dx_w = x1_w - x0_w;
	m_l = dx_w.norm();

	if (m_L == 0.0) {
		m_L = m_l;
	}

	Matrix3d R0, R1;
	R0 = E0.block<3, 3>(0, 0);
	R1 = E1.block<3, 3>(0, 0);
	
	Matrix3x6d G0, G1;
	G0 = SE3::gamma(m_r0);
	G1 = SE3::gamma(m_r1);
	

	Vector3d v0_w, v1_w;
	v0_w.noalias() = R0 * G0 * phi0;
	v1_w.noalias() = R1 * G1 * phi1;

	double v = (dx_w / m_l).transpose() * (v1_w - v0_w);

	Vector6d fx_0, fx_1;
	fx_0.noalias() = -G0.transpose() * R0.transpose() * dx_w;
	fx_1.noalias() = G1.transpose() * R1.transpose() * dx_w;

	Vector12d fx, fn;
	fx << fx_0, fx_1;

	fn = (1.0 / m_l) *fx;
	double fs = m_K * (m_l - m_L) / m_L - m_damping * v;

	f = -fs * fn;

	// Kn0
	Vector3d ddxinvdx0, ddxinvdx1;
	ddxinvdx0 = dx_w / (m_l * m_l * m_l);
	ddxinvdx1 = -ddxinvdx0;

	Vector6d ddxinvdE0, ddxinvdE1;
	ddxinvdE0.noalias() = ddxinvdx0.transpose() * R0 * G0;
	ddxinvdE1.noalias() = ddxinvdx1.transpose() * R1 * G1;
	Vector12d ddxinvdE;
	ddxinvdE << ddxinvdE0, ddxinvdE1;

	Matrix12d Kn0, Kn1, Kn;
	Kn0 = fx * ddxinvdE.transpose();
	//cout << Kn0 << endl;
	// Kn1
	Kn1.setZero();
	Vector3d p0, p1;
	p0 = E0.block<3, 1>(0, 3);
	p1 = E1.block<3, 1>(0, 3);

	Matrix3d x0b, x1b;
	x0b = SE3::bracket3(m_r0);
	x1b = SE3::bracket3(m_r1);

	Matrix3d I3;
	I3.setIdentity();

	Matrix3d R1R0, R0R1;
	R1R0.noalias() = R1.transpose() * R0;
	R0R1 = R1R0.transpose();

	Kn1.block<3, 3>(3, 0).noalias() = SE3::bracket3(R0.transpose() * (p0 - x1_w));
	Kn1.block<3, 3>(0, 0).noalias() = x0b * Kn1.block<3, 3>(3, 0);
	Kn1.block<3, 3>(9, 0).noalias() = R1R0 * x0b;//
	Kn1.block<3, 3>(6, 0).noalias() = x1b * Kn1.block<3, 3>(9, 0);//
	Kn1.block<3, 3>(3, 3) = I3;
	Kn1.block<3, 3>(0, 3) = x0b;//
	Kn1.block<3, 3>(9, 3) = -R1R0;
	Kn1.block<3, 3>(6, 3).noalias() = x1b * Kn1.block<3, 3>(9, 3);
	Kn1.block<3, 3>(3, 6).noalias() = R0R1 * x1b;
	Kn1.block<3, 3>(0, 6).noalias() = x0b * Kn1.block<3, 3>(3, 6);
	Kn1.block<3, 3>(9, 6).noalias() = SE3::bracket3(R1.transpose() *(p1 - x0_w));
	Kn1.block<3, 3>(6, 6).noalias() = x1b * Kn1.block<3, 3>(9, 6);
	Kn1.block<3, 3>(3, 9) = -R0R1;
	Kn1.block<3, 3>(0, 9).noalias() = x0b * Kn1.block<3, 3>(3, 9);
	Kn1.block<3, 3>(9, 9) = I3;
	Kn1.block<3, 3>(6, 9) = x1b;
	Kn1 /= m_l;
	// Stiffness term for vector part
	Kn = Kn0 + Kn1;

	// Stiffness scalar part
	Vector3d dfsdx0 = - m_K / m_L * dx_w/ m_l;
	Vector6d dfsdE0 = dfsdx0.transpose() * R0 * G0;
	Vector6d dfsdE1 = -dfsdx0.transpose() * R1 * G1;

	Vector12d dfsdE;
	dfsdE << dfsdE0, dfsdE1;
	K.noalias() = fn * dfsdE.transpose() + fs * Kn;

	Matrix12d K_sym = -0.5 * (K + K.transpose());

	K = K_sym; // symmetrize
	
	// Damping scalar part
	Vector3d dir_w = dx_w / m_l;
	Vector6d dfmdphi0, dfmdphi1;
	Vector3d dfmdv0 = m_damping * dir_w;
	dfmdphi0.noalias() = dfmdv0.transpose() * R0 * G0;
	dfmdphi1.noalias() = -dfmdv0.transpose() * R1 * G1;
	Vector12d dfmdphi;
	dfmdphi << dfmdphi0, dfmdphi1;
	D.noalias() = -fn * dfmdphi.transpose();

}

void SpringDamper::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, shared_ptr<MatrixStack> P) const {
	// Draw nodes
	prog->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	int n_nodes = (int)m_nodes.size();
	for (int i = 0; i < n_nodes; i++) {
		m_nodes[i]->draw(MV, prog);
	}

	MV->popMatrix();
	prog->unbind();

	// Draw line segments

	progSimple->bind();
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));

	glLineWidth(5);
	glBegin(GL_LINES);
	glColor3f(0.0f, 0.0f, 0.0f);

	for (int i = 0; i < n_nodes - 1; i++) {

		Vector3f x0 = m_nodes[i]->x.cast<float>();
		Vector3f x1 = m_nodes[i + 1]->x.cast<float>();

		glVertex3f(x0(0), x0(1), x0(2));
		glVertex3f(x1(0), x1(1), x1(2));
	}

	glEnd();
	progSimple->unbind();
	
}


