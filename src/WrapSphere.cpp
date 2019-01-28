#include "rmpch.h"
#include "WrapSphere.h"
#include "CompSphere.h"

WrapSphere::WrapSphere()
{
	m_type = sphere;
}

WrapSphere::WrapSphere(const std::shared_ptr<Node> &P,
	const std::shared_ptr<Node> &S,
	const std::shared_ptr<CompSphere> &compSphere, int num_points)
	: WrapObst(P, S, num_points)
{
	m_type = sphere;
	m_compSphere = compSphere;
	m_point_O = compSphere->getOrigin();
	m_radius = compSphere->getRadius();
}

void WrapSphere::load(const std::string &RESOURCE_DIR) {
	m_point_O->load(RESOURCE_DIR);
	m_point_P->load(RESOURCE_DIR);
	m_point_S->load(RESOURCE_DIR);
}

void WrapSphere::init() {
	m_point_O->init();
	m_point_P->init();
	m_point_S->init();
}

void WrapSphere::compute()
{
	Eigen::Vector3d OS = m_point_S->x - m_point_O->x;
	OS = OS / OS.norm();
	Eigen::Vector3d OP = m_point_P->x - m_point_O->x;
	OP = OP / OP.norm();
	Eigen::Vector3d N = OP.cross(OS);
	N = N / N.norm();

	if (N.dot(Eigen::Vector3d(0.0, 0.0, 1.0)) < 0) {
		N = -N;
	}

	this->M << OS.transpose(), N.cross(OS).transpose(), N.transpose();
	//std::cout << this->M << std::endl;

	Eigen::Vector3d p = this->M * (m_point_P->x - m_point_O->x);
	Eigen::Vector3d s = this->M * (m_point_S->x - m_point_O->x);

	double denom_q = p(0)*p(0) + p(1)*p(1);
	double denom_t = s(0)*s(0) + s(1)*s(1);
	double R = this->m_radius;

	this->m_status = wrap;

	if ((denom_q - R*R < 0.0) || (denom_t - R*R < 0.0))
	{
		this->m_status = inside_radius;
	}

	double root_q = sqrt(denom_q - R*R);
	double root_t = sqrt(denom_t - R*R);

	Eigen::Vector3d q(0.0, 0.0, 0.0);
	Eigen::Vector3d t(0.0, 0.0, 0.0);
	q(0) = (p(0) * R*R + R * p(1) * root_q) / denom_q;
	q(1) = (p(1) * R*R - R * p(0) * root_q) / denom_q;
	t(0) = (s(0) * R*R - R * s(1) * root_t) / denom_t;
	t(1) = (s(1) * R*R + R * s(0) * root_t) / denom_t;

	if (R * (q(0) * t(1) - q(1) * t(0)) > 0.0)
	{
		this->m_status = no_wrap;
	}

	m_point_q->x = q;
	m_point_t->x = t;

	//std::cout << q << std::endl << t << std::endl;

	Eigen::Vector3d Q = this->M.transpose() * q + m_point_O->x;
	Eigen::Vector3d T = this->M.transpose() * t + m_point_O->x;

	//  std::cout << Q.transpose() << std::endl << T.transpose() << std::endl;

	m_path_length = R * acos(1.0 - 0.5 *
		((Q(0) - T(0)) * (Q(0) - T(0))
			+ (Q(1) - T(1)) * (Q(1) - T(1))) / (R*R));
}

Eigen::MatrixXd WrapSphere::getPoints(int num_points)
{
	double theta_q = atan(m_point_q->x(1) / m_point_q->x(0));
	if (m_point_q->x(0) < 0.0) {
		theta_q += M_PI;
	}
		

	double theta_t = atan(m_point_t->x(1) / m_point_t->x(0));
	if (m_point_t->x(0) < 0.0) {
		theta_t += M_PI;
	}
		
	Eigen::MatrixXd points(3, num_points + 1);

	double theta_s, theta_e;

	if (theta_q < theta_t)
	{
		theta_s = theta_q; theta_e = theta_t;
	}
	else
	{
		theta_s = theta_t; theta_e = theta_q;
	}

	if (theta_e - theta_s > theta_s + 2 * M_PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * M_PI;
	}

	int col = 0;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{
		if (col == num_points + 1) {
			break;
		}

		Eigen::Vector3d point = this->m_radius * this->M.transpose() *
			Eigen::Vector3d(cos(i), sin(i), 0.0) + m_point_O->x;
		points.col(col++) = point;
	}

	return points;
}

void WrapSphere::update_() {
	m_point_P->update();
	m_point_S->update();
	m_compSphere->update();

	compute();

	if (m_status == wrap) {
		m_arc_points = getPoints(m_num_points);
	}
}

void WrapSphere::draw_(std::shared_ptr<MatrixStack> MV, const std::shared_ptr<Program> prog, const std::shared_ptr<Program> prog2, std::shared_ptr<MatrixStack> P) const {

	prog->bind();

	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniform3f(prog->getUniform("lightPos1"), 1.0f, 1.0f, 1.0f);
	glUniform1f(prog->getUniform("intensity_1"), 0.8f);
	glUniform3f(prog->getUniform("lightPos2"), -1.0f, 1.0f, 1.0f);
	glUniform1f(prog->getUniform("intensity_2"), 0.2f);
	glUniform1f(prog->getUniform("s"), 200.0f);
	glUniform3f(prog->getUniform("ka"), 0.2f, 0.2f, 0.2f);
	glUniform3f(prog->getUniform("kd"), 0.0f, 0.0f, 1.0f);
	glUniform3f(prog->getUniform("ks"), 0.0f, 1.0f, 0.0f);
	MV->pushMatrix();

	// Draw P, S points
	this->m_point_P->draw(MV, prog);
	this->m_point_S->draw(MV, prog);
	this->m_point_O->draw(MV, prog);

	MV->popMatrix();
	prog->unbind();

	// Draw wrapping
	prog2->bind();
	glUniformMatrix4fv(prog2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glColor3f(0.0f, 0.0f, 0.0f);
	glLineWidth(4);
	glBegin(GL_LINE_STRIP);
	glVertex3f(float(m_point_S->x(0)), float(m_point_S->x(1)), float(m_point_S->x(2)));

	if (m_status == wrap) {
		for (int i = 0; i < m_arc_points.cols(); i++) {
			Eigen::Vector3f p = m_arc_points.block<3, 1>(0, i).cast<float>();
			glVertex3f(p(0), p(1), p(2));
		}
	}

	glVertex3f(m_point_P->x(0), m_point_P->x(1), m_point_P->x(2));
	glEnd();
	MV->popMatrix();
	prog2->unbind();
}