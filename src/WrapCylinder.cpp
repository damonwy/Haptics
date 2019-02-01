#include "rmpch.h"
#include "WrapCylinder.h"
#include "Body.h"
#include "Node.h"
#include "Vector.h"
#include "CompCylinder.h"

using namespace std;
using namespace Eigen;

WrapCylinder::WrapCylinder()
{
	m_type = cylinder;
}

WrapCylinder::WrapCylinder(const shared_ptr<Node> &P, const shared_ptr<Node> &S, const shared_ptr<CompCylinder> compCylinder, const int num_points): 
	WrapObst(P, S, num_points), m_compCylinder(compCylinder)
{
	m_type = cylinder;
	m_radius = compCylinder->getRadius();
	m_vec_z = compCylinder->getZAxis();
	m_point_O = compCylinder->getOrigin();
	m_arc_points.resize(3, m_num_points + 1);
}

void WrapCylinder::load(const std::string &RESOURCE_DIR) {
	m_point_O->load(RESOURCE_DIR);
	m_point_P->load(RESOURCE_DIR);
	m_point_S->load(RESOURCE_DIR);
}

void WrapCylinder::init() {
	m_point_O->init();
	m_point_P->init();
	m_point_S->init();
}

void WrapCylinder::compute()
{
	Vector3d OP = m_point_P->x - m_point_O->x;
	OP = OP / OP.norm();
	Vector3d vec_Z = m_vec_z->dir / m_vec_z->dir.norm();
	Vector3d vec_X = vec_Z.cross(OP);
	vec_X = vec_X / vec_X.norm();
	Vector3d vec_Y = vec_Z.cross(vec_X);
	vec_Y = vec_Y / vec_Y.norm();

	this->M << vec_X.transpose(), vec_Y.transpose(), vec_Z.transpose();

	Vector3d p = this->M * (m_point_P->x - m_point_O->x);
	Vector3d s = this->M * (m_point_S->x - m_point_O->x);

	double denom_q = p(0)*p(0) + p(1)*p(1);
	double denom_t = s(0)*s(0) + s(1)*s(1);
	double R = m_radius;

	m_status = wrap;

	if ((denom_q - R*R < 0.0) || (denom_t - R*R < 0.0))
	{
		m_status = inside_radius;
	}

	double root_q = sqrt(denom_q - R*R);
	double root_t = sqrt(denom_t - R*R);

	Vector3d q(0.0, 0.0, 0.0);
	Vector3d t(0.0, 0.0, 0.0);
	q(0) = (p(0) * R*R + R * p(1) * root_q) / denom_q;
	q(1) = (p(1) * R*R - R * p(0) * root_q) / denom_q;
	t(0) = (s(0) * R*R - R * s(1) * root_t) / denom_t;
	t(1) = (s(1) * R*R + R * s(0) * root_t) / denom_t;

	if (R * (q(0) * t(1) - q(1) * t(0)) > 0.0)
	{
		m_status = no_wrap;
	}

	std::complex<double> qt_i = 1.0 - 0.5 *
		((q(0) - t(0)) * (q(0) - t(0))
			+ (q(1) - t(1)) * (q(1) - t(1))) / (R*R);
	double qt_xy = abs(R * acos(qt_i));// changed
	m_path_length = qt_xy;

	double pq_xy = sqrt((p(0) - q(0)) * (p(0) - q(0)) +
		(p(1) - q(1)) * (p(1) - q(1)));
	double ts_xy = sqrt((t(0) - s(0)) * (t(0) - s(0)) +
		(t(1) - s(1)) * (t(1) - s(1)));
	q(2) = p(2) + (s(2) - p(2)) * pq_xy / (pq_xy + qt_xy + ts_xy);
	t(2) = s(2) - (s(2) - p(2)) * ts_xy / (pq_xy + qt_xy + ts_xy);

	m_point_q->x = q;
	m_point_t->x = t;

	Vector3d Q = this->M.transpose() * q + m_point_O->x;
	Vector3d T = this->M.transpose() * t + m_point_O->x;

	// std::cout << Q.transpose() << std::endl << T.transpose() << std::endl;
}

MatrixXd WrapCylinder::getPoints(int num_points) const
{
	double theta_q = atan(m_point_q->x(1) / m_point_q->x(0));
	if (m_point_q->x(0) < 0.0)
		theta_q += M_PI;

	double theta_t = atan(m_point_t->x(1) / m_point_t->x(0));
	if (m_point_t->x(0) < 0.0) {
		theta_t += M_PI;
	}

	MatrixXd points(3, num_points + 1);

	double z_s, z_e;
	//theta_s, theta_e
	double theta_s, theta_e;

	if (theta_q < theta_t)
	{
		theta_s = theta_q; theta_e = theta_t;
		z_s = m_point_q->x(2); z_e = m_point_t->x(2);
	}
	else
	{
		theta_s = theta_t; theta_e = theta_q;
		z_s = m_point_t->x(2); z_e = m_point_q->x(2);
	}

	if (theta_e - theta_s > theta_s + 2 * M_PI - theta_e)
	{
		double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * M_PI;
		tmp = z_s; z_s = z_e; z_e = tmp;
	}

	int col = 0;
	double z_i = z_s, dz = (z_e - z_s) / num_points;
	for (double i = theta_s; i <= theta_e + 0.001;
		i += (theta_e - theta_s) / num_points)
	{
		if (col == num_points + 1) {
			break;
		}

		Eigen::Vector3d point = this->M.transpose() *
			Eigen::Vector3d(m_radius * cos(i), m_radius * sin(i), z_i) +
			m_point_O->x;
		z_i += dz;
		points.col(col++) = point;
	}

	return points;
}

void WrapCylinder::update_() {
	m_point_P->update();
	m_point_S->update();
	m_compCylinder->update();
	compute();
	
	if (m_status == wrap) {
		m_arc_points = getPoints(m_num_points);
	}
}

void WrapCylinder::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const {

	prog->bind();
	
	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniform3f(prog->getUniform("lightPos1"), 1.0f, 1.0f, 1.0f);
	glUniform1f(prog->getUniform("intensity_1"), 0.8f);
	glUniform3f(prog->getUniform("lightPos2"), -1.0f, 1.0f, 1.0f);
	glUniform1f(prog->getUniform("intensity_2"), 0.2f);
	glUniform1f(prog->getUniform("s"), 200.0f);
	glUniform3f(prog->getUniform("ka"), 0.2f, 0.7f, 0.2f);
	glUniform3f(prog->getUniform("kd"), 0.0f, 0.0f, 1.0f);
	glUniform3f(prog->getUniform("ks"), 0.0f, 1.0f, 0.0f);
	MV->pushMatrix();
		
	// Draw P, S points
	this->m_point_P->draw(MV, prog);
	this->m_point_S->draw(MV, prog);
	this->m_point_O->draw(MV, prog);

	prog->unbind();

	// Draw wrapping
	prog2->bind();
	glUniformMatrix4fv(prog2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glColor3f(0.0, 0.0, 0.0);
	glLineWidth(4);
	glBegin(GL_LINE_STRIP);
	glVertex3f(float(m_point_S->x(0)), float(m_point_S->x(1)), float(m_point_S->x(2)));

	if (m_status == wrap) {
		for (int i = 0; i < m_arc_points.cols(); i++) {
			Vector3f p = m_arc_points.block<3, 1>(0, i).cast<float>();
			glVertex3f(p(0), p(1), p(2));
		}
	}

	glVertex3f(static_cast<GLfloat>(m_point_P->x(0)), static_cast<GLfloat>(m_point_P->x(1)), static_cast<GLfloat>(m_point_P->x(2)));
	glEnd();

	// Draw z axis
	m_vec_z->draw(MV, P, prog2);

	MV->popMatrix();
	prog2->unbind();

}