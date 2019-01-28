#include "rmpch.h"
#include "WrapDoubleCylinder.h"
#include "Body.h"
#include "Node.h"
#include "Vector.h"
#include "CompDoubleCylinder.h"

using namespace std;
using namespace Eigen;

WrapDoubleCylinder::WrapDoubleCylinder()
{
	m_type = double_cylinder;
	m_point_g = make_shared<Node>();
	m_point_g->x0.setZero();
	m_point_h = make_shared<Node>();
	m_point_h->x0.setZero();

}

WrapDoubleCylinder::WrapDoubleCylinder(const shared_ptr<Node> &P,
	const shared_ptr<Node> &S,
	const shared_ptr<CompDoubleCylinder> compDoubleCylinder,
	const int num_points)
	: WrapObst(P, S, num_points), m_compDoubleCylinder(compDoubleCylinder)
{
	m_type = double_cylinder;
	m_arc_points.resize(3, 3 * m_num_points + 1);
	m_point_g = std::make_shared<Node>();
	m_point_g->x0.setZero();
	m_point_h = std::make_shared<Node>();
	m_point_h->x0.setZero();
	m_radius_U = compDoubleCylinder->getRadiusA();
	m_radius_V = compDoubleCylinder->getRadiusB();
	m_point_U = compDoubleCylinder->getOriginA();
	m_point_V = compDoubleCylinder->getOriginB();
	m_vec_z_U = compDoubleCylinder->getZAxisA();
	m_vec_z_V = compDoubleCylinder->getZAxisB();

}

void WrapDoubleCylinder::load(const string &RESOURCE_DIR) {
	m_point_O->load(RESOURCE_DIR);
	m_point_P->load(RESOURCE_DIR);
	m_point_S->load(RESOURCE_DIR);
	m_point_g->load(RESOURCE_DIR);
	m_point_h->load(RESOURCE_DIR);
}

void WrapDoubleCylinder::init() {
	m_point_O->init();
	m_point_P->init();
	m_point_S->init();
	m_point_g->init();
	m_point_h->init();
}

void WrapDoubleCylinder::compute()
{
	// compute Matrix U and V
	Vector3d OP = m_point_P->x - m_point_U->x;
	OP = OP / OP.norm();
	Vector3d vec_Z_U = m_vec_z_U->dir / m_vec_z_U->dir.norm();
	Vector3d vec_X_U = vec_Z_U.cross(OP);
	vec_X_U = vec_X_U / vec_X_U.norm();
	Vector3d vec_Y_U = vec_Z_U.cross(vec_X_U);
	vec_Y_U = vec_Y_U / vec_Y_U.norm();

	Vector3d OS = m_point_S->x - m_point_V->x;
	OS = OS / OS.norm();
	Vector3d vec_Z_V = m_vec_z_V->dir / m_vec_z_V->dir.norm();
	Vector3d vec_X_V = vec_Z_V.cross(OS);
	vec_X_V = vec_X_V / vec_X_V.norm();
	Vector3d vec_Y_V = vec_Z_V.cross(vec_X_V);
	vec_Y_V = vec_Y_V / vec_Y_V.norm();

	this->M_U.resize(3, 3);
	this->M_U.row(0) = vec_X_U.transpose();
	this->M_U.row(1) = vec_Y_U.transpose();
	this->M_U.row(2) = vec_Z_U.transpose();

	this->M_V.resize(3, 3);
	this->M_V.row(0) = vec_X_V.transpose();
	this->M_V.row(1) = vec_Y_V.transpose();
	this->M_V.row(2) = vec_Z_V.transpose();

	// step 1: compute H and T
	Vector3d pv = this->M_V * (m_point_P->x - m_point_V->x);
	Vector3d sv = this->M_V * (m_point_S->x - m_point_V->x);

	double denom_h = pv(0)*pv(0) + pv(1)*pv(1);
	double denom_t = sv(0)*sv(0) + sv(1)*sv(1);
	double Rv = m_radius_V;

	double root_h = sqrt(denom_h - Rv*Rv);
	double root_t = sqrt(denom_t - Rv*Rv);

	Vector3d h(0.0, 0.0, 0.0);
	Vector3d t(0.0, 0.0, 0.0);
	h(0) = (pv(0) * Rv*Rv + Rv * pv(1) * root_h) / denom_h;
	h(1) = (pv(1) * Rv*Rv - Rv * pv(0) * root_h) / denom_h;
	t(0) = (sv(0) * Rv*Rv - Rv * sv(1) * root_t) / denom_t;
	t(1) = (sv(1) * Rv*Rv + Rv * sv(0) * root_t) / denom_t;

	if (Rv * (h(0) * t(1) - h(1) * t(0)) > 0.0)
	{
		status_V = no_wrap;
		h(0) = sv(0);
		h(1) = sv(1);
	}
	else {
		status_V = wrap;
	}

	m_status = wrap;

	complex<double> ht_i = 1.0 - 0.5 *
		((h(0) - t(0)) * (h(0) - t(0))
			+ (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
	complex<double> ph_i = 1.0 - 0.5 *
		((pv(0) - h(0)) * (pv(0) - h(0))
			+ (pv(1) - h(1)) * (pv(1) - h(1))) / (Rv*Rv);
	complex<double> ts_i = 1.0 - 0.5 *
		((t(0) - sv(0)) * (t(0) - sv(0))
			+ (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

	double ht_xy = abs(Rv * acos(ht_i));
	double ph_xy = abs(Rv * acos(ph_i));
	double ts_xy = abs(Rv * acos(ts_i));

	h(2) = pv(2) + (sv(2) - pv(2)) * ph_xy / (ph_xy + ht_xy + ts_xy);
	t(2) = sv(2) - (sv(2) - pv(2)) * ts_xy / (ph_xy + ht_xy + ts_xy);

	Vector3d H = this->M_V.transpose() * h + m_point_V->x;
	Vector3d T = this->M_V.transpose() * t + m_point_V->x;
	Vector3d H0 = H;

	Vector3d q(0.0, 0.0, 0.0);
	Vector3d g(0.0, 0.0, 0.0);
	Vector3d Q, G;

	double len = 0.0;
	Vector3d pu = this->M_U * (m_point_P->x - m_point_U->x);

	for (int i = 0; i < 30; i++)
	{
		len = 0.0;

		// step 2: compute Q and G
		Eigen::Vector3d hu = this->M_U * (H - m_point_U->x);

		double denom_q = pu(0)*pu(0) + pu(1)*pu(1);
		double denom_g = hu(0)*hu(0) + hu(1)*hu(1);
		double Ru = -m_radius_U;

		double root_q = sqrt(denom_q - Ru*Ru);
		double root_g = sqrt(denom_g - Ru*Ru);

		q(0) = (pu(0) * Ru*Ru + Ru * pu(1) * root_q) / denom_q;
		q(1) = (pu(1) * Ru*Ru - Ru * pu(0) * root_q) / denom_q;
		g(0) = (hu(0) * Ru*Ru - Ru * hu(1) * root_g) / denom_g;
		g(1) = (hu(1) * Ru*Ru + Ru * hu(0) * root_g) / denom_g;

		if (Ru * (q(0) * g(1) - q(1) * g(0)) > 0.0)
		{
			status_U = no_wrap;
			g(0) = pu(0);
			g(1) = pu(1);
		}
		else {
			status_U = wrap;
		}

		complex<double> qg_i = 1.0 - 0.5 *
			((q(0) - g(0)) * (q(0) - g(0))
				+ (q(1) - g(1)) * (q(1) - g(1))) / (Ru*Ru);
		complex<double> pq_i = 1.0 - 0.5 *
			((pu(0) - q(0)) * (pu(0) - q(0))
				+ (pu(1) - q(1)) * (pu(1) - q(1))) / (Ru*Ru);
		complex<double> gh_i = 1.0 - 0.5 *
			((g(0) - hu(0)) * (g(0) - hu(0))
				+ (g(1) - hu(1)) * (g(1) - hu(1))) / (Ru*Ru);

		double qg_xy = abs(Rv * acos(qg_i));
		double pq_xy = abs(Rv * acos(pq_i));
		double gh_xy = abs(Rv * acos(gh_i));
		len += qg_xy;

		q(2) = pu(2) + (hu(2) - pu(2)) * pq_xy / (pq_xy + qg_xy + gh_xy);
		g(2) = hu(2) - (hu(2) - pu(2)) * gh_xy / (pq_xy + qg_xy + gh_xy);

		Q = this->M_U.transpose() * q + m_point_U->x;
		G = this->M_U.transpose() * g + m_point_U->x;

		// step 3: compute H based on G and T
		Vector3d gv = this->M_V * (G - m_point_V->x);

		double denom_h = gv(0)*gv(0) + gv(1)*gv(1);
		double root_h = sqrt(denom_h - Rv*Rv);

		h = Vector3d(0.0, 0.0, 0.0);
		h(0) = (gv(0) * Rv*Rv + Rv * gv(1) * root_h) / denom_h;
		h(1) = (gv(1) * Rv*Rv - Rv * gv(0) * root_h) / denom_h;

		complex<double> ht_i = 1.0 - 0.5 *
			((h(0) - t(0)) * (h(0) - t(0))
				+ (h(1) - t(1)) * (h(1) - t(1))) / (Rv*Rv);
		gh_i = 1.0 - 0.5 *
			((gv(0) - h(0)) * (gv(0) - h(0))
				+ (gv(1) - h(1)) * (gv(1) - h(1))) / (Rv*Rv);
		complex<double> ts_i = 1.0 - 0.5 *
			((t(0) - sv(0)) * (t(0) - sv(0))
				+ (t(1) - sv(1)) * (t(1) - sv(1))) / (Rv*Rv);

		double ht_xy = abs(Rv * acos(ht_i));
		gh_xy = abs(Rv * acos(gh_i));
		double ts_xy = abs(Rv * acos(ts_i));
		len += ht_xy;

		h(2) = gv(2) + (sv(2) - gv(2)) * gh_xy / (gh_xy + ht_xy + ts_xy);
		t(2) = sv(2) - (sv(2) - gv(2)) * ts_xy / (gh_xy + ht_xy + ts_xy);

		if (Rv * (h(0) * t(1) - h(1) * t(0)) > 0.0)
		{
			status_V = no_wrap;
			h(0) = sv(0);
			h(1) = sv(1);
		}
		else {
			status_V = wrap;
		}

		H = this->M_V.transpose() * h + m_point_V->x;
		T = this->M_V.transpose() * t + m_point_V->x;

		len += (G - H).norm();

		double dist = (H - H0).norm();
		if (dist == 0) break;

		H0 = H;
	}

	if (status_V == no_wrap)
	{
		t = sv;
		T = this->M_V.transpose() * t + m_point_V->x;
	}
	else
	{
		status_V = wrap;
	}

	if (status_U == no_wrap)
	{
		q = pu;
		Q = this->M_U.transpose() * q + m_point_U->x;
	}
	else {
		status_U = wrap;
	}

	m_path_length = len;
	m_point_q->x = q;
	m_point_g->x = g;
	m_point_h->x = h;
	m_point_t->x = t;
	/*
	std::cout << Q.transpose() << std::endl << G.transpose() << std::endl
	<< H.transpose() << std::endl << T.transpose() << std::endl;
	*/
}


Eigen::MatrixXd WrapDoubleCylinder::getPoints(int num_points) const
{
	int col = 0;

	double theta_s, theta_e, theta_q, theta_g, theta_h, theta_t;
	double z_i, dz, z_s, z_e;

	MatrixXd points;

	if (status_U == wrap && status_V == wrap)
		points = Eigen::MatrixXd(3, 3 * num_points + 1);
	else if (status_U == wrap || status_V == wrap)
		points = Eigen::MatrixXd(3, 2 * num_points + 1);
	else
		points = Eigen::MatrixXd(3, 1 * num_points + 1);

	theta_q = atan(m_point_q->x(1) / m_point_q->x(0));

	if (m_point_q->x(0) < 0.0)
		theta_q += M_PI;

	theta_g = atan(m_point_g->x(1) / m_point_g->x(0));
	if (m_point_g->x(0) < 0.0)
		theta_g += M_PI;

	theta_h = atan(m_point_h->x(1) / m_point_h->x(0));
	if (m_point_h->x(0) < 0.0)
		theta_h += M_PI;

	theta_t = atan(m_point_t->x(1) / m_point_t->x(0));
	if (m_point_t->x(0) < 0.0)
		theta_t += M_PI;

	// q to g
	if (status_U == wrap)
	{
		if (theta_q < theta_g)
		{
			theta_s = theta_q; theta_e = theta_g;
			z_s = m_point_q->x(2); z_e = m_point_g->x(2);

		}
		else
		{
			theta_s = theta_g; theta_e = theta_q;
			z_s = m_point_g->x(2); z_e = m_point_q->x(2);
		}

		if (theta_e - theta_s > theta_s + 2 * M_PI - theta_e)
		{
			double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * M_PI;
			tmp = z_s; z_s = z_e; z_e = tmp;
		}

		z_i = z_s;
		dz = (z_e - z_s) / num_points;
		for (double i = theta_s; i <= theta_e + 0.001;
			i += (theta_e - theta_s) / num_points)
		{
			if (col == num_points + 1) {
				break;
			}
			Vector3d point = this->M_U.transpose() *
				Vector3d(m_radius_U * cos(i),
					m_radius_U * sin(i), z_i) +
				m_point_U->x;
			z_i += dz;
			points.col(col++) = point;
		}

	}
	else
	{
		points.col(col++) = m_point_P->x;
	}

	// g to h
	Vector3d G = this->M_U.transpose() * m_point_g->x + m_point_U->x;
	Vector3d H = this->M_V.transpose() * m_point_h->x + m_point_V->x;
	Vector3d diff = H - G;

	for (int i = 1; i < num_points; i++)
	{
		Vector3d point = G + diff / num_points * i;
		points.col(col++) = point;
	}

	if (status_V == wrap)
	{
		// h to t
		if (theta_h < theta_t)
		{
			theta_s = theta_h; theta_e = theta_t;
			z_s = m_point_h->x(2); z_e = m_point_t->x(2);
		}
		else
		{
			theta_s = theta_t; theta_e = theta_h;
			z_s = m_point_t->x(2); z_e = m_point_h->x(2);
		}

		if (theta_e - theta_s > theta_s + 2 * M_PI - theta_e)
		{
			double tmp = theta_s; theta_s = theta_e; theta_e = tmp + 2 * M_PI;
			tmp = z_s; z_s = z_e; z_e = tmp;
		}

		z_i = z_s;
		dz = (z_e - z_s) / num_points;
		int flag = 0;
		col = col + num_points;
		for (double i = theta_s; i <= theta_e + 0.001;
			i += (theta_e - theta_s) / num_points)
		{
			if (flag == num_points + 1) {
				break;
			}

			Vector3d point = this->M_V.transpose() *
				Vector3d(m_radius_V * cos(i),
					m_radius_V * sin(i), z_i) +
				m_point_V->x;
			z_i += dz;
			points.col(col--) = point;
			flag++;
		}
	}
	else
	{
		points.col(col++) = m_point_V->x;
	}

	return points;
}

void WrapDoubleCylinder::update_() {
	m_point_P->update();
	m_point_S->update();
	m_compDoubleCylinder->update();
	
	compute();

	if (status_U == wrap && status_V == wrap) {
		m_arc_points.resize(3, 3 * m_num_points + 1);
	}
	else if (status_U == wrap || status_V == wrap) {
		m_arc_points.resize(3, 2 * m_num_points + 1);
	}
	else {
		m_arc_points.resize(3, 1 * m_num_points + 1);
	}

	if (status_U == wrap || status_V == wrap) {
		m_arc_points = getPoints(m_num_points);
	}
}

void WrapDoubleCylinder::draw_(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> prog2, shared_ptr<MatrixStack> P) const {
	prog->bind();

	glUniformMatrix4fv(prog->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	MV->pushMatrix();
	glUniform3f(prog->getUniform("lightPos1"), 1.0f, 1.0f, 1.0f);
	glUniform1f(prog->getUniform("intensity_1"), 0.8f);
	glUniform3f(prog->getUniform("lightPos2"), -1.0f, 1.0f, 1.0f);
	glUniform1f(prog->getUniform("intensity_2"), 0.2f);
	glUniform1f(prog->getUniform("s"), 200.0f);
	glUniform3f(prog->getUniform("ka"), 0.4f, 0.3f, 0.5f);
	glUniform3f(prog->getUniform("kd"), 0.0f, 0.0f, 1.0f);
	glUniform3f(prog->getUniform("ks"), 0.0f, 1.0f, 0.0f);


	// Draw P, S, U, V points
	this->m_point_P->draw(MV, prog);
	this->m_point_S->draw(MV, prog);
	this->m_point_U->draw(MV, prog);
	this->m_point_V->draw(MV, prog);

	prog->unbind();

	// Draw wrapping
	prog2->bind();
	glUniformMatrix4fv(prog2->getUniform("P"), 1, GL_FALSE, glm::value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog2->getUniform("MV"), 1, GL_FALSE, glm::value_ptr(MV->topMatrix()));
	MV->pushMatrix();
	glColor3f(0.3f, 0.4f, 0.5f);
	glLineWidth(4);

	glBegin(GL_LINE_STRIP);
	glVertex3f(m_point_P->x(0), m_point_P->x(1), m_point_P->x(2));

	if (status_U == wrap || status_V == wrap) {
		for (int i = 0; i < m_arc_points.cols(); i++) {
			Vector3f p = m_arc_points.block<3, 1>(0, i).cast<float>();
			glVertex3f(p(0), p(1), p(2));
		}
	}

	glVertex3f(m_point_S->x(0), m_point_S->x(1), m_point_S->x(2));
	glEnd();

	// Draw z axis
	m_vec_z_U->draw(MV, P, prog2);
	m_vec_z_V->draw(MV, P, prog2);

	MV->popMatrix();
	prog2->unbind();
}

