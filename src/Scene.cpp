#include "rmpch.h"
#include "Scene.h"
#include "Node.h"
#include "Joint.h"
#include "Vector.h"
#include "JsonEigen.h"
#include "World.h"
#include "Solver.h"
#include "SolverDense.h"
#include "SolverSparse.h"
//#include "Spring.h"
#include "Deformable.h"
#include "DeformableSpring.h"
//#include <boost/numeric/odeint.hpp>
#include "EigenOdeintHelper.h"

using namespace std;
using namespace Eigen;
using json = nlohmann::json;
//using namespace boost::numeric::odeint;

#include <unsupported/Eigen/MatrixFunctions> // TODO: avoid using this later, write a func instead

Scene::Scene() :
	t(0.0),
	h(1e-2),
    drawHz(10),
	grav(0.0, 0.0, 0.0)
{
}

Scene::~Scene()
{
}

void Scene::load(const string &RESOURCE_DIR)
{	
	//read a JSON file
	ifstream i(RESOURCE_DIR + "input.json");
	i >> js;
	i.close();

	// Units: meters, kilograms, seconds
	h = js["h"];
	Eigen::from_json(js["grav"], grav);
	drawHz = js["drawHz"];

	m_world = make_shared<World>(MUSCLE_INERTIA);//_INVERTIBLE
	m_world->load(RESOURCE_DIR);

	//m_solver = make_shared<SolverDense>(m_world, REDMAX_EULER);
	m_solver = make_shared<SolverSparse>(m_world, REDMAX_EULER, LU);	
}


void Scene::init()
{
	count = 0;
	m_world->init();
	VectorXd y0, y1;
	y0.resize(2 * m_world->nr);
	y0.setZero();
	
	//y1 = m_solver->dynamics(y0);
	y.resize(2 * m_world->nr);
	y.setZero();
	m_world->getJoint0()->reparam();
	m_world->getJoint0()->gatherDofs(y, m_world->nr);
	m_world->getDeformable0()->gatherDofs(y, m_world->nr);
	//m_solution = m_solver->solve();
	//vec_to_file(m_solution->t, "t");
	//mat_to_file(m_solution->y, "y");

	//tk = m_solution->t(0);
	drawH = 1.0 / drawHz;
	search_idx = 0;
}

void Scene::reset()
{
	
}


void Scene::solve() {
}

int torend = 0;
void Scene::step()
{
	VectorXd ys;
	//runge_kutta_cash_karp54<Eigen::VectorXd> stepper;
	//vector<Eigen::VectorXd> x_vec;
	//vector<double> times;

	//integrate_adaptive(stepper, , y, 0.0, h, h, Observer(x_vec, times));
	//cout << x_vec[1] << endl;	
	y = m_solver->dynamics(y);
	


	//m_world->getJoint0()->reparam();
	//m_world->getJoint0()->gatherDofs(y, m_world->nr);
	m_world->update();
	m_world->incrementTime();
	torend++;

	count++;
	if (count == 99) {
		cout << count << endl;
	}
	
	saveEnergyData(100);
	//if(tk < m_solution->t(n_steps-1)) {
	//	m_solution->searchTime(tk, search_idx, output_idx, s);
	//	search_idx = output_idx;
	//	ys = (1 - s)* m_solution->y.row(output_idx) + s * m_solution->y.row(output_idx + 1);

	//	m_world->getJoint0()->scatterDofs(ys, m_world->nr);
	//	m_world->getSpring0()->scatterDofs(ys, m_world->nr);
	//	m_world->getSoftBody0()->scatterDofs(ys, m_world->nr);
	//	tk = tk + drawH;
	//}
	//else {
	//	// reset
	//	tk = m_solution->t(0);
	//}	

}

void Scene::saveEnergyData(int num_steps) {
	// Save data in matlab
	if (count % 1 == 0) {
		m_energy_vector.push_back(m_solver->m_energy);
		m_time_vector.push_back(m_world->getTime());
	}

	if (count == num_steps) {
		cout << "finished!" << endl;
		VectorXd Kv;
		VectorXd Vv;
		VectorXd Tv;

		Kv.resize(m_energy_vector.size());
		Vv.resize(m_energy_vector.size());
		Tv.resize(m_energy_vector.size());

		for (int i = 0; i < m_energy_vector.size(); i++) {
			Kv(i) = m_energy_vector[i].K;
			Vv(i) = m_energy_vector[i].V;
			Tv(i) = m_time_vector[i];
		}
		
		vec_to_file(Kv, "Kcpp");
		vec_to_file(Vv, "Vcpp");
		vec_to_file(Tv, "Tcpp");

	}
}


void Scene::draw(shared_ptr<MatrixStack> MV, const shared_ptr<Program> prog, const shared_ptr<Program> progSimple, const shared_ptr<Program> progSoft, shared_ptr<MatrixStack> P) const
{
	m_world->draw(MV, prog, progSimple, progSoft, P);
	
}
