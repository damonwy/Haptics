#include "rmpch.h"
#include "Solver.h"
#include "World.h"

using namespace std;
using namespace Eigen;

Solver::Solver()
{
	m_solutions = make_shared<Solution>();
}

Solver::Solver(shared_ptr<World> world, Integrator integrator) :
	m_world(world),
	m_integrator(integrator), step(0)
{
	m_solutions = make_shared<Solution>();
}

void Solver::reset() {
	nr = m_world->nr;
	nm = m_world->nm;
	// constraints
	nem = m_world->nem;
	ner = m_world->ner;
	ne = nem + ner;
	step = 0;
}
