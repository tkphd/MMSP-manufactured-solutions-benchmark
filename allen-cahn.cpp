// allen-cahn.cpp
// Algorithms for 2D and 3D Allen-Cahn model
// Questions/comments to gruberja@gmail.com (Jason Gruber)

#ifndef ALLENCAHN_UPDATE
#define ALLENCAHN_UPDATE

#include<cmath>
#include "mpi.h"

#include "MMSP.hpp"
#include "allen-cahn.hpp"
#include "manufactured.h"

const double runtime = 8.0;

const double A1 = 0.0075;
const double B1 = 8.0 * M_PI;
const double A2 = 0.03;
const double B2 = 22.0 * M_PI;

using std::cos;
using std::sin;
using std::tanh;

namespace MMSP
{

template<int dim, typename T>
double analyze(grid<dim, T>& Grid, const double elapsed, const double kappa, const double C2)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	double l2 = 0.0;
	const double dV = dx(Grid, 0) * dx(Grid, 1);
	for (int n = 0; n < nodes(Grid); n++) {
		vector<int> r = position(Grid, n);
		double x = dx(Grid, 0) * r[0];
		double y = dx(Grid, 1) * r[1];
		const double dEta = eta(A1, A2, B1, B2, C2, kappa, elapsed, x, y) - Grid(n);
		l2 += dEta * dEta * dV;
	}

	#ifdef MPI_VERSION
	double myL2(l2);
	MPI::COMM_WORLD.Allreduce(&l2, &myL2, 1, MPI_DOUBLE, MPI_SUM);
	#endif

	return std::sqrt(l2);
}

template<int dim, typename T>
void bc(grid<dim, T>& Grid)
{
	vector<int> r(dim, 0);

	if (x0(Grid, 1) == g0(Grid, 1)) {
		// lower boundary (y = 0)
		r[1] = x0(Grid, 1) - 1;
		for (r[0] = x0(Grid, 0); r[0] < x1(Grid, 0); r[0]++)
			Grid(r) = 1.0;
	}
	if (x1(Grid, 1) == g1(Grid, 1)) {
		// upper boundary (y = ½)
		r[1] = x1(Grid, 1);
		for (r[0] = x0(Grid, 0); r[0] < x1(Grid, 0); r[0]++)
			Grid(r) = 0.0;
	}
}

void generate(int dim, const char* filename, const int L, const double kappa, const double C2)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif
	if (dim == 2) {
		GRID2D initGrid(0, 0, 2 * L, 0, L);

		for (int d = 0; d < dim; d++)
			dx(initGrid, d) = 1.0 / (2 * L);

		vector<int> r(2, 0);

		for (int n = 0; n < nodes(initGrid); n++) {
			r = position(initGrid, n);
			double x = dx(initGrid, 0) * r[0];
			double y = dx(initGrid, 1) * r[1];
			double t = 0.0;

			initGrid(n) = eta0(A2, B2, kappa, x, y);
		}

		bc(initGrid);

		output(initGrid, filename);
	} else {
		std::cerr << "Error: PFHub Benchmark 7 is 2-D, only." << std::endl;
	}
}

template <int dim, typename T>
double update(grid<dim, T>& Grid, const unsigned steps, const double dt, const double kappa, const double C2)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	ghostswap(Grid);

	const double LinStab = dx(Grid, 0) * dx(Grid, 1) / kappa;
	assert(LinStab < 1.0);

	grid<dim, T> newGrid(Grid);

	double elapsed = 0.0;

	for (int step = 0; step < steps; step++) {
		if (rank == 0)
			print_progress(step, steps);

		for (int n = 0; n < nodes(Grid); n++) {
			vector<int> r = position(Grid, n);
			double x = dx(Grid, 0) * r[0];
			double y = dx(Grid, 1) * r[1];
			T phi = Grid(n);
			newGrid(n) = phi + dt * (- 4 * phi * (phi - 1) * (phi - 0.5) + kappa * laplacian(Grid, n)
			                         + S(A1, A2, B1, B2, C2, kappa, elapsed, x, y));
		}

		swap(Grid, newGrid);
		ghostswap(Grid);
		bc(Grid);

		elapsed += dt;
	}
	return elapsed;
}

} // namespace MMSP

#endif

#include"main.cpp"
