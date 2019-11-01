// MMSP.main.hpp

#ifndef MMSP_MAIN
#define MMSP_MAIN
#include <cassert>
#include <cctype>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <gsl/gsl_fit.h>
#include <iostream>
#include <string.h>
#include <vector>

void spatial(const char* part, const double kappa, const double C2,
			 const double lnX0, const double lnX1, const double dlnX,
			 const double dt)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	if (rank == 0)
		std::cout << "Running Part " << char(*part - 32) << std::endl;

	std::string prefix = std::string("bm7") + std::string(part);

	std::ofstream ofs;
	ofs.open((prefix + std::string(".csv")).c_str());

	if (rank == 0) {
		ofs << "NX,"
		    << "dx,"
		    << "dt,"
		    << "L2,"
		    << "fin,"
		    << "memory,"
		    << "runtime"
		    << std::endl;
	}

	std::vector<double> E, R;
	double lnX(lnX0);
	const int increment = floor(1.0 + 8.0 / dt);

	while (lnX > lnX1) {
		const unsigned NX = std::exp(lnX);
		auto start = std::chrono::steady_clock::now();

		char filename[512];
		sprintf(filename, "bm7a-%03d-ini.dat", NX);
		assert(std::string(filename).length() < 512);

		MMSP::generate(2, filename, NX, kappa, C2);

		// file open error check
		std::ifstream input(filename);
		if (!input) {
			std::cerr << "File input error: could not open " << filename << ".\n\n";
			MMSP::Abort(-1);
		} else {
			input.close();
		}

		// construct grid object
		GRID2D grid(filename);

		const double elapsed = MMSP::update(grid, increment, dt, kappa, C2);

		const double l2 = MMSP::analyze(grid, elapsed, kappa, C2);

		E.push_back(l2);
		R.push_back(dx(grid, 0));

		// write grid output to file
		sprintf(filename, "bm7a-%03u-fin.dat", NX);
		assert(std::string(filename).length() < 512);
		MMSP::output(grid, filename);

		const double gridsize = 2 * (glength(grid, 0) + 2 * ghosts(grid))
		                        * (glength(grid, 1) + 2 * ghosts(grid))
		                        * sizeof(double);

		auto end = std::chrono::steady_clock::now();

		if (rank == 0) {
			ofs << NX << ','
			    << MMSP::dx(grid, 0) << ','
			    << dt << ','
			    << l2 << ','
			    << elapsed << ','
			    << gridsize << ','
			    << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0
			    << std::endl;
		}

		lnX *= dlnX;
	}

	ofs.close();
	ofs.open((prefix + std::string(".log")).c_str());

	if (rank == 0) {
		double p, b, cov00, cov01, cov11, sumsq;
		gsl_fit_linear(R.data(), 1, E.data(), 1, R.size(), &p, &b, &cov00, &cov01, &cov11, &sumsq);

		ofs << "log(E) = p log(R) + b:\n"
		    << "    p = " << p << '\n'
		    << "    b = " << b << '\n'
		    << "    R² = " << sumsq << '\n'
		    << "    Covariance: "
		    << cov00 << '\t' << cov01 << '\t' << cov11
		    << std::endl;
	}
}

int main(int argc, char* argv[])
{
	MMSP::Init(argc, argv);

	const char* part(argv[1]);

	// check argument list
	if (argc != 2) {
		std::cout << "Usage: " << PROGRAM << "part" << std::endl;
		std::cout << "e.g., for Part a: " << PROGRAM << " a" << std::endl;
		MMSP::Abort(-1);
	}

	if (*part == char(97)) {
		const double kappa = 4.0e-4;
		const double C2 = 0.0625 * M_PI;
		const double lnX0 = 5.5219;
		const double lnX1 = 3.1638;
		const double dlnX = 0.975;
		const double dt = 2.0e-4;
		spatial(part, kappa, C2, lnX0, lnX1, dlnX, dt);
	} else if (*part == char(98)) {
		const double kappa = 1.5625e-6;
		const double C2 = 0.0625 * M_PI;
		const double lnX0 = 5.5219;
		const double lnX1 = 3.1638;
		const double dlnX = 0.975;
		const double dt = 2.0e-4;
		spatial(part, kappa, C2, lnX0, lnX1, dlnX, dt);
	}  else if (*part == char(99)) {
		const double kappa = 0.0004;
		const double C2 = 0.5;
		const double lnX0 = 5.5219;
		const double lnX1 = 3.1638;
		const double dlnX = 0.975;
		const double dt = 2.0e-4;
		spatial(part, kappa, C2, lnX0, lnX1, dlnX, dt);
	} else {
		std::cerr << "Error: Undefined argument \"" << argv[1] << "\": is it lower-case?" << std::endl;
	}

	MMSP::Finalize();
}

#endif
