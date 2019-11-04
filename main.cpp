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
		std::cout << "Running Part " << char(*part - 32) << "(x)" <<  std::endl;

	std::string prefix = std::string("bm7") + std::string(part) + std::string("-x");

	std::ofstream ofs;
	ofs.open((prefix + std::string(".csv")).c_str());

	if (rank == 0) {
		ofs << "NX,"
		    << "dt,"
		    << "dx,"
		    << "L2,"
		    << "lgR,"
		    << "lgE,"
		    << "fin,"
		    << "memory,"
		    << "runtime"
		    << std::endl;
	}

	std::vector<double> E, R;
	double lnX(lnX0);
	const int increment = floor(1.0 + runtime / dt);

	while (lnX > lnX1) {
		const unsigned NX = std::exp(lnX);
		auto start = std::chrono::steady_clock::now();

		char filename[512];
		sprintf(filename, "%s-%03d-ini.dat", prefix.c_str(), NX);
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
		const double lgE = std::log(l2);
		const double lgR = std::log(dx(grid, 0));

		E.push_back(lgE);
		R.push_back(lgR);

		// write grid output to file
		sprintf(filename, "%s-%03d-fin.dat", prefix.c_str(), NX);
		assert(std::string(filename).length() < 512);
		MMSP::output(grid, filename);

		const double gridsize = 2 * (glength(grid, 0) + 2 * ghosts(grid))
		                        * (glength(grid, 1) + 2 * ghosts(grid))
		                        * sizeof(double);

		auto end = std::chrono::steady_clock::now();

		if (rank == 0) {
			ofs << NX << ','
			    << dt << ','
			    << MMSP::dx(grid, 0) << ','
			    << l2 << ','
			    << lgR << ','
			    << lgE << ','
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

void temporal(const char* part, const double kappa, const double C2,
              const double lnX, const double lnT0, const double lnT1, const double dlnT)
{
	int rank = 0;
	#ifdef MPI_VERSION
	rank = MPI::COMM_WORLD.Get_rank();
	#endif

	if (rank == 0)
		std::cout << "Running Part " << char(*part - 32) << "(t)" << std::endl;

	std::string prefix = std::string("bm7") + std::string(part) + std::string("-t");

	std::ofstream ofs;
	ofs.open((prefix + std::string(".csv")).c_str());

	if (rank == 0) {
		ofs << "NX,"
		    << "dx,"
		    << "dt,"
		    << "L2,"
		    << "lgR,"
		    << "lgE,"
		    << "fin,"
		    << "memory,"
		    << "runtime"
		    << std::endl;
	}

	std::vector<double> E, R;
	double lnT(lnT0);
	const unsigned NX = std::exp(lnX);

	while (lnT > lnT1) {
		auto start = std::chrono::steady_clock::now();

		const double dt = 1.0 / exp(lnT);
		const unsigned NT = 1.0 / dt;
		const int increment = floor(1.0 + runtime / dt);

		char filename[512];
		sprintf(filename, "%s-%05u-ini.dat", prefix.c_str(), NT);
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
		const double lgE = std::log(l2);
		const double lgR = std::log(dt);

		E.push_back(lgE);
		R.push_back(lgR);

		// write grid output to file
		sprintf(filename, "%s-%05u-fin.dat", prefix.c_str(), NT);
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
			    << lgR << ','
			    << lgE << ','
			    << elapsed << ','
			    << gridsize << ','
			    << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() / 1000.0
			    << std::endl;
		}

		lnT *= dlnT;
	}

	ofs.close();
	ofs.open((prefix + std::string(".log")).c_str());

	if (rank == 0) {
		double b, m, cov00, cov01, cov11, sumsq;
		gsl_fit_linear(R.data(), 1, E.data(), 1, R.size(), &b, &m, &cov00, &cov01, &cov11, &sumsq);

		ofs << "log(E) = p log(R) + b:\n"
		    << "    p = " << m << '\n'
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

	// check argument list
	if (argc != 3) {
		std::cout << "Usage: " << PROGRAM << "part type" << std::endl;
		std::cout << "e.g., for Part a spatial: " << PROGRAM << " a x" << std::endl;
		std::cout << "                temporal: " << PROGRAM << " a t" << std::endl;
		MMSP::Abort(-1);
	}

	const char* part(argv[1]);
	const char* disc(argv[2]);

	if (*part == char(97) /* "a" */) {
		const double kappa = 4.0e-4;
		const double C2 = 0.0625 * M_PI;

		if (*disc == char(120) /* "x" */) {
			const double lnX0 = 5.7041;
			const double lnX1 = 5.25;
			const double dlnX = 0.985;
			const double dt = 2.0e-4;
			spatial(part, kappa, C2, lnX0, lnX1, dlnX, dt);
		}
		if (*disc == char(116)  /* "t" */) {
			const double lnT0 = 9.2104;
			const double lnT1 = 8.0;
			const double dlnT = 0.985;
			const double dx = 0.015625;
			temporal(part, kappa, C2, lnT0, lnT1, dlnT, dx);
		}
	} else if (*part == char(98) /* "b" */) {
		/*
		const double kappa = 1.5625e-6;
		const double C2 = 0.0625 * M_PI;
		*/
	}  else if (*part == char(99) /* "c" */) {
		/*
		const double kappa = 0.0004;
		const double C2 = 0.5;
		*/
	} else {
		std::cerr << "Error: Undefined argument \"" << argv[1] << "\": is it lower-case?" << std::endl;
	}

	MMSP::Finalize();
}

#endif
