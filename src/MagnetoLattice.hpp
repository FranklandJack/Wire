#ifndef MagnetoLattice_hpp
#define MagnetoLattice_hpp
#include <cmath>
#include <random>
#include <array>
#include <iostream>

/**
 *\file
 *\class MagnetoLattice
 *\brief Lattice for evolving electrostatic potential according to the Magneto equation
 * with a selection of algorithms.
 *
 * 2D lattice consisting of an array of floating points which represent the values of the potential
 * at some time t. The lattice can be evolved through computer time with a selection of algorithms to
 * converge on the correct function.
 */
class MagnetoLattice
{
private:
	/// range of x-values.
	int m_xRange;

	/// Range of y-values.
	int m_yRange;

	/// Permittivity constant.
	double m_permativity;

	/// Lattice space discretisation step size.
	double m_dx;

	/// Charge density for the potential.
	std::vector<double> m_chargeDensity;

	/// The actual potential on the lattice.
	std::vector<double> m_potential;

public:
	/**
	 *\brief Constructs a MagnetoLattice of the specified size, permittivity and step size.
	 *
	 * The actual potential will consist of a cube where each dimension is one larger than specified, and
	 * we will introduce a halo of sites where \phi = 0 to impose the boundary conditions.
	 *
	 *\param xRange range of x values in lattice.
	 *\param yRange range of y values in lattice.
	 *\param zRange range of z values in lattice.
	 *\param permittivity permittivity in Magneto equation.
	 *\param dx spatial discretisation step size.
	 *
	 */
	MagnetoLattice(int xRange, int yRange, double permativity, double dx);

	/**
	 *\brief Initialises the non-boundary entries in the lattice with a value and some uniformly distributed
	 * noise of a magnitude specified by the user.
	 *\param initialValue floating point representing the initial value at each lattice site.
	 *\param noise magnitude of random noise.
	 *\param generator for generating the random numbers according to the uniform distribution.
	 */
	void initialise(double initialValue, double noise, std::default_random_engine &generator);

	/**
	 *\brief Operator overload for accessing the correct elements in the lattice and taking into account the
	 * halo of the boundary condition.
	 *\param i x index.
	 *\param j y index.
	 *\param k z index.
	 */
	double& operator()(int i, int j);

	/**
	 *\brief Operator overload for accessing the correct elements in the lattice and taking into account the
	 * halo of the boundary condition (constant version).
	 *\param i x index.
	 *\param j y index.
	 *\param k z index.
	 */
	const double& operator()(int i, int j) const;

	/**
	 *\brief gets the charge density at a site.
	 *\param i x index of site.
	 *\param j y index of site.
	 *\param k z index of site.
	 *\return charge density at that site.
	 */
	double getChargeDensity(int i, int j) const;

	/**
	 *\brief sets the charge density at a site with a specified value.
	 *\param i x index of site.
	 *\param j y index of site.
	 *\param k z index of site.
	 *\param charge charge density at that point.
	 */
	void setChargeDensity(int i, int j, double charge);

	/**
	 *\brief updates the potential on the lattice according to the Jacobi algorithm
	 *\param latticeCurrent lattice to be used to do the update based on.
	 *\param latticeUpdated lattice to be updated based on current lattice.
	 *\param returns floating point representing how close to lattices are in terms of convergence.
	 */
	friend double jacobiUpdate(MagnetoLattice &currentLattice, MagnetoLattice &updatedLattice);

	/**
	 *\brief updates the potential based on the GaussSeidel algorithm.
	 *\param lattice reference to lattice to be updated.
	 *\return floating point representing how close old lattice was to updated one.
	 */
	friend double gaussSeidelUpdate(MagnetoLattice &lattice);

	/*
	 *\brief updates the potential based on the GS algorithm with successive over-relaxation.
	 *
	 * Updates are computed according to the SOR algorithm x(n+1) = wxf(x(n)) + (1-w)x(n) where
	 * x(n) is our original update rule.
	 *
	 *\param sorParameter floating point value representing the SOR-parameter omega.
	 *\param lattice Magneto lattice to be updated via the SOR algorithm.
	 *\return floating point representing how close old lattice was to updated one.
	 */
	friend double sorUpdate(double sorParameter, MagnetoLattice &lattice);

	/**
	 *\brief Calculates the next value of the potential at that site based on the Jacobi update.
	 *\param i x index.
	 *\param j y index.
	 *\param k z index.
	 *\return floating point representing the next value at that site.
	 */
	double nextValueJacobi(int i, int j) const;

	/**
	 *\brief Calculates the magnetic field based on the magnetostatic equation B = \grad x A.
	 *\param i x coordinate of point.
	 *\param j y coordinate of point.
	 *\param k z coordinate of point.
	 *\return array of three doubles representing the x,y and z components of B.
	 */
	std::array<double,2> magneticField(int i, int j) const;


	/**
	 *\brief operator overload to print various pieces of data into a single file for efficiency
	 *
	 * Will print to file in the form x y z r \phi E where x y z are coordinates r is their magnitude \phi is
	 * the solution at that point and E is the electric field.
	 *
	 *\param out output stream reference to stream to.
	 *\param lattice Magneto lattice to print from.
	 *\return ostream reference so output can be chained.
	 */
	friend std::ostream& operator<<(std::ostream &out, const MagnetoLattice &lattice);

	/**
	 *\brief sets wire distribution.
	 */
	 void setWireDist();


};

#endif /* MagnetoLattice_hpp */
