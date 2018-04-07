#ifndef MagnetoInputParameters_hpp
#define MagnetoInputParameters_hpp
#include <iostream>
#include <iomanip>
/**
 *\file
 *\class MagnetoInputParameters
 *\brief Class for easily handling input parameters of Magneto differential equation solver.
 *
 * This class essentially just holds some values and has an operator to easily output
 * them to a stream.
 */
class MagnetoInputParameters
{

public:

    /**
     *\enum to hold the three solution methods we will consider when solving the Magneto equation
     */
    enum SolutionMethod
    {
        Jacobi,
        GaussSeidel,
        SOR
    };

    /// Solution method.
    SolutionMethod solutionMethod;

	/// Spatial discretisation step.
    double spaceStep;

    /// Permittivity in the Magneto equation.
    double permittivity;

    /// \phi_0 initial value of potential.
    double initialValue;

    /// Maximum magnitude of initial noise.
    double noise;

    /// Precision of the final answer in terms of convergence.
    double precision;

    /// Range of x-values lattice domain.
    int xRange;

    /// Range of y-values lattice domain.
    int yRange;

    /// Name of output directory to save any output into.
    std::string outputName;

    /// Successive over relaxation parameter.
    double sorParameter;

    /**
	 *\brief operator<< overload for outputting the results.
	 *\param out std::ostream reference that is the stream being outputted to.
	 *\param params constant MagnetoInputParameters instance to be output.
	 *\return std::ostream reference so the operator can be chained.
	 *
	 * Results will be output in a formatted table for easy viewing in the command line or a file.
	 */
    friend std::ostream& operator<<(std::ostream& out, const MagnetoInputParameters& params);

};
#endif /* MagnetoInputParameters_hpp */
