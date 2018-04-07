#include <iostream> // For file IO.
#include <boost/filesystem.hpp> // For constructing directories for file IO.
#include <boost/program_options.hpp> // For command line arguments.
#include <fstream> // For file output.
#include <chrono> // For timing.
#include <ctime>  // For timing.
#include <random> // For generating random numbers.
#include <cmath> // For any maths functions.
#include <iomanip> // For manipulating output.
#include <string> // For naming output directory.
#include <fstream> // For file output.
#include <algorithm> // For swapping the lattices.
#include "Timer.hpp" // For custom timer.
#include "makeDirectory.hpp" // For making directories.
#include "MagnetoInputParameters.hpp" // For neatly packaging together input parameters.
#include "MagnetoLattice.hpp"


int main(int argc, char const *argv[])
{
/*************************************************************************************************************************
************************************************* Preparations **********************************************************
*************************************************************************************************************************/
    // Start the clock so execution time can be calculated.
    Timer timer;

    // Seed the pseudo random number generator using the system clock.
    unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());

    // Create a generator that can be fed to any distribution to produce pseudo random numbers according to that distribution.
    std::default_random_engine generator(seed);

/*************************************************************************************************************************
******************************************************** Input **********************************************************
*************************************************************************************************************************/
    // Input parameters for the simulation.

    // Method for relaxation algorithm.
    MagnetoInputParameters::SolutionMethod solutionMethod;

    // Spatial discretisation step.
    double spaceStep;

    // Permittivity in Magneto equation.
    double permittivity;

    // \phi_0 initial value on lattice away from boundary.
    double initialValue;

    // Maximum magnitude of initial noise.
    double noise;

    // Precision of convergence.
    double precision;

    // Number of x values in cubic lattice domain.
    int xRange;

    // Number of y values in cubic lattice domain.
    int yRange;

    // Name of output directory to save any output into.
    std::string outputName;

    // SOR update parameter.
    double sorParameter;

    // Set up optional command line argument.
    boost::program_options::options_description desc("Options for Magneto simulation");

    desc.add_options()
        ("spatial-discretisation,x", boost::program_options::value<double>(&spaceStep)->default_value(1), "Spatial discretisation step size.")
        ("permittivity,p", boost::program_options::value<double>(&permittivity)->default_value(1), "M parameter from Cahn-Hilliard equation.")
        ("initial-value,v", boost::program_options::value<double>(&initialValue)->default_value(0), "Initial value of order parameter.")
        ("noise,n",boost::program_options::value<double>(&noise)->default_value(0.0), "Maximum magnitude of initial noise.")
        ("precision,d", boost::program_options::value<double>(&precision)->default_value(0.001),"Precision of convergence.")
        ("x-range,r", boost::program_options::value<int>(&xRange)->default_value(100),"Total number of x points in domain of simulation domain.")
        ("y-range,c", boost::program_options::value<int>(&yRange)->default_value(100),"Total number of y points in domain of simulation domain.")
        ("output,o",boost::program_options::value<std::string>(&outputName)->default_value(getTimeStamp()), "Name of output directory to save output files into.")
        ("sor-parameter,w",boost::program_options::value<double>(&sorParameter)->default_value(1),"Parameter for the successive over-relaxation algorithm.")
        ("Jacobi","Use Jacobi relaxation method")
        ("Gauss-Seidel","Use Gauss-Seidel relaxation method (will take precedence over Gauss-Seidel")
        ("SOR","Use successive over relaxation method with Gauss-Seidel algorithm, will take overall precedence")
        ("help,h","Display help message.");


    boost::program_options::variables_map vm;
    boost::program_options::store(boost::program_options::parse_command_line(argc,argv,desc), vm);
    boost::program_options::notify(vm);

    // If the user asks for help display it then exit.
    if(vm.count("help"))
    {
        std::cout << desc << "\n";
        return 1;
    }

    // If the user asks for specific algorithm use it.
    if(vm.count("SOR"))
    {
        solutionMethod = MagnetoInputParameters::SOR;
    }
    else if(vm.count("Gauss-Seidel"))
    {
        solutionMethod = MagnetoInputParameters::GaussSeidel;
    }
    else
    {
        solutionMethod = MagnetoInputParameters::Jacobi;
    }


    // Construct an input parameter object, this just makes printing a lot cleaner.
    MagnetoInputParameters inputParameters
    {
        solutionMethod,
        spaceStep,
        permittivity,
        initialValue,
        noise,
        precision,
        xRange,
        yRange,
        outputName,
        sorParameter
    };


/*************************************************************************************************************************
************************************************* Create Output Files ***************************************************
*************************************************************************************************************************/

    // Create an output directory from either the default time stamp or the user defined string.
    makeDirectory(outputName);

    // Create output file for the input parameters.
    std::fstream inputParameterOutput(outputName+"/input.txt", std::ios::out);

    // Create single output file for lattice so we can print the potential and anything else at the same time
    // which is quicker than iterating through 3D lattice multiple times.
    std::fstream magnetoOutput(outputName+"/magnetoOutput.dat", std::ios::out);

    // Output file to hold any statistical results e.g. number of iterations until convergence.
    std::fstream outputResults(outputName+"/results.txt", std::ios::out);

    // Print input parameters to command line.
    std::cout << inputParameters << '\n';

    // Print input parameters to file.
    inputParameterOutput << inputParameters << '\n';

/*************************************************************************************************************************
************************************************* The Simulation ********************************************************
*************************************************************************************************************************/

// There will be three possible algorithm choices and the structure of the simulation will depend on which is used so here
// the program takes three possible branches.

// Create a lattice to hold the current state of the potential.
    MagnetoLattice currentLattice(xRange, yRange, permittivity, spaceStep);

// Initialise the lattice with some value and random noise.
    currentLattice.initialise(initialValue, noise, generator);

// Set the current density on the wire.
    currentLattice.setWireDist();

// Create a lattice to hold the updated state of the lattice.
    MagnetoLattice  updatedLattice = currentLattice;

// Create a counter to calculate the number of iterations required for convergence.
    int counter = 0;

// Create a variable to hold how ``converged'' the lattice is relative to the user defined precision.
    double convergence;

switch(solutionMethod)
{
    // The case the user specifies to use the Jacobi update rule.
    case MagnetoInputParameters::Jacobi:
        {
            while(true)
            {
                // Count the number of times we have to do an update before convergence.
                ++counter;
                // Update the lattice based on its current state.
                convergence = jacobiUpdate(currentLattice, updatedLattice);
                // Swap the lattices so the current lattice becomes the updated one and the updated one the one will update into.
                std::swap(currentLattice, updatedLattice);

                if(0==counter%1000)
                {
                    std::cout << counter << ' ' << convergence << '\n';
                }

                // Check to see if the lattice has converged and if it has stop updating the lattice.
                if(convergence < precision)
                {
                    break;
                }


            }

        }

        break;

    // The case the user specifies to use the Gauss Seidel algorithm.
    case MagnetoInputParameters::GaussSeidel:
            {
                while(true)
                {
                    // Count the number of times we have to do an update before convergence.
                    ++counter;

                    // Update the lattice based on the GaussSeidal algorithm
                    convergence = gaussSeidelUpdate(currentLattice);

                    if(0==counter%1000)
                    {
                        std::cout << counter << ' ' << convergence << '\n';
                    }

                    // Check to see if the lattice has converged and if it has stop updating the lattice.
                    if(convergence < precision)
                    {
                        break;
                    }

                }
            }

            break;

    // The case the user specifies to use the successive over-relaxation algorithm.
    case MagnetoInputParameters::SOR:
            {
                while(true)
                {
                    // Count the number of times we have to do an update before convergence.
                    ++counter;

                    // Update the lattice based on the GaussSeidal algorithm
                    convergence = sorUpdate(sorParameter, currentLattice);

                    if(0==counter%1000)
                    {
                        std::cout << counter << ' ' << convergence << '\n';
                    }

                    // Check to see if the lattice has converged and if it has stop updating the lattice.
                    if(convergence < precision)
                    {
                        break;
                    }

                }
            }

            break;

    // Default statement to stop compiler throwing warning-doesn't actually do anything.
    default:
        break;

}

/*************************************************************************************************************************
***********************************************  Output/Clean Up ********************************************************
*************************************************************************************************************************/

    // Save the potential and field to a file.
    magnetoOutput << currentLattice;

    // Report how many iterations the program took and how long the program took to execute in time and save that data to file.
    double runTime = timer.elapsed();
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Number-of-iterations-until-convergence: " << std::right << counter << std::endl;
    std::cout << std::setw(30) << std::setfill(' ') << std::left << "Time-take-to-execute(s): " << std::right << runTime << std::endl << std::endl;

    outputResults << std::setw(30) << std::setfill(' ') << std::left << "Number-of-iterations-until-convergence: " << std::right << counter << std::endl;
    outputResults << std::setw(30) << std::setfill(' ') << std::left << "Time-take-to-execute(s): " << std::right << runTime << std::endl << std::endl;


    return 0;
}
