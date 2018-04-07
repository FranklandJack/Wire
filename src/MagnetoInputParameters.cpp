#include "MagnetoInputParameters.hpp"

std::ostream& operator<<(std::ostream& out, const MagnetoInputParameters& params)
{
	int outputColumnWidth = 30;
    out << "Input-Parameters..." << '\n';
    switch(params.solutionMethod)
    {
        case MagnetoInputParameters::Jacobi:
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Solution-method: " << std::right << "Jacobi" <<'\n';
            break;

        case MagnetoInputParameters::GaussSeidel:
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Solution-method: " << std::right << "Gauss-Seidel" <<'\n';
            break;

        case MagnetoInputParameters::SOR:
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Solution-method: " << std::right << "SOR" <<'\n';
            out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "SOR-parameter: " << std::right << params.sorParameter <<'\n';
            break;

        default:
            break;

    }
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Spatial-discretistation: " << std::right << params.spaceStep << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Permittivity: " << std::right << params.permittivity <<'\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Initial-value: " << std::right << params.initialValue << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Initial-noise: " << std::right << params.noise << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Convergence-precision: " << std::right << params.precision<< '\n';
	out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Domain-x-range: " << std::right << params.xRange<< '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Domain-y-range: " << std::right << params.yRange << '\n';
    out << std::setw(outputColumnWidth) << std::setfill(' ') << std::left << "Output-directory: " << std::right << params.outputName << '\n';
    return out;
}
