#include "MagnetoLattice.hpp"

MagnetoLattice::MagnetoLattice(int xRange, int yRange, double permativity, double dx): m_xRange(xRange),
																								   m_yRange(yRange),
																								   m_permativity(permativity),
																								   m_dx(dx),
																								   m_chargeDensity(xRange * yRange, 0.0),
																								   m_potential(xRange * yRange, 0.0)
{

}

void MagnetoLattice::initialise(double initialValue, double noise, std::default_random_engine &generator)
{
	// Create the uniform distribution for generating the random numbers.
	std::uniform_real_distribution<double> noiseDistribution(-noise, noise);
	// Only update from 1 to range-1 since boundaries should be fixed by initial conditions.
		for(int j = 1; j < m_yRange-1; ++j)
		{
			for(int i = 1; i < m_xRange-1; ++i )
			{
				(*this)(i,j) = initialValue + noiseDistribution(generator);
			}
		}


}

double& MagnetoLattice::operator()(int i, int j)
{
	return m_potential[i + j*m_xRange];
}

const double& MagnetoLattice::operator()(int i, int j) const
{
	return m_potential[i + j*m_xRange];
}

double MagnetoLattice::getChargeDensity(int i, int j) const
{
	return m_chargeDensity[i + j*m_xRange];
}


void MagnetoLattice::setChargeDensity(int i, int j, double charge)
{
	m_chargeDensity[i + j*m_xRange] = charge;
}

double MagnetoLattice::nextValueJacobi(int i, int j) const
{
	return (((*this)(i+1,j) + (*this)(i-1,j)
		         + (*this)(i,j+1) + (*this)(i,j-1)
		          + (pow(m_dx,2)/m_permativity) * getChargeDensity(i,j))/4.0);

}


double jacobiUpdate(MagnetoLattice &currentLattice, MagnetoLattice &updatedLattice)
{
	double convergenceMeasure = 0;

	// Only need to check from 1 to range-1 since bounaries are fixed.
	for(int i = 1; i < currentLattice.m_xRange-1; ++i)
	{
		for(int j = 1; j < currentLattice.m_yRange-1; ++j)
		{
				updatedLattice(i,j) = currentLattice.nextValueJacobi(i,j);

				convergenceMeasure += abs(updatedLattice(i,j)-currentLattice(i,j));
				std::cout << convergenceMeasure << '\n';

		}
	}

	return convergenceMeasure;

}

double gaussSeidelUpdate(MagnetoLattice &lattice)
{
	double convergenceMeasure = 0;

	double updatedValue = 0;
	double currentValue = 0;

	// Only update from 1 to range-1 since boundaries should be fixed by initial conditions.
	for(int i = 1; i < lattice.m_xRange-1; ++i)
	{
		for(int j = 1; j < lattice.m_yRange-1; ++j)
		{


				updatedValue = lattice.nextValueJacobi(i,j);

				currentValue = lattice(i,j);

				convergenceMeasure += abs(updatedValue-currentValue);

				lattice(i,j) = updatedValue;


		}
	}

	return convergenceMeasure;

}

double sorUpdate(double sorParameter, MagnetoLattice &lattice)
{
	// Local variable to measure how much the lattice changes on this update.
	double convergenceMeasure = 0;

	// Local variable to hold the standard GS update.
	double updatedGSValue = 0;

	// Local variable to hold current value at given lattice site.
	double currentValue = 0;

	// Local variable to hold the final update value under the SOR update.
	double updatedSORValue = 0;

	// Only update from 1 to range-1 since boundaries should be fixed by initial conditions.
	for(int i = 1; i < lattice.m_xRange-1; ++i)
	{
		for(int j = 1; j < lattice.m_yRange-1; ++j)
		{
				currentValue = lattice(i,j);

				updatedGSValue = lattice.nextValueJacobi(i,j);

				updatedSORValue = (1-sorParameter) * currentValue + sorParameter * updatedGSValue;

				lattice(i,j) = updatedSORValue;

				convergenceMeasure += abs(updatedSORValue-currentValue);


		}
	}

	return convergenceMeasure;

}







std::array<double,2> MagnetoLattice::magneticField(int i, int j) const
{
	std::array<double,2> magneticField = {((*this)(i,j+1)-(*this)(i,j-1))/(2*m_dx),
							-((*this)(i+1,j)-(*this)(i-1,j))/(2*m_dx)};

	return magneticField;

}



std::ostream& operator<<(std::ostream &out, const MagnetoLattice &lattice)
{
	// Calculate the centre of the lattice utilising integer division.
	double xCentre = lattice.m_xRange/2;
	double yCentre = lattice.m_yRange/2;

	std::array<double,2> magneticFieldTemp;

		for(int j = 0; j < lattice.m_yRange; ++j)
		{
			for(int i = 0; i < lattice.m_xRange; ++i )
			{
				// Calculate distance of current point from centre of lattice.
				double xDistance = xCentre - i;
				double yDistance = yCentre - j;
				double radialDistance = sqrt(xDistance*xDistance+yDistance*yDistance);

				double fieldStrength;


				// Make sure we aren't calculating field at edge of lattice since this will try and access
				// memory the program doesn't own via the indexing during the calculation - field is just
				// set to zero at these points.
				if(i==0 || j==0 || i==lattice.m_xRange-1 || j==lattice.m_yRange-1)
				{
					magneticFieldTemp = std::array<double,2>{0,0};
				}
				else
				{
					magneticFieldTemp = lattice.magneticField(i,j);
				}

				fieldStrength = sqrt(magneticFieldTemp[0]*magneticFieldTemp[0]+magneticFieldTemp[1]*magneticFieldTemp[1]);

				out << i << ' ' << j << ' ' <<
				radialDistance << ' ' << lattice(i,j) <<
				' ' << magneticFieldTemp[0] << ' ' << magneticFieldTemp[1] << ' ' << fieldStrength << '\n';
			}

			out << '\n';
		}



	return out;

}



 void MagnetoLattice::setWireDist()
 {
	 // Utilise integer division to find the centre of the box.
	 int xCentre = m_xRange/2;
	 int yCentre = m_yRange/2;

	 // Set current magnitude.
	 double current = 1;

	setChargeDensity(xCentre, yCentre, current);

 }
