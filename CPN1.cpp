#include "stdafx.h"
#include <math.h>
#include <omp.h>
#include <chrono>
#include <iostream>
#include <random>
#include <fstream>
#include <string>
#include <random>
/*#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
*/

//Simulation Parameters
#define vd 1	//Dimension of Vector. vd>=1
#define ld 2	//Dimension of Lattice. ld>=1
#define n 200	//No of lattice points along each side
#define latlen 40000	//Total number of lattice points = n^ld .You'll get segmentation fault if you set it wrong.
#define thetamax 3.14/4.0	//Maximum angle while generating a random angle
#define tempmin 1 //min temperature in units of 0.01
#define tempmax 1 //max temperature in units of 0.01
#define tempsteps 5 //temperature steps in units of 0.01
#define maxthreads 1	//Maximum number of parallel threads to use
#define measurelen 5 //Number of sweeps after which measurement has to be taken
#define measurenum 0 //Number of measurements to take
#define equibsweep 10000 //Number of sweeps to run to equilibriate
#define folder "filest"  //Name of the output directory. Create the folder before running
#define j 1		//coupling parameter in the hamiltonian

//Some definitions to make timing easier. Change 'steady_clock' to 'monotonic_clock' for older c++ standards
#define startt starttime = chrono::steady_clock::now()  //start counting
#define endt endtime = chrono::steady_clock::now()   //End counting
#define showt chrono::duration_cast<chrono::milliseconds>(endtime - starttime).count()  //show time elapsed


using namespace std;

//Globally defining random number generators
uniform_real_distribution<double> angle(-thetamax, +thetamax);
uniform_real_distribution<double> metro(0, 1);
uniform_int_distribution<unsigned int> vect(0, 2*vd - 1);
uniform_int_distribution<unsigned int> posi(0, n - 1);
mt19937 rng1, rng2;
#pragma omp threadprivate(vect,posi,angle,metro,rng1,rng2)

/*
//Mersenne twister Psuedo Random Number Generator with seed = 1 . u means the datatype is long
boost::mt19937 rng1(1u);
//threadprivate makes an independent copy of variable rng1 for each thread so that they dont have to share the same memory
#pragma omp threadprivate(rng1)
//Some random variables we'll eventually use
boost::variate_generator<boost::mt19937&, boost::uniform_int<> > vect(rng1, boost::uniform_int<>(0, 2 * vd - 1));
boost::variate_generator<boost::mt19937&, boost::uniform_int<> > posi(rng1, boost::uniform_int<>(0, n - 1));
boost::variate_generator<boost::mt19937&, boost::uniform_real<> > angle(rng1, boost::uniform_real<>(-tmax, +tmax));
boost::uniform_01<boost::mt19937> metro(rng1);
//Making a copy of the random variables for each thread
#pragma omp threadprivate(vect,posi,angle,metro)
*/
unsigned int l(unsigned int* pos)
{
	unsigned int len = 0;
	for (unsigned int i = 0;i < ld;i++)
		len += pos[i] * (unsigned int)pow(n, i);
	return len;

}
double dot(double *vec1, double *vec2)
{
	double temp = 0.0;
	for (unsigned int i = 0;i < 2 * vd;i++)
		temp += vec1[i] * vec2[i];
	return temp;
}

void vectorzero(double* vec)
{
	for (unsigned int i = 0;i < 2 * vd;i++)
		vec[i] = 0;
}
void vectorcpy(double* tempvec, double* vec)
{
	for (unsigned int i = 0; i < 2 * vd; i++)
		tempvec[i] = vec[i];
}
void vectoradd(double* vec1, double* vec2, double val)
{
	for (unsigned int i = 0; i < 2 * vd; i++)
		vec1[i] += val*vec2[i];
}
void vectoradd(double* vec1, double* vec2, double val2, double* vec3, double val3)
{
	for (unsigned int i = 0; i < 2 * vd; i++)
		vec1[i] += val2*vec2[i] + val3*vec3[i];
}
void vectorrandrot(double* tempvec, double* vec)
{

	double theta = angle(rng1);
	unsigned int t1 = vect(rng1), t2 = vect(rng1);

	while (t1 == t2)
	{
		t2 = vect(rng1);
	}
	for (unsigned int i = 0;i < 2 * vd;i++)
	{
		if (i == t1) 	tempvec[i] = sin(theta)*vec[t2] + cos(theta)*vec[t1];
		else if (i == t2) tempvec[i] = -sin(theta)*vec[t1] + cos(theta)*vec[t2];
		else tempvec[i] = vec[i];
	}

}
void vectorrand(double* tempvec)
{
	double abs = 0.0;
	for (unsigned int i = 0;i < 2 * vd;i++)
	{
		tempvec[i] = angle(rng1);
		abs += tempvec[i] * tempvec[i];
	}
	abs = sqrt(abs);
	for (unsigned int i = 0;i < 2 * vd;i++)
	{
		tempvec[i] = tempvec[i] / abs;
	}
}
double vectorabs(double* tempvec)
{
	double abs = 0.0;
	for (unsigned int i = 0;i < 2 * vd;i++)
	{
		abs += tempvec[i] * tempvec[i];
	}
	return sqrt(abs);
}

void poscpy(unsigned int* temppos, unsigned int* pos)
{
	for (unsigned int i = 0; i < ld; i++)
		temppos[i] = pos[i];
}
void posrand(unsigned int* pos)
{
	for (unsigned int i = 0;i < ld;i++)
	{
		pos[i] = posi(rng1);
	}
}

void nearestneighbour(double **lattice, unsigned int pos[ld], double npossum[2 * vd])
{
	unsigned int temppos[ld];
	vectorzero(npossum);
	poscpy(temppos, pos);
	for (unsigned int i = 0;i < ld;i++)
	{
		if (pos[i] == n - 1)
		{
			temppos[i] --;vectoradd(npossum, lattice[l(temppos)], 1.0);
			temppos[i] = 0;vectoradd(npossum, lattice[l(temppos)], 1.0);
		}
		else if (pos[i] == 0)
		{
			temppos[i] ++;vectoradd(npossum, lattice[l(temppos)], 1.0);
			temppos[i] = n - 1;vectoradd(npossum, lattice[l(temppos)], 1.0);
		}
		else
		{
			temppos[i] ++;vectoradd(npossum, lattice[l(temppos)], 1.0);
			temppos[i] -= 2;vectoradd(npossum, lattice[l(temppos)], 1.0);
		}
		temppos[i] = pos[i];
	}
}

void latticeinirand(double** templat)
{
	for (unsigned int i = 0;i < latlen;i++)
	{
		vectorrand(templat[i]);
	}
}
void latticeini1(double** templat)
{
	for (unsigned int i = 0;i < latlen;i++)
	{
		for (unsigned int k = 0;k < 2 * vd;k++)
			templat[i][k] = 1.0;
	}
}

////////////////////////////////////////////////////////////
#define diffcalci \
tem_real = 0;tem_imag = 0;\
for (unsigned int i = 0;i < vd;i++)\
{\
	tem_real += vec[i] * nearvec[i] + vec[vd + i] * nearvec[vd + i];\
	tem_imag += vec[i] * nearvec[vd + i] - vec[vd + i] * nearvec[i];\
}\
	tempenergy+= tem_imag*tem_imag + tem_real*tem_real;while(0)
//////////////////////////////////////////////////////////////

double latticeenergy(double **lattice)
{
	double tempenergy = 0.0;
	double npossum[2 * vd] = {};
	double *nearvec, *vec;
	unsigned int pos[ld] = {};
	unsigned int temppos[ld] = {};
	unsigned int kld = 0;
	double tem_real, tem_imag;
	while (pos[ld - 1] != n)
	{
		kld = 0;
		if (pos[kld] < n)
		{
			poscpy(temppos, pos);
			vec = lattice[l(pos)];
			for (unsigned int k = 0;k < ld;k++)
			{
				if (pos[k] == n - 1)
				{
					temppos[k] --;nearvec = lattice[l(temppos)];diffcalci;
					temppos[k] = 0;nearvec = lattice[l(temppos)];diffcalci;
				}
				else if (pos[k] == 0)
				{
					temppos[k] ++;nearvec = lattice[l(temppos)];diffcalci;
					temppos[k] = n - 1;nearvec = lattice[l(temppos)];diffcalci;
				}
				else
				{
					temppos[k] ++;nearvec = lattice[l(temppos)];diffcalci;
					temppos[k] -= 2;nearvec = lattice[l(temppos)];diffcalci;
				}
				temppos[k] = pos[k];
			}
			pos[kld]++;
		}
		else
			while (pos[kld] == n && pos[ld - 1] != n)
			{
				pos[kld] = 0;
				kld++;
				pos[kld]++;
			}
	}
	return j*tempenergy / 2.0;
}
void latticemag(double **lattice, double* mag)
{
	vectorzero(mag);
	for (unsigned int i = 0;i < latlen;i++)
		vectoradd(mag, lattice[i], 1.0);
}
double latticemagabs(double **lattice)
{
	double mag[2 * vd] = {};
	for (unsigned int i = 0;i < latlen;i++)
		vectoradd(mag, lattice[i], 1.0);
	return vectorabs(mag);
}
void latticeexport(double** lattice, double t)
{
	ofstream lwrite;
	lwrite.open("./"folder"/cpnlat_ld" + to_string((long long)ld) + "_vd" + to_string((long long)vd) + "_n" + to_string((long long)n) + "_t" + to_string((long double)t) + ".csv");
	for (unsigned int i1 = 0;i1 < latlen;i1++)
	{
		lwrite << lattice[i1][0];
		unsigned int i2 = 1;
		while (i2 < 2 * vd)
		{
			lwrite << "," << lattice[i1][i2];
			i2++;
		}
		lwrite << "\n";
	}
	lwrite.close();
}
void latticecopy(double** newlattice, double** oldlattice)
{
	for (unsigned int i = 0;i < latlen;i++)
		for (unsigned int k = 0;k < 2 * vd;k++)
			newlattice[i][k] = oldlattice[i][k];
}

////////////////////////////////////////////////////////////
#define diffcalc \
tem_real = 0;tem_imag = 0;tem1_real = 0;tem1_imag = 0;\
for (unsigned int i = 0;i < vd;i++)\
{\
	tem_real += vec[i] * nearvec[i] + vec[vd + i] * nearvec[vd + i];\
	tem_imag += vec[i] * nearvec[vd + i] - vec[vd + i] * nearvec[i];\
	tem1_real += tempvec[i] * nearvec[i] + tempvec[vd + i] * nearvec[vd + i];\
	tem1_imag += tempvec[i] * nearvec[vd + i] - tempvec[vd + i] * nearvec[i];\
}\
ediff += tem1_imag*tem1_imag + tem1_real*tem1_real - tem_imag*tem_imag - tem_real*tem_real;while(0)
//////////////////////////////////////////////////////////////

double runmcmc(double **lattice, unsigned int sweeps, double t, int measure)
{
	unsigned int pos[ld], temppos[ld];
	double acceptance = 0.0;
	double tempvec[2 * vd], ediff, energy, mag[2 * vd];
	ofstream fileenergy;
	//ofstream filemag;
	double tem_real, tem_imag, tem1_real, tem1_imag;
	double *vec, *nearvec;


	//Initializing energy and magnetization files
	if (measure > 1)
	{
		energy = latticeenergy(lattice);
		latticemag(lattice, mag);

		fileenergy.open("./" folder  "/cpnenergy_ld" + to_string((long long)ld) + "_vd" + to_string((long long)vd) + "_n" + to_string((long long)n) + "_t" + to_string((long double)t) + ".csv");
		//filemag.open("./" folder "/cpnmag_ld" + to_string((long long)ld) + "_vd" + to_string((long long)vd) + "_n" + to_string((long long)n) + "_t" + to_string((long double)t) + ".csv");
	}

	for (unsigned int i1 = 0;i1 < sweeps;i1++)
	{
		for (unsigned int i2 = 0;i2< latlen;i2++)
		{
			posrand(pos);


			vec = lattice[l(pos)];
			vectorrandrot(tempvec, vec);

			ediff = 0;
			poscpy(temppos, pos);
			for (unsigned int k = 0;k < ld;k++)
			{
				if (pos[k] == n - 1)
				{
					temppos[k] --;nearvec = lattice[l(temppos)];diffcalc;
					temppos[k] = 0;nearvec = lattice[l(temppos)];diffcalc;
				}
				else if (pos[k] == 0)
				{
					temppos[k] ++;nearvec = lattice[l(temppos)];diffcalc;
					temppos[k] = n - 1;nearvec = lattice[l(temppos)];diffcalc;
				}
				else
				{
					temppos[k] ++;nearvec = lattice[l(temppos)];diffcalc;
					temppos[k] -= 2;nearvec = lattice[l(temppos)];diffcalc;
				}
				temppos[k] = pos[k];
			}
			ediff *= j;



			if (metro(rng1) < exp(-ediff / t))
			{
				if (measure > 1)
				{
					energy += ediff;
					//vectoradd(mag, vec, -1.0, tempvec, +1.0);
				}
				vectorcpy(vec, tempvec);
				acceptance += 1.0;
			}
		}

		if (measure > 1 && i1 % measurelen == 0)
		{
			fileenergy << energy << "\n";
			/*filemag << mag[0];
			for (unsigned int i = 1;i < 2*vd;i++)
			filemag << "," << mag[i];
			filemag << "\n";*/
		}
	}
	//filemag.close();
	fileenergy.close();
	return acceptance / (sweeps*latlen);
}



int main()
{
	omp_set_num_threads(maxthreads);
	auto starttimehead = chrono::steady_clock::now();
#pragma omp parallel 
		{
			
			auto starttime = chrono::steady_clock::now();
			auto endtime = chrono::steady_clock::now();
			double acc;
			//Dynamically allocating array to store the lattice (each thread has its own separate copy)
			double **lattice = new double*[latlen];
			for (unsigned int i = 0; i < latlen; ++i)
			{
				lattice[i] = new double[2 * vd];
			}

#pragma omp for  
			for (int t = tempmin;t <= tempmax;t += tempsteps)
			{

				startt;
				latticeinirand(lattice);
				acc = runmcmc(lattice, equibsweep, 0.001*(double)t, 0);
				endt;
				printf("Equilibration:: t=%f__acc=%f__Time=%I64d__Thread=%d\n", 0.001*(double)t, acc, showt, omp_get_thread_num());
				
				latticeexport(lattice, 0.01*(double)t);

				startt;
				acc = runmcmc(lattice,measurelen*measurenum, 0.001*(double)t, 2);
				endt;
				printf("Measurement:: t=%f__acc=%f__Time=%I64d__Thread=%d\n", 0.001*(double)t, acc, showt, omp_get_thread_num());
		
			}

			//Deleting dynamically allocated array storing the lattice
			for (unsigned int i = 0; i < latlen; ++i)
			{
				delete[] lattice[i];
			}
			delete[] lattice;
		}
		//latticeinirand(lattice);
		//runmcmc(lattice, 300, 0.01, 0);
		//latticeexport(lattice,0.01);
		//endt;
		auto endtimehead = chrono::steady_clock::now();
		cout << "Total Time:" << chrono::duration_cast<chrono::milliseconds>(endtimehead - starttimehead).count() << endl;
	cout << "The End";
	return 0;
}







