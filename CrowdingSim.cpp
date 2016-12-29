// CrowdingSim.cpp : Defines the entry point for the console application.
//

//#include "stdafx.h"



//  Paul Mlynarczyk
//	Modified by Charles Chin
//  This program simulates particles reacting and diffusing in a biochemical reaction network defined by a simple 2 particle bursting system and mRNA productiion
//  The reaction volume is 3-dimensional
//  Modifed to have transcription in a crowded environment

//
//  TODO
//	change reaction scheme to 2 particles diffusing, making mRNA when they meet
//	Add crowder functionality

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <fstream>
#include <ctime>
#include <algorithm>
#include <vector>
#include <array>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <functional>
#include <numeric>
//#include <Windows.h>
#include <ctime>
#include <sstream>

using namespace std;

mt19937::result_type seed = (unsigned)time(0);
auto real_rand = std::bind(std::uniform_real_distribution<double>(0, 1), mt19937(seed));

class Reactor
{
public:
	static int const xmax = 16;  //reactor dimensions
	static int const ymax = 16;
	static int const zmax = 16;
	double h = 0.05;	//side length
	double crowdingFrac = 0; //crowding fraction by volume
	double t0;  //initial time
	double dt;  //sampling time
	double tfinal;  //final time
	double km, gm, k_hop, k_hopC;  //rate constants
	int nG0, nTF0;  //number of total particles and initial number of X particles
	double a[xmax][ymax][zmax][5]; //propensity functions by voxel
	double amn[xmax][ymax][zmax];  //total propensity by voxel
	float a0;  //total propensity
	int Gnv[xmax][ymax][zmax]; //A concentration by voxel
	int TFnv[xmax][ymax][zmax]; //X concentration by voxel
	int mRNAnv[xmax][ymax][zmax]; //mRNA population total
	int Cnv[xmax][ymax][zmax]; //crowders by voxel
	int TFCurrent[3]; //TF tracking for encounters
	int GCurrent[3]; //G tracking for encounters
	double encT;
	double encD;

	ofstream myfileG;   //declare output files
	ofstream myfileTF;
	ofstream myfileC;
	ofstream myfilemRNA;
	ofstream myfilemRNAsum;
	ofstream myfileEncT;
	ofstream myfileEncD;
	ofstream myfileTime;
	void setRates(double, double, double, double);
	int sumTotalX(void);
	void populateParticles(int, int, char);
	void calcPropensities(int, int, int);
	void setTimes(double, double, double);
	void writeOut(void);
}  R;

void Reactor::setRates(double kmin, double gmin, double k_hop_in, double k_hopC_in)
{
	double volume = (xmax*ymax*zmax)*pow(h, 3);
	km = kmin*(volume / pow(h, 3));
	gm = gmin;
	k_hop = (6 * k_hop_in) / pow(h, 2);
	k_hopC = (6 * k_hopC_in) / pow(h, 2);
}

void Reactor::populateParticles(int nXinit, int nptot, char pos)
{
	for (int b = 0; b<xmax; b++) {   //initialize populations to 0
		for (int c = 0; c<ymax; c++) {
			for (int d = 0; d<zmax; d++) {
				Gnv[b][c][d] = 0;
				TFnv[b][c][d] = 0;
				mRNAnv[b][c][d] = 0;
				Cnv[b][c][d] = 0;
			}
		}
	}

	nTF0 = nXinit;
	nG0 = nptot;
	//int nA0 = np - nX0;

	if (crowdingFrac > 0) { // populate crowders
		int latticeSites = xmax*ymax*zmax;
		double crowdTot = (crowdingFrac / 100)*latticeSites;
		int i = 0;
		while (i < crowdTot) {
			int rx = floor(real_rand()*xmax);
			int ry = floor(real_rand()*ymax);
			int rz = floor(real_rand()*zmax);
			if (Cnv[rx][ry][rz] == 0) {
				Cnv[rx][ry][rz] += 1;
				i++;
			}
		}
	}

	int i = 0;
	while (i < nG0) { //populate initial A particles
		int rx = floor(real_rand()*xmax);
		int ry = floor(real_rand()*ymax);
		int rz = floor(real_rand()*zmax);
		if (Cnv[rx][ry][rz] == 0) {
			Gnv[rx][ry][rz] = Gnv[rx][ry][rz] + 1;   //increment [A]
			GCurrent[1] = rx;
			GCurrent[2] = ry;
			GCurrent[3] = rz;

			i++;
		}

	}

	i = 0;
	while (i < nTF0) { //populate initial X particle
		int rx = floor(real_rand()*xmax);
		int ry = floor(real_rand()*ymax);
		int rz = floor(real_rand()*zmax);
		if (Cnv[rx][ry][rz] == 0) {
			TFnv[rx][ry][rz] = TFnv[rx][ry][rz] + 1;   //increment [X]
			TFCurrent[1] = rx;
			TFCurrent[2] = ry;
			TFCurrent[3] = rz;
			i++;
		}

	}


}


void Reactor::calcPropensities(int x, int y, int z)
{

	a[x][y][z][0] = km*Gnv[x][y][z] * TFnv[x][y][z];  //compute propensities in voxel (x,y,z)
	a[x][y][z][1] = gm*mRNAnv[x][y][z];
	a[x][y][z][2] = k_hop*(Gnv[x][y][z] + TFnv[x][y][z]);
	a[x][y][z][3] = k_hopC*Cnv[x][y][z];
	amn[x][y][z] = a[x][y][z][0] + a[x][y][z][1] + a[x][y][z][2] + a[x][y][z][3];
}


void Reactor::setTimes(double tstart, double tsample, double tend)
{
	t0 = tstart;
	dt = tsample;
	tfinal = tend;
}


int Reactor::sumTotalX()
{
	int totalX = 0;   //declare and initialize variable to store total number of X particles
	for (int b = 0; b<xmax; b++) {     //loop over all voxels
		for (int c = 0; c<ymax; c++) {
			for (int d = 0; d<zmax; d++) {
				totalX += TFnv[b][c][d];   //increment total by population in each voxel
			}
		}
	}
	return totalX;
}

void Reactor::writeOut()
{
	for (int g = 0; g<xmax; g++) {
		for (int h = 0; h<ymax; h++) {
			for (int i = 0; i<zmax; i++) {
				myfileG << Gnv[g][h][i] << "    ";
			}
			myfileG << endl;
		}
		myfileG << endl;
	}
	myfileG << endl;

	for (int g = 0; g<xmax; g++) {
		for (int h = 0; h<ymax; h++) {
			for (int i = 0; i<zmax; i++) {
				myfileTF << TFnv[g][h][i] << "    ";
			}
			myfileTF << endl;
		}
		myfileTF << endl;
	}
	myfileTF << endl;

	for (int g = 0; g<xmax; g++) {
		for (int h = 0; h<ymax; h++) {
			for (int i = 0; i<zmax; i++) {
				myfileC << Cnv[g][h][i] << "    ";
			}
			myfileC << endl;
		}
		myfileC << endl;
	}
	myfileC << endl;

	int mRNAsum = 0;
	for (int g = 0; g<xmax; g++) {
		for (int h = 0; h<ymax; h++) {
			for (int i = 0; i<zmax; i++) {
				myfilemRNA << mRNAnv[g][h][i] << "    ";
				mRNAsum += mRNAnv[g][h][i];
			}
			myfilemRNA << endl;
		}
		myfilemRNA << endl;
	}
	myfilemRNA << endl;



	myfilemRNAsum << mRNAsum << "    ";
	myfilemRNAsum << endl;
}

int main(int argc, char* const argv[])
{
	clock_t start = clock();
	int SGE_TASK_ID = atoi(argv[1]);
	mt19937::result_type seed = ((unsigned)time(NULL)*SGE_TASK_ID);
	auto real_rand = std::bind(std::uniform_real_distribution<double>(0, 1), mt19937(seed));
	R.a0 = 0;  //initial total propensity
	R.setTimes(0, 1, 1000);
	R.setRates(10, 1, 60, 0.6);
	cout << R.km << endl;
	//R.setSize(100,100,1);  //set reaction volume grid dimensions
	R.populateParticles(1, 1, 's');
	int tpoints = (R.tfinal - R.t0) / R.dt + 1;  //number of time points
	double t[1000];  //vector of time points //changed for compiler
	t[0] = R.t0;
	for (int i = 0; i<tpoints - 1; i++)
	{
		t[i + 1] = t[i] + R.dt;
	}

	int meet = 0; //initialize meeting var for encounter tracking
	int dir = 0;  //initialize moving direction of chosen particle
	int m = 0;   //initialize x,y,z position of target voxel
	int n = 0;
	int p = 0;
	int v[6][3] = { { 1,0,0 },{ -1,0,0 },{ 0,1,0 },{ 0,-1,0 },{ 0,0,1 },{ 0,0,-1 } };  //movement position change matrix
	int vj[2][1] = { 1,-1 }; //stoich coeff matrix

	stringstream ssmRNA;
	stringstream encTimes;
	stringstream encDur;
	stringstream elapsedTime;
	(R.myfileG).open("G_output.txt");   //open output files
	(R.myfileTF).open("TF_output.txt");
	(R.myfileC).open("C_output.txt");
	(R.myfilemRNA).open("mRNA_output.txt");

	ssmRNA << "/data/home/cchin/mRNA_Run" << R.xmax << "x" << R.ymax << "x" << R.zmax << "D6" << "Crowd" << R.crowdingFrac << "Run" << SGE_TASK_ID << ".txt";
	string filenameX = ssmRNA.str();
	(R.myfilemRNAsum).open(filenameX);
	encTimes << "/data/home/cchin/encT_Run" << R.xmax << "x" << R.ymax << "x" << R.zmax << "D6" << "Crowd" << R.crowdingFrac << "Run" << SGE_TASK_ID << ".txt";
	string filenameeT = encTimes.str();
	(R.myfileEncT).open(filenameeT);
	encDur << "/data/home/cchin/encD_Run" << R.xmax << "x" << R.ymax << "x" << R.zmax << "D6" << "Crowd" << R.crowdingFrac << "Run" << SGE_TASK_ID << ".txt";
	string filenameeD = encDur.str();
	(R.myfileEncD).open(filenameeD);
	elapsedTime << "/data/home/cchin/Time_Run" << R.xmax << "x" << R.ymax << "x" << R.zmax << "D6" << "Crowd" << R.crowdingFrac << "Run" << SGE_TASK_ID << ".txt";
	string filenametime = elapsedTime.str();
	(R.myfileTime).open(filenametime);
	//(R.myfilemRNAsum).open("mRNAsum_output.txt");


	int ct = 1;  //current time point
	double tn = R.t0;  //current time

	for (int j = 0; j<R.xmax; j++) {   //calculate initial propensities for all voxels
		for (int k = 0; k<R.ymax; k++) {
			for (int l = 0; l<R.zmax; l++) {
				R.calcPropensities(j, k, l);
				R.a0 += R.amn[j][k][l];   //increment total propensity by voxel propensity
			}
		}
	}
	R.writeOut();

	while (tn < R.tfinal) {

		double rtime = real_rand();
		double rvox = real_rand();
		double ra = real_rand();
		double rmove = real_rand();

		double tau = (1.0 / (R.a0))*log(1.0 / rtime);  //time-step increment
		tn += tau;   //increment current time by tau
		double ajk = 0.0;
		int flag = 0;
		int j, k, l;
		for (j = 0; j<R.xmax; j++) {   //select voxel (j,k,l)
			for (k = 0; k<R.ymax; k++) {
				for (l = 0; l<R.zmax; l++) {
					ajk += R.amn[j][k][l];
					if (ajk > rvox*(R.a0)) {
						flag = 1;
						break;
					}
				}
				if (flag == 1) {
					break;
				}
			}
			if (flag == 1) {
				break;
			}
		}
		int q = 0;
		double ajkl = R.a[j][k][l][q];
		while (ajkl <= ra*(R.amn[j][k][l])) {
			q += 1;
			ajkl += R.a[j][k][l][q];
		}
		if (q == 0) {  //reaction chosen
			R.mRNAnv[j][k][l] += 1;// vj[q][0];
								   //R.Anv[j][k][l] += vj[q][0];
								   //R.Xnv[j][k][l] += vj[q][1];
		}
		else if (q == 1) {
			R.mRNAnv[j][k][l] -= 1;
		}
		else if (q == 2) {  //particle move chosen
			int reject = 0;
			dir = floor(6 * rmove);
			m = j + v[dir][0];  //x-index of target voxel
			n = k + v[dir][1];  //y-index of target voxel
			p = l + v[dir][2];  //z-index of target voxel

			if (R.Cnv[m][n][p] == 1) { // if crowder at target voxel
				reject = 1;
			}
			if (m > R.xmax - 1)  reject = 1;  //hard wall boundary conditions
			if (m < 0) reject = 1;
			if (n > R.ymax - 1) reject = 1;
			if (n < 0) reject = 1;
			if (p > R.zmax - 1) reject = 1;
			if (p < 0) reject = 1;

			if (reject == 0) {
				double probA = (double)R.Gnv[j][k][l] / (R.Gnv[j][k][l] + R.TFnv[j][k][l]);
				if (real_rand() < probA) {  //if A particle
					R.Gnv[j][k][l] -= 1;  //decrement conc of G in vox (j,k)
					R.Gnv[m][n][p] += 1;  //increment conc of G in vox (m,n)
					R.GCurrent[1] = m;
					R.GCurrent[2] = n;
					R.GCurrent[3] = p;
				}
				else {   //if X particle
					R.TFnv[j][k][l] -= 1;  //decrement conc of TF in vox (j,k)
					R.TFnv[m][n][p] += 1;  //increment conc of TF in vox (m,n)
					R.TFCurrent[1] = m;
					R.TFCurrent[2] = n;
					R.TFCurrent[3] = p;
				}

				//calculate encounter properties
				if (R.GCurrent[1] == R.TFCurrent[1]) {
					if (R.GCurrent[2] == R.TFCurrent[2]) {
						if (R.GCurrent[3] == R.TFCurrent[3]) {
							if (meet == 0) {
								meet = 1;
								R.encT = tn;
								//write out encounter time
								R.myfileEncT << R.encT << endl;
							}
						}
					}
				}
				else if (meet == 1) {
					meet = 0;
					//write out encounter duration
					R.encD = tn - R.encT;
					R.myfileEncD << R.encD << endl;
				}

				R.a0 -= R.amn[m][n][p];
				R.calcPropensities(m, n, p);  //update propensities of modified voxel
				R.a0 += R.amn[m][n][p];
			}

		}
		else {	//crowder move chosen
			int reject = 0;
			dir = floor(6 * rmove);
			m = j + v[dir][0];  //x-index of target voxel
			n = k + v[dir][1];  //y-index of target voxel
			p = l + v[dir][2];  //z-index of target voxel

			if (R.Cnv[m][n][p] == 1) { // if crowder at target voxel
				reject = 1;
			}
			if (m > R.xmax - 1)  reject = 1;  //hard wall boundary conditions
			if (m < 0) reject = 1;
			if (n > R.ymax - 1) reject = 1;
			if (n < 0) reject = 1;
			if (p > R.zmax - 1) reject = 1;
			if (p < 0) reject = 1;

			if (reject == 0) {
				R.Cnv[j][k][l] -= 1;  //decrement conc of C in vox (j,k)
				R.Cnv[m][n][p] += 1;  //increment conc of C in vox (m,n)

				R.a0 -= R.amn[m][n][p];
				R.calcPropensities(m, n, p);  //update propensities of modified voxel
				R.a0 += R.amn[m][n][p];
			}
		}
		R.a0 -= R.amn[j][k][l];
		R.calcPropensities(j, k, l);  //update propensities of modified voxel
		R.a0 += R.amn[j][k][l];

		while (tn >= ct*(R.dt)) {  //sampled time
			R.writeOut();
			cout << tn << endl;
			ct += 1;
		}
	}
	cout << endl;

	clock_t end = clock();
	double time = (double)(end - start) / CLOCKS_PER_SEC;
	R.myfileTime << time << endl;

	(R.myfileG).close();
	(R.myfileTF).close();
	(R.myfileC).close();
	(R.myfilemRNA).close();
	(R.myfilemRNAsum).close();
	(R.myfileEncT).close();
	(R.myfileEncD).close();
	(R.myfileTime).close();
	cout << "Done" << endl;
	//OutputDebugStringA("done");

	cout << "Time Elapsed: " << time << endl;
	return 0;   //indicate successful termination
}
