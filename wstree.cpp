/***************************************************

	wstree.cpp

	A 3D watershed transform implementation with
	emphasis on the heirarchy of voids

	by Devon Powell

***************************************************/

#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdint>
#include <cfloat>
#include "HDF_IO.hh"

// custom macros and typedefs
using namespace std;
typedef uint32_t uint;
#define uint_max UINT32_MAX

// forward declarations
void read_hdf5(string filename, string fieldname);
void write_hdf5(string filename, string fieldname);
void argsort();
void watershed();

// global data arrays
uint nx, ny, nz, ntot;
vector<double> field;
vector<uint> inds_sorted;
uint nzones;
vector<uint> zones; 

int main(int argc, char **argv) {

	// we want four arguments
	if(argc != 5) {
		printf("-------------------------------------\n");
		printf(" Usage:\n");
		printf("   ./wstree input.hdf5 input_field output.hdf5 output_field\n");
		printf(" Example:\n");
		printf("   ./wstree data/dset128.hdf5 RHO output/ws128.hdf5 WS\n");
		printf("-------------------------------------\n");
		return 0;
	}

	printf("-------------------------------------\n");
	read_hdf5(argv[1], argv[2]);
	printf("-------------------------------------\n");
	argsort();
	printf("-------------------------------------\n");
	watershed();
	printf("-------------------------------------\n");
	write_hdf5(argv[3], argv[4]);
	printf("-------------------------------------\n");

	return 0;    
}


void watershed() {

	printf(" Running watershed transform...");

	nzones = 0;
	zones.assign(ntot, uint_max); // unassigned zones use uint_max

	for(uint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		uint ind_flat = inds_sorted[ind_uns];


		// get 3D indices from ind_flat
		uint ix = ind_flat/(ny*nz);
		uint iy = (ind_flat - ix*ny*nz)/nz;
		uint iz = ind_flat - ix*ny*nz - iy*nz;

		// iterate over the 27 neighboring cells
		double f0 = field[ind_flat];
		//double dmin = f0; 
		double grad_max = 0.0; 
		uint zmin = zones[ind_flat];// = uint_max;
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					uint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f0 - field[tmp_flat])/sqrt(ox*ox + oy*oy + oz*oz);

					//if(zones[tmp_flat] < zmin) { // This disambiguation has an inherent bias towards deeper voids!
					//if(field[tmp_flat] < dmin) { // finds the neighboring cell with the lowest density 
					if(grad > grad_max) { // finds the largest gradient to a neighboring cell 
					
						//dmin = field[tmp_flat];
						grad_max = grad;
						
						zmin = zones[tmp_flat];
						zones[ind_flat] = zmin;
					}
				}
			}
		}
		if(zmin == uint_max) {
			zones[ind_flat] = nzones++;
		}
	}

	printf(" done.\n");
	printf("   Found %u distinct zones.\n", nzones);

}

// read in a multidimensional hdf5 file
void read_hdf5(string filename, string fieldname) {
	printf(" Reading field %s from file %s...", fieldname.c_str(), filename.c_str());
	vector<int> dims;
	HDFGetDatasetExtent(filename, fieldname, dims);
	nx = dims[0]; ny = dims[1]; nz = dims[2];
	ntot = nx*ny*nz;
	HDFReadDataset(filename, fieldname, field);
	printf(" done.\n");
	printf("   Read %d data values.\n", (int)field.size());
	printf("   nx = %d, ny = %d, nz = %d for %d total.\n", nx, ny, nz, ntot);
	return;
}

void write_hdf5(string filename, string fieldname) {
	printf(" Writing field %s to file %s...", fieldname.c_str(), filename.c_str());
	HDFCreateFile(filename);
	uint dims[3] = {nx, ny, nz};
	HDFWriteDataset3D(filename, fieldname, dims, zones);
	printf(" done.\n");
}

bool indcmp(int a, int b) { // helper function for argsort
	return field[a] < field[b];
}
void argsort() {
	printf(" Sorting array indices...");
	inds_sorted.resize(ntot);
	uint itmp = 0;
	generate(inds_sorted.begin(), inds_sorted.end(), [&] { return itmp++; });
	sort(inds_sorted.begin(), inds_sorted.end(), indcmp);
	printf(" done.\n");
}

