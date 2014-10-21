/***************************************************

	wstree.cpp

	A 3D watershed transform implementation with
	emphasis on the heirarchy of voids

	by Devon Powell

***************************************************/

#include <cstdio>
//#include <cmath>
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
void write_hdf5();
void argsort();
void watershed();

// global data arrays
uint nx, ny, nz, ntot;
vector<double> field;
vector<uint> inds_sorted;
uint nzones;
vector<uint> zones; 

int main(int argc, char **argv) {

	// we want two arguments
	/*if(argc != 3) {
		printf("Requires one input file and one output file!\n");
		exit(0);
	}

	printf("\n"); 

	*/

	//readInput(argv[1]);
	printf("-------------------------------------\n");
	read_hdf5("data/dset128.hdf5", "RHO");
	printf("-------------------------------------\n");
	argsort();
	printf("-------------------------------------\n");
	watershed();
	printf("-------------------------------------\n");

	//write_hdf5();

	return 0;    
}


void watershed() {

	printf("Running watershed transform...");

	nzones = 0;
	zones.assign(ntot, uint_max); // unassigned zones use uint_max

	for(uint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		uint ind_flat = inds_sorted[ind_uns];
		//printf("%u\t%f\n", ind_flat, field[ind_flat]);	

		// get 3D indices from ind_flat
		uint i_x = ind_flat%(ny*nz);
		uint i_y = (ind_flat - i_x*ny*nz)%nz;
		uint i_z = ind_flat - i_x*ny*nz - i_y*nz;

		// iterate over the 27 neighboring cells
		double dmin = DBL_MAX; 
		uint zmin = uint_max;
		for(int o_x = -1; o_x <= 1; ++o_x) {
			for(int o_y = -1; o_y <= 1; ++o_y) {
				for(int o_z = -1; o_z <= 1; ++o_z) {

					// get neighboring flat indices, accounting for periodicity
					uint tmp_flat = ((i_x + o_x + nx)%nx)*ny*nz + ((i_y + o_y + ny)%ny)*nz + (i_z + o_z + nz)%nz;

					// a neighboring zone has been assigned
					if(zones[tmp_flat] < zmin) {
						// TODO: This disambiguation has an inherent bias towards deeper voids!
						zmin = zones[tmp_flat];
						zones[ind_flat] = zmin;
					}
				}
			}
		}
		if(zmin == uint_max) {
			zones[ind_flat] = ++nzones;
		}
		
	}

	printf(" done.\n");
	printf("  Found %u distinct zones.\n", nzones);

}

void write_hdf5() {
/*
	char* outfilename = "output/test.hdf5";

	printf("Writing to file: %s...\n", outfilename);
	vector<float> pos_v(pos_h, pos_h + npart);
	vector<float> vel_v(vel_h, vel_h + npart);
	vector<float> phi_v(phi_h, phi_h + ngrid);
	vector<float> rho_v(rho_h, rho_h + ngrid);

	HDFCreateFile(outfilename);
	HDFWriteDataset(outfilename, "POS", pos_v);
	HDFWriteDataset(outfilename, "VEL", vel_v);
	HDFWriteDataset(outfilename, "PHI", phi_v);
	HDFWriteDataset(outfilename, "RHO", rho_v);

	printf("...done.\n");
*/
}

// read in a multidimensional hdf5 file
void read_hdf5(string filename, string fieldname) {
	printf("Reading field %s from file %s...", fieldname.c_str(), filename.c_str());
	vector<int> dims;
	HDFGetDatasetExtent(filename, fieldname, dims);
	nx = dims[0]; ny = dims[1]; nz = dims[2];
	ntot = nx*ny*nz;
	HDFReadDataset(filename, fieldname, field);
	printf(" done.\n");
	printf("  Read %d data values.\n", field.size());
	printf("  nx = %d, ny = %d, nz = %d for %d total.\n", nx, ny, nz, ntot);
	return;
}


// helper function for argsort
bool indcmp(int a, int b) {
	return field[a] < field[b];
}
void argsort() {
	printf("Sorting array indices...");
	inds_sorted.resize(ntot);
	int i_tmp = 0;
	generate(inds_sorted.begin(), inds_sorted.end(), [&] { return i_tmp++; });
	sort(inds_sorted.begin(), inds_sorted.end(), indcmp);
	printf(" done.\n");
}

