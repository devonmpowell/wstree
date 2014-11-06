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
#include <cmath>
#include <unordered_map>
#include <set>
#include "HDF_IO.hh"

// custom macros and typedefs
using namespace std;
typedef uint32_t luint;
#define LUINT_MAX UINT32_MAX

// a POD struct for zone info
typedef struct {
	double vol;
	double mass;
	double min;
	unordered_map<luint, double> neighbors; // key is the neighbor zone ID; value is minimum barrier between neighbors
} Zone;

// a POD struct for void info
typedef struct {
	double vol;
	double mass;
	double min;
	set<luint> subzones; // the void's constituent zones
} Void;


// forward declarations
void read_hdf5(string filename, string field_cubename);
void write_hdf5(string filename, string field_cubename);
void argsort();
void watershed();
void neighbors();
void merge_zobov();
void write_tree(string filename);

// global data arrays
// volumetric data
luint nx, ny, nz, ntot;
vector<double> field_cube;
vector<luint> inds_sorted;
vector<luint> zone_cube; 

// zones (raw zones)
luint nzones;
vector<Zone> zones;

// voids (merged zones)
luint nvoids;
vector<Void> nvoids;

int main(int argc, char **argv) {

	// we want four arguments
	if(argc != 5) {
		printf("-------------------------------------\n");
		printf(" Usage:\n");
		printf("   ./wstree input.hdf5 input_field_cube output.hdf5 output_field_cube\n");
		printf(" Example:\n");
		printf("   ./wstree data/dset128.hdf5 RHO output/ws128.hdf5 WS\n");
		printf("-------------------------------------\n");
		return 0;
	}

	printf("-------------------------------------\n");
	printf(" Reading field_cube %s from file %s...", argv[2], argv[1]);
	read_hdf5(argv[1], argv[2]);
	printf(" done.\n");
	printf("   Read %d data values.\n", (int)field_cube.size());
	printf("   nx = %d, ny = %d, nz = %d for %d total.\n", nx, ny, nz, ntot);
	printf("-------------------------------------\n");
	printf(" Sorting array indices...");
	argsort();
	printf(" done.\n");
	printf("-------------------------------------\n");
	printf(" Running watershed transform...");
	watershed();
	printf(" done.\n");
	printf("   Found %u distinct zones.\n", nzones);
	printf("-------------------------------------\n");
	printf(" Writing field_cube %s to file %s...", argv[4], argv[3]);
	write_hdf5(argv[3], argv[4]);
	printf(" done.\n");
	printf("-------------------------------------\n");

	printf("  Finding neighbors...\n");
	neighbors();

	double nav = 0.0;
	for(luint z = 0; z < nzones; ++z) {
		nav += zones[z].neighbors.size();
	}
	printf("Average number of zone neighbors: %f\n", nav/nzones);

	printf("  Merging zones...\n")
	merge_zobov();


	//printf(" Building heirarchy...");
	//tree();

	//for(luint z = 0; z < 20; ++z) {
		//printf(" Zone %u:\n", z);
		//printf("   vol = %f\tmass = %f\tfavg = %f\tfmin = %f\n", zones[z].vol, zones[z].mass, zones[z].favg, zones[z].fmin);
	//}
	printf("nzones = %u\n", nzones);

	write_tree("output/tree128.dat");

	printf(" done.\n");
	printf("-------------------------------------\n");

	return 0;    
}

void write_tree(string filename) {

	FILE* file = fopen(filename.c_str(), "w");
	printf("\tWriting to file: %s...", filename.c_str());

	if(file != NULL) {

		fprintf(file, "###########################################################################\n");
		fprintf(file, "#    TREE											                     #\n");
		fprintf(file, "###########################################################################\n");

/*typedef struct {
//	luint zone_id;
	luint subzone_0, subzone_1;
	double vol;
	double mass;
	double fmin;
//	double favg;

	double barrier;
//	luint barrier_neighbor;

} Zone;
*/


		fprintf(file, "# sub0\tsub1\tvol\tmass\tfmin\tbar\n");

		for(luint z = 0; z < nzones; ++z) {
			Zone ws = zones[z];
			//fprintf(file, "%011u\t%011u\t%.11f\t%.11f\t%.11f\t%.11f\n",
			fprintf(file, "%lf\t%lf\t%lf\n",
					ws.vol, ws.mass, ws.min);//ws.barrier);
	
		}

		fclose(file);
	}
	else {
		printf("\n\tFailed to open file.\n");
	}

}

void mhlp(luint vcur, luint vhome) {


	

}

void merge_zobov() {

	// Zobov-style void merging
	nvoids = nzones;
	voids.resize(nvoids);

	// initialize void information
		
	for(luint v = 0; v < nvoids; ++v) {
		voids[v].vol = zones[v].vol;
		voids[v].mass = zones[v].mass;
		voids[v].min = zones[v].min;
		voids[v].subzones.insert(v);
	}

	// traverse backwards, taking advantage of the fact that voids are ordered by
	for(luint v = nvoids - 1; v >= 0; --v) mhlp(v, v);



//	for(luint z = 0; z < nzones; ++z)
//		zones[z].favg = zones[z].mass/zones[z].vol;

/*	
	// find barrier saddle points and merge zones
	for(luint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		luint ind_flat = inds_sorted[ind_uns];
		luint z = zones[ind_flat];
		double f0 = field_cube[ind_flat];

		// get 3D indices from ind_flat
		luint ix = ind_flat/(ny*nz);
		luint iy = (ind_flat - ix*ny*nz)/nz;
		luint iz = ind_flat - ix*ny*nz - iy*nz;

		// inspect the 26 neighboring cells
		double fnmin = DBL_MAX;
		luint znmin = LUINT_MAX;
		double grad_min = DBL_MAX; 
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					luint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					luint z_neighbor = zones[tmp_flat];
					double f_neighbor = field_cube[tmp_flat];

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f_neighbor - f0)/sqrt(ox*ox + oy*oy + oz*oz);

					// check for a neighboring watershed with the shallowest gradient
					if(z != z_neighbor && grad < grad_min) {	

						grad_min = grad;

						fnmin = f_neighbor;
						znmin = z_neighbor;
					}
				}
			}
		}
		// merge the two zones
		if(znmin != LUINT_MAX) {

			// create a new watershed from the union of the two
			Zone ws_new;
			luint z0, z1;
			if(zones[z].fmin < zones[znmin].fmin) {
				z0 = z; z1 = znmin;
			}
			else {
				z0 = znmin; z1 = z;
			}
			zones[z0].barrier = fnmin;
			zones[z1].barrier = fnmin;

			ws_new.subzone0 = z0;
			ws_new.subzone1 = z1;
			ws_new.vol = zones[z0].vol + zones[z1].vol;
			ws_new.mass = zones[z0].mass + zones[z1].mass;
			ws_new.fmin = zones[z0].fmin;
			ws_new.barrier = DBL_MAX;//fnmin;

			// flood the two subzones with the new zone
			// TODO: make the loop bounds smarter
			for(luint i = 0; i < ntot; ++i) {
				if(zones[i] == z0 || zones[i] == z1)
					zones[i] = nzones;
			}

			// add the newest parent zone
			zones.push_back(ws_new);
			++nzones;
		}
	}*/

	return;
}

void neighbors() {

	// fill in zone information
	zones.resize(nzones);
	for(luint z = 0; z < nzones; ++z) {
		zones[z].vol = 0.0;
		zones[z].mass = 0.0;
		zones[z].min = DBL_MAX;
	}
	for(luint i = 0; i < ntot; ++i) {
		luint z = zone_cube[i];
		double f = field_cube[i];
		// TODO: use physical units for volume and mass
		zones[z].vol += 1.0;
		zones[z].mass += f;
		if(f < zones[z].min)
			zones[z].min = f;
	}

	// find barrier saddle points and merge zones
	for(luint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		luint ind_flat = inds_sorted[ind_uns];
		luint z = zone_cube[ind_flat];
		double f0 = field_cube[ind_flat];

		// get 3D indices from ind_flat
		luint ix = ind_flat/(ny*nz);
		luint iy = (ind_flat - ix*ny*nz)/nz;
		luint iz = ind_flat - ix*ny*nz - iy*nz;

		// inspect the 26 neighboring cells
		double barmin = DBL_MAX;
		luint znmin = LUINT_MAX;
		double grad_min = DBL_MAX; 
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					luint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					luint z_neighbor = zone_cube[tmp_flat];
					double f_neighbor = field_cube[tmp_flat];

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f_neighbor - f0)/sqrt(ox*ox + oy*oy + oz*oz);

					// check for a neighboring watershed with the shallowest gradient
					if(z != z_neighbor && grad < grad_min) { // Add a tolerance here?	
						grad_min = grad; // TODO: Is it possible to mess up here?
						barmin = f_neighbor;
						znmin = z_neighbor;
					}
				}
			}
		}
		// have we found a new neighbor zone?
		if(znmin != LUINT_MAX && !zones[z].neighbors.count(znmin)) {
			// add the mutual neighbors
			zones[z].neighbors[znmin] = barmin;
			zones[znmin].neighbors[z] = barmin;
		}
	}

	return;
}



void watershed() {

	nzones = 0;
	zone_cube.assign(ntot, LUINT_MAX); // unassigned zones use LUINT_MAX

	for(luint ind_uns = 0; ind_uns < ntot; ++ind_uns) {

		// get the flattened array index
		luint ind_flat = inds_sorted[ind_uns];

		// get 3D indices from ind_flat
		luint ix = ind_flat/(ny*nz);
		luint iy = (ind_flat - ix*ny*nz)/nz;
		luint iz = ind_flat - ix*ny*nz - iy*nz;

		// iterate over the 26 neighboring cells
		double f0 = field_cube[ind_flat];
		double grad_max = 0.0; 
		luint zmin = zone_cube[ind_flat];
		for(int ox = -1; ox <= 1; ++ox) {
			for(int oy = -1; oy <= 1; ++oy) {
				for(int oz = -1; oz <= 1; ++oz) {

					if(ox == 0 && oy == 0 && oz == 0) continue;

					// get neighboring flat indices, accounting for periodicity
					luint tmp_flat = ((ix + ox + nx)%nx)*ny*nz + ((iy + oy + ny)%ny)*nz + ((iz + oz + nz)%nz);

					// divide by the pixel distance to isotropize the neighbor stencil
					double grad = (f0 - field_cube[tmp_flat])/sqrt(ox*ox + oy*oy + oz*oz);

					if(grad > grad_max) { // finds the largest gradient to a neighboring cell 
					
						grad_max = grad;
						zmin = zone_cube[tmp_flat];
						zone_cube[ind_flat] = zmin;
					}
				}
			}
		}
		if(zmin == LUINT_MAX) {
			zone_cube[ind_flat] = nzones++;
		}
	}
	return;
}

// read in a multidimensional hdf5 file
void read_hdf5(string filename, string field_cubename) {
	vector<int> dims;
	HDFGetDatasetExtent(filename, field_cubename, dims);
	nx = dims[0]; ny = dims[1]; nz = dims[2];
	ntot = nx*ny*nz;
	HDFReadDataset(filename, field_cubename, field_cube);
	return;
}

void write_hdf5(string filename, string field_cubename) {
	HDFCreateFile(filename);
	luint dims[3] = {nx, ny, nz};
	HDFWriteDataset3D(filename, field_cubename, dims, zone_cube);
	return;
}

void argsort() {
	inds_sorted.resize(ntot);
	luint itmp = 0;
	generate(inds_sorted.begin(), inds_sorted.end(), [&] { return itmp++; });
	sort(inds_sorted.begin(), inds_sorted.end(), [&](luint a, luint b) { return field_cube[a] < field_cube[b]; });
	return;
}

