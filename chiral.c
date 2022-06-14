/****************************************************************************/
/* This file deals with making a chiral nanocrystal */

#include "nc.h"

/****************************************************************************/
// Changes surface and passivation atom types in order to make the nanocrystal chiral
// Input:
// 		atom *atoms = a pointer to the first atom structure of a passivated nanocrystal
//
//
//
// Output:
// 		Changes the atom types that atoms point such that atoms has chiral centers now
// 		Returns the number of chiral centers added to the nanocrystal

int addChiralCenters(atom *atoms, int nAtoms, int nMaxBonds, chiralNC chiralParams) {
	int i, j, nChiralCenters = 0;
	int handedness = 0; // 1 = right handed and -1 = left handed
	vector crossProduct;

	for (i = 0; i < nAtoms; i++) {
		if (! newIsAtomASurfaceAtom(atoms, i, nAtoms)) {
			continue; // do not make interior atoms chiral
		}
		else {
			nChiralCenters++;
		}
	}

	return nChiralCenters;
}

/****************************************************************************/
