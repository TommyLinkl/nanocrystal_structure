/****************************************************************************/
/* This file deals with removing the outer layer of a nanocrystals */

#include "nc.h"

/****************************************************************************/
// Removes the outer layer of semiconductor atoms from atoms and stores the 
// pointers to only the new nanocrystal atoms in newAtoms

int cutOneLayer(atom *newAtoms, atom *origAtoms, int nOrigAtoms, int nMaxBonds) {
	atom *tmpAtoms;
	int i, nNewAtoms = 0;

	// Dynamically allocate memory to be used only within this function
	tmpAtoms = (atom *) calloc(2*nOrigAtoms, sizeof(atom));
	for (i = 0; i < 2*nOrigAtoms; i++) {
		tmpAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
		tmpAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
	}

	// Create new section for cutting procedure
	writeSeparation(stdout);
	printf("Cutting one layer of the nanocrystal\n\n");

	for (i = 0; i < nOrigAtoms; i++) {
		if (isAtomASurfaceAtom(origAtoms, i, nMaxBonds) && newIsAtomASurfaceAtom(origAtoms, i, nOrigAtoms)) { 
			continue;
		}
		else {
			deepCopySingleAtom(tmpAtoms, nNewAtoms, origAtoms, i, nMaxBonds);
			nNewAtoms++;
		}
	}
	calcNearestNeighborLists(tmpAtoms, nMaxBonds, nNewAtoms);
	deepCopyAllAtoms(newAtoms, tmpAtoms, nNewAtoms, nMaxBonds);
	//nNewAtoms = cutDanglingAtoms(newAtoms, tmpAtoms, nMaxBonds, nMaxBonds-2, nNewAtoms);

	// Print new size of nanocrystal
	printf("\nNumber of original atoms = %d\n", nOrigAtoms);
	printf("Number of semiconductor atoms removed during cutting this layer = %d\n", nOrigAtoms-nNewAtoms);
	printf("Number of new atoms after cutting one layer = %d\n", nNewAtoms);

	// Free dynamically allocated memory
	for (i = 0; i < nOrigAtoms; i++) {
		free(tmpAtoms[i].neighborList);
		free(tmpAtoms[i].neighborPos);
	}
	free(tmpAtoms);

	return nNewAtoms;
}

/****************************************************************************/
// Fills the newAtoms pointers with only the semiconductor atoms from origAtoms that 
// have between (inclusive) nMinBonds and nMaxBonds to other semiconductor atoms

int cutDanglingAtoms(atom *newAtoms, atom *origAtoms, int nMaxBonds, int nMinBonds, int nOrigAtoms) {
	int i, iteration, nTmpAtoms, nNewAtoms = 0;
	atom *tmpAtoms;

	// Make sure the nearest neighbor list for origAtoms is up to date
	calcNearestNeighborLists(origAtoms, nMaxBonds, nOrigAtoms);

	// Dynamically allocate memory to be used only within this function
	tmpAtoms = (atom *) calloc(nOrigAtoms, sizeof(atom));
	for (i = 0; i < nOrigAtoms; i++) {
		tmpAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
		tmpAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
	}	
	deepCopyAllAtoms(tmpAtoms, origAtoms, nOrigAtoms, nMaxBonds);

	// Self-consistent loop
	printf("\nRemoving atoms with missing bonds:\n"); 
	nTmpAtoms = nOrigAtoms;
	iteration = 0;
	// TODO: remove doNotCutCdSurfaceFlag - convert to params.doNotCutCdSurfaceFlag
	int doNotCutCdSurfaceFlag = 0;
	while (nNewAtoms != nTmpAtoms) {
		iteration++;
		printf("\nIteration %d\n", iteration);
		if (iteration > 1) nTmpAtoms = nNewAtoms;
		printf("Starting number of atoms = %d\n", nTmpAtoms);
		nNewAtoms = 0;
		for (i = 0; i < nOrigAtoms; i++) {
			if (hasMinNumBonds(tmpAtoms, i, nMaxBonds, nMinBonds) || (! strcmp(tmpAtoms[i].symbol, "Cd") && doNotCutCdSurfaceFlag)) {
				deepCopySingleAtom(newAtoms, nNewAtoms, tmpAtoms, i, nMaxBonds);
				nNewAtoms++;
			} 
			else { // used to remove atom 
				tmpAtoms[i].pos.x = tmpAtoms[i].pos.y = 0.0;
				tmpAtoms[i].pos.z = 1000000.0; 
				tmpAtoms[i].pos.mag = 1000000.0; 
			}
		}
		calcNearestNeighborLists(tmpAtoms, nMaxBonds, nOrigAtoms);
		printf("Final number of atoms = %d\n", nNewAtoms);
	}

	// Free dynamically allocated memory
	for (i = 0; i < nOrigAtoms; i++) {
		free(tmpAtoms[i].neighborList);
		free(tmpAtoms[i].neighborPos);
	}
	free(tmpAtoms);

	return nNewAtoms;
}

/****************************************************************************/
// Return 1 (True) if the atom corresponding to atom at atomIndex has at least 
// nMinBonds to other semiconductor atoms; else return 0 (False)

int hasMinNumBonds(atom *atoms, int atomIndex, int nMaxBonds, int nMinBonds) {
	int i, nBonds = 0;
	
	// Count the number of bonds -> neighborPos.mag must be non-zero and the
	// symbol of the neighboring atom must not be a passivation symbol to count as a bond
	for (i = 0; i < nMaxBonds; i++) {
		if (atoms[atomIndex].neighborPos[i].mag > EPS && 
			! isAPassivationSymbol(atoms[atoms[atomIndex].neighborList[i]].symbol)) nBonds++;
	}

	if (nBonds >= nMinBonds) return 1; // true -> the atom has at least nMinBonds
	else return 0; // false -> the atom does not have nMinBonds
}

/****************************************************************************/
// Removes atoms from atoms array that have less than nMinBonds
// and return the new number of atoms

int removeAtomsWithMissingBonds(atom *atoms, int nMaxBonds, int nMinBonds, int nOrigAtoms) {
	int i, j, nBonds, nNewAtoms = nOrigAtoms;

	for (i = 0; i < nOrigAtoms; i++) {
		nBonds = 0;
		for (j = 0; j < nMaxBonds; j++) if (atoms[i].neighborPos[j].mag > EPS) nBonds++; 
		if (nBonds >= nMinBonds) continue;
		else {
			// TODO: remove this atom from the atoms list
			nNewAtoms--;
		}

	}

	printf("The number of new atoms = %d\n", nNewAtoms);

	return nNewAtoms;
}

/****************************************************************************/
// Return 1 (True) if the atom corresponding to atom at atomIndex has either:
// less than nMaxBonds number of bonds or is bonded to P1, P2, P3 or P4
// else returns 0 (False)

int isAtomASurfaceAtom(atom *atoms, int atomIndex, int nMaxBonds) {
	int i, isSurfaceAtom = 0;
	
	for (i = 0; i < nMaxBonds; i++) {
		if (atoms[atomIndex].neighborPos[i].mag > EPS && 
			! isAPassivationSymbol(atoms[atoms[atomIndex].neighborList[i]].symbol)) {
			isSurfaceAtom = 0;
		}
		else {
			isSurfaceAtom = 1;
			break;
		}
	}

	return isSurfaceAtom;
}

/****************************************************************************/
// Return 1 (True) if when one moves the atom corresponding to atomIndex 
// in all 14 directions by 3.0 A that in at least one of the directions 
// the new atom position would have zero bonds else returns 0 (False)
// Note: this should return 1 for atoms in an edge dislocation cite

int newIsAtomASurfaceAtom(atom *atoms, int atomIndex, int nAtoms) {
	int i, j, nTestPositions = 14;
	double testShift = 4.0;
	double testBondLength = 2.5; // standard for Cd, Cdz, Se and Sez -> smaller than testShift
	int nNeighbors[nTestPositions];
	vector diff, testPos[nTestPositions];
	
	// Return 1 (True) if atomIndex is a passivation atom
	if (isAPassivationSymbol(atoms[atomIndex].symbol)) return 1;

	// Test position vectors
	for (i = 0; i < nTestPositions; i++) testPos[i] = retZeroVector();
	testPos[0].x = testShift;
	testPos[1].x = -testShift;
	testPos[2].y = testShift;
	testPos[3].y = -testShift;
	testPos[4].z = testShift;
	testPos[5].z = -testShift;
	testPos[6].x = testShift; testPos[6].y = testShift; testPos[6].z = testShift;
	testPos[7].x = testShift; testPos[7].y = testShift; testPos[7].z = -testShift;
	testPos[8].x = testShift; testPos[8].y = -testShift; testPos[8].z = testShift;
	testPos[9].x = testShift; testPos[9].y = -testShift; testPos[9].z = -testShift;
	testPos[10].x = -testShift; testPos[10].y = testShift; testPos[10].z = testShift;
	testPos[11].x = -testShift; testPos[11].y = testShift; testPos[11].z = -testShift;
	testPos[12].x = -testShift; testPos[12].y = -testShift; testPos[12].z = testShift;
	testPos[13].x = -testShift; testPos[13].y = -testShift; testPos[13].z = -testShift;
	for (i = 6; i < nTestPositions; i++) testPos[i] = retScaledVector(testPos[i], 1.0/sqrt(testShift));
	for (i = 0; i < nTestPositions; i++) testPos[i] = retAddedVectors(atoms[atomIndex].pos, testPos[i]);

	// Calculate how many nearest neighbors each of the test positions would have
	for (i = 0; i < nTestPositions; i++) {
		nNeighbors[i] = 0;
		for (j = 0; j < nAtoms; j++) {
			if (isAPassivationSymbol(atoms[j].symbol)) continue;
			else { // test if would be a neighbor
				diff = retSubtractedVectors(testPos[i], atoms[j].pos);
				if (diff.mag < testBondLength) nNeighbors[i]++;
			}
		}
		if (! nNeighbors[i]) return 1; // if nNneighbors[i] = 0 then it is a surface atom
	}

	return 0; // if nNeighbors[i] > 0 for all i then atom is not a surface atom -> return 0 (False)
}

/****************************************************************************/
// Return 1 (True) if the atomSymbol signifies a passivation ligand (P1, P2, P3 or P4)
// and return 0 (False) otherwise

int isAPassivationSymbol(char *atomSymbol) {
	// check if a the symbol is P1, P2, P3 or P4 -> return 1/True if so
	if (! strcmp(atomSymbol, "P1") || ! strcmp(atomSymbol, "P2") ||
		! strcmp(atomSymbol, "P3") || ! strcmp(atomSymbol, "P4")) return 1;
	else return 0; // not a passivation ligand symbol -> return 0/False
}

/****************************************************************************/
