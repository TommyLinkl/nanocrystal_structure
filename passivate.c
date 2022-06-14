/****************************************************************************/
/* This file deals with passivating a nanocrystal */

#include "nc.h"

/****************************************************************************/
// Replaces the outer semiconductor layer with the passivation ligands

int passivateNanocrystal(atom *passivatedAtoms, atom *origAtoms, int nOrigAtoms, int nMaxBonds) {
    FILE *pf;
    int i, j, k, tmpIndex, nAtomTypes;
    int nSemiCondAtoms, nSurfaceAtoms, nRemovedAtoms, nTotalAtoms; 
    vector tmpVector, passPosition;
    atom *tmpAtoms;

    // Dynamically allocate memory to be used only within this function
    tmpAtoms = (atom *) calloc(2*nOrigAtoms, sizeof(atom));
    for (i = 0; i < 2*nOrigAtoms; i++) {
        tmpAtoms[i].neighborList = (int *) calloc(nMaxBonds, sizeof(int));
        tmpAtoms[i].neighborPos = (vector *) calloc(nMaxBonds, sizeof(vector));
    }

    // Print out conf without passivation ligands
    calcNearestNeighborLists(origAtoms, nMaxBonds, nOrigAtoms); // make sure NN list is up to date
    pf = fopen("prepassivated_conf.xyz", "w");
    writeXYZFile(origAtoms, nOrigAtoms, 1.0, pf);
    fclose(pf);	

    // Write output
    writeSeparation(stdout);
    printf("Passivating the nanocrystal\n\n");

    // Remove the surface atoms
    nSemiCondAtoms = nSurfaceAtoms = 0;
    for (i = 0; i < nOrigAtoms; i++) {
        if (isAtomASurfaceAtom(origAtoms, i, nMaxBonds) && newIsAtomASurfaceAtom(origAtoms, i, nOrigAtoms)) { 
            nSurfaceAtoms++;
        }  
        else {
            deepCopySingleAtom(tmpAtoms, nSemiCondAtoms, origAtoms, i, nMaxBonds);
            nSemiCondAtoms++;
        }
    }
    printf("Number of original atoms = %d\n", nOrigAtoms);
    printf("Number of surface atoms removed = %d\n", nSurfaceAtoms);
    printf("Number of atoms excluding passivation (possible dangling atoms) = %d\n\n", nSemiCondAtoms);

    // Remove atoms with missing bonds and store final in passivatedAtoms
    printf("Finalizing the semiconductor atoms in the passivated nanocrystal:\n\n");
    calcNearestNeighborLists(tmpAtoms, nMaxBonds, nSemiCondAtoms);
    nSemiCondAtoms = cutDanglingAtoms(passivatedAtoms, tmpAtoms, nMaxBonds, nMaxBonds-2, nSemiCondAtoms);
    calcNearestNeighborLists(passivatedAtoms, nMaxBonds, nSemiCondAtoms); // probably not needed

    // Print out conf without passivation ligands
    pf = fopen("unpassivated_conf.xyz", "w");
    writeXYZFile(passivatedAtoms, nSemiCondAtoms, 1.0, pf);
    fclose(pf);

    // Build list of atoms in tmpAtoms that have been removed from origAtoms in making passivatedAtoms so far
    nRemovedAtoms = 0;
    vector diff;
    for (i = 0; i < nOrigAtoms; i++) {
        for (j = 0; j < nSemiCondAtoms; j++) {
            diff = retSubtractedVectors(origAtoms[i].pos, passivatedAtoms[j].pos);
            if (diff.mag < EPS) break;
        }
        if (j == nSemiCondAtoms) { // the origAtoms[i] atom is not in passivatedAtoms
            deepCopySingleAtom(tmpAtoms, nRemovedAtoms, origAtoms, i, nMaxBonds);
            nRemovedAtoms++;
        }
    }
    printf("\nTotal number of semiconductor atoms removed during passivation = %d\n", nRemovedAtoms);

    // Passivate the outer layer of passivatedAtoms
    printf("\nAdding in the passivation atoms:\n\n");
    nTotalAtoms = nSemiCondAtoms;
    // TODO: remove allCdSurfaceFlag - convert to params.allCdSurfaceFlag
    int allCdSurfaceFlag = 0;
    int nBonds;
    for (i = 0; i < nSemiCondAtoms; i++) {
        if (! isAtomASurfaceAtom(passivatedAtoms, i, nMaxBonds) || 
            ! newIsAtomASurfaceAtom(passivatedAtoms, i, nSemiCondAtoms)) {
                continue; // do not passivate interior atoms
        }
        else {
            tmpIndex = 0;
            nBonds = 0;
	    // Get number of missing bonds - needed for Se-P2 distance
	    for (j = 0; j < nMaxBonds; j++) if (passivatedAtoms[i].neighborPos[j].mag > EPS) nBonds++;
	    for (j = 0; j < nMaxBonds; j++) {
		if (passivatedAtoms[i].neighborPos[j].mag > EPS) continue; // this neighbor is a semiconductor atom
		else { // add passivation atom to passivatedAtoms
		    for (k = tmpIndex; k < nRemovedAtoms; k++) { // find its neighbor in tmpAtoms
			if (areAtomsNearestNeighbors(passivatedAtoms[i], tmpAtoms[k])) {
			    tmpIndex = k+1;
			    tmpVector = retSubtractedVectors(passivatedAtoms[i].pos, tmpAtoms[k].pos);
			    // TODO: have separate function for for determining if cation
			    // if (atomIsCation(passivatedAtoms[i].symbol)) { // True of symbol = Cd, Zn, Cdz, Cd1, ...
			    if (! strcmp(passivatedAtoms[i].symbol, "Cd") || ! strcmp(passivatedAtoms[i].symbol, "Cdz") ||
				! strcmp(passivatedAtoms[i].symbol, "Zn") || ! strcmp(passivatedAtoms[i].symbol, "In")  ||
				! strcmp(passivatedAtoms[i].symbol, "Ga") ) { 
				    passPosition = retScaledVector(tmpVector, 0.45); // 45% shorter 
				    strcpy(passivatedAtoms[nTotalAtoms].symbol, "P1");
			    } 
			    // else if (atomIsAnion(passivatedAtoms[i].symbol)) { // True of symbol = Se, S, Te, Sez, Se1, ...
			    else if (! strcmp(passivatedAtoms[i].symbol, "Se") || ! strcmp(passivatedAtoms[i].symbol, "S") ||
				! strcmp(passivatedAtoms[i].symbol, "Te") || ! strcmp(passivatedAtoms[i].symbol, "Sez") ||
				! strcmp(passivatedAtoms[i].symbol, "As") || ! strcmp(passivatedAtoms[i].symbol, "P") ) {
				    if (allCdSurfaceFlag) {
					passPosition = retZeroVector();
				    }
				    else if (nBonds == 3) {
					passPosition = retScaledVector(tmpVector, 0.7); // 75% shorter 
				    }
				    else if (nBonds == 2) { 
					passPosition = retScaledVector(tmpVector, 0.75); // 75% shorter 
				    }
				    else { // should not reach here
					passPosition = retScaledVector(tmpVector, 0.80); // 80% shorter 
				    }
				    if (allCdSurfaceFlag) {
					strcpy(passivatedAtoms[nTotalAtoms].symbol, "Cd");
				    }
				    else {
					strcpy(passivatedAtoms[nTotalAtoms].symbol, "P2");
				    }
			    }
                            else {
                                printf("Unknown atom, %s, that needs to be passivated!!\n", passivatedAtoms[i].symbol);
                                printf("The program is exiting due to this error (can be fixed in passivate.c)\n");
                                fflush(stdout);
                                exit(EXIT_FAILURE);
                            }

			    passivatedAtoms[nTotalAtoms].pos = retAddedVectors(tmpAtoms[k].pos, passPosition); 
			    // TODO: save passivation ligand position / index to passivatedAtoms neighbor list
			    passivatedAtoms[i].neighborList[j] = nTotalAtoms;
			    nTotalAtoms++;
			    break;
			}
		    }
		}
	    }
	}
    }
    nAtomTypes = assignAtomTypes(passivatedAtoms, nTotalAtoms); // updates the atom type labels
    printf("Number of atoms excluding passivation = %d\n", nSemiCondAtoms);
    printf("Number of passivation atoms = %d\n", nTotalAtoms-nSemiCondAtoms);
    printf("Number of atoms including passivation = %d\n", nTotalAtoms); fflush(stdout);
    
    // Print the final, passivated nanocrystal
    pf = fopen("passivated_conf.dat", "w");
    writeConfFile(passivatedAtoms, nTotalAtoms, ANGTOAU, pf);
    fclose(pf);
    pf = fopen("passivated_conf_H.xyz", "w");
    writeXYZFile(passivatedAtoms, nTotalAtoms, 1.0, pf);
    fclose(pf);
	
    // Temporary - making chiral nanocrystal
    int chiralFlag = 0;
    if (chiralFlag) {
        long alreadyChangedFlag, chiralType;
        vector diffCdPassivation, crossProduct;
        double dotProduct;
        for (i = 0; i < nSemiCondAtoms; i++) {
            if (! isAtomASurfaceAtom(passivatedAtoms, i, nMaxBonds) || 
                ! newIsAtomASurfaceAtom(passivatedAtoms, i, nSemiCondAtoms)) {
                continue; // do not passivate interior atoms
            }
            else {
                tmpIndex = 0;
                nBonds = 0;
            	// Get number of missing bonds - needed for Se-P2 distance
            	for (j = 0; j < nMaxBonds; j++) if (passivatedAtoms[i].neighborPos[j].mag > EPS) nBonds++;
            	if (nBonds == 2 && ! strcmp(passivatedAtoms[i].symbol, "Cd")) {
            	    alreadyChangedFlag = 0;
            	    for (j = 0; j  < nMaxBonds; j++) {
            		if (! strcmp(passivatedAtoms[passivatedAtoms[i].neighborList[j]].symbol, "Se1")) {
            		    alreadyChangedFlag += 1;
            		}
            	    }
            	    if (alreadyChangedFlag == 2) {
            		continue;
            	    }
            	    else if (alreadyChangedFlag == 1) {
            		if (! strcmp(passivatedAtoms[passivatedAtoms[i].neighborList[0]].symbol, "Se")) {
            		    crossProduct = retCrossProduct(passivatedAtoms[i].neighborPos[0], passivatedAtoms[i].neighborPos[1]);
            		}
            		else {
            		    crossProduct = retCrossProduct(passivatedAtoms[i].neighborPos[1], passivatedAtoms[i].neighborPos[0]);
            		}
            		diffCdPassivation = retSubtractedVectors(passivatedAtoms[i].pos, passivatedAtoms[passivatedAtoms[i].neighborList[2]].pos);
            		dotProduct = retDotProduct(crossProduct, diffCdPassivation);
            		if (dotProduct > 0.0) {
            		    strcpy(passivatedAtoms[passivatedAtoms[i].neighborList[3]].symbol, "P3");
            		}
            		else {
            		    strcpy(passivatedAtoms[passivatedAtoms[i].neighborList[2]].symbol, "P3");
            		}
            	    }
            	    else { // make Cd a chiral center -> have freedom to change both Se atoms and ligands as desired
            		crossProduct = retCrossProduct(passivatedAtoms[i].neighborPos[0], passivatedAtoms[i].neighborPos[1]);
            		diffCdPassivation = retSubtractedVectors(passivatedAtoms[i].pos, passivatedAtoms[passivatedAtoms[i].neighborList[2]].pos);
            		dotProduct = retDotProduct(crossProduct, diffCdPassivation);
            		strcpy(passivatedAtoms[passivatedAtoms[i].neighborList[1]].symbol, "Se1");
            		if (dotProduct > 0.0) {
            		    strcpy(passivatedAtoms[passivatedAtoms[i].neighborList[3]].symbol, "P3");
            		}
            		else {
            		    strcpy(passivatedAtoms[passivatedAtoms[i].neighborList[2]].symbol, "P3");
            		}
            	    }
            	}
            }
	}

	// Print the chiral passivated nanocrystal
	pf = fopen("chiral_passivated_conf.dat", "w");
	writeConfFile(passivatedAtoms, nTotalAtoms, ANGTOAU, pf);
	fclose(pf);
	pf = fopen("chiral_passivated_conf_H.xyz", "w");
	writeXYZFile(passivatedAtoms, nTotalAtoms, 1.0, pf);
	fclose(pf);
    }

    // Free dynamically allocated memory
    for (i = 0; i < nOrigAtoms; i++) {
    	free(tmpAtoms[i].neighborList);
    	free(tmpAtoms[i].neighborPos);
    }
    free(tmpAtoms);
    
    return nTotalAtoms;
} 

/****************************************************************************/
// Returns 1 (True) if a passivation ligand is one of the symbols in atoms
// and returns 0 (False) if all atoms are semiconductor atoms

int isCrystalPassivated(atom *atoms, int nAtoms) {
	int i;

	for (i = 0; i < nAtoms; i++) // return 1 (True) if a passivation ligand is found
		if (isAPassivationSymbol(atoms[i].symbol)) return 1;  

	return 0;
}

/****************************************************************************/
