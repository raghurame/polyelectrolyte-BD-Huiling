#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

typedef struct position
{
	int molType;
	double x, y, z;
} BEAD_POSITIONS;

typedef struct bonds
{
	int id, atom1, atom2, bondType;
} BONDS;

double computeDistance (double x1, double y1, double z1, double x2, double y2, double z2)
{
	double distance;

	distance = sqrt (
		(x2 - x1) * (x2 - x1) +
		(y2 - y1) * (y2 - y1) +
		(z2 - z1) * (z2 - z1)
		);

	return distance;
}

BEAD_POSITIONS **generatePositions (BEAD_POSITIONS **beads, int nBeads, float bondLength, BEAD_POSITIONS initPosition, int currentChain, int currentMolType)
{
	double checkingDistance;

	beads[currentChain][0].x = initPosition.x;
	beads[currentChain][0].y = initPosition.y;
	beads[currentChain][0].z = initPosition.z;
	beads[currentChain][0].molType = currentMolType;

	// srand(time(0));
	float randomAngle1, randomAngle2;

	for (int i = 1; i < nBeads; ++i)
	{
		randomAngle1 = rand() % 360;
		randomAngle2 = rand() % 360;
		randomAngle1 = randomAngle1 / 57.2958;
		randomAngle2 = randomAngle2 / 57.2958;

		beads[currentChain][i].x = beads[currentChain][i - 1].x + cos (randomAngle2) * (double)bondLength * cos (randomAngle1);
		beads[currentChain][i].y = beads[currentChain][i - 1].y + cos (randomAngle2) * (double)bondLength * sin (randomAngle1);
		beads[currentChain][i].z = beads[currentChain][i - 1].z + (double)bondLength * sin (randomAngle2);

		checkingDistance = computeDistance (beads[currentChain][i - 1].x, beads[currentChain][i - 1].y, beads[currentChain][i - 1].z, beads[currentChain][i].x, beads[currentChain][i].y, beads[currentChain][i].z);

		printf("%lf\n", checkingDistance);
		usleep (100000);

		beads[currentChain][i].molType = currentMolType;
	}

	return beads;
}

BEAD_POSITIONS generateRandomInitialPositions (BEAD_POSITIONS initPosition, float boxLength)
{
	initPosition.x = fmod ((double)rand (), (double)boxLength);
	initPosition.y = fmod ((double)rand (), (double)boxLength);
	initPosition.z = fmod ((double)rand (), (double)boxLength);

	return initPosition;
}

BEAD_POSITIONS generatePeriodicInitialPositions (BEAD_POSITIONS initPosition, float boxLength, bool *periodicChainPlacement)
{
	(*periodicChainPlacement) = true;

	initPosition.x = (double)boxLength / 2;
	initPosition.y = (double)boxLength / 2;
	initPosition.z = (double)boxLength / 2;

	return initPosition;
}

bool checkChainBoundary (BEAD_POSITIONS **beads, int currentChain, int nBeads, float boxLength, bool chainFinalized)
{
	double x_max, y_max, z_max;
	double x_min, y_min, z_min;

	for (int i = 0; i < nBeads; ++i)
	{
		if (i == 0)
		{
			x_max = beads[currentChain][i].x;
			x_min = beads[currentChain][i].x;
			y_max = beads[currentChain][i].y;
			y_min = beads[currentChain][i].y;
			z_max = beads[currentChain][i].z;
			z_min = beads[currentChain][i].z;
		}
		else
		{
			if (beads[currentChain][i].x < x_min) {
				x_min = beads[currentChain][i].x; }
			else if (beads[currentChain][i].x > x_max) {
				x_max = beads[currentChain][i].x; }

			if (beads[currentChain][i].y < y_min) {
				y_min = beads[currentChain][i].y; }
			else if (beads[currentChain][i].y > y_max) {
				y_max = beads[currentChain][i].y; }

			if (beads[currentChain][i].z < z_min) {
				z_min = beads[currentChain][i].z; }
			else if (beads[currentChain][i].z > z_max) {
				z_max = beads[currentChain][i].z; }
		}
	}

	if (x_min < 0 || y_min < 0 || z_min < 0 || x_max > (double)boxLength || y_max > (double)boxLength || z_max > (double)boxLength)
	{
		chainFinalized = false;
	}
	else
	{
		chainFinalized = true;
	}

	return chainFinalized;
}

void printDataHeader (FILE *output, int nAtoms, int nBonds, float boxLength)
{
	printf("Printing header information in the data file...\n");
	FILE *getCurrentDate;
	getCurrentDate = popen ("date", "r");
	char lineString[3000], firstLine[3000];

	fgets (lineString, 3000, getCurrentDate);
	snprintf (firstLine, 3000, "BD polyelectrolyte simulations - created on %s", lineString);

	fprintf(output, "%s\n\n%d atoms\n%d atom types\n%d bonds\n%d bond types\n\n0 %.3f xlo xhi\n0 %.3f ylo yhi\n0 %.3f zlo zhi\n\nMasses\n\n1 1\n\nAtoms\n\n", firstLine, nAtoms, 1, nBonds, 1, boxLength, boxLength, boxLength);
}

void printDataAtoms (FILE *output, BEAD_POSITIONS **beads1, int nBeads1, BEAD_POSITIONS **beads2, int nBeads2, int nChains)
{
	FILE *getCurrentDate;
	getCurrentDate = popen ("date", "r");

	char lineString[3000], firstLine[3000];
	fgets (lineString, 3000, getCurrentDate);

	FILE *outputXYZ;
	outputXYZ = fopen ("chain.xyz", "w");
	int nBeads = nBeads1 + nBeads2;

	fprintf(outputXYZ, "%d\n", nBeads * nChains);
	fprintf(outputXYZ, "XYZ file created on %s", lineString);

	int sino = 1, currentAtomID = 0, currentPolycation = 0, currentPolyanion = 0;

	printf("Printing atomic coordinates..\n");
	for (int i = 0; i < nChains; ++i)
	{
		if (i % 2 == 0)
		{
			for (int j = 0; j < nBeads1; ++j)
			{
				currentAtomID++;
				fprintf(output, "%d %d 1 %lf %lf %lf\n", currentAtomID, beads1[currentPolycation][j].molType, beads1[currentPolycation][j].x, beads1[currentPolycation][j].y, beads1[currentPolycation][j].z);
				fprintf(outputXYZ, "C %lf %lf %lf\n", beads1[currentPolycation][j].x, beads1[currentPolycation][j].y, beads1[currentPolycation][j].z);
			}
			currentPolycation++;
		}
		else
		{
			for (int j = 0; j < nBeads2; ++j)
			{
				currentAtomID++;
				fprintf(output, "%d %d 2 %lf %lf %lf\n", currentAtomID, beads2[currentPolyanion][j].molType, beads2[currentPolyanion][j].x, beads2[currentPolyanion][j].y, beads2[currentPolyanion][j].z);
				fprintf(outputXYZ, "N %lf %lf %lf\n", beads2[currentPolyanion][j].molType, beads2[currentPolyanion][j].x, beads2[currentPolyanion][j].y, beads2[currentPolyanion][j].z);
			}
			currentPolyanion++;
		}
	}
}

BONDS *generateBonds (int nChains, BEAD_POSITIONS **beads1, int nBeads1, BEAD_POSITIONS **beads2, int nBeads2, int nBonds, BONDS *polymerBonds, int nAtoms)
{
	int currentBondID = 0, bondType = 1, currentPolyanion = 0, currentPolycation = 0, currentAtomID = 0;

	for (int i = 0; i < nChains; ++i)
	{
		if (i % 2 == 0)
		{
			for (int j = 0; j < nBeads1; ++j)
			{
				currentAtomID++;
				if (j > 0)
				{
					currentBondID++;
					polymerBonds[currentBondID - 1].id = currentBondID;
					polymerBonds[currentBondID - 1].bondType = 1;
					polymerBonds[currentBondID - 1].atom1 = currentAtomID - 1;
					polymerBonds[currentBondID - 1].atom2 = currentAtomID;
				}
			}

			currentPolycation++;
		}
		else
		{
			for (int j = 0; j < nBeads2; ++j)
			{
				currentAtomID++;
				if (j > 0)
				{
					currentBondID++;
					polymerBonds[currentBondID - 1].id = currentBondID;
					polymerBonds[currentBondID - 1].bondType = 2;
					polymerBonds[currentBondID - 1].atom1 = currentAtomID - 1;
					polymerBonds[currentBondID - 1].atom2 = currentAtomID;
				}
			}

			currentPolyanion++;
		}
	}

	return polymerBonds;
}

void printDataBonds (FILE *output, BONDS *polymerBonds, int nBonds)
{
	fprintf(output, "\nBonds\n\n");
	printf("Printing bonds...\n");

	for (int i = 0; i < nBonds; ++i)
	{
		fprintf(output, "%d %d %d %d\n", polymerBonds[i].id, polymerBonds[i].bondType, polymerBonds[i].atom1, polymerBonds[i].atom2);
	}
}

float computeEndToEndDistance (float endToEndDistance, BEAD_POSITIONS **beads, int currentChain, int nBeads)
{
	endToEndDistance = 0;
	float currentDistance;

	for (int i = 0; i < nBeads; ++i)
	{
		for (int j = 0; j < nBeads; ++j)
		{
			if (i != j)
			{
				currentDistance = computeDistance (beads[currentChain][i].x, beads[currentChain][i].y, beads[currentChain][i].z, beads[currentChain][j].x, beads[currentChain][j].y, beads[currentChain][j].z);

/*				currentDistance = sqrt (
					(beads[currentChain][i].x - beads[currentChain][j].x) * (beads[currentChain][i].x - beads[currentChain][j].x) +
					(beads[currentChain][i].y - beads[currentChain][j].y) * (beads[currentChain][i].y - beads[currentChain][j].y) +
					(beads[currentChain][i].z - beads[currentChain][j].z) * (beads[currentChain][i].z - beads[currentChain][j].z)
					);
*/

				if (currentDistance > endToEndDistance)
				{
					endToEndDistance = currentDistance;
				}
			}
		}
	}

	return endToEndDistance;
}

BEAD_POSITIONS computeCenterOfMass (BEAD_POSITIONS com, BEAD_POSITIONS **beads, int nBeads, int currentChain)
{
	com.x = 0; com.y = 0; com.z = 0;

	for (int i = 0; i < nBeads; ++i)
	{
		com.x += beads[currentChain][i].x;
		com.y += beads[currentChain][i].y;
		com.z += beads[currentChain][i].z;
	}

	com.x /= nBeads;
	com.y /= nBeads;
	com.z /= nBeads;

	return com;
}

void packPolymers (BEAD_POSITIONS ***beads1, int nBeads1, BEAD_POSITIONS ***beads2, int nBeads2, int nChains, float maxEndToEndDistance)
{
	int currentPolycation = 0, currentPolyanion = 0;
	BEAD_POSITIONS com, *lattice1, *lattice2;

	for (int i = 0; i < nChains; ++i)
	{
		com = computeCenterOfMass (com, (*beads1), nBeads1, i);
		com = computeCenterOfMass (com, (*beads2), nBeads2, i);
	}
}

int main(int argc, char const *argv[])
{
	if (argc != 7)
	{
		printf("REQUIRED ARGUMENTS:\n~~~~~~~~~~~~~~~~~~~\n\n {~} argv[0] = program\n {~} argv[1] = Number of beads of type 1\n {~} argv[2] = Number of beads of type 2\n {~} argv[3] = bond length of type 1\n {~} argv[4] = bond length of type 2\n {~} argv[5] = Number of chains to add\n {~} argv[6] = Simulation box length\n\n");
		exit (1);
	}

	FILE *output;
	output = fopen ("chain.data", "w");

	int nBeads1 = atoi (argv[1]), nBeads2 = atoi (argv[2]), nBeads = nBeads1 + nBeads2;
	float bondLength1 = atof (argv[3]), bondLength2 = atof (argv[4]);
	int nChains = atoi (argv[5]), nPolycation = ceil (nChains / 2), nPolyanion = ceil (nChains / 2);
	float boxLength = atof (argv[6]);

	BEAD_POSITIONS **beads1, **beads2, initPosition;
	beads1 = (BEAD_POSITIONS **) malloc (nPolycation * sizeof (BEAD_POSITIONS *));
	beads2 = (BEAD_POSITIONS **) malloc (nPolyanion * sizeof (BEAD_POSITIONS *));

	for (int i = 0; i < nChains; ++i)
	{
		beads1[i] = (BEAD_POSITIONS *) malloc (nBeads1 * sizeof (BEAD_POSITIONS));
		beads2[i] = (BEAD_POSITIONS *) malloc (nBeads2 * sizeof (BEAD_POSITIONS));
	}

	BONDS *polymerBonds;
	int nBonds1 = (nBeads1 - 1) * nPolycation, nBonds2 = (nBeads2 - 1) * nPolyanion, nBonds = nBonds1 + nBonds2;
	polymerBonds = (BONDS *) malloc (nBonds * sizeof (BONDS));

	srand(time(0));

	bool chainFinalized, periodicChainPlacement = false;
	float maxEndToEndDistance = 0, currentEndToEndDistance;

	int currentPolycation = 0, currentPolyanion = 0;

	for (int i = 0; i < nChains; )
	{
		// initPosition = generateRandomInitialPositions (initPosition, boxLength);
		initPosition = generatePeriodicInitialPositions (initPosition, boxLength, &periodicChainPlacement);

		if (i % 2 == 0) {
			beads1 = generatePositions (beads1, nBeads1, bondLength1, initPosition, currentPolycation, 1);
			chainFinalized = checkChainBoundary (beads1, currentPolycation, nBeads1, boxLength, chainFinalized);
			currentEndToEndDistance = computeEndToEndDistance (currentEndToEndDistance, beads1, currentPolycation, nBeads1);
			if (currentEndToEndDistance > maxEndToEndDistance) {
				maxEndToEndDistance = currentEndToEndDistance; }
			currentPolycation++; }
		else {
			beads2 = generatePositions (beads2, nBeads2, bondLength2, initPosition, currentPolyanion, 2);
			chainFinalized = checkChainBoundary (beads2, currentPolyanion, nBeads2, boxLength, chainFinalized);
			currentEndToEndDistance = computeEndToEndDistance (currentEndToEndDistance, beads2, currentPolyanion, nBeads2);
			if (currentEndToEndDistance > maxEndToEndDistance) {
				maxEndToEndDistance = currentEndToEndDistance; }
			currentPolyanion++; }

		if (chainFinalized) {
			++i; }
	}

	packPolymers (&beads1, nBeads1, &beads2, nBeads2, nChains, maxEndToEndDistance);

	currentPolycation = 0; currentPolyanion = 0;
	printf("max end to end distance: %f\n", maxEndToEndDistance);

	int nAtoms = (nBeads1 * nPolycation) + (nBeads2 * nPolyanion);
	printDataHeader (output, nAtoms, nBonds, boxLength);
	printDataAtoms (output, beads1, nBeads1, beads2, nBeads2, nChains);
	polymerBonds = generateBonds (nChains, beads1, nBeads1, beads2, nBeads2, nBonds, polymerBonds, nAtoms);
	printDataBonds (output, polymerBonds, nBonds);

	fclose (output);
	return 0;
}
