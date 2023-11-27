#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <stdbool.h>

typedef struct position
{
	float x, y, z;
} BEAD_POSITIONS;

typedef struct bonds
{
	int id, atom1, atom2, bondType;
} BONDS;

BEAD_POSITIONS **generatePositions (BEAD_POSITIONS **beads, int nBeads, float bondLength, BEAD_POSITIONS initPosition, int currentChain)
{
	beads[currentChain][0].x = initPosition.x;
	beads[currentChain][0].y = initPosition.y;
	beads[currentChain][0].z = initPosition.z;

	srand(time(0));
	float randomAngle1, randomAngle2;

	for (int i = 1; i < nBeads; ++i)
	{
		randomAngle1 = rand() % 360;
		randomAngle1 = randomAngle1 / 57.2958;
		randomAngle2 = rand() % 360;
		randomAngle2 = randomAngle2 / 57.2958;

		beads[currentChain][i].x = beads[currentChain][i - 1].x + bondLength * sin (randomAngle1);
		beads[currentChain][i].y = beads[currentChain][i - 1].y + bondLength * cos (randomAngle1);
		beads[currentChain][i].z = beads[currentChain][i - 1].z + bondLength * sin (randomAngle2);

	}

	return beads;
}

BEAD_POSITIONS generateInitialPositions (BEAD_POSITIONS initPosition, float boxLength)
{
	initPosition.x = fmod ((float)rand (), (float)boxLength);
	initPosition.y = fmod ((float)rand (), (float)boxLength);
	initPosition.z = fmod ((float)rand (), (float)boxLength);

	return initPosition;
}

bool checkChainBoundary (BEAD_POSITIONS **beads, int currentChain, int nBeads, float boxLength, bool chainFinalized)
{
	float x_max, y_max, z_max;
	float x_min, y_min, z_min;

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

	if (x_min < 0 || y_min < 0 || z_min < 0 || x_max > boxLength || y_max > boxLength || z_max > boxLength)
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

void printDataAtoms (FILE *output, BEAD_POSITIONS **beads, int nChains, int nBeads)
{
	int sino = 1;

	printf("Printing atomic coordinates..\n");
	for (int i = 0; i < nChains; ++i)
	{
		for (int j = 0; j < nBeads; ++j)
		{
			if ((i + 1) % 2 == 0)
			{
				fprintf(output, "%d %d 1 1.0 %f %f %f\n", sino, (i + 1), beads[i][j].x, beads[i][j].y, beads[i][j].z);
				sino++;
			}
			else
			{
				fprintf(output, "%d %d 1 -1.0 %f %f %f\n", sino, (i + 1), beads[i][j].x, beads[i][j].y, beads[i][j].z);
				sino++;
			}
		}
	}
}

BONDS *generateBonds (int nChains, int nBeads, int nBonds, BONDS *polymerBonds)
{
	int currentBondID = 0, bondType = 1;

	for (int i = 0; i < nBonds; ++i)
	{
		if ((i + 1) % nBeads != 0)
		{
			polymerBonds[currentBondID].id = currentBondID + 1;
			polymerBonds[currentBondID].bondType = 1;
			polymerBonds[currentBondID].atom1 = i + 1;
			polymerBonds[currentBondID].atom2 = i + 2;

			currentBondID++;
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

int main(int argc, char const *argv[])
{
	FILE *output;
	output = fopen ("chain.data", "w");

	int nBeads = atoi (argv[1]);
	float bondLength = atof (argv[2]);
	int nChains = atoi (argv[3]);
	float boxLength = atof (argv[4]);

	BEAD_POSITIONS **beads, initPosition;
	beads = (BEAD_POSITIONS **) malloc (nChains * sizeof (BEAD_POSITIONS *));

	for (int i = 0; i < nChains; ++i)
	{
		beads[i] = (BEAD_POSITIONS *) malloc (nBeads * sizeof (BEAD_POSITIONS));
	}

	BONDS *polymerBonds;
	int nBonds = (nBeads - 1) * nChains;
	polymerBonds = (BONDS *) malloc (nBonds * sizeof (BONDS));

	srand(time(0));

	bool chainFinalized;

	for (int i = 0; i < nChains; )
	{
		initPosition = generateInitialPositions (initPosition, boxLength);
		beads = generatePositions (beads, nBeads, bondLength, initPosition, i);
		chainFinalized = checkChainBoundary (beads, i, nBeads, boxLength, chainFinalized);

		if (chainFinalized) {
			++i; }
	}

	int nAtoms = nBeads * nChains;
	printDataHeader (output, nAtoms, nBonds, boxLength);
	printDataAtoms (output, beads, nChains, nBeads);
	polymerBonds = generateBonds (nChains, nBeads, nBonds, polymerBonds);
	printDataBonds (output, polymerBonds, nBonds);

	fclose (output);
	return 0;
}
