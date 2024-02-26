#include "processingModule.h"


/**
 * Mise a jour des lignes de start a end
 * @param :
 *      data : anciennes valeurs
 *      newData : nouvelles valeurs
 *      nbCols : nombre de donnees en Y
 */
void updateNLigs(double *old_data, double *new_data,
                 int start, int end, int nbCols,
                 double *masque)
{
    int ix, iy, ixm, iym;
    double *ptr_masque;
    
    for (ix=start; ix<=end; ix++) {
        for (iy=1; iy<nbCols-1; iy++) {
            ptr_masque = masque;
            *(new_data+ix*nbCols+iy) = 0.0;
            
            for (ixm=-1; ixm<=1; ixm++)
                for (iym=-1; iym<=1; iym++) {
                    *(new_data+ix*nbCols+iy) += *ptr_masque * *(old_data + (ix+ixm)*nbCols + iy+iym);
                    ptr_masque++;
                }
            
        }
    }
}


/**
 * Initialisation de donnees : un domaine de nbLigs x nbCols points
 */
double * initHeat2D(int nbLigs, int nbCols)
{
	double *data = NULL;
	int    i, j;

	data = (double *) malloc(nbLigs * nbCols * sizeof(double));

	if (data) {
		for (i=0; i<nbLigs; i++)
			for (j=0; j<nbCols; j++)
				*(data + i*nbCols + j) = ((double) (i * (nbLigs-i-1) * j * (nbCols-j-1)) / (pow(nbLigs/2.0,2.0) * pow(nbCols/2.0, 2.0))) * MAXVALUE;
	}

	return data;
}


/**
 * Savegarde de donnees dans un fichier
 * @param :
 *      data : pointeur de nbLigs x nbCols double
 *      fname : nom du fichier de sauvegarde
 */
void saveData(double *data, int nbLigs, int nbCols, char *fname)
{
	FILE *fp = fopen(fname, "w");
    int   i, j;

	if (fp) {
		fprintf(fp, "%d %d\n", nbLigs, nbCols);

		for (i=0; i<nbLigs; i++)
			for (j=0; j<nbCols; j++)
				fprintf(fp, "%6.1f%c", *(data+i*nbCols+j),(j==nbCols-1)?'\n':' ');
	
		fclose(fp);
	}

	
}


/**
 * Ajout de donnees dans un fichier
 * @param :
 *      data : pointeur de nbLigs x nbCols double
 *      fname : nom du fichier de sauvegarde
 */
void appendData(double *data, int nbLigs, int nbCols, char *fname)
{
    FILE *fp = fopen(fname, "a");
    int   i, j;
    
    if (fp) {
        fprintf(fp, "\n");
        
        for (i=0; i<nbLigs; i++)
            for (j=0; j<nbCols; j++)
                fprintf(fp, "%6.1f%c", *(data+i*nbCols+j),(j==nbCols-1)?'\n':' ');
        
        fclose(fp);
    }
    
    
}

