//
// heat2D_Jacobi.c
//
// Diffusion de chaleur par la methode de Jacobi
//
// Utilisation : heat2D_Jacobi nbLigs nbCols
// les arguments sur la ligne de commande sont optionnels, 
// leurs valeurs sont =80 par defaut 
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXVALUE 511  // la valeur maximale des donnees, utilisee dans initData()
#define NBSTEP   5    // le nombre d'iterations

/* Masque pour la mise a jour des donnees */
double MASK[] = {0.0, 0.1, 0.0, 0.1, 0.6, 0.1, 0.0, 0.1, 0.0};

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
double * initData(int nbLigs, int nbCols)
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
 * Savegarde des donnees dans un fichier
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
	    printf("\n");
	
		fclose(fp);
	}

	
}


/**
 * Ajout des donnees a un fichier
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

/**
 * Echange de 2 pointers sur double
 */
void exchangeData(double **data1, double **data2)
{
    double *tmp = *data1;
    
    *data1 = *data2;
    *data2 = tmp;
}


int main(int argc, char **argv)
{
    double *data1=NULL, *data2=NULL;
    int     nbLigs=80, nbCols=80;
    int     step=0;
    
    if (argc>=3) {
        nbLigs = atoi(argv[1]);
        nbCols = atoi(argv[2]);
    }
    
    printf("Taille du probleme : %d x %d\n", nbLigs, nbCols);
    printf("Pour la modifier : %s nbLigs nbCols\n\n", *argv);

    /* Initialisation des donnees */
    data1 = initData(nbLigs, nbCols);
    data2 = (double *) calloc(nbLigs*nbCols, sizeof(double));
    
    if (data1 && data2) {
        saveData(data1, nbLigs, nbCols, "saveData");
        
	/* Calcul de la diffusion de la chaleur */
        for (step=0; step<NBSTEP; step++) {
            updateNLigs(data1, data2, 1, nbLigs-2, nbCols, MASK);
            appendData(data2, nbLigs, nbCols, "saveData");
            exchangeData(&data1, & data2);
        }
     }
    
    if (data1) free(data1);
    if (data2) free(data2);
    
   return EXIT_SUCCESS;
}
