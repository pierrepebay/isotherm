//
//  imgProcessing.h
//  
//
//  Created by Jian-Jin LI on 01/02/2017.
//
//

#ifndef imgProcessing_h
#define imgProcessing_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#define MAXVALUE 511   /* la valeur maximale des donnees, utilisee dans initData() */

/**
 * Initialisation de donnees : un domaine de nbLigs x nbCols points
 */
double * initHeat2D(int nbLigs, int nbCols);

/**
 * Savegarde de donnees dans un fichier
 * @param :
 *      data : pointeur de nbLigs x nbCols double
 *      fname : nom du fichier de sauvegarde
 */
void saveData(double *data, int nbLigs, int nbCols, char *fname);

/**
 * Ajout de donnees dans un fichier
 * @param :
 *      data : pointeur de nbLigs x nbCols double
 *      fname : nom du fichier de sauvegarde
 */
void appendData(double *data, int nbLigs, int nbCols, char *fname);


/**
 * Mise a jour des lignes de start a end
 * @param :
 *      data : anciennes valeurs
 *      newData : nouvelles valeurs
 *      nbCols : nombre de donnees en Y
 */
void updateNLigs(double *old_data, double *new_data,
                 int start, int end, int nbCols,
                 double *masque);

#endif /* dataProcessing_h */
