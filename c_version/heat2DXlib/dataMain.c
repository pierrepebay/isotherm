//
//  dataMain.c
//  
// Diffusion de chaleur par la methode de Jacobi
// Taille du probleme : nbLigs x nbCols; par defaut, 80 x 80
//                      peut etre modifiee sur la ligne de commande
//
//  Created by Jian-Jin LI on 01/02/2017.
//
//

#include "processingModule.h"

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/keysym.h>
#include <X11/Xatom.h>

#define NBCOLORS 256

static double MASK[] = {0.0, 0.1, 0.0, 0.1, 0.6, 0.1, 0.0, 0.1, 0.0};
static XColor greyScale[256];

/**
 * Affichage graphique des donnees de nbLigs x nbCols points
 * @param : 
 *      zoom : zoom en X et en Y de chaque point
 */
void drawData(Display *dpy, Window win, GC gContext,
              double *data, int nbLigs, int nbCols, int zoomX, int zoomY)
{
    int lig, col;
    int color;
    
    for (lig=0; lig<nbLigs; lig++) {
        for (col=0; col<nbCols; col++) {
            color = (int) (*(data + lig*nbCols + col) / MAXVALUE * 255);
            XSetForeground(dpy, gContext, greyScale[color].pixel);
            XFillRectangle(dpy, win, gContext, lig*zoomX, col*zoomY, zoomX, zoomY);
        }
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
    double      *data1=NULL, *data2=NULL;
    int          nbLigs=80, nbCols=80;
    int          step=0;
    
    Display      *dpy;           /* display                                          */
    Window        win;           /* fenetre graphique                                */
    XSizeHints    indic;
   
    int           ecran;
    unsigned long noir, blanc;   /* couleurs noir et blanche                         */
    
    unsigned long masque_valeur;
    XGCValues     valeurs_gc;
    GC            cg_cp;         /* contexte graphique pour le dessin definitif      */
    
    Colormap      cmap;          /* table de couleurs par defaut                     */
    
    XEvent        ev;            /* evenement                                        */
    
    char          caractere[10];
    KeySym        touche;
    
    Atom          protocoles[1]; /* Pour la fermeture de la fenetre avec le bouton x */
    
    int           encore, i;
    int           colorStep=1, zoom_x=4, zoom_y=4;
    
    if (argc>=3) { /* Taille du probleme */
        nbLigs = atoi(argv[1]);
        nbCols = atoi(argv[2]);
    }
    
    printf("Taille du probleme : %d x %d\n", nbLigs, nbCols);
    printf("Pour la modifier : %s nbLigs nbCols\n\n", *argv);
    
    data1 = initHeat2D(nbLigs, nbCols);
    data2 = (double *) calloc(nbLigs*nbCols, sizeof(double));
    
    if (data1 && data2) {
    
        /* --------------------------------------------------------- */
        /* Connexion au serveur X                                    */
        /* --------------------------------------------------------- */
        dpy = XOpenDisplay(NULL);
        if (!dpy) {
            fprintf(stderr, "Impossible d'ouvrir le display.\n");
            fprintf(stderr, "Assurez-vous que la variable DISPLAY est initialisee.\n");
            fprintf(stderr, "Elle doit contenir, par exemple : \"termx:0\"\n");
            exit(1);
        }
        
        /* --------------------------------------------------------- */
        /* Consultation de certaines caracteristiques du serveur X.  */
        /* --------------------------------------------------------- */
        ecran = DefaultScreen(dpy);
        blanc = WhitePixel(dpy, ecran);
        noir = BlackPixel(dpy, ecran);
        
        /* --------------------------------------------------------- */
        /* Allocation des couleurs                                   */
        /* --------------------------------------------------------- */
        cmap = DefaultColormap(dpy, ecran);
        colorStep = 256 / NBCOLORS;
        for (i=0; i<NBCOLORS; i+=colorStep) {
            greyScale[i].red = i << 8;
            greyScale[i].green = 0;
            greyScale[i].blue = 0;
            greyScale[i].flags = DoRed | DoGreen | DoBlue;
            XAllocColor( dpy, cmap, &(greyScale[i]));
        }
        
        /* --------------------------------------------------------- */
        /* Creation de la fenetre top-niveau du nom de "Jacobi"      */
        /* --------------------------------------------------------- */
        indic.x = 200; indic.y = 300;
        indic.width = nbLigs*zoom_x; indic.height = nbCols*zoom_y;
        indic.flags = PPosition | PSize;
        win = XCreateSimpleWindow(dpy, DefaultRootWindow(dpy),
                                     indic.x, indic.y, indic.width, indic.height,
                                     6, noir, blanc);
        /* Definition des proprietes standards de la fenetre top-niveau
         a l'intention du gestionnaire de fenetres */
        XStoreName(dpy, win, "Heat 2D");
        
        /* Pour la fermeture de la fenetre avec le bouton x */
        protocoles[0] = XInternAtom(dpy, "WM_DELETE_WINDOW", False);
        XSetWMProtocols(dpy, win, protocoles, 1);
        
        /* Selection des evenements interessants pour la fenetre top-niveau */
        XSelectInput(dpy, win, ButtonPressMask | KeyPressMask | ExposureMask | StructureNotifyMask);
        
        /* --------------------------------------------------------- */
        /* Creation d'un contexte graphique dessin */
        /* --------------------------------------------------------- */
        valeurs_gc.foreground = noir;
        valeurs_gc.background = blanc;
        valeurs_gc.function = GXcopy;
        masque_valeur = GCForeground | GCBackground | GCFunction;
        cg_cp = XCreateGC(dpy, win, masque_valeur, &valeurs_gc);
        
        XMapRaised(dpy, win);
        
        /* --------------------------------------------------------- */
        /* Boucle Principale */
        /* --------------------------------------------------------- */
        encore = True;
        while (encore) {
            XNextEvent(dpy, &ev);
            switch (ev.type) {
                case Expose :
                    while ( XCheckTypedEvent(dpy, Expose, &ev) ); // purge des ev. Expose
                
                    if (ev.xexpose.count == 0) {
                        printf("Je redessine : Window size=(%d, %d)\n", indic.width, indic.height);
                        drawData(dpy, win, cg_cp, data1, nbLigs, nbCols, zoom_x, zoom_y);
                        saveData(data1, nbLigs, nbCols, "saveData");
        
                    }
                    break;
                
                case ButtonPress: /* Evenement souris */
                    updateNLigs(data1, data2, 1, nbLigs-2, nbCols, MASK);
                    step++;
                    printf("Iteration %d\n", step);
                    
                    exchangeData(&data1, &data2);
                    drawData(dpy, win, cg_cp, data1, nbLigs, nbCols, zoom_x, zoom_y);
                    appendData(data1, nbLigs, nbCols, "saveData");
                    break;
                
                case ConfigureNotify :
                    if (ev.xconfigure.window == win) {
                        while ( XCheckTypedEvent(dpy, ConfigureNotify, &ev) ); // purge des ev. Configure
                    
                        /* recuperation de la nouvelle position et la nouvelle taille de la fenetre racine */
                        indic.x = ev.xconfigure.x;
                        indic.y = ev.xconfigure.y;
                    
                        if( (indic.width != ev.xconfigure.width) ||
                           (indic.height != ev.xconfigure.height)) {
                            indic.width = ev.xconfigure.width;
                            indic.height = ev.xconfigure.height;
                            printf("configreNotify: Window size=(%d, %d)\n", indic.width, indic.height);
                            
                            zoom_x = indic.width/nbCols;
                            zoom_y = indic.height/nbLigs;
                            drawData(dpy, win, cg_cp, data1, nbLigs, nbCols, zoom_x, zoom_y);
                        }
                    }
                    break;
                
                case ClientMessage : /* Pour la fermeture de la fenetre avec le bouton x */
                    printf("Fermeture de la fenetre\n");
                    encore = False;
                    break;
                
                case KeyPress :
                    i = XLookupString(&(ev.xkey), caractere, 10, &touche, 0);
                    printf("Touche Fin\n");
                        if (touche == XK_End) encore = False;
                    break;
                
            }
    }
    
    if (data1) free(data1);
    if (data2) free(data2);

    XDestroyWindow(dpy,win);
    XFreeGC(dpy, cg_cp);
    XCloseDisplay(dpy);
    }
    
   return EXIT_SUCCESS;
}

