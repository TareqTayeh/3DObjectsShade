/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #3
 * Date: 9 November 2017
 * Program purpose: The purpose of this program is to design and implementat a program that allows a user to shade
 *                  objects (spheres, cones, and torii) in various colours. The program should scale the parametric
 *                  objects within a scene (using 3D homogeneous transformation matrices) and render them using
 *                  a simple Lambertian shading model that assumes diffuse light reflection from objects.
 * File: dialog.cpp
 */


#include "dialog.h"
#include "ui_dialog.h"
#include "matrix.h"
#include <vector>

#define Ex 15.0
#define Ey 15.0
#define Ez 15.0

#define Gx 0.0
#define Gy 0.0
#define Gz 0.0

#define UPx 0.0
#define UPy 0.0
#define UPz 1.0

#define NP 5.0
#define FP 50.0

#define THETA 90.0
#define ASPECT 1.0

#define W  512
#define H  512

//Light Source
#define Lx 256.0
#define Ly 256.0 + 400.0
#define Lz 256.0 + 300.0

//Vectors to store lambertians values for each shape
std::vector<double> coneLambertianVector,sphereLambertianVector,torusLambertianVector;

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog)
{
    ui->setupUi(this);
}

Dialog::~Dialog()
{
    delete ui;
}

/* Author: Dr. Steven S. Beauchemin
 * Date of creation: 12 October 2017
 *
 * Function:build_camera_matrix()
 * Parameter:
 *  In: dmatrix_t *E - Camera's "Eye" matrix
 *      dmatrix_t *G - Camera's "Gaze" matrix
 *  Out: None
 * Returns: dmatrix_t - A matrix of the computed built camera
 * Desc: Computes the camera matrix
 */
dmatrix_t *build_camera_matrix(dmatrix_t *E, dmatrix_t *G) {

    dmatrix_t N ; /* Viewing axis */

    N = *dmat_normalize(dmat_sub(E,G)) ;
    N.l = 3 ;

    dmatrix_t UP ;
    dmat_alloc(&UP,4,1) ;
    UP.l = 3 ;

    UP.m[1][1] = UPx ;
    UP.m[2][1] = UPy ;
    UP.m[3][1] = UPz ;
    UP.m[4][1] = 1.0 ;

    dmatrix_t U ;

    U = *dmat_normalize(dcross_product(&UP,&N)) ;

    dmatrix_t V ;
    V = *dcross_product(&N,&U) ;

    dmatrix_t Mv ; /* Build matrix M_v */
    dmat_alloc(&Mv,4,4) ;

    Mv.m[1][1] = U.m[1][1] ;
    Mv.m[1][2] = U.m[2][1] ;
    Mv.m[1][3] = U.m[3][1] ;
    Mv.m[1][4] = -1.0*((*E).m[1][1]*U.m[1][1] + (*E).m[2][1]*U.m[2][1] + (*E).m[3][1]*U.m[3][1]) ;

    Mv.m[2][1] = V.m[1][1] ;
    Mv.m[2][2] = V.m[2][1] ;
    Mv.m[2][3] = V.m[3][1] ;
    Mv.m[2][4] = -1.0*((*E).m[1][1]*V.m[1][1] + (*E).m[2][1]*V.m[2][1] + (*E).m[3][1]*V.m[3][1]) ;

    Mv.m[3][1] = N.m[1][1] ;
    Mv.m[3][2] = N.m[2][1] ;
    Mv.m[3][3] = N.m[3][1] ;
    Mv.m[3][4] = -1.0*((*E).m[1][1]*N.m[1][1] + (*E).m[2][1]*N.m[2][1] + (*E).m[3][1]*N.m[3][1]) ;

    Mv.m[4][1] = 0.0 ;
    Mv.m[4][2] = 0.0 ;
    Mv.m[4][3] = 0.0 ;
    Mv.m[4][4] = 1.0 ;

    dmatrix_t Mp ; /* Build matrix Mp */
    dmat_alloc(&Mp,4,4) ;
    Mp = *dmat_identity(&Mp) ;

    float a = -1.0*(FP + NP)/(FP - NP) ;
    float b = -2.0*(FP*NP)/(FP - NP) ;

    Mp.m[1][1] = NP ;
    Mp.m[2][2] = NP ;
    Mp.m[3][3] = a ;
    Mp.m[3][4] = b ;
    Mp.m[4][3] = -1.0 ;
    Mp.m[4][4] = 0.0 ;

    /* Build matrices T_1 and S_1 */

    /* Work out coordinates of near plane corners */

    float top = NP*tan(M_PI/180.0*THETA/2.0) ;
    float right = ASPECT*top ;
    float bottom = -top ;
    float left = -right ;

    dmatrix_t T1 ;
    dmat_alloc(&T1,4,4) ;

    T1 = *dmat_identity(&T1) ;
    T1.m[1][4] = -(right + left)/2.0 ;
    T1.m[2][4] = -(top + bottom)/2.0 ;

    dmatrix_t S1 ;
    dmat_alloc(&S1,4,4) ;

    S1 = *dmat_identity(&S1) ;
    S1.m[1][1] = 2.0/(right - left) ;
    S1.m[2][2] = 2.0/(top - bottom) ;

    /* Build matrices T2, S2, and W2 */

    dmatrix_t T2 ;
    dmatrix_t S2 ;
    dmatrix_t W2 ;

    dmat_alloc(&T2,4,4) ;
    dmat_alloc(&S2,4,4) ;
    dmat_alloc(&W2,4,4) ;

    T2 = *dmat_identity(&T2) ;
    S2 = *dmat_identity(&S2) ;
    W2 = *dmat_identity(&W2) ;

    T2.m[1][4] = 1.0 ;
    T2.m[2][4] = 1.0 ;

    S2.m[1][1] = W/2.0 ;
    S2.m[2][2] = H/2.0 ;

    W2.m[2][2] = -1.0 ;
    W2.m[2][4] = (double)H ;

    return dmat_mult(&W2,dmat_mult(&S2,dmat_mult(&T2,dmat_mult(&S1,dmat_mult(&T1,dmat_mult(&Mp,&Mv)))))) ;
}

/* Author: Dr. Steven S. Beauchemin
 * Date of creation: 12 October 2017
 *
 * Function:perspective_projection()
 * Parameter:
 *  In: dmatrix_t *P - A matrix pointer
 *  Out: None
 * Returns: dmatrix_t - A matrix of the computed perspective projection
 * Desc: Computes the perspective projection of a matrix and returns it
 */
dmatrix_t *perspective_projection(dmatrix_t *P) {

    (*P).m[1][1] /= (*P).m[4][1] ;
    (*P).m[2][1] /= (*P).m[4][1] ;
    (*P).m[3][1] /= (*P).m[4][1] ;
    (*P).m[4][1] /= (*P).m[4][1] ;

    return P ;
}

/* Author: Dr. Steven S. Beauchemin
 * Date of creation: 12 October 2017
 *
 * Function: camera_setup()
 * Parameter:
 *  In: None
 *  Out: None
 * Returns: dmatrix_t - A matrix of the camera
 * Desc: Sets up the camera matrix by using its E (eye) and G (gaze), and calling the build_camera_matrix() function
 */
dmatrix_t camera_setup(){
    dmatrix_t E ; /* The centre of projection for the camera */

    dmat_alloc(&E,4,1) ;

    E.m[1][1] = Ex ;
    E.m[2][1] = Ey ;
    E.m[3][1] = Ez ;
    E.m[4][1] = 1.0 ;

    dmatrix_t G ; /* Point gazed at by camera */

    dmat_alloc(&G,4,1) ;

    G.m[1][1] = Gx ;
    G.m[2][1] = Gy ;
    G.m[3][1] = Gz ;
    G.m[4][1] = 1.0 ;

    dmatrix_t C ; /* The camera matrix */

    dmat_alloc(&C,4,4) ;
    C = *build_camera_matrix(&E,&G) ;

    return C;
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: NormalVector()
 * Parameter:
 *  In: dmatrix_t *P1 - A matrix pointer P1
 *  In: dmatrix_t *P2 - A matrix pointer P2
 *  In: dmatrix_t *P3 - A matrix pointer P3
 *  In: dmatrix_t *P4 - A matrix pointer P4
 *  Out: None
 * Returns: std::vector<double> - A computed normal vector
 * Desc: Computes the normal vector of face with the corresponding points
 */
std::vector<double> NormalVector(dmatrix_t *P1, dmatrix_t *P2, dmatrix_t *P3, dmatrix_t *P4){
    //vector a = (P2 - P1)
    std::vector<double> a;
    a.push_back(P2->m[1][1] - P1->m[1][1]);
    a.push_back(P2->m[2][1] - P1->m[2][1]);
    a.push_back(P2->m[3][1] - P1->m[3][1]);
    double denominator = sqrt(pow(a.at(0),2)+pow(a.at(1),2)+pow(a.at(2),2));
    a.at(0) = a.at(0) / denominator;
    a.at(1) = a.at(1) / denominator;
    a.at(2) = a.at(2) / denominator;

    //vector b = (P3 - P2)
    std::vector<double> b;
    b.push_back(P3->m[1][1] - P2->m[1][1]);
    b.push_back(P3->m[2][1] - P2->m[2][1]);
    b.push_back(P3->m[3][1] - P2->m[3][1]);
    //At certain poles P3 = P2 so use b = (P4 - P1) instead
    if(b.at(0) == 0 && b.at(1) == 0 && b.at(2) == 0){
        b.push_back(P4->m[1][1] - P1->m[1][1]);
        b.push_back(P4->m[2][1] - P1->m[2][1]);
        b.push_back(P4->m[3][1] - P1->m[3][1]);
    }
    denominator = sqrt(pow(b.at(0),2)+pow(b.at(1),2)+pow(b.at(2),2));
    b.at(0) = b.at(0) / denominator;
    b.at(1) = b.at(1) / denominator;
    b.at(2) = b.at(2) / denominator;

    //vector normal = a x b
    std::vector<double> normal;
    normal.push_back((a.at(1)*b.at(2)) - (b.at(1)*a.at(2)));
    normal.push_back((a.at(2)*b.at(0)) - (b.at(2)*a.at(0)));
    normal.push_back((a.at(0)*b.at(1)) - (b.at(0)*a.at(1)));
    denominator = sqrt(pow(normal.at(0),2)+pow(normal.at(1),2)+pow(normal.at(2),2));
    normal.at(0) = normal.at(0) / denominator;
    normal.at(1) = normal.at(1) / denominator;
    normal.at(2) = normal.at(2) / denominator;
    return normal;
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: LightVectorToFace()
 * Parameter:
 *  In: dmatrix_t *P1 - A matrix pointer P1
 *  Out: None
 * Returns: std::vector<double> - A computed light vector
 * Desc: Computes the light vector of light with a point from face
 */
std::vector<double> LightVectorToFace(dmatrix_t *P1){
    //vector light = (point from face - light)
    std::vector<double> light;
    light.push_back(P1->m[1][1] - Lx);
    light.push_back(P1->m[2][1] - Ly);
    light.push_back(P1->m[3][1] - Lz);
    double denominator = sqrt(pow(light.at(0),2)+pow(light.at(1),2)+pow(light.at(2),2));
    light.at(0) = (light.at(0) / denominator);
    light.at(1) = (light.at(1) / denominator);
    light.at(2) = (light.at(2) / denominator);
    return light;
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: PolygonFilling()
 * Parameter:
 *  In: double polyX[] - y coordinates of the points in the face
 *  In: double polyX[] - x coordinates of the points in the face
 *  In: int lightIntensity - the lightintensity for the pen using the RGB model
 *  In: int shape - each shape corresponds to a number
 *  Out: None
 * Returns: None
 * Desc: 2D Polygon Scan Line Filling Algorithm
 */
void PolygonFilling(double polyX[], double polyY[], int lightIntensity, int shape){
    QPainter painter(this);
    int  nodes, nodeX[4], pixelX, pixelY, i, j, swap;

    //Loop through the rows of the image.
    for (pixelY=*std::min_element(polyY,polyY+4); pixelY<*std::max_element(polyY,polyY+4); pixelY++) {

        //Build a list of nodes.
        nodes=0; j=3;
        for (i=0; i<3; i++) {
            if (polyY[i]<(double) pixelY && polyY[j]>=(double) pixelY
                    ||  polyY[j]<(double) pixelY && polyY[i]>=(double) pixelY) {
                nodeX[nodes++]=(int) (polyX[i]+(pixelY-polyY[i])/(polyY[j]-polyY[i])
                                      *(polyX[j]-polyX[i])); }
            j=i;
        }

        //Sort the nodes, via a simple “Bubble” sort.
        i=0;
        while (i<nodes-1) {
            if (nodeX[i]>nodeX[i+1]) {
                swap=nodeX[i]; nodeX[i]=nodeX[i+1]; nodeX[i+1]=swap; if (i) i--;
            }
            else {
                i++;
            }
        }

        //Fill the pixels between node pairs.
        for (i=0; i<nodes; i+=2) {
            if   (nodeX[i]>=*std::max_element(polyX,polyX+4)) break;
            if   (nodeX[i+1]> *std::min_element(polyX,polyX+4) ) {
                if (nodeX[i]< *std::min_element(polyX,polyX+4) ) nodeX[i]=*std::min_element(polyX,polyX+4) ;
                if (nodeX[i+1]> *std::max_element(polyX,polyX+4)) nodeX[i+1]=*std::max_element(polyX,polyX+4);
                for (pixelX=nodeX[i]; pixelX<nodeX[i+1]; pixelX++){
                    if(shape == 1){
                        painter.setPen(QColor(1,1,lightIntensity)); // Torus
                    }
                    else if(shape == 2){
                        painter.setPen(QColor(lightIntensity,1,1)); // Sphere
                    }
                    else if(shape == 3){
                        painter.setPen(QColor(1,lightIntensity,1)); // Cone
                    }
                   painter.drawPoint(pixelX,pixelY);
                }
            }
        }
    }
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: DrawTorus()
 * Parameter:
 *  In: double outerRadius - The Torus outer radius
 *  In: double innerRadius - The Torus inner radius
 *  In: double translation - The Torus translation
 *  Out: None
 * Returns: std::vector<dmatrix_t> - A matrix with all the torus points from the parametric equation
 * Desc: Algorithm to draw the front face of a Torus. Computes points from the parametric equations and
 * stores them in a matrix
 */
std::vector<dmatrix_t> DrawTorus(double outerRadius, double innerRadius, double translation){
    //Variables and Initializations for the Torus' Algorithm
    double theta = 0, phi = 0; //Angles
    float DTOR = 0.01745329252; //(2 * PI / 36)
    dmatrix_t P1,P2,P3,P4; //Points
    std::vector<dmatrix_t> torusPVector;

    //Algorithm to Draw a Torus
    int du = 10, dv = 10; //u and v increments

    for (int u = 0;u < 360;u += du) {
        for (int v = 0;v < 360;v += dv) {

            //Segment 1
            theta = (u) * DTOR;
            phi   = (v) * DTOR;
            dmat_alloc(&P1,4,1) ;
            P1.m[1][1] = translation + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x1
            P1.m[2][1] = translation + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y1
            P1.m[3][1] = translation + innerRadius * sin(phi); //z1
            P1.m[4][1] = 1.0 ;

            //Segment 2
            theta = (u+du) * DTOR;
            phi   = (v) * DTOR;
            dmat_alloc(&P2,4,1) ;
            P2.m[1][1] = translation + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x2
            P2.m[2][1] = translation + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y2
            P2.m[3][1] = translation + innerRadius * sin(phi); //z2
            P2.m[4][1] = 1.0 ;

            //Segment 3
            theta = (u+du) * DTOR;
            phi   = (v+dv) * DTOR;
            dmat_alloc(&P3,4,1) ;
            P3.m[1][1] = translation + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x3
            P3.m[2][1] = translation + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y3
            P3.m[3][1] = translation + innerRadius * sin(phi); //z3
            P3.m[4][1] = 1.0 ;

            //Segment 4
            theta = (u) * DTOR;
            phi   = (v+dv) * DTOR;
            dmat_alloc(&P4,4,1) ;
            P4.m[1][1] = translation + cos(theta) * ( outerRadius + innerRadius * cos(phi) ); //x4
            P4.m[2][1] = translation + sin(theta) * ( outerRadius + innerRadius * cos(phi) ); //y4
            P4.m[3][1] = translation + innerRadius * sin(phi); //z4
            P4.m[4][1] = 1.0 ;

            //Calculate normal
            std::vector<double> normal = NormalVector(&P1,&P2,&P3,&P4);

            //Calculate light vector to face
            std::vector<double> light = LightVectorToFace(&P1);

            //Remove backfaces
            double toDraw = (P1.m[1][1] * normal.at(0)) + (P1.m[2][1] * normal.at(1)) + (P1.m[3][1] * normal.at(2));
            if(toDraw > 0){
                //Push light vector dot product normal vector to torusLambertianVector
                torusLambertianVector.push_back((light.at(0) * normal.at(0)) + (light.at(1) * normal.at(1)) + (light.at(2) * normal.at(2)));
                //Push points to torusPVector
                torusPVector.push_back(P1);
                torusPVector.push_back(P2);
                torusPVector.push_back(P3);
                torusPVector.push_back(P4);
            }
        }
    }

    return torusPVector;
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: DrawSphere()
 * Parameter:
 *  In: double sphereRadius - The radius of the Sphere
 *  In: double translation - The translation of the Sphere
 *  Out: None
 * Returns: std::vector<dmatrix_t> - A matrix with all the sphere points from the parametric equation
 * Desc: Algorithm to draw the front face of a Sphere. Computes points from the parametric equations and stores
 * them in a matrix
 */
std::vector<dmatrix_t> DrawSphere(double sphereRadius, double translation){
    //Variables and Initializations for the Sphere's Algorithms
    double theta = 0, phi = 0; //Angles
    float DTOR = 0.01745329252; //(2 * PI / 36)
    dmatrix_t P1,P2,P3,P4; //Points
    std::vector<dmatrix_t> spherePVector;

    //Algorithm to Draw a Sphere
    int di = 10, dj = 10; //i and j increments

    for (int i = 0; i < 360; i += di){
        for (int j = 0; j < 180; j += dj){

            //Segment 1
            theta = i * DTOR;
            phi = j * DTOR;
            dmat_alloc(&P1,4,1) ;
            P1.m[1][1] = translation + cos(theta) * sin(phi) * sphereRadius; //x1
            P1.m[2][1] = translation + sin(theta) * sin(phi) * sphereRadius; //y1
            P1.m[3][1] = translation + cos(phi) * sphereRadius; //z1
            P1.m[4][1] = 1.0 ;

            //Segment 2
            theta = (i + di) * DTOR;
            phi = j * DTOR;
            dmat_alloc(&P2,4,1) ;
            P2.m[1][1] = translation + cos(theta) * sin(phi) * sphereRadius; //x2
            P2.m[2][1] = translation + sin(theta) * sin(phi) * sphereRadius; //y2
            P2.m[3][1] = translation + cos(phi) * sphereRadius; //z2
            P2.m[4][1] = 1.0 ;

            //Segment 3
            theta = (i + di) * DTOR;
            phi = (j + dj) * DTOR;
            dmat_alloc(&P3,4,1) ;
            P3.m[1][1] = translation + cos(theta) * sin(phi) * sphereRadius; //x3
            P3.m[2][1] = translation + sin(theta) * sin(phi) * sphereRadius; //y3
            P3.m[3][1] = translation + cos(phi) * sphereRadius; //z3
            P3.m[4][1] = 1.0 ;

            //Segment 4
            theta = i * DTOR;
            phi = (j + dj) * DTOR;
            dmat_alloc(&P4,4,1) ;
            P4.m[1][1] = translation + cos(theta) * sin(phi) * sphereRadius; //x4
            P4.m[2][1] = translation + sin(theta) * sin(phi) * sphereRadius; //y4
            P4.m[3][1] = translation + cos(phi) * sphereRadius; //z4
            P4.m[4][1] = 1.0 ;

            //Calculate normal
            std::vector<double> normal = NormalVector(&P1,&P2,&P3,&P4);

            //Calculate light vector to face
            std::vector<double> light = LightVectorToFace(&P1);

            //Remove Backfaces
            double toDraw = (P1.m[1][1] * normal.at(0)) + (P1.m[2][1] * normal.at(1)) + (P1.m[3][1] * normal.at(2));
            if(toDraw > 0){
                //Push light vector dot product normal vector to sphereLambertianVector
                sphereLambertianVector.push_back((light.at(0) * normal.at(0)) + (light.at(1) * normal.at(1)) + (light.at(2) * normal.at(2)));
                //Push points to spherePVector
                spherePVector.push_back(P1);
                spherePVector.push_back(P2);
                spherePVector.push_back(P3);
                spherePVector.push_back(P4);
            }
        }
    }

    return spherePVector;
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: DrawCone()
 * Parameter:
 *  In: double height - The height of the Cone
 *  In: double coneRadius - The radius of the Cone
 *  In: double translation - The translation of the Cone
 *  Out: None
 * Returns: std::vector<dmatrix_t> - A matrix with all the cone points from the parametric equation
 * Desc: Algorithm to draw the front face of a Cone. Computes points from the parametric equations and stores
 * them in a matrix
 */
std::vector<dmatrix_t> DrawCone(double height, double coneRadius, double translation){
    //Variables and Initializations for the Cone's Algorithms
    double theta = 0; //Angle
    double heightStep = 0; // Step of the Height
    float DTOR = 0.01745329252; //(2 * PI / 36)
    dmatrix_t P1,P2,P3,P4; //Points
    std::vector<dmatrix_t> conePVector;

    //Algorithm to Draw a Cone
    double di = 10, dj = 10; //i and j increments

    for (double i = 0; i < 360; i += di){
        for (double j = height; j > 0; j -= dj){

            //Segment 1
            theta = i * DTOR;
            heightStep = j;
            dmat_alloc(&P1,4,1) ;
            P1.m[1][1] = translation + (((height - heightStep) / height) * coneRadius * cos(theta)); //x1
            P1.m[2][1] = translation + (((height - heightStep) / height) * coneRadius * sin(theta)); //y1
            P1.m[3][1] = heightStep; //z1
            P1.m[4][1] = 1.0 ;

            //Segment 2
            theta = (i + di) * DTOR;
            heightStep = j;
            dmat_alloc(&P2,4,1) ;
            P2.m[1][1] = translation + (((height - heightStep) / height) * coneRadius * cos(theta)); //x2
            P2.m[2][1] = translation + (((height - heightStep) / height) * coneRadius * sin(theta)); //y2
            P2.m[3][1] = heightStep; //z2
            P2.m[4][1] = 1.0 ;

            //Segment 3
            theta = (i + di) * DTOR;
            heightStep = (j - dj);
            dmat_alloc(&P3,4,1) ;
            P3.m[1][1] = translation + (((height - heightStep) / height) * coneRadius * cos(theta)); //x3
            P3.m[2][1] = translation + (((height - heightStep) / height) * coneRadius * sin(theta)); //y3
            P3.m[3][1] = heightStep; //z3
            P3.m[4][1] = 1.0 ;

            //Segment 4
            theta = i * DTOR;
            heightStep = (j - dj);
            dmat_alloc(&P4,4,1) ;
            P4.m[1][1] = translation + (((height - heightStep) / height) * coneRadius * cos(theta)); //x4
            P4.m[2][1] = translation + (((height - heightStep) / height) * coneRadius * sin(theta)); //y4
            P4.m[3][1] = heightStep; //z4
            P4.m[4][1] = 1.0 ;

            //Calculate normal
            std::vector<double> normal = NormalVector(&P1,&P2,&P3,&P4);

            //Calculate light vector to face
            std::vector<double> light = LightVectorToFace(&P1);

            double toDraw = (P1.m[1][1] * normal.at(0)) + (P1.m[2][1] * normal.at(1)) + (P1.m[3][1] * normal.at(2));
            if(toDraw > 0){
                //Push light vector dot product normal vector to coneLambertianVector
                coneLambertianVector.push_back((light.at(0) * normal.at(0)) + (light.at(1) * normal.at(1)) + (light.at(2) * normal.at(2)));
                //Push points to conePVector
                conePVector.push_back(P1);
                conePVector.push_back(P2);
                conePVector.push_back(P3);
                conePVector.push_back(P4);
            }
        }
    }

    return conePVector;
}

/* Author: Tareq Tayeh
 * Date of creation: 12 October 2017
 *
 * Function: keyPressEvent()
 * Parameter:
 *  In: QKeyEvent *event - Any key event when the dialog is open
 *  Out: None
 * Returns: None
 * Desc: Quits the program when the key "q" is pressed
 */
void Dialog::keyPressEvent(QKeyEvent *event)
{
    switch(event->key())
    {
    case Qt::Key_Q:
        close();
        break;
    default:
        QDialog::keyPressEvent(event);
    }
}

/* Author: Tareq Tayeh
 * Date of creation: 8 November 2017
 *
 * Function: paintEvent()
 * Parameter:
 *  In: QPaintEvent - A paint event when the dialog is open
 *  Out: None
 * Returns: None
 * Desc: Dialog's painter function where the shape algorithms are called, perspective projections and light
 *       intensities are calculated, and the polygon filling algorithm function is called where the dots are "painted"
 */
void Dialog::paintEvent(QPaintEvent *e)
{
    //Camera's Matrix and Setup
    printf("Camera Matrix:\n") ;
    dmatrix_t C = camera_setup();
    write_dmatrix(&C) ;

    //Variables and Initializations for the Shapes
    dmatrix_t *PP1,*PP2, *PP3, *PP4; //Perspective Projections
    std::vector<dmatrix_t> torusPVector, spherePVector, conePVector; //Shape's points vector
    qreal pointx1, pointy1, pointx2, pointy2, pointx3, pointy3, pointx4, pointy4; //Shape's points

    //Draws the Torus
    int lightIntensity = 0, x = 0; //x indicates element position for vectors
    torusPVector = DrawTorus(220,55,256);
    for(int a = 0; a < torusPVector.size() - 3; a+=4){
        PP1 = perspective_projection(dmat_mult(&C,&(torusPVector.at(a))));
        PP2 = perspective_projection(dmat_mult(&C,&(torusPVector.at(a+1))));
        PP3 = perspective_projection(dmat_mult(&C,&(torusPVector.at(a+2))));
        PP4 = perspective_projection(dmat_mult(&C,&(torusPVector.at(a+3))));
        pointx1 = PP1->m[1][1];
        pointy1 = PP1->m[2][1];
        pointx2 = PP2->m[1][1];
        pointy2 = PP2->m[2][1];
        pointx3 = PP3->m[1][1];
        pointy3 = PP3->m[2][1];
        pointx4 = PP4->m[1][1];
        pointy4 = PP4->m[2][1];
        double polyY[4] = {pointy1,pointy2,pointy3,pointy4};
        double polyX[4] = {pointx1,pointx2,pointx3,pointx4};
        lightIntensity = torusLambertianVector.at(x) > 0 ? 255 * torusLambertianVector.at(x) : 0;
        x++;
        PolygonFilling(polyX,polyY,lightIntensity,1);
    }

    //Draws the Sphere
    lightIntensity = 0, x = 0;
    spherePVector = DrawSphere(100,256);
    for(int a = 0; a < spherePVector.size() - 3; a+=4){
        PP1 = perspective_projection(dmat_mult(&C,&(spherePVector.at(a))));
        PP2 = perspective_projection(dmat_mult(&C,&(spherePVector.at(a+1))));
        PP3 = perspective_projection(dmat_mult(&C,&(spherePVector.at(a+2))));
        PP4 = perspective_projection(dmat_mult(&C,&(spherePVector.at(a+3))));
        pointx1 = PP1->m[1][1];
        pointy1 = PP1->m[2][1];
        pointx2 = PP2->m[1][1];
        pointy2 = PP2->m[2][1];
        pointx3 = PP3->m[1][1];
        pointy3 = PP3->m[2][1];
        pointx4 = PP4->m[1][1];
        pointy4 = PP4->m[2][1];
        double polyY[4] = {pointy1,pointy2,pointy3,pointy4};
        double polyX[4] = {pointx1,pointx2,pointx3,pointx4};
        lightIntensity = sphereLambertianVector.at(x) > 0 ? 255 * sphereLambertianVector.at(x) : 0;
        x++;
        PolygonFilling(polyX,polyY,lightIntensity,2);
    }

    //Draw the Cone
    lightIntensity = 0, x = 0;
    conePVector = DrawCone(120,100,256);
    for(int a = 0; a < conePVector.size() - 3; a+=4){
        PP1 = perspective_projection(dmat_mult(&C,&(conePVector.at(a))));
        PP2 = perspective_projection(dmat_mult(&C,&(conePVector.at(a+1))));
        PP3 = perspective_projection(dmat_mult(&C,&(conePVector.at(a+2))));
        PP4 = perspective_projection(dmat_mult(&C,&(conePVector.at(a+3))));
        pointx1 = PP1->m[1][1];
        pointy1 = PP1->m[2][1];
        pointx2 = PP2->m[1][1];
        pointy2 = PP2->m[2][1];
        pointx3 = PP3->m[1][1];
        pointy3 = PP3->m[2][1];
        pointx4 = PP4->m[1][1];
        pointy4 = PP4->m[2][1];
        double polyY[4] = {pointy1,pointy2,pointy3,pointy4};
        double polyX[4] = {pointx1,pointx2,pointx3,pointx4};
        lightIntensity = coneLambertianVector.at(x) > 0 ? 255 * coneLambertianVector.at(x) : 0;
        x++;
        PolygonFilling(polyX,polyY,lightIntensity,3);
    }

    //Freeing Memory
    delete (PP1->m);
    delete (PP2->m);
    delete (PP3->m);
    delete (PP4->m);
    delete PP1;
    delete PP2;
    delete PP3;
    delete PP4;
}

/* Author: Tareq Tayeh
 * Date of creation: 21 September 2017
 *
 * Function: Bresenham()
 * Parameter:
 *  In: int x1 - coordinate x1 as an integer
 *      int y1 - coordinate y1 as an integer
 *      int x2 - coordinate x2 as an integer
 *      int y2 - coordinate y2 as an integer
 *  Out: None
 * Returns: None
 * Desc: Bresenham's line algorithm, and where the dots are "painted" on the dialog
 */
void Dialog::Bresenham(int x1, int y1, int x2, int y2)
{
  QPainter painter(this);
  int dx, dy, error, ystep, y, maxX;

  const bool steep = (fabs(y2 - y1) > fabs(x2 - x1));
  if(steep)
  {
    std::swap(x1, y1);
    std::swap(x2, y2);
  }

  if(x1 > x2)
  {
    std::swap(x1, x2);
    std::swap(y1, y2);
  }

  dx = x2 - x1;
  dy = fabs(y2 - y1);

  error = dx / 2.0f;
  ystep = (y1 < y2) ? 1 : -1;
  y = y1;

  maxX = x2;

  for(int x=x1; x<maxX; x++)
  {
    if(steep)
    {
        painter.drawPoint(y,x);
    }
    else
    {
        painter.drawPoint(x,y);
    }

    error -= dy;
    if(error < 0)
    {
        y += ystep;
        error += dx;
    }
  }
}
