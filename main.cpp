/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #3
 * Date: 9 November 2017
 * Program purpose: The purpose of this program is to design and implementat a program that allows a user to shade
 *                  objects (spheres, cones, and torii) in various colours. The program should scale the parametric
 *                  objects within a scene (using 3D homogeneous transformation matrices) and render them using
 *                  a simple Lambertian shading model that assumes diffuse light reflection from objects.
 * File: main.cpp
 */


#include "dialog.h"
#include <QApplication>

/* Author: Tareq Tayeh
 * Date of creation: 6 November 2017
 *
 * Function: main()
 * Purpose: Instantiates an instance of Dialog, and sets its size, background and shows it.
 */
int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Dialog w;

    w.setFixedSize(512, 512); //Opening a window of size 512 by 512 pixels
    w.setStyleSheet("QWidget { background-color: white; }"); //Setting the window background color to white
    w.show(); //Show dialog

    return a.exec();
}
