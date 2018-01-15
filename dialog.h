/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #3
 * Date: 9 November 2017
 * Program purpose: The purpose of this program is to design and implementat a program that allows a user to shade
 *                  objects (spheres, cones, and torii) in various colours. The program should scale the parametric
 *                  objects within a scene (using 3D homogeneous transformation matrices) and render them using
 *                  a simple Lambertian shading model that assumes diffuse light reflection from objects.
 * File: dialog.h
 */

#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <QKeyEvent>
#include <QtGui>
#include <QtCore>
#include <iostream>
#include <QVector3D>

namespace Ui {
class Dialog;
}

class Dialog : public QDialog
{
    Q_OBJECT

public:
    explicit Dialog(QWidget *parent = 0);
    ~Dialog();
    void keyPressEvent(QKeyEvent *event);

protected:
    void paintEvent(QPaintEvent *e);
    void Bresenham( int x1, int y1, int x2, int y2);

private:
    Ui::Dialog *ui;
};

#endif // DIALOG_H
