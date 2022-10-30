#include "Render.h"
#include <Windows.h>
#include <GL\GL.h>
#include <GL\GLU.h>
#include <corecrt_math.h>
#include <iostream>

#define pairDD pair<double, double>
#define _height 9
#define fi 0.01
#define oboroty 0.7

using namespace std;


pairDD findCircumCenter(pairDD P, pairDD Q, pairDD R);
pairDD lineLineIntersection(double a1, double b1, double c1, double a2, double b2, double c2);
void perpendicularBisectorFromLine(pairDD P, pairDD Q, double& a, double& b, double& c);
void lineFromPoints(pairDD P, pairDD Q, double& a, double& b, double& c);
void roundFlor(double step, double pointA[], double pointB[]);
void roundWall(double step, double pointA[], double pointB[], double pointA1[], double pointB1[], double hight);
void roundInWall(double step, double pointA[], double pointB[], double pointC[], double pointA1[], double pointB1[], double pointC1[]);
void func();
void twisted(double point[], double Fi, double height, double oborots);
void twisted(double point[], double Fi);
void roundInFloor(double step, double pointA[], double pointB[], double pointC[], double E[], double H[]);

void Render(double delta_time)
{
	func();
}

void func()
{
	double A[] = { 4, -3,  0 };
	double B[] = { 7, 1, 0 };
	double C[] = { 3, 5, 0 };
	double D[] = { 0, 1, 0 };
	double E[] = { -4, 3, 0 };
	double F[] = { -6, -2, 0 };
	double G[] = { -2, -6, 0 };
	double H[] = { 0, -1, 0 };

	double AA[] = { 4, -3,  0 };
	double BB[] = { 7, 1, 0 };
	double CC[] = { 3, 5, 0 };
	double DD[] = { 0, 1, 0 };
	double EE[] = { -4, 3, 0 };
	double FF[] = { -6, -2, 0 };
	double GG[] = { -2, -6, 0 };
	double HH[] = { 0, -1, 0 };

	double A1[] = { 4, -3,   _height };
	double B1[] = { 7, 1, _height };
	double C1[] = { 3, 5, _height };
	double D1[] = { 0, 1, _height };
	double E1[] = { -4, 3, _height };
	double F1[] = { -6, -2, _height };
	double G1[] = { -2, -6, _height };
	double H1[] = { 0, -1, _height };

	glBegin(GL_TRIANGLES);
	glColor3d(0.2, 0.7, 0);
	glVertex3dv(G);
	glVertex3dv(F);
	glVertex3dv(H);
	glEnd();



	glBegin(GL_TRIANGLES);
	glColor3d(0.9, 0.9, 0.9);
	glVertex3dv(F);
	glVertex3dv(E);
	glVertex3dv(D);
	glEnd();


	glBegin(GL_TRIANGLES);
	glColor3d(0.4, 0.7, 0.9);
	glVertex3dv(F);
	glVertex3dv(D);
	glVertex3dv(H);
	glEnd();




	glBegin(GL_TRIANGLES);
	glColor3d(0.7, 0.5, 0.2);
	glVertex3dv(A);
	glVertex3dv(D);
	glVertex3dv(H);
	glEnd();

	roundFlor(0.001, F, G);
	double pointC[] = { 4, 2 ,0 };
	double pointC1[] = { 4, 2, 0 };
	double pointC2[] = { 4, 2, 0 };
	roundInFloor(0.01, B, C, pointC, A, D);

	double step = abs(_height / (oboroty / fi));


	for (double i = 0; i < _height - step; i += step)
	{
		twisted(AA, fi, _height, oboroty);
		twisted(BB, fi, _height, oboroty);
		twisted(CC, fi, _height, oboroty);
		twisted(DD, fi, _height, oboroty);
		twisted(EE, fi, _height, oboroty);
		twisted(FF, fi, _height, oboroty);
		twisted(GG, fi, _height, oboroty);
		twisted(HH, fi, _height, oboroty);
		twisted(pointC2, fi, _height, oboroty);


		glBegin(GL_QUADS);
		glColor3d(0.1, 0.1, 0.9);
		glVertex3dv(E);
		glVertex3dv(EE);
		glVertex3dv(DD);
		glVertex3dv(D);
		glEnd();

		glBegin(GL_QUADS);
		glColor3d(0.2, 0.2, 0.9);
		glVertex3dv(C);
		glVertex3dv(CC);
		glVertex3dv(DD);
		glVertex3dv(D);
		glEnd();



		glBegin(GL_QUADS);
		glColor3d(0.4, 0.4, 0.9);
		glVertex3dv(B);
		glVertex3dv(BB);
		glVertex3dv(AA);
		glVertex3dv(A);
		glEnd();

		glBegin(GL_QUADS);
		glColor3d(0.5, 0.5, 0.9);
		glVertex3dv(A);
		glVertex3dv(AA);
		glVertex3dv(HH);
		glVertex3dv(H);
		glEnd();

		glBegin(GL_QUADS);
		glColor3d(0.6, 0.6, 0.9);
		glVertex3dv(G);
		glVertex3dv(GG);
		glVertex3dv(HH);
		glVertex3dv(H);
		glEnd();



		glBegin(GL_QUADS);
		glColor3d(0.8, 0.8, 0.9);
		glVertex3dv(F);
		glVertex3dv(FF);
		glVertex3dv(EE);
		glVertex3dv(E);
		glEnd();

		glBegin(GL_QUADS);
		glColor3d(0.9, 0.9, 0.9);
		glVertex3dv(E);
		glVertex3dv(EE);
		glVertex3dv(DD);
		glVertex3dv(D);
		glEnd();

		roundWall(0.01, F, G, FF, GG, _height);
		roundInWall(0.01, B, C, pointC1, BB, CC, pointC2);


	}



	twisted(A1, oboroty);
	twisted(B1, oboroty);
	twisted(C1, oboroty);
	twisted(D1, oboroty);
	twisted(E1, oboroty);
	twisted(F1, oboroty);
	twisted(G1, oboroty);
	twisted(H1, oboroty);
	twisted(pointC, oboroty);



	glBegin(GL_TRIANGLES);
	glColor3d(0.2, 0.7, 0);
	glVertex3dv(G1);
	glVertex3dv(F1);
	glVertex3dv(H1);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0.4, 0.7, 0.9);
	glVertex3dv(F1);
	glVertex3dv(D1);
	glVertex3dv(H1);
	glEnd();

	glBegin(GL_TRIANGLES);
	glColor3d(0.9, 0.9, 0.9);
	glVertex3dv(F1);
	glVertex3dv(E1);
	glVertex3dv(D1);
	glEnd();


	glBegin(GL_TRIANGLES);
	glColor3d(0.7, 0.5, 0.2);
	glVertex3dv(A1);
	glVertex3dv(D1);
	glVertex3dv(H1);
	glEnd();

	roundFlor(0.001, F1, G1);
	pointC[2] = _height;
	roundInFloor(0.01, B1, C1, pointC, A1, D1);
}




// Function to find the line given two points
void lineFromPoints(pairDD P, pairDD Q, double& a, double& b, double& c)
{
	a = Q.second - P.second;
	b = P.first - Q.first;
	c = a * (P.first) + b * (P.second);
}

// Function which converts the input line to its
// perpendicular bisector. It also inputs the points
// whose mid-point lies on the bisector
void perpendicularBisectorFromLine(pairDD P, pairDD Q, double& a, double& b, double& c)
{
	pairDD mid_point = make_pair((P.first + Q.first) / 2,
		(P.second + Q.second) / 2);

	// c = -bx + ay
	c = -b * (mid_point.first) + a * (mid_point.second);

	double temp = a;
	a = -b;
	b = temp;
}

// Returns the intersection point of two lines
pairDD lineLineIntersection(double a1, double b1, double c1, double a2, double b2, double c2)
{
	double determinant = a1 * b2 - a2 * b1;
	if (determinant == 0)
	{
		// The lines are parallel. This is simplified
		// by returning a pair of FLT_MAX
		return make_pair(FLT_MAX, FLT_MAX);
	}

	else
	{
		double x = (b2 * c1 - b1 * c2) / determinant;
		double y = (a1 * c2 - a2 * c1) / determinant;
		return make_pair(x, y);
	}
}

pairDD findCircumCenter(pairDD P, pairDD Q, pairDD R)
{
	// Line PQ is represented as ax + by = c
	double a, b, c;
	lineFromPoints(P, Q, a, b, c);

	// Line QR is represented as ex + fy = g
	double e, f, g;
	lineFromPoints(Q, R, e, f, g);

	// Converting lines PQ and QR to perpendicular
	// vbisectors. After this, L = ax + by = c
	// M = ex + fy = g
	perpendicularBisectorFromLine(P, Q, a, b, c);
	perpendicularBisectorFromLine(Q, R, e, f, g);

	// The point of intersection of L and M gives
	// the circumcenter
	pairDD circumcenter =
		lineLineIntersection(a, b, c, e, f, g);




	return circumcenter;
}

void roundFlor(double step, double pointA[], double pointB[]) //крышка
{
	glColor3d(0, 0, 0);

	double vect_AB[] = { pointA[0] - pointB[0], pointA[1] - pointB[1] };
	double centre[] = { (pointA[0] + pointB[0]) / 2,  (pointA[1] + pointB[1]) / 2 };

	double length = sqrt(vect_AB[0] * vect_AB[0] + vect_AB[1] * vect_AB[1]);

	double radius = length / 2;

	double Fi = acos(vect_AB[0] / (length + 1));


	glBegin(GL_TRIANGLE_FAN);

	if (Fi > 2.5)
	{
		Fi = -Fi;
		for (double i = Fi - 0.65; i <= 2.5 + Fi; i += step)
		{

			double point[] = { radius * cos(i) + centre[0], radius * sin(i) + centre[1], pointA[2] };
			glVertex3dv(point);

		}
	}
	else
	{
		for (double i = Fi + 0.137; i <= 3.283 + Fi; i += step)
		{

			double point[] = { radius * cos(i) + centre[0], radius * sin(i) + centre[1], pointA[2] };
			glVertex3dv(point);

		}
	}



	glEnd();
}

void roundWall(double step, double pointA[], double pointB[], double pointA1[], double pointB1[], double hight)
{


	double vect_AB[] = { pointA[0] - pointB[0], pointA[1] - pointB[1] };
	double centre[] = { (pointA[0] + pointB[0]) / 2,  (pointA[1] + pointB[1]) / 2 };

	double vect_A1B1[] = { pointA1[0] - pointB1[0], pointA1[1] - pointB1[1] };
	double centre1[] = { (pointA1[0] + pointB1[0]) / 2,  (pointA1[1] + pointB1[1]) / 2 };

	double length = sqrt(vect_AB[0] * vect_AB[0] + vect_AB[1] * vect_AB[1]);

	double radius = length / 2;

	double Fi = acos(vect_AB[0] / (length));
	double Fi1 = acos(vect_A1B1[0] / (length));
	double delta = Fi1 - Fi;
	glBegin(GL_QUADS);
	glColor3d(0.8, 0.1, 0.3);

	for (double i = -Fi - 1.57; i < 1.57 - Fi; i += step)
	{

		double point[] = { radius * cos(i) + centre[0], radius * sin(i) + centre[1], pointA[2] };
		double point1[] = { radius * cos(i + delta) + centre1[0], radius * sin(i + delta) + centre1[1], pointA1[2] };
		double point2[] = { radius * cos(i + step) + centre[0], radius * sin(i + step) + centre[1], pointA[2] };
		double point3[] = { radius * cos(i + delta + step) + centre1[0], radius * sin(i + delta + step) + centre1[1], pointA1[2] };
		glVertex3dv(point);
		glVertex3dv(point2);
		glVertex3dv(point3);
		glVertex3dv(point1);




	}


	glEnd();


}

void roundInFloor(double step, double pointA[], double pointB[], double pointC[], double E[], double H[]) // закругленность в крышки
{
	pairDD P = make_pair(pointA[0], pointA[1]);
	pairDD Q = make_pair(pointB[0], pointB[1]);
	pairDD R = make_pair(pointC[0], pointC[1]);
	pairDD pointO = findCircumCenter(P, Q, R);



	double c = sqrt(pow((pointB[0] - pointA[0]), 2) + pow((pointB[1] - pointA[1]), 2));//длины
	double a = sqrt(pow((pointC[0] - pointB[0]), 2) + pow((pointC[1] - pointB[1]), 2));//сторон
	double b = sqrt(pow((pointA[0] - pointC[0]), 2) + pow((pointA[1] - pointC[1]), 2));//треугольника
	double s = 0.5 * abs(pointA[0] * (pointB[1] - pointC[1]) - pointA[1] * (pointB[0] - pointC[0]) + pointB[0] * pointC[1] - pointB[1] * pointC[0]);//площадь треуг

	double radius = (a * b * c) / (4 * s);//радиус описанной окр

	double Fi1 = -0.3 - acos((-pointA[0]) / (sqrt(pow(pointA[0], 2) + pow(pointA[1], 2))));

	double Fi2 = 0.7 - acos((-pointB[0]) / (sqrt(pow(pointB[0], 2) + pow(pointB[1], 2))));




	glBegin(GL_TRIANGLE_FAN);
	glColor3d(0.5, 0.5, 0.5);
	glVertex3dv(E);


	for (double i = Fi1; i < Fi2; i += step)
	{
		double point[] = { radius * cos(i) + pointO.first, radius * sin(i) + pointO.second, pointA[2] };

		glVertex3dv(point);




	}
	glVertex3dv(H);
	glEnd();



	glBegin(GL_TRIANGLE_FAN);

	glVertex3dv(H);


	for (double i = Fi1; i < Fi2; i += step)
	{
		double point[] = { radius * cos(i) + pointO.first, radius * sin(i) + pointO.second, pointA[2] };

		glVertex3dv(point);




	}
	glVertex3dv(E);
	glEnd();


}


void roundInWall(double step, double pointA[], double pointB[], double pointC[], double pointA1[], double pointB1[], double pointC1[]) //вогнутость стенки
{
	pairDD P = make_pair(pointA[0], pointA[1]);
	pairDD Q = make_pair(pointB[0], pointB[1]);
	pairDD R = make_pair(pointC[0], pointC[1]);
	pairDD pointO = findCircumCenter(P, Q, R);

	pairDD P1 = make_pair(pointA1[0], pointA1[1]);
	pairDD Q1 = make_pair(pointB1[0], pointB1[1]);
	pairDD R1 = make_pair(pointC1[0], pointC1[1]);
	pairDD pointO1 = findCircumCenter(P1, Q1, R1);


	double c = sqrt(pow((pointB[0] - pointA[0]), 2) + pow((pointB[1] - pointA[1]), 2));//длины
	double a = sqrt(pow((pointC[0] - pointB[0]), 2) + pow((pointC[1] - pointB[1]), 2));//сторон
	double b = sqrt(pow((pointA[0] - pointC[0]), 2) + pow((pointA[1] - pointC[1]), 2));//треугольника
	double s = 0.5 * abs(pointA[0] * (pointB[1] - pointC[1]) - pointA[1] * (pointB[0] - pointC[0]) + pointB[0] * pointC[1] - pointB[1] * pointC[0]);//площадь треуг

	double radius = (a * b * c) / (4 * s);//радиус описанной окр

	double Fi1 = -0.3 - acos((-pointA[0]) / (sqrt(pow(pointA[0], 2) + pow(pointA[1], 2))));

	double Fi2 = 0.7 - acos((-pointB[0]) / (sqrt(pow(pointB[0], 2) + pow(pointB[1], 2))));

	double Fi3 = -0.3 - acos((-pointA1[0]) / (sqrt(pow(pointA1[0], 2) + pow(pointA1[1], 2))));

	double Fi4 = 0.7 - acos((-pointB1[0]) / (sqrt(pow(pointB1[0], 2) + pow(pointB1[1], 2))));


	double delta1 = Fi3 - Fi1;
	double delta2 = Fi4 - Fi2;

	glBegin(GL_POLYGON);



	for (double i = Fi1; i < Fi2; i += step)
	{
		double point[] = { radius * cos(i) + pointO.first, radius * sin(i) + pointO.second, pointA[2] };
		double point1[] = { radius * cos(i + delta1) + pointO1.first, radius * sin(i + delta1) + pointO1.second, pointA1[2] };
		double point2[] = { radius * cos(i + step) + pointO.first, radius * sin(i + step) + pointO.second, pointA[2] };
		double point3[] = { radius * cos(i + step + delta2) + pointO1.first, radius * sin(i + step + delta2) + pointO1.second, pointA1[2] };
		glVertex3dv(point);
		glVertex3dv(point2);
		glVertex3dv(point3);
		glVertex3dv(point1);



	}
	glEnd();
}
void twisted(double point[], double Fi, double height, double oborots)
{
	double deltaZ = abs(height / (oborots / Fi));
	double x = point[0], y = point[1];
	point[0] = (x * cos(Fi) - y * sin(Fi));
	point[1] = (x * sin(Fi) + y * cos(Fi));
	point[2] += deltaZ;

}

void twisted(double point[], double Fi) // две перегруженные функции
{

	double x = point[0], y = point[1];
	point[0] = (x * cos(Fi) - y * sin(Fi));
	point[1] = (x * sin(Fi) + y * cos(Fi));


}