#include <iostream>
#include <cmath>
#include <fstream>

//====================PROTOTYPES====================
void mattimes(double t, double a[][3], double m[][3]);
void matsub(double a[][3], double b[][3], double c[][3]);
void matvecmult(double m[][3], double *v, double *prod);
void probTwo(void);
void tests(void);
double convolve(double *f1, int len1, double *f2, int len2, double* y);
double f(int);
double h(int);

double f1[] = {0,1,2,3,2,1};
double len1 = 6;
double f2[] = {-2,-2,-2,-2,-2,-2,-2};
double len2 = 7;
double f3[] = {1,-1,1,-1};
double len3 = 4;
double f4[] = {0,0,0,-3,-3};
double len4 = 5;

//====================F & H====================
double f(double t) {
	return sin(4 * 3.14159 * t);
}

double h(double t) {
	return (-2 / 153)*exp(-2 * t) + 0.16392*exp(-0.5*t)*cos(6 * t - 85.4261);
}
using namespace std;
//====================MATRIX FUNCTIONS====================
void mattimes(double t, double a[][3], double m[][3])
{
int i, j;
for(i = 0; i < 3; i++)
  for(j = 0; j < 3; j++)
    m[i][j] = t*a[i][j];
}

void matsub(double a[][3], double b[][3], double c[][3])
{
  int i, j;
  for(i = 0; i < 3; i++)
    for(j = 0; j < 3; j++)
      c[i][j] = a[i][j]-b[i][j];
}
void matvecmult(double m[][3], double *v, double *prod)
{
  double sum;
  int i, j;
  for(i = 0; i < 3; i++)
  {
    sum = 0;
    for(j = 0; j < 3; j++)
    {
      sum += m[i][j] * v[j];
    }
    prod[i] = sum;
  }
}


//====================OLD CODE====================
void probTwo(double * y0)
{
  double deltat, lbound, ubound, a0, a1, a2, x[3], n, SAMPLE_RATE;
  double mat1[3][3];
  double mat2[3][3];
  double I[3][3] = { { 1, 0, 0}, {0, 1, 0}, {0, 0, 1} };
   
  cout << "Enter number of iterations: ";
  cin >> n;
  cout << "Enter lower bound: ";
  cin >> lbound;
  cout << "Enter upper bound: ";
  cin >> ubound;

  x[0] = -2;
  x[1] = 3;
  x[2] = -1.7;

  a0 = 72.5;//replace these
  a1 = 38.25;
  a2 = 3; //s^2
  
  double a[3][3] = { { 0, -1, 0}, {0, 0, -1}, {a0, a1, a2} };

  /**********INPUT VERIFY**********/
  cout << "Conditions entered for problem 3:" << endl;
  cout << "n = " << n << '\n' << "Lower Bound = " << lbound << '\n' << "Upper Bound = " << ubound << '\n' << "Initial Conditions y0: " << x[0] << '\n' << "Initial Conditions 'y0: " << x[1] << '\n' << "Initial Conditions ''y0: " << x[2] << '\n' << endl;

  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      cout << a[i][j] << " ";
    }
    cout << '\n';
  }

  /**********ITERATIONS**********/
  y0[0] = x[0];


  float t = 0;
  deltat =0.001;
  SAMPLE_RATE=n/deltat;
  
  for(int i = 1; i < SAMPLE_RATE; i++) 
  {
    t=i*deltat;
    mattimes(deltat, a, mat1);
    matsub(I, mat1, mat2);
    matvecmult(mat2, x, x);
    y0[i] = x[0];
  }
}

//====================CONVOLUTION====================
double convolve(double *f1, int len1, double *f2, int len2, double* y)
{
  // int totLen = len1 + len2 - 1;
	// int val, val2;
	// double temp;
	// for (int n = 0; n < totLen; n++) {
	// 	temp = 0;
	// 	if (n >= len2) { val = n - len2 + 1; }
	// 	else { val = 0; }
	// 	if (n >= len1) { val2 = len1 - 1; }
	// 	else { val2 = n; }
	// 	for (int i = val; i <= val2; i++) {
	// 		temp += (f1[i] * f2[n - i]);
	// 	}
	// 	y[n] = temp;
	// }
	// return totLen;
  int max_lenght=len1+len2-1;
  double temp_f1, temp_f2;
  int i,j;
  for(i=0;i<max_lenght;i++)
  {
    y[i]=0;
    for(j=0;j<=i;j++)
    {
      temp_f1=j<len1?f1[j]:0;
      temp_f2=(i-j)<len2?f2[i-j]:0;
      y[i]+=temp_f1*temp_f2;
    }
  }
  return max_lenght;
}

//====================MAIN====================
int main() 
{
  int SIZE = 10000;
  double ft[SIZE];
  double ht[SIZE];
  double fh[SIZE*2-1];
  double y0[SIZE];
  double y[SIZE*2-1];
  float T = 0.001;
  ofstream outfile("Problem3_2");
	ofstream outfilef("f");
	ofstream outfileh("h");
	ofstream outfilefh("fh");
	ofstream outfiley0("y0");

  tests();//comment out to skip test convolutions

  //create f(t)
	for (int i = 0; i < SIZE; i ++) {
		ft[i] = f(i*0.001);
		outfilef << (i*0.001) << " " << ft[i] << endl;
	}

	//create h(t)
	for (int i = 0; i < SIZE; i++) {
		ht[i] = h(i*0.001);
		outfileh << (i*0.001) << " " << ht[i] << endl;
	}

	//convolute f and h to create fh
	convolve(ft, SIZE, ht, SIZE, fh);
	for (int i = 0; i < SIZE*2-1; i++) {
		outfilefh << (i*0.001) << " " << (0.001*fh[i]) << endl;
	}
	
	//create y0
	probTwo(y0);
	for (int i = 0; i < SIZE; i++) {
		outfiley0 << (i*0.001) << " " << y0[i] << endl;
	}

	for (int t = 0; t < SIZE; t++) {
		y[t] = (0.001 * fh[t]) + y0[t];
		outfile << (t*0.001) << "   " << y[t] << endl;

	}
	for (int t = 0; t < SIZE-1; t++) {
		y[t+SIZE] = (0.001*fh[t+SIZE]);

		outfile << t*0.001 << " " << (y[t+SIZE])<<endl;
	
	}
    return 0;
}

void tests() {
  
  ofstream outfilea("con1");
  ofstream outfileb("con2");
  ofstream outfilec("con3");
  ofstream outfiled("con4");
  ofstream outfilee("con5");

	double leny;
	double y[100] = { 0 };
	double f1[] = { 0,1,2,3,2,1 };
	int len1 = 6;
	double f2[] = { -2,-2,-2,-2,-2,-2,-2 };
	int len2 = 7;
	double f3[] = { 1,-1,1,-1 };
	int len3 = 4;
	double f4[] = { 0,0,0,-3,-3 };
	int len4 = 5;

	leny = convolve(f1, len1, f1, len1, y);
	cout << "(a) f1*f1 = ";
	for (int i = 0; i < leny; i++) {
		cout << y[i] << " ";
    outfilea << y[i] << " ";
	}
	cout << endl;
	leny = convolve(f1, len1, f2, len2, y);
	cout << "(b) f1*f2 = ";
	for (int i = 0; i < leny; i++) {
		cout << y[i] << " ";
    outfileb << y[i] << " ";

	}
	cout << endl;
	leny = convolve(f1, len1, f3, len3, y);
	cout << "(c) f1*f3 = ";
	for (int i = 0; i < leny; i++) {
		cout << y[i] << " ";
    outfilec << y[i] << " ";

	}
	cout << endl;
	leny = convolve(f2, len2, f3, len3, y);
	cout << "(d) f2*f3 = ";
	for (int i = 0; i < leny; i++) {
		cout << y[i] << " ";
    outfiled << y[i] << " ";

	}
	cout << endl;
	leny = convolve(f1, len1, f4, len4, y);
	cout << "(e) f1*f4 = ";
	for (int i = 0; i < leny; i++) {
		cout << y[i] << " ";
    outfilee << y[i] << " ";

	}
	cout << '\n' << "Convolution tests complete" << endl << endl;
}