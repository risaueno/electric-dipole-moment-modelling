#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

/* Define parameters for coax cable
---------------------------------------------------------------------*/
const int res = 10;                             //Pixels per 1cm. Increase for higher resolution
const int bsize = 2;                            //Inner bar thickness
const int vsize = 2;                            //Vacuum region thickness
const int tsize = 1;                            //Outer tube thickness
const int Vbar = 10;                            //Inner bar voltage
const int side = res*(bsize+2*vsize+2*tsize);   //Total pixels per side

/* Function prototypes
---------------------------------------------------------------------*/
void setvoltage (double[][side], double, int);  //Set initial potential matrix
void truefalse (bool, int);                     //Set T/F in bool matrix
void update(int);                               //Update potential in vacuum (input method 1 (update whole mesh at a time) or method 2 (update element by element)
double sum();                                   //Sum all elements in potential matrix
double Vfield(int, ofstream&);                  //Create potential field and returns no. of iterations performed
void Efield();                                  //Create E-field
void printmatrix(double[][side], ofstream&);    //Output matrix
void printmatrixb(ofstream&);                   //Output bool matrix
void extract1D(double[][side], ofstream&);      //Output 1D cross-section of potential

/* Initialise matrices
---------------------------------------------------------------------*/
double a[side][side] = {0};         //Potential matrix
double a0[side][side] = {0};        //Copy of potential matrix
bool b[side][side] = {false};       //Bool matrix
double Ex[side][side] = {0};        //Electric field i component
double Ey[side][side] = {0};        //Electric field j component
double E[side][side] = {0};         //Electric field magnitude

int main()
{
    //CREATE .TXT FILES
    ofstream txtfile0("Bool matrix.txt");
    ofstream txtfile1("Iteration steps.txt");
    ofstream txtfile2("Potential matrix.txt");
    ofstream txtfile3("1D cross-section.txt");
    ofstream txtfile4("E-field matrix.txt");

    //CREATE INITIAL POTENTIAL MATRIX
    setvoltage(a,Vbar,(vsize+tsize));       //Set inner bar voltage
    setvoltage(a0,Vbar,(vsize+tsize));      //Copy of initial matrix to use in update method 1

    //CREATE BOOL MATRIX
    truefalse(true, tsize);                 //Vacuum area = true
    truefalse(false, (tsize+vsize));        //Fixed potential = false
    printmatrixb(txtfile0);                 //Print matrix

    //UDATE MATRIX UNTIL CONVERGENCE REACHED
    double it = Vfield(2,txtfile1);         //Create potential field and return no. of iterations performed (note: update method 2 is more efficient)
    printmatrix(a,txtfile2);                //Print matrix
    extract1D(a,txtfile3);                  //Extract 1D cross-section
    cout << "Resolution = " << res << endl;
    cout << "No. of iterations to converge = " << it  << endl;

    //CREATE E FIELD MATRIX
    Efield();                   //Create E-field
    printmatrix(E,txtfile4);    //Print matrix

    return 0;
}

/* Functions
------------------------------------------------------------------*/

//SET INITIAL POTENTIAL MESH
void setvoltage(double element[][side], double potential, int outer)
{
    for (int j=(outer*res); j<(side-outer*res); j++)
    {
        for (int i=(outer*res); i<(side-outer*res); i++)
        {
            element[i][j] = potential;
        }
    }
}

//SET BOOL MESH
void truefalse(bool TorF, int outer)
{
    for (int j=(outer*res); j<(side-outer*res); j++)
    {
        for (int i=(outer*res); i<(side-outer*res); i++)
        {
            b[i][j] = TorF;
        }
    }
}

//UPDATE POTENTIAL INSIDE VACUUM
void update(int method)
{
    if (method == 1)
    {
       for (int i=0; i<side; i++)
        {
            for (int j=0; j<side; j++)
            {
                if (b[i][j] == true) //Only update potential in vacuum
                {
                    //Use elements from initial matrix a0 to find new matrix a
                    a[i][j]=0.25*(a0[i+1][j]+a0[i-1][j]+a0[i][j+1]+a0[i][j-1]);
                }
            }
        }
        for (int i=0; i<side; i++)  //Update a0
        {
            for (int j=0; j<side; j++)
            {
                a0[i][j] = a[i][j];
            }
        }
    }
    if (method == 2) //Updates element by element
    {
       for (int i=0; i<side; i++)
        {
            for (int j=0; j<side; j++)
            {
                if (b[i][j] == true)
                {
                    a[i][j]=0.25*(a[i+1][j]+a[i-1][j]+a[i][j+1]+a[i][j-1]);
                }
            }
        }
    }
}

//SUM ALL ELEMENTS IN POTENTIAL MATRIX
double sum()
{
    double S=0;
    for (int i=0; i<side; i++)
    {
        for (int j=0; j<side; j++)
        {
            if (b[i][j]==true)
            {
                S += abs(a[i][j]);
            }
        }
    }
    return S;
}

//BUILD POTENTIAL MATRIX
double Vfield(int method, ofstream& print)
{
    double Si = 0;                          //Initial sum of matrix
    double Sf = 0;                          //Final sum of matrix
    int kconv = 0;                          //Number of iterations needed
    double diff = 1;                        //Sf-Si (initially set to 1 to start while-loop)
    double tol = 0.00000000000001;          //Tolerance set to 1E-14
    while(abs(diff) > abs(tol*Si) || Si==0 || Sf==0)
    {
        Si = sum();
        update(method);
        Sf = sum();
        diff=abs(Sf-Si);
        kconv++;
        print << kconv << '\t' << Si << '\t' << Sf << '\t' << diff/Si << endl; //Output sum difference for each iteration
    }
    return kconv;
}

//BUILD ELECTRIC FIELD MATRIX
void Efield()
{
    for (int i=0; i<side; i++)
    {
        for (int j=0; j<side; j++)
        {
            if (b[i][j] == true)     //Only update potential in vacuum
            {
                Ex[i][j] = ((a[i-1][j]-a[i+1][j])*res)/2;               //i component
                Ey[i][j] = ((a[i][j-1]-a[i][j+1])*res)/2;               //j component
                E[i][j]= sqrt((Ex[i][j])*(Ex[i][j])+(Ey[i][j])*(Ey[i][j])); //Magnitude
            }
        }
    }
}

//PRINT MATRIX FUNCTION (POTENTIAL)
void printmatrix(double element[][side], ofstream& print)
{
    for (int i=0; i<side; i++)
    {
        for (int j=0; j<side; j++)
        {
                print << element[i][j] << '\t'; //Print element at i,j
        }
        print << endl;
    }
}

//PRINT MATRIX FUNCTION (BOOL)
void printmatrixb(ofstream& print)
{
    for (int i=0; i<side; i++)
    {
        for (int j=0; j<side; j++)
        {
                print << b[i][j] << '\t';
        }
        print << endl;
    }
}

//EXTRACT 1D CROSS-SECTION IN MIDDLE
void extract1D(double element[][side], ofstream& print)
{
    int j=ceil(res*(tsize+vsize+bsize/2)); //Take a line in midpoint
    for (int i=0; i<side; i++)
    {
            print << element[i][j] << endl;
    }
}
