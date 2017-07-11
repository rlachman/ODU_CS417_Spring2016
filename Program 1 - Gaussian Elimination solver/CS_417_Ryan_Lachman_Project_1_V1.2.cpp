#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
using namespace std;

typedef std::chrono::high_resolution_clock myclock;

double **A;
double *y;
int n;

void generate_row_dominant()
{
     myclock::time_point beginning = myclock::now();

//int n;
cout<<" Enter the random size of a matrix: "; cin>>n;

myclock::duration d = myclock::now() - beginning;
unsigned seed2 = d.count();
//cout<<seed2<<endl;
minstd_rand0 generator (seed2);
uniform_int_distribution<int> distribution(-9999999,9999999);


A = new double *[n];
for(int i=0; i<n; i++)
{
    A[i] = new double [n];
}

cout << endl << "The generated random matrix A after enforcing diagonal dominance is: " << endl;
cout << "******************************************************************** " << endl;
for(int r=0; r<n; r++)
{
    for(int c=0; c<n; c++)
    {
        if(r == c)
        A[r][c]=double(distribution(generator))/100.0;
        else
        A[r][c] = double(distribution(generator))/10000.0;
        cout<<A[r][c]<<"  ";
    }
    cout << endl;
}

cout << endl;

}

void generate_random_solution()
{
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 = d.count();
   // cout<<seed2<<endl;
    minstd_rand0 generator (seed2);
    uniform_int_distribution<int> distribution(-9999999,9999999);

    y = new double [n];

    cout << "The generated random solution is: " << endl;
    cout << "********************************* " << endl;
    for(int i=0;i<n;i++)
    {
        y[i] = double(distribution(generator))/1000000.0;
        cout << y[i] << " ";
    }

    cout << endl;
    cout << endl;
}

void multiply(double* mul)
{

    for(int i=0;i<n;i++)
    {
        mul[i] = 0;
    }
    for(int r=0;r<n;r++)
    {
        for(int c = 0;c < n; c++)
        {
            mul[r]+= (double)A[r][c]*y[c];
        }
    }

    cout << "b matrix generated after multiplication of A*Y: " << endl;
    cout << "*********************************************** " << endl;
    for(int r=0;r<n;r++)
    {
            cout << mul[r]<<" ";
    }
    cout << endl;
    cout << endl;

}

void upper_row_echelon(double* b)
{
    int curr_row = 0;
    int max1;
    int index;
    double temp,temp2;

    for(int i = 0; i<n; i++)
    {
        for(int c = i; c<n; c++)
        {
            if(c == i)
            {
                max1 = A[i][c];
                index = i;
            }
            else
            {
                if(A[i][c] > max1)
                {
                    max1 = A[i][c];
                    index = c;
                }
            }
        }

        if(index != i)
        {

            for(int c = i; c< n; c++)      // Exchange Rows
            {
                // Exchange Rows of A matrix.

                temp = A[curr_row][c];
                A[curr_row][c] = A[index][c];
                A[index][c] = temp;

            }
            // Exchange Rows of b matrix

                temp = b[curr_row];
                b[curr_row] = b[index];
                b[index] = temp;


        }

            //reduce the row to 1.0

            temp = A[curr_row][curr_row];

            for(int c= curr_row; c < n; c++)
            {
                A[curr_row][c] = A[curr_row][c]/temp;
            }
            b[curr_row] = b[curr_row]/temp;

            // Reduce the below entries of curr_row to zero.

            for(int i=curr_row+1;i<n;i++)
            {
                temp = A[i][curr_row];

                for(int c = curr_row; c < n ; c++)
                {
                    temp2 = A[curr_row][c];
                    A[i][c] = A[i][c] - (temp*temp2);
                }

                b[i] = b[i] - (b[curr_row]*temp);
            }
        curr_row++;
    }


    // Print A matrix
    cout << "Matrix generated in upper row echelon: " <<endl;
    cout << "************************************* " << endl;
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            cout << A[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl;
    // backtrack and solve
}

void back_solve(double* sol,double* b)
{
        sol[n-1] = b[n-1];

//        sol[n-2] = (b[n-2] - a[n-2][n-1]*sol[n-1])/a[n-2][n-2];

//        sol[n-3] = (b[n-3] - a[n-3][n-1]*sol[n-1] - a[n-3][n-2]*sol[n-2])/a[n-3][n-3];


        double temp;
        for(int i = 1;i<n;i++)
        {
            temp = b[n-i-1];

            for(int j=0; j<i;j++)
            {
                temp = temp - (A[n-i-1][n-j-1]*sol[n-j-1]);
            }

            temp = temp/A[n-i-1][n-i-1];
            sol[n-i-1] = temp;
        }

        cout <<"Solution Matrix obtained after backtracking: " << endl;
        cout <<"******************************************** " << endl;
        for(int i=0;i<n;i++)
        {
            cout << sol[i] << " ";
        }
        cout << endl;
        cout << endl;

}

void compute_b1(double* sol,double* b1)
{

    for(int i=0;i<n;i++)
    {
        b1[i] = 0;
        for(int j = 0; j<n; j++)
        {
            b1[i] = b1[i] + sol[j]*A[i][j];
        }
    }
}

void compute_e(double* e,double* b,double* b1)
{
    double norm = 0;

    cout << "The error matrix generated is: " << endl;
    cout << "****************************** " << endl;
    for(int i=0;i<n;i++)
    {
        e[i] = b[i] - b1[i];
        if(e[i] < 0)
        {
            e[i] = (-1)*e[i];
        }
        norm = norm + (e[i]*e[i]);
        cout << e[i] <<" ";
    }
	norm = (double)sqrt(norm);
    cout << endl;
    cout << endl;
    cout <<"Euclidian Norm of the error matrix is: " << endl;
    cout <<"************************************** " << endl;
    cout << norm <<endl;
}
int main()
{

generate_row_dominant();            // Generate Diagonally Row Dominant Matrix [ A matrix]
generate_random_solution();         // Generate Random Solution [y matrix ]
double* b;
b = new double [n];

multiply(b);                        // Multiply matrix A with matrix y to generate matrix b
upper_row_echelon(b);               // Convert the matrix to upper row echelon form
double* sol;
sol = new double[n];
back_solve(sol,b);                  // Back solve for x.

double* b1;
b1 = new double [n];
compute_b1(sol,b1);                 // Compute A*x

double* e;
e = new double [n];
compute_e(e,b,b1);                  // Compute Error matrix i.e. e = ~b-b

return 0;

}
