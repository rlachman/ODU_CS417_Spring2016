#include <iostream>
#include <random>
#include <ctime>
#include <chrono>
#include <limits>
using namespace std;

typedef std::chrono::high_resolution_clock myclock;

double** A;
double** A_copy;
double** A_copy2;

void generate_row_dominant(int n)
{
     myclock::time_point beginning = myclock::now();

//int n;


myclock::duration d = myclock::now() - beginning;
unsigned seed2 = d.count();
//cout<<seed2<<endl;
minstd_rand0 generator (seed2);
uniform_int_distribution<int> distribution(-9999999,9999999);




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

for(int r = 0;r < n;r++)
{
    for(int c = 0;c < n; c++)
    {
        A_copy[r][c] = A[r][c];
        A_copy2[r][c] = A[r][c];
    }
}

cout << endl;

}

void generate_random_solution(double* y, int n)
{
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;
    unsigned seed2 = d.count();
   // cout<<seed2<<endl;
    minstd_rand0 generator (seed2);
    uniform_int_distribution<int> distribution(-9999999,9999999);

    cout << "The generated random solution Y is: " << endl;
    cout << "*********************************** " << endl;
    for(int i=0;i<n;i++)
    {
        y[i] = double(distribution(generator))/1000000.0;
        cout << y[i] << " ";
    }

    cout << endl;
    cout << endl;
}

void multiply(double* mul,int n,double* y)
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

void upper_row_echelon(double* b,int n)
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

void back_solve(double* sol,double* b,int n)
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

void compute_b1(double* sol,double* b1,int n)
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

void compute_e(double* e,double* b,double* b1,int n)
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

// New Assignment ( LU Decomposition )

void LUDecomposer(double* b,int n)
{
    int curr_row = 0,counter = 0,i,j;
    double temp,temp2,temp3;
    double** U = new double*[n];
    double** L = new double*[n];
    U = new double *[n];
    L = new double *[n];
    for(int i=0; i<n; i++)
    {
        U[i] = new double [n];
        L[i] = new double [n];
    }

    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            U[i][j] = 0;
            L[i][j] = 0;
        }
    }

    // Set the diagonal elements of U to be 1.0
    for(int i=0;i<n;i++)
        U[i][i] = 1.0;

    // Copy the first column of A into the first column of L.
    for(int j=0;j<n;j++)
        L[j][0] = A_copy[j][0];

    // Divide the first row of A by A[0][0] and copy it to the first row of U.
    for(int i=0;i<n;i++)
    {
        A_copy2[0][i] = A_copy2[0][i] / A_copy2[0][0];
        U[0][i] = A_copy2[0][i];
    }

    // Compute The row echlon form of the matrix which will give U matrix
     for(i=0;i<n-1;i++)
     {
            counter++;
            temp3 = A_copy[curr_row][curr_row];
            for(int c= curr_row; c < n; c++)
            {
                A_copy[curr_row][c] = A_copy[curr_row][c]/temp3;
            }

            for(j = curr_row+1;j<n;j++)
            {

                temp = A_copy[j][curr_row];
                for(int c = curr_row; c < n ; c++)
                {
                    temp2 = A_copy[curr_row][c];
                    A_copy[j][c] = A_copy[j][c] - (temp*temp2);


                }
            }

        // Compute L matrix using multipliers at each stage.
        curr_row++;
        int c = curr_row;
        for(j = curr_row; j<n;j++)
        {
            L[j][c] = A_copy[j][c];
        }
     }


    // The matrix A in row echlon gives us the U matrix.
     for(int i=0;i<n;i++)
     {

         for(int j=0;j<n;j++)
         {
             if(i == j)
                U[i][i] = 1;
            else
                U[i][j] = A_copy[i][j];
         }

     }

    // Print L matrix
    cout << "L matrix is " << endl;
    cout << "*****************************" << endl;
    cout <<endl;
     for(int i=0;i<n;i++)
     {
         for(int j=0;j<n;j++)
         {
             cout << L[i][j] << " ";
         }
         cout << endl;
     }
    cout << endl;


    // Print U matrix
    cout << "U matrix is " << endl;
    cout << "*****************************" << endl;
    cout <<endl;
     for(int i=0;i<n;i++)
     {

         for(int j=0;j<n;j++)
         {
             cout << U[i][j] << " ";
         }

         cout << endl;
     }
     cout << endl;


    // Compute y matrix. Y matrix is denoted here by sol[] array.

    // Using forward solve to compute y matrix.
    double sol[n];
    sol[0] = b[0]/L[0][0];

        for(int i = 1;i<n;i++)
        {
            temp = b[i];

            for(int j=0; j<i;j++)
            {
                temp = temp - (L[i][j]*sol[j]);
            }

            temp = temp/L[i][i];
            sol[i] = temp;
        }



        cout <<"Solution Matrix obtained after Forward solve Ly=b: " << endl;
        cout <<"******************************************** " << endl;
        for(int i=0;i<n;i++)
        {
            cout << sol[i] << " ";
        }
        cout << endl;
        cout << endl;

        // Compute x matrix.

        // Using Back solve to find x matrix.
        cout << "Vector x generated by back solve Ux = y" << endl;
        double x[n];
        x[n-1] = sol[n-1];

        for(int i = 1;i<n;i++)
        {
            temp = sol[n-i-1];

            for(int j=0; j<i;j++)
            {
                temp = temp - (U[n-i-1][n-j-1]*sol[n-j-1]);
            }

            temp = temp/U[n-i-1][n-i-1];
            x[n-i-1] = temp;
        }

        cout <<"******************************************** " << endl;
        for(int i=0;i<n;i++)
        {
            cout << x[i] << " ";
        }

        // Computing A*x to generate b1 matrix
        double b1[n];
        for(int i=0;i<n;i++)
        {
            b1[i] = 0;
            for(int j = 0; j<n; j++)
            {
                b1[i] = b1[i] + x[j]*A[i][j];
            }
        }


        cout << endl;
        cout << endl;

        double ans = 0;
        double e[n];

        // Computing error matrix and norm of the error.
        for(int i=0;i<n;i++)
        {
            if(b[i] > b1[i])
            {
                e[i] = b[i] - b1[i];
                ans = ans + e[i]*e[i];
            }
            else
            {
                e[i] = b1[i] - b[i];
                ans = ans + e[i]*e[i];
            }
        }
        ans = sqrt(ans);
        cout << "Norm  ||b-Ax|| is " << ans << endl;
        cout << endl;

}

double error(int n,double* sol,double* y,double* b)
{
        double b1[n];
        double e[n];
        double ans;
        for(int i=0;i<n;i++)
        {
            b1[i] = 0;
            for(int j = 0; j<n; j++)
            {
                b1[i] = b1[i] + sol[j]*A[i][j];
            }
        }

        for(int i=0;i<n;i++)
        {
            if(b[i] > b1[i])
            {
                e[i] = b[i] - b1[i];
                ans = ans + e[i]*e[i];
            }
            else
            {
                e[i] = b1[i] - b[i];
                ans = ans + e[i]*e[i];
            }
        }
        ans = sqrt(ans);

        return ans;


}

void Jacobi_Iteration(int n,double* b,double* y)
{
    cout << "**************Solution by Jacobi_Iteration**************" << endl;
    cout << endl;
    double sol[n];
    double c[n];
    int i = 0;
    // Set initial solution to be zero.

    for(int i=0;i<n;i++)
        sol[i] = 0;

    // Compute solutiobs for each iteration and print them.

    int count = 1;
    double error_check =  INT_MAX;
    while(error_check > 0.0000001)
    {

        for(int i=0;i<n;i++)
        {
            c[i] = b[i];

            for(int j=0;j<n;j++)
            {
                if(i!=j)
                {
                    c[i] = c[i] - A[i][j]*sol[j];
                }
            }
        }

        for(int i=0;i<n;i++)
        {
            sol[i] = c[i] / A[i][i];
        }

        cout << "Solution after iteration " << count << ":" << endl;
        for(int i=0;i<n;i++)
        {
            cout << sol[i] << " ";
        }
        cout << endl;
        error_check = error(n,sol,y,b);
        cout << endl << "Error after iteration " << count << ":" << endl << error_check << endl;
        cout << endl;
        count++;
    }

    cout << endl;
    cout << "Iteration complete error is < or = 10^(-7)";
    cout << endl;
    cout << endl;
    cout << endl;
}

void Gauss_Siedel(int n,double* b,double* y)
{
    cout << "**************Solution by Gauss_Siedel**************" << endl;
    cout << endl;
    double sol[n];
    int i = 0;
    double c[n];
    double error_check = INT_MAX;
    // Set initial solution to be zero.

    for(int i=0;i<n;i++)
        sol[i] = 0;

    // Compute solutiobs for each iteration and print them.

    int count = 1;

    while(error_check > 0.0000001)
    {
        for(int i=0;i<n;i++)
        {
            c[i] = b[i]/A[i][i];
            for(int j=0;j<n;j++)
            {
                    if(i != j)
                    {
                        c[i] = c[i] - ((A[i][j]/A[i][i])*sol[j]);
                    }
            }
            sol[i] = c[i];
        }

        cout << "Solution after iteration " << count << ":" << endl;
        for(int i=0;i<n;i++)
        {
            cout << sol[i] << " ";
        }
        cout << endl;
        error_check = error(n,sol,y,b);
        cout << endl << "Error after iteration " << count << ":" << endl << error_check << endl;
        cout << endl;
        count++;
    }

    cout << endl;
    cout << "Iteration complete error is < or = 10^(-7)";
    cout << endl;
    cout << endl;
    cout << endl;

}

void SOR(int n,double* b,double* y)
{
    cout << "**************Solution by SOR**************" << endl;
    double lambda;
    int num_itr;
    cout << endl << "Enter the value of omega: ";
    cin >> lambda;
    cout << endl << endl;

    int count = 1;
    double sol[n];
    double c[n];
    int i = 0;
    double error_check = INT_MAX;

    for(int i=0;i<n;i++)
        sol[i] = 0;

    while(error_check > 0.0000001)
    {
        for(int i=0;i<n;i++)
        {

            c[i] = lambda * (b[i]/A[i][i]);

            for(int j=0;j<n;j++)
            {
                if(i!=j)
                {
                    c[i] = c[i] - (lambda*((A[i][j]/A[i][i])*sol[j]));
                }
            }
            sol[i] = c[i] + (1-lambda)*sol[i];
        }
        cout << "Solution after iteration " << count << ":" << endl;
        for(int i=0;i<n;i++)
        {
            cout << sol[i] << " ";
        }
        cout << endl << endl;
        error_check = error(n,sol,y,b);
        cout <<"Error after iteration " << count << ":" << endl << error_check << endl;
        cout << endl;
        count++;
    }

    cout << endl;
    cout << "Iteration complete error is < or = 10^(-7)";
    cout << endl;
    cout << endl;
    cout << endl;
}

int main()
{
double *y;
double* b;
int n;

cout<<" Enter the random size of a matrix: ";
cin>>n;

A = new double *[n];
A_copy = new double *[n];
A_copy2 = new double * [n];
y = new double [n];
b = new double [n];

for(int i=0; i<n; i++)
{
    A[i] = new double [n];
    A_copy[i] = new double [n];
    A_copy2[i] = new double [n];
}


generate_row_dominant(n);            // Generate Diagonally Row Dominant Matrix [ A matrix]
generate_random_solution(y,n);         // Generate Random Solution [y matrix ]

multiply(b,n,y);

/*                        // Multiply matrix A with matrix y to generate matrix b
upper_row_echelon(b,n);               // Convert the matrix to upper row echelon form
double* sol;
sol = new double[n];
back_solve(sol,b,n);                  // Back solve for x.

double* b1;
b1 = new double [n];
compute_b1(sol,b1,n);                 // Compute A*x

double* e;
e = new double [n];
compute_e(e,b,b1,n);                  // Compute Error matrix i.e. e = ~b-b

*/


//LUDecomposer(b,n);                // Compute LU Decomposer
Jacobi_Iteration(n,b,y);              // Jacobi Iteration
Gauss_Siedel(n,b,y);                  // Gauss Siedel Iteration
SOR(n,b,y);
return 0;

}
