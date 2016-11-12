#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**< Define some macros */
#define K 10


/**< Define some static constants */
const double epsilon = 1.e-12;
const int M=2000;

/**< Statements of miscellaneous functions */
int vector_add(double * vector1, double * vector2, double * sum, int n, int flag);
double norm_row(double * vector, int n);
int matrix_multiplication(double **A, double **B, double **C, int m, int n, int s);
int matrix_transpose(double ** A, double ** C, int m, int n);
int solve_linear_equation(double **A, double *b, double *solution, int size);
double gaussian_algorithm(double a[4][4],double b[4],double x[4]);
int get_Fun(double *F, double * X, double x, double y);
int get_Jacobi_Fun(double **JocabiF, double * X, int size);
int product_vector_and_scalar(double *X, double * output, double c, int size);
int newton_algorithm(double * F, double * Fb, double ** JacobiF, double *X, double *dX, double xi, double yj, int size);
double newton(double x, double y,double tuvw[4]);
int get_t_u(double *tij, double *uij, double xi, double yj);
int get_T_U(double ** T, double **U, int m, int n);
double L_u(double y,int j,int r);
double L_t(double x,int i,int k);
double order2_interp(double tt,double uu);
double FindMax(double Arr[4]);


////////////////////////////////////////////////
int get_Z(double **Z);
int find_opt_k(void);
double surf(double x, double y, double **Z, int k, double **optC);
int get_C(double **C, double **Z, int k);
int get_GtGG(double **GtGG, double ** Z, int k);
int get_G(double G[][10], int k);
int get_A(double **A, double ** Z, int k);
int product_matrix_and_vector(double** A, double*x, double *output, int m, int n);
int get_Uj(double ** Z, double *U, int j, int m, int n);
int get_B(double B[][10], int k);
void matrix_inverse(double A[][10],int n);
double S(double x,double y,double C[][10],int k,double Z[][21]);
int surf_with_optk(int optk);

int test1(void);
int test2(void);
int test3(void);
int find_opt_k2(void);

int main()
{
    printf("By Janshan Zhou (BY1613123)\n\n");
    //test1();
//    test2();
//    test3();
    find_opt_k2();
    return 0;
}



/**< Miscellaneous functions involved in the numerical computing */
int test1(void)
{
    /*Test some things.*/
    double v[3] = {1.,2.,3.}, u[3] = {4.,5.,6.};
    double w[3]={0.0};
    double a[][3] = {{1.,2.,3.},{4.,5.,6.}};
    double b[][2] = {{1.,2.},{3.,4.},{5.,6.}};
    double c[2][2] = {0.0};
    double d[3][2] = {0.};

    double norm_s = norm_row(v, 3);
    printf("The row norm of v is %.3f\n",norm_s);

    vector_add(v,u,w,3,-1);
    int i;
    for(i=0;i<3;i++)
    {
        printf("w[%d] = %.3f+%.3f= %.3f\n",i+1,*(v+i),*(u+i),*(w+i));
    }
    matrix_multiplication(a,b,c,2,3,2);
    int j;
    for(i=0;i<2;i++)
    {
        for(j=0;j<2;j++)
        {
            printf("c[%d][%d]=%.3f\t",i+1,j+1,*((double*)c+2*i+j));
        }
        printf("\n");
    }

    matrix_transpose(a,d,2,3);
    for(i=0;i<3;i++)
    {
        for(j=0;j<2;j++)
        {
            printf("d[%d][%d]=%.3f\t",i+1,j+1,*((double*)d+2*i+j));
        }
        printf("\n");
    }
    printf("\n\n");
    double A[][3]={{2.,3.,1.},{4.,2.,3.},{7.,1.,-1.0}};
    double bb[3] = {4.,17.,1.};
    double solution[3]={0.0};
    solve_linear_equation(A,bb,solution,3);
    for(i=0;i<3;i++)
    {
        printf("x[%d]=%.3f\n",i+1,*(solution+i));
    }

    return 0;
}

int test2(void)
{
    double T[10+1][20+1]={{0.0}};
    double U[10+1][20+1]={{0.0}};
    double Z[10+1][20+1]={{0.0}};
    int m=10+1, n=20+1;
    get_T_U(T, U, m, n);
    int i,j;
//    for(i=0;i<m;i++)
//    {
//        for(j=0;j<n;j++)
//        {
//            printf("T[%d][%d]=%.3f\t",i,j,*((double*)T+n*i+j));
//        }
//        printf("\n");
//    }
//
//    printf("\n");
//    for(i=0;i<m;i++)
//    {
//        for(j=0;j<n;j++)
//        {
//            printf("U[%d][%d]=%.3f\t",i,j,*((double*)U+n*i+j));
//        }
//        printf("\n");
//    }

    double tt, uu;
    //solve for f(xi,yj)
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            tt = *((double*)T+n*i+j);
            uu = *((double*)U+n*i+j);
            *((double*)Z+n*i+j)=order2_interp(tt,uu);
        }
    }
    //print (xi,yj, f(xi,yj))
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            printf("(x%d, y%d, f(x%d, y%d)) = (%.3f, %.3f, %.11e)\n",i,j,i,j,0.08*i,0.5+0.05*j,*((double*)Z+n*i+j));
        }
        printf("\n");
    }
    return 0;
}

int test3(void)
{
    find_opt_k();
    return 0;
}

int get_T_U(double ** T, double **U, int m, int n)
{
    /*Solve for all the (t,u) given all (xi,yj).
    T is in m*n and U is in m*n.*/
    int i,j;
    double xi,yj;
    double tij, uij;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            xi = 0.08*i;
            yj = 0.5+0.05*j;
            get_t_u(&tij, &uij, xi, yj);
            *((double*)T+n*i+j) = tij;
            *((double*)U+n*i+j) = uij;
            //printf("t[%d][%d]=%.3f\t",i,j,tij);
        }
        //printf("\n");
    }
    return 0;
}

int get_t_u(double *tij, double *uij, double xi, double yj)
{
    /*This function calculates (t,u) given a pair (xi,yj).*/
    double F[4]={0.0};
    double Fb[4]={0.0};
    double JacobiF[4][4]={0.0};
    double X[4]={1.0,1.0,1.0,1.0};
    double dX[4]={1.0,1.0,1.0,1.0};
    int size=4;
    //solve for X=[t,u,v,w]^T
    newton_algorithm(F,Fb,JacobiF,X,dX,xi,yj,size);
    *(tij) = *(X+0);
    *(uij) = *(X+1);
    return 0;
}


int newton_algorithm(double * F, double * Fb, double ** JacobiF, double *X, double *dX, double xi, double yj, int size)
{
    /*The Newton algorithm is used to solve the nonlinear equation F(X)=0*/
    int i;
    double norm_dX, norm_X;

    for(i=0;i<=M;i++)
    {
        get_Fun(F, X, xi, yj);//get F
        get_Jacobi_Fun(JacobiF, X, size);// get JacobiF
        product_vector_and_scalar(F,Fb,-1.0,size);//Fb=-F
        solve_linear_equation(JacobiF,Fb,dX,size);//solve JacobiF*dx=Fb for dX
        norm_dX = norm_row(dX,size);
        norm_X = norm_row(X,size);
        if((norm_dX/norm_X)<=epsilon)
        {
            break;
        }
        vector_add(X,dX,X,size,1);//X=X+dX
    }
    if(i==M)
    {
        printf("The Newton algorithm cannot converge!\n\n");
        exit(-1);
    }
    return 0;
}

int get_Fun(double *F, double * X, double x, double y)
{
    /*This function calculates the vectorized function F(X),
    where X=[t,u,v,w]^T given (xi,yj)*/
    double t=*(X+0);
    double u=*(X+1);
    double v=*(X+2);
    double w=*(X+3);

    *(F+0) = 0.5*cos(t)+u+v+w-x-2.67;
    *(F+1) = t+0.5*sin(u)+v+w-y-1.07;
    *(F+2) = 0.5*t+u+cos(v)+w-x-3.74;
    *(F+3) = t+0.5*u+v+sin(w)-y-0.79;
    return 0;
}


int get_Jacobi_Fun(double **JocabiF, double * X, int size)
{
    /*This function obtains the Jacobi matrix of the vectorized function F(X)
     w.r.t. X*/
    int i,j;
    double t=*(X+0);
    double u=*(X+1);
    double v=*(X+2);
    double w=*(X+3);

    for(i=0;i<size;i++)
    {
        for(j=0;j<size;j++)
        {
            *((double*)JocabiF+size*i+j) = 1.0;
        }
    }
    *((double*)JocabiF+size*0+0) = -0.5*sin(t);
    *((double*)JocabiF+size*1+1) = 0.5*cos(u);
    *((double*)JocabiF+size*2+0) = 0.5;
    *((double*)JocabiF+size*2+2) = -sin(v);
    *((double*)JocabiF+size*3+1) = 0.5;
    *((double*)JocabiF+size*3+3) = cos(w);
    return 0;
}

int product_vector_and_scalar(double *X, double * output, double c, int size)
{
    int i;
    for(i=0;i<size;i++)
    {
        *(output+i)=(*(X+i))*c;
    }
    return 0;
}

int vector_add(double * vector1, double * vector2, double * sum, int n, int flag)
{
    /*This function calculates the sum of two given column vectors in size n*/
    int i;
    for(i=0;i<n;i++)
    {
        *(sum+i) = *(vector1+i)+(flag)*(*(vector2+i));
    }
    return 0;
}


double norm_row(double * vector, int n)
{
    /*This function calculates the row-oriented norm of a given vector*/
    double s = fabs(*(vector + 0)), temp_s;
    int i;
    for(i=0;i<n;i++)
    {
        temp_s=fabs(*(vector+i));
        if(temp_s>s)
        {
            s = temp_s;
        }
    }
    return s;
}


int matrix_multiplication(double **A, double **B, double **C, int m, int n, int s)
{
    /*This function calculates the product of a matrix A in m*n and a matrix B in n*s.
    The product is denoted by a matrix C in m*s*/
    double st = 0.0;
    int i,j,k;
    for(i=0;i<m;i++)
    {
        for(j=0;j<s;j++)
        {
            st = 0.0;
            for(k=0;k<n;k++)
            {
                st += (*((double*)A+n*i+k))*(*((double*)B+s*k+j));
            }
            *((double*)C+s*i+j) = st;
        }
    }
    return 0;
}


int matrix_transpose(double ** A, double ** C, int m, int n)
{
    /*This function transposes the given matrix A in m*n.
    The output matrix C is is n*m.*/
    int i, j;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            *((double*)C+m*j+i) = *((double*)A+n*i+j);
        }
    }
    return 0;
}


int solve_linear_equation(double **A, double *b, double *solution, int size)
{
    /*This function solves a group of linear equations by using Gaussian elimination algorithm.
    The coefficient matrix A is in m*m, i.e., a square matrix;
    the constant column vector b is in m*1, and
    the output solution is in m*1.*/
    double Aik, S, temp;
    int i, j, k;
    double max;//the maximum value in a given column column
    int col;//the row index associated with the maximum column element

    //Do elimination
    for(k=0;k<size;k++)
    {
        max = fabs(*(double*)A+size*k+k);
        col = k;
        //see for the row index associated with the maximum column element
        for(i=k;i<size;i++)
        {
            if(max<fabs(*((double*)A+size*i+k)))
            {
                max = fabs(*((double*)A+size*i+k));
                col = i;
            }
        }
        //exchange the row of the col-th and the k-th rows
        for(j=k;j<size;j++)
        {
            temp = *((double*)A+size*col+j);
            *((double*)A+size*col+j) = *((double*)A+size*k+j);
            *((double*)A+size*k+j) = temp;
        }
        //exchange the row of the col-th and the k-th b- entries
        temp = *(b+col);
        *(b+col) = *(b+k);
        *(b+k) = temp;
        if(!(*((double*)A+size*k+k)))
        {
            printf("Cannot solve the linear equations!");
            exit(-1);
        }

        for(i=k+1;i<size;i++)
        {
            Aik=1.0*(*((double*)A+size*i+k))/(*((double*)A+size*k+k));
            for(j=k;j<size;j++)
            {
                *((double*)A+size*i+j)=*((double*)A+size*i+j)-Aik*(*((double*)A+size*k+j));
            }
            *(b+i) = *(b+i) - Aik*(*(b+k));
        }
    }

    //Do re-iteration
    *(solution+size-1) = (*(b+size-1))/(*((double*)A+size*(size-1)+size-1));
    for(k=size-2;k>=0;k--)
    {
        S=*(b+k);
        for(j=k+1;j<size;j++)
        {
            S=S-(*((double*)A+size*k+j))*(*(solution+j));
        }
        *(solution+k)=S/(*((double*)A+size*k+k));
    }

    return 0;

}


double L_t(double x,int i,int k)
{
    /*Lagrange-type l_k(xi)*/
	int m;
	double value=1,t[6]={0,0.2,0.4,0.6,0.8,1.0 };
	for (m=i-1;m<=i+1;m++)
	if(m!=k)value=value*(x-t[m])/(t[k]-t[m]);
	return value;
}

double L_u(double y,int j,int r)
{
     /*Lagrange-type Lr(yj)*/
	int m;
	double value=1;
	double u[6]={0,0.4,0.8,1.2,1.6,2};
	for (m=j-1;m<=j+1;m++)
	if(m!=r)value=value*(y-u[m])/(u[r]-u[m]);
	return value;
 }

double order2_interp(double tt,double uu)
 {
     /*Do the order-2 interpolation over t and u*/
	int i,j,k,r;
	double zvalue;
	double z[6][6] = {{-0.5,-0.34,0.14,0.94,2.06,3.5},
	{-0.42,-0.5,-0.26,0.3,1.18,2.38},
	{-0.18,-0.5,-0.5,-0.18,0.46,1.42},
	{0.22,-0.34,-0.58,-0.5,-0.1,0.62},
	{0.78,-0.02,-0.5,-0.66,-0.5,-0.02},
	{1.5,0.46,-0.26,-0.66,-0.74,-0.5}};
	//determine the row index of the central node
	switch((int)(10*tt))
	{
		case 0:
		case 1:
		case 2:i=1;break;
		case 3:
		case 4:i=2;break;
		case 5:
		case 6:i=3;break;
		case 7:
		case 8:
		case 9:i=4;break;
		default:
		    {printf("Some errors in t!\n");
		    exit(-1);}
	}
	//Determine the column index of the central node
	switch((int)(5*uu))
	{
		case 0:
		case 1:
		case 2:j=1;break;
		case 3:
		case 4:j=2;break;
		case 5:
		case 6:j=3;break;
		case 7:
		case 8:
		case 9:j=4;break;
		default:
		    {printf("Some errors in t!\n");
		    exit(-1);}
	}
	zvalue=0.0;
	for(k=i-1;k<=i+1;k++)
    {
        for (r=j-1;r<=j+1;r++)
        {
            zvalue += z[k][r]*L_t(tt,i,k)*L_u(uu,j,r);
        }
    }
	return zvalue;
 }




 //////////////////////////
 int get_Uj(double ** Z, double *U, int j, int m, int n)
 {
     int i;
     for(i=0;i<(m+1);i++)
     {
         *(U+i) = *((double*)Z+(n+1)*i+j);
     }
     return 0;
 }

int get_B(double B[][10], int k)
{
	double x[11],sum;
	int i,j,m;
	for (i=0;i<11;i++)
	x[i]=0.08*i;
	for (i=0;i<11;i++)
	{
		B[i][0]=1;
		for (j=1;j<=k;j++)
		{
			sum=1;
			for (m=1;m<=j;m++)
				sum=sum*x[i];
			B[i][j]=sum;
		}
	}
	return 0;
 }


 int product_matrix_and_vector(double** A, double*x, double *output, int m, int n)
 {
     /*A is in m*n, and x is in n*1, output is in m*1*/
     int i, j, k;
     double s;
     for(i=0;i<m;i++)
     {
         s = 0.0;
         for(j=0;j<n;j++)
         {
             s += (*((double*)A+n*i+j))*(*(x+j));
         }
         *(output+i)=s;
     }
     return 0;
 }

int get_A(double **A, double ** Z, int k)
{
    /*A is in (K+1)*(n+1)*/
    int m=10, n=20;

    double B[10+1][K+1]={0.0};
    double Bt[K+1][10+1] = {0.0};
    double uj[10+1]={0.0};
    double Btuj[K+1]={0.0};
    double BtB[K+1][K+1]={0.0};
    double alpha[K+1]={0.0};

    get_B(B,k);
    matrix_transpose(B,Bt,m+1,K+1);
    matrix_multiplication(Bt,B,BtB,K+1,m+1,K+1);
    int r,j;
    for(j=0;j<(n+1);j++)
    {
        get_Uj(Z,uj,j,m,n);
        product_matrix_and_vector(Bt,uj,Btuj,K+1,10+1);
        solve_linear_equation(BtB,Btuj,alpha,k);
        for(r=0;r<(k+1);r++)
        {
            *((double*)A+(n+1)*r+j) = *(alpha+r);
        }
    }
    return 0;
}

int get_G(double G[][10], int k)
{
	int i,j,m;
	double y[21],sum;
	for (i=0;i<21;i++)
	y[i]=0.5+0.05*i;
	for (i=0;i<21;i++)
	{
		G[i][0]=1;
		for (j=1;j<=k;j++)
		{
			sum=1;
			for (m=1;m<=j;m++)
				sum=sum*y[i];
			G[i][j]=sum;
		}
	}
	return 0;
}


int get_GtGG(double **GtGG, double ** Z, int k)
{
    /*GtGG is in (K+1)*(n+1)*/
    int m=10, n=20;
    double G[20+1][K+1]={0.0};
    double Gt[K+1][20+1]={0.0};
    double GtG[K+1][K+1]={0.0};
    double gamma[K+1]={0.0};
    double b[K+1]={0.0};
    int j,s;
    get_G(G,k);
    matrix_transpose(G,Gt,n+1,K+1);
    matrix_multiplication(Gt,G,GtG,K+1,n+1,K+1);
    for(j=0;j<(n+1);j++)
    {
        for(s=0;s<(k+1);s++)
        {
            *(b+s)=*((double*)Gt+(n+1)*s+j);
        }
        solve_linear_equation(GtG,b,gamma,k);
        for(s=0;s<(k+1);s++)
        {
            *((double*)GtGG+(n+1)*s+j)=*(gamma+s);
        }
    }

    return 0;
}

int get_C(double **C, double **Z, int k)
{
    int n=20;
    double A[K+1][20+1]={0.0};
    double GtGG[K+1][20+1]={0.0};
    double GGtG[20+1][K+1]={0.0};
    get_A(A,Z,k);
    get_GtGG(GtGG,Z,k);
    matrix_transpose(GtGG,GGtG,K+1,n+1);
    matrix_multiplication(A,GGtG,C,K+1,n+1,K+1);
}

double surf(double x, double y, double **Z, int k, double **optC)
{
    double zvalue=0.0;
    //double C[K+1][K+1]={0.0};
    int m=10, n=20;
    get_Z(Z);
    zvalue = S(x,y,optC,k,Z);
    int i,j;
//    for(i=0;i<(k+1);i++)
//    {
//        for(j=0;j<(k+1);j++)
//        {
//            *((double*)optC+(K+1)*i+j) = *((double*)C+(K+1)*i+j);
//        }
//    }
    return zvalue;
}

int surf_with_optk(int optk)
{
    double xi,yj;
    double tij,uij;
    double X[4]={0.0};
    int m=8,n=5;
    double Z[10+1][20+1]={{0.0}};
    double fvalue[8][5]={{0.0}};
    double optC[K+1][K+1]={0.0};
    int i,j;

    for(i=1;i<=m;i++)
    {
        for(j=1;j<=n;j++)
        {
            xi=0.1*i;
            yj=0.5+0.2*j;
            //newton(xi, yj, X);
            get_t_u(&tij,&uij,xi,yj);
            //*((double*)Z+n*(i-1)+j-1)=order2_interp(X[0], X[1]);
            *((double*)Z+(20+1)*(i-1)+j-1)=order2_interp(tij, uij);
            //printf("%.11e\n",*((double*)Z+(20+1)*(i-1)+j-1));
        }
    }

    printf("**********************************************\n");
    for(i=1;i<=m;i++)
    {
        for(j=1;j<=n;j++)
        {
            fvalue[i-1][j-1] = *((double*)Z+(20+1)*(i-1)+j-1);
        }
    }
    double pvalue;

    for(i=1;i<=m;i++)
    {
        for(j=1;j<=n;j++)
        {
            xi=0.1*i;
            yj=0.5+0.2*j;
            pvalue = surf(xi,yj,Z,optk,optC);
            printf("(x*%d,y*%d,f(x*%d,y*%d),p(x*%d,y*%d))=(%.3f,%.3f,%.11e,%.11e)\n",
                   i,j,i,j,i,j,xi,yj,
                   fvalue[i-1][j-1],pvalue);
        }
    }

   return 0;
}

int find_opt_k(void)
{

    printf("*******************************************\n");

    double Z[10+1][20+1]={{0.0}};
    double C[K][K]={{0.0}};
    get_Z(Z);
    double sigma[K]={0.0};
    int k;
    int i,j;
    double xi,yj;
    double ep=0.0;
    int m=10,n=20;
    double pvalue,fvalue;
    int optk;
    int r,s;
    for(k=1;k<K;k++)
    {
        ep = 0.0;
        for(i=0;i<(m+1);i++)
        {
            for(j=0;j<(n+1);j++)
            {
                xi=0.08*i;
                yj=0.5+0.05*j;
                pvalue = surf(xi,yj,Z,k,C);
                fvalue = *((double*)Z+(n+1)*i+j);
                ep+=(pvalue-fvalue)*(pvalue-fvalue);
            }
        }
        *(sigma+k-1)=ep;
        printf("(%d,sigma(%d))=(%d,%.11e)\n",k,k,k,ep);
        if(ep<=1.e-7)
        {
            printf("*********************************************\n");
            optk=k;
            for(r=0;r<(optk+1);r++)
            {
                for(s=0;s<(optk+1);s++)
                {
                    printf("c[%d,%d]=%.11e\n",r,s,*((double*)C+(K)*r+s));
                }
            }
            printf("\n");
            surf_with_optk(optk);

            break;
        }
    }

    return 0;
}

int get_Z(double **Z)
{
    double T[10+1][20+1]={{0.0}};
    double U[10+1][20+1]={{0.0}};
    //double Z[10+1][20+1]={{0.0}};
    int m=10+1, n=20+1;
    get_T_U(T, U, m, n);
    int i,j;

    double tt, uu;
    //solve for f(xi,yj)
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            tt = *((double*)T+n*i+j);
            uu = *((double*)U+n*i+j);
            *((double*)Z+n*i+j)=order2_interp(tt,uu);
        }
    }
//    //print (xi,yj, f(xi,yj))
//    for(i=0;i<m;i++)
//    {
//        for(j=0;j<n;j++)
//        {
//            printf("(x%d, y%d, f(x%d, y%d)) = (%.3f, %.3f, %.11e)\n",i,j,i,j,0.08*i,0.5+0.05*j,*((double*)Z+n*i+j));
//        }
//    }
    return 0;
}


double S(double x,double y,double C[][10],int k,double Z[][21])
{
	double B[11][10],G[21][10];
	double BTB[10][10],GTG[10][10];
	double U[11][21],UU[21][21];
	int i,j,m;
	double sum;
	get_B(B,k);
	get_G(G,k);
	for (i=0;i<11;i++)
		for (j=0;j<21;j++)
			U[i][j]=Z[i][j];
	for (i=0;i<=k;i++)
		for (j=0;j<=k;j++)
		{
			sum=0;
			for(m=0;m<11;m++)
				sum=sum+B[m][i]*B[m][j];
			BTB[i][j]=sum;
		}
	for (i=0;i<=k;i++)
		for (j=0;j<=k;j++)
		{
			sum=0;
			for(m=0;m<21;m++)
				sum=sum+G[m][i]*G[m][j];
			GTG[i][j]=sum;
		}
	matrix_inverse(BTB,k+1);
	matrix_inverse(GTG,k+1);
	for (i=0;i<=k;i++)
	{
		for (j=0;j<21;j++)
		{
			sum=0;
			for(m=0;m<11;m++)
				sum=sum+B[m][i]*U[m][j];
			UU[i][j]=sum;
		}
	}
	for (i=0;i<=k;i++)
	{
		for (j=0;j<=k;j++)
		{
			sum=0;
			for(m=0;m<21;m++)
				sum=sum+UU[i][m]*G[m][j];
			C[i][j]=sum;
		}
	}
	for (i=0;i<=k;i++)
	{
		for (j=0;j<=k;j++)
		{
			sum=0;
			for(m=0;m<=k;m++)
				sum=sum+BTB[i][m]*C[m][j];
			UU[i][j]=sum;
		}
	}
	for (i=0;i<=k;i++)
	{
		for (j=0;j<=k;j++)
		{
			sum=0;
			for(m=0;m<=k;m++)
				sum=sum+UU[i][m]*GTG[m][j];
			C[i][j]=sum;
		}
	}
	double result=0,r,s;
	for (i=0;i<=k;i++)
	{
		r=1;
		for (m=1;m<=i;m++)
			r=r*x;
			for (j=0;j<=k;j++)
			{
				s=1;
				for (m=1;m<=j;m++)
					s=s*y;
				result=result+C[i][j]*r*s;
			}
	}
	return result;
}


void matrix_inverse(double A[][10],int n)
 {
     /*solve for the inverse of a matrix by using Gaussian elimination algorithm.*/
	double b[10][20],c[10][10],temp;
	int i,j,k;
	for (i=0;i<n;i++)
	{
		for (j=0;j<2*n;j++)
		{
			if(j<n)
				b[i][j]=A[i][j];
			else
				b[i][j]=0;
		}
	}
	for (i=0;i<n;i++)
	b[i][n+i]=1;
	for (k=0;k<n;k++)
	{
		temp=b[k][k];
		i=k;
		while (b[k][k]==0)
		{
			b[k][k]=b[i+1][k];
			i++;
		}
		if (i>k)
		{
			b[i][k]=temp;
			for (j=0;j<k;j++)
			{
				temp=b[k][j];
				b[k][j]=b[i][j];
				b[i][j]=temp;
			}
		for (j=k+1;j<2*n;j++)
		{
			temp=b[k][j];
			b[k][j]=b[i][j];
			b[i][j]=temp;
		}
	}
	for (i=k+1;i<n;i++)
	for (j=2*n-1;j>=k;j--)
		b[i][j]=b[i][j]-b[i][k]*b[k][j]/b[k][k];
	for (j=2*n-1;j>=k;j--)
		b[k][j]=b[k][j]/b[k][k];
 }
 for (k=n-1;k>0;k--)
 {
	for (i=0;i<k;i++)
	for(j=2*n-1;j>=k;j--)
		b[i][j]=b[i][j]-b[i][k]*b[k][j];
 }
 for (i=0;i<n;i++)
 for (j=0;j<n;j++)
	c[i][j]=b[i][n+j];
 for (i=0;i<n;i++)
	for (j=0;j<n;j++)
		A[i][j]=c[i][j];
 }


double newton(double x, double y,double tuvw[4])
 {
     /*Another form of Newton algorithm*/
	 double A[4][4],b[4],var[4]={1,2,1,2},D_var[4],MaxVar,MaxD_Var;
	 int i;
	 for(;;)
	 {
		A[0][0]=-0.5*sin(var[0]);
		A[0][1]=1;
		A[0][2]=1;
		A[0][3]=1;
		A[1][0]=1;
		A[1][1]=0.5*cos(var[1]);
		A[1][2]=1;
		A[1][3]=1;
		A[2][0]=0.5;
		A[2][1]=1;
		A[2][2]=-sin(var[2]);
		A[2][3]=1;
		A[3][0]=1;
		A[3][1]=0.5;
		A[3][2]=1;
		A[3][3]=cos(var[3]);
		b[0]=-(0.5*cos(var[0])+var[1]+var[2]+var[3]-x-2.67);
		b[1]=-(var[0]+0.5*sin(var[1])+var[2]+var[3]-y-1.07);
		b[2]=-(0.5*var[0]+var[1]+cos(var[2])+var[3]-x-3.74);
		b[3]=-(var[0]+0.5*var[1]+var[2]+sin(var[3])-y-0.79);
		//solve_linear_equation(A,b,D_var,4);
		//gauss(A,b,D_var);
		gaussian_algorithm(A,b,D_var);
		MaxVar=FindMax(var);
		//MaxVar=norm_row(var,4);
		MaxD_Var=FindMax(D_var);
		//MaxD_Var=norm_row(D_var,4);
		if ((MaxD_Var/MaxVar)<(1e-12))
		{
			for (i=0;i<4;i++)
				tuvw[i]=var[i];
			break;
		}
		for (i=0;i<4;i++)
			var[i]=var[i]+D_var[i];
	}
	return 0;
}


double gaussian_algorithm(double a[4][4],double b[4],double x[4])
 {
     /*Another form of Gaussian algorithm*/
	int k,i,j,ik;
	double mik,temp,maxa,sum,tempb;
	for (k=0;k<3;k++)
	{
		temp=fabs(a[k][k]);
		for(i=k;i<4;i++)
		{
			if (fabs(a[i][k])>temp)
			{
				ik=i;
				maxa=fabs(a[i][k]);
			}
		}
		if(ik!=k)
		{
			for (j=k;j<4;j++)
			{
				temp=a[k][j];
				a[k][j]=a[ik][j];
				a[ik][j]=temp;
			}
			tempb=b[k];
			b[k]=b[ik];
			b[ik]=tempb;
		}
		for(i=k+1;i<4;i++)
		{
			mik=a[i][k]/a[k][k];
			for (j=k+1;j<4;j++)
			a[i][j]=a[i][j]-mik*a[k][j];
			b[i]=b[i]-mik*b[k];
		}
	}
	x[3]=b[3]/a[3][3];
	for (k=2;k>=0;k--)
	{
		sum=0;
		for (j=k+1;j<4;j++)
		{
			sum+=a[k][j]*x[j];
			x[k]=(b[k]-sum)/a[k][k];
		}
	}
	return 0;
 }

double FindMax(double Arr[4])
 {
	int i;
	double MaxVal;
	MaxVal=fabs(Arr[0]);
	for (i=0;i<4;i++)
	{
		if(fabs(Arr[i])>MaxVal)
		MaxVal=fabs(Arr[i]);
	}
	return (MaxVal);
 }


int find_opt_k2(void)
{

    printf("*******************************************\n");

    double Z[10+1][20+1]={{0.0}};
    double C[K][K]={{0.0}};
    get_Z(Z);
    double sigma[K]={0.0};
    int k;
    int i,j;
    double xi,yj;
    double ep=0.0;
    int m=10,n=20;
    double pvalue,fvalue;
    int optk;
    int r,s;
    for(k=1;k<K;k++)
    {
        ep = 0.0;
        for(i=0;i<(m+1);i++)
        {
            for(j=0;j<(n+1);j++)
            {
                xi=0.08*i;
                yj=0.5+0.05*j;
                pvalue = surf(xi,yj,Z,k,C);
                fvalue = *((double*)Z+(n+1)*i+j);
                ep+=(pvalue-fvalue)*(pvalue-fvalue);
            }
        }
        *(sigma+k-1)=ep;
        printf("(%d,sigma(%d))=(%d,%.11e)\n",k,k,k,ep);
    }

    return 0;
}

