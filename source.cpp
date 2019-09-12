#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {

	//Firstly, I want to read A and b from two input files. 
	//By dint of dynamic memory,without using unnecessary memory, I readed the files.
    ifstream infile1("A.txt");
    ifstream infile2("b.txt");
    ifstream infile3("b.txt");
    ofstream out("output.txt");//also, the output file that will show the outputs is created as output.txt .
	string line="";
    int n=0;
	
    while(getline(infile3,line)){ // by counting the lines in b.txt, program determines the size of matrix
        n++;
	}
	if(n==0) {        // if 'b.txt' has no entry, it express in the output.txt.
		cout<<"Error: b.txt has no entry!!";
		out<<"Error: b.txt has no entry!!"<<endl;
		out.close();
		return 0;
	}
	double **a = new double*[n]; // To increase precision, instead of float, I prefer to use double; while saving A as matrix, b as vector.
	double **con = new double*[n];
	double *b=new double[n];
	double *x=new double[n];	
	for(int i=0;i<n;i++){
	infile2>>b[i];
	x[i]=0;
}

	for(int i=0; i<n; i++){
  	a[i]=new double[n];
  	con[i]=new double[n];
    }
    
    for(int i=0;i<n;i++){
  		for(int j=0;j<n;j++){
  			infile1>>a[i][j] ;
  			con[i][j]=a[i][j];
			  }
		  }

        
    // Here, we start gaussian elimination with partial pivoting
    for(int k=0;k<n;k++){
        int p=k;
        for(int e=k;e<n;e++){
            if(fabs(a[e][k])>=fabs(a[p][k]))   //Initially, program finds the number that has maximum absolute value in the column 
                p=e;
        }

        if(p!=k){
            for(int i=0;i<n;i++){
                double temp=a[p][i]; //the columns are changed by program, the one that has absolute maximum value goes up.
                a[p][i]=a[k][i];
                a[k][i]=temp;
            }
            
            double temp2=b[p];
            b[p]=b[k];   // Vector b is changed simultaneously with A matrix.
            b[k]=temp2;
        }

        if(a[k][k]==0)  //program sees 0,which means we dont need to make it 0;therefore, we pass the equations below.
            continue;

        double m[n];
        for(int i=k+1;i<n;i++){  // Here, we calculate the matrix m which includes the numbers that are needed to make columns except those including pivot equal to 0.
            m[i]=a[i][k]/a[k][k];          
        }

        for(int j=k;j<n;j++){
            for(int i=k+1;i<n;i++){  // By means of m, the columns are made equal to 0.
                a[i][j]-=m[i]*a[k][j];
                
            }
        }

        for(int i=k+1;i<n;i++) //Vector b is changed simultaneously with A matrix.
            b[i]-=m[i]*b[k];

    }

    for(int j=n-1;j>=0;j--){ //when program comes to here, we have A as upper triangular matrix
        if(fabs(a[j][j])<pow(10,-8)){ //machine precision is made 10^-8.
        	cout<<"Error, A is singular.";
			out<<"Error, A is singular."<<endl;  // If one of the pivots is smaller than 10^-8, matrix is assumed as singular and program saves this expression in the output.txt .
        	if(n==2){
        	out<<"Condition number at 1: Infinity"<<endl;//According to description, we must calculate the condition numbers of every 2x2 matrix.
        	out<<"Condition number at infinity: Infinity"<<endl;//If the matrix is singular and 2x2, program saves these expressions in the output.txt . 
        }
			out.close();
			return 0; //Afterwards, quits.
		}

        x[j]=b[j]/a[j][j]; // by starting the lowermost pivot, program calculates the vector x.

        for(int i=0;i<=j-1;i++) //Then, every x value found above is subtracted from vector b;which possibles to calculate other x values.
            b[i]-=a[i][j]*x[j];   
    }
    
    for(int i=0;i<n;i++)// the function of that for loop is to save vector x in the output.txt .
       if(x[i]>-pow(10,-8) &&x[i]<pow(10,-8)) {// the machion precision is adjusted for vector x.
       	   cout<< "0"<<"\n";
		   out<<"0"<<endl;
	   } 
	   else {
	   out<<x[i]<<endl;
	   cout<<x[i]<<"\n";
	   };
    if(n==2){
        double det=0;            // program has to calculate the condition numbers of 2x2 matrixes.therefore, after checking whether the matrix is 2x2,or not;
        det=a[0][0]*a[1][1]; // program starts calculating condition numbers.Firstly,it calculates the determinant of the matrix A.
        	
        double con1i=fabs(con[0][0])+fabs(con[0][1]); // selects the row having the maximum absolute value in A
        if(con1i<fabs(con[1][0])+fabs(con[1][1]))
        	con1i=fabs(con[1][0])+fabs(con[1][1]);
        	
        double con1o=fabs(con[0][0])+fabs(con[1][0]);  //selects the column having the maximum absolute value in A
        if(con1o<fabs(con[0][1])+fabs(con[1][1]))
        	con1o=fabs(con[0][1])+fabs(con[1][1]);
        	
        double a,b,c,d;
        a=con[0][0];
        b=con[0][1];
        c=con[1][0];
        d=con[1][1];
        	
        con[0][0]=d/det;
        con[0][1]=(-1*b)/det;  // forming the inverse matrix of A.
        con[1][0]=(-1*c)/det;
        con[1][1]=a/det;
        	
        double con2i=fabs(con[0][0])+fabs(con[0][1]);
        if(con2i<fabs(con[1][0])+fabs(con[1][1]))      // selects the row having the maximum absolute value in the inverse of A 
        	con2i=fabs(con[1][0])+fabs(con[1][1]);
        	
        double con2o=fabs(con[0][0])+fabs(con[1][0]); // selects the column having the maximum absolute value in the inverse of A
        if(con2o<fabs(con[0][1])+fabs(con[1][1]))
        	con2o=fabs(con[0][1])+fabs(con[1][1]);
        		
        out<<"Condition number at 1: "<< (con1o*con2o)<<endl;                    //calculates the cond. numb. at 1 and saves in output.txt
        out<<"Condition number at infinity: "<< (con1i*con2i)<<endl;            //calculates the cond. numb. at infinity and saves in output.txt
        cout<<"Condition number at 1: "<< (con1o*con2o)<<"\n";
        cout<<"Condition number at infinity: "<< (con1i*con2i)<<"\n";
		}

	out.close();
    return 0;
}
