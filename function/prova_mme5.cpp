
#include<iostream>
#include <RcppArmadillo.h>
#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]


using namespace std;
using namespace arma;
using namespace Eigen;

// class of sparse matrix for rccp
// copy form https://www.r-bloggers.com/2021/04/constructing-a-sparse-matrix-class-in-rcpp/


/*
//conversion 
SparseMatrixXd cast_eigen(sp_mat arma_A) {

  SparseMatrixXd eigen_B = Map<SparseMatrixXd>(arma_A.memptr(),
                                                        arma_A.n_rows,
                                                        arma_A.n_cols);
  return eigen_B;
}

// 
mat cast_arma(MatrixXd eigen_A) {

  mat arma_B = mat(eigen_A.data(), eigen_A.rows(), eigen_A.cols(),
                               false, false);

  return arma_B;
}
*/


//MatrixXd PCG_solve(SparseMatrixXd  A) {
 // int size = A.cols(); 
  //MatrixXd Caa(size,size);
  //MatrixXd b(size, size)  ; // dichiara b
  //MatrixXd x(size, size)  ; // dichiara b
  //b=MatrixXd::Identity(size, size);
  //x=b;
  //ConjugateGradient<SparseMatrixXd>, Lower|Upper> cg;
  //cg.compute(A);
  //x = cg.solve(b);
  //cout << "#iterations:     " << cg.iterations() << endl;
  //cout << "estimated error: " << cg.error() << endl;
  // update b, and solve again
  //Caa = cg.solve(b);
  //return Caa;       
//}


// [[Rcpp::export]]

sp_mat random_sparse() {
for (i = 0; i < 1000; i++)
{
    if (rand() % a == 0)
    {
        array[i] = rand() % 100;
    }
    else
    {
        array[i] = 0;
    }
      return lhs;
}
}
// [[Rcpp::export]]
sp_mat  fast_five_traitsMM( sp_mat   A, sp_mat  Z1, sp_mat  Z2,sp_mat VARGi,sp_mat  e) {
 
  int size = std::sqrt(Z1.size());
  //cout << size  <<endl;
  sp_mat  lhs11(size*3, size);
  sp_mat  lhs1(size*5, size);
  sp_mat  lhs22(size*3, size);
  sp_mat  lhs2(size*5, size);
  sp_mat  lhs33(size*3, size);
  sp_mat  lhs3(size*5, size);
  sp_mat  lhs44(size*3, size);
  sp_mat  lhs4(size*5, size);
  sp_mat  lhs55(size*3, size);
  sp_mat  lhs5(size*5, size);
  sp_mat  lhs_tmp(size*3, size*5);
  sp_mat  lhs(size*5, size*5);
  //mat caa (size*5, size*5);
  //SparseMatrixXd lhs_e(size*5, size*5);
 //MatrixXd I(size*5,size*5)  ; 
  //I = MatrixXd::Identity(size*5, size*5);
  //MatrixXd Caa(size*5, size*5);
  //Rcpp::NumericMatrix lhs2 = lhs1;  Rcpp::NumericMatrix lhs2 = lhs1; Rcpp::NumericMatrix lhs3 = lhs1; 
   //Rcpp::NumericMatrix lhs5 = lhs1; 
  
  sp_mat Z3 = Z1;
  sp_mat Z4 = Z1;
  //sex limit traits
  sp_mat Z5 = Z2;
  
  sp_mat  Z1_Z1A = (trans(Z1)*Z1)*e(0,0) + A*VARGi(0,0);
  sp_mat  Z1_Z2A = (trans(Z1)*Z2)*e(0,1) + A*VARGi(0,1);
  sp_mat  Z1_Z3A = (trans(Z1)*Z3)*e(0,2) + A*VARGi(0,2);
  sp_mat  Z1_Z4A = (trans(Z1)*Z4)*e(0,3) + A*VARGi(0,3);
  sp_mat  Z1_Z5A = (trans(Z1)*Z5)*e(0,4) + A*VARGi(0,4);
  
  sp_mat  Z2_Z1A = (trans(Z2)*Z1)*e(1,0) + A*VARGi(1,0);
  sp_mat  Z2_Z2A = (trans(Z2)*Z2)*e(1,1) + A*VARGi(1,1);
  sp_mat  Z2_Z3A = (trans(Z2)*Z3)*e(1,2) + A*VARGi(1,2);
  sp_mat  Z2_Z4A = (trans(Z2)*Z4)*e(1,3) + A*VARGi(1,3);
  sp_mat  Z2_Z5A = (trans(Z2)*Z5)*e(1,4) + A*VARGi(1,4);
  
  sp_mat  Z3_Z1A = (trans(Z3)*Z1)*e(2,0) + A*VARGi(2,0);
  sp_mat  Z3_Z2A = (trans(Z3)*Z2)*e(2,1) + A*VARGi(2,1);
  sp_mat  Z3_Z3A = (trans(Z3)*Z3)*e(2,2) + A*VARGi(2,2);
  sp_mat  Z3_Z4A = (trans(Z3)*Z4)*e(2,3) + A*VARGi(2,3);
  sp_mat  Z3_Z5A = (trans(Z3)*Z5)*e(2,4) + A*VARGi(2,4);

  sp_mat  Z4_Z1A = (trans(Z4)*Z1)*e(3,0) + A*VARGi(3,0);
  sp_mat  Z4_Z2A = (trans(Z4)*Z2)*e(3,1) + A*VARGi(3,1);
  sp_mat  Z4_Z3A = (trans(Z4)*Z3)*e(3,2) + A*VARGi(3,2);
  sp_mat  Z4_Z4A = (trans(Z4)*Z4)*e(3,3) + A*VARGi(3,3);
  sp_mat  Z4_Z5A = (trans(Z4)*Z5)*e(3,4) + A*VARGi(3,4);
 
  sp_mat  Z5_Z1A = (trans(Z5)*Z1)*e(4,0) + A*VARGi(4,0);
  sp_mat  Z5_Z2A = (trans(Z5)*Z2)*e(4,1) + A*VARGi(4,1);
  sp_mat  Z5_Z3A = (trans(Z5)*Z3)*e(4,2) + A*VARGi(4,2);
  sp_mat  Z5_Z4A = (trans(Z5)*Z4)*e(4,3) + A*VARGi(4,3);
  sp_mat  Z5_Z5A = (trans(Z5)*Z5)*e(4,4) + A*VARGi(4,4);

  lhs11 = join_cols( Z1_Z1A, Z2_Z1A,Z3_Z1A);
  lhs1 = join_cols(lhs11,Z4_Z1A,Z5_Z1A);
  //
  lhs22 = join_cols( Z2_Z1A, Z2_Z2A,Z2_Z3A);
  lhs2 = join_cols(lhs22,Z2_Z4A,Z2_Z5A);
  //
  lhs33 = join_cols( Z3_Z1A, Z3_Z2A,Z3_Z3A);
  lhs3 = join_cols(lhs33,Z3_Z4A,Z3_Z5A);
  //
  lhs44 = join_cols( Z4_Z1A, Z4_Z2A,Z4_Z3A);
  lhs4 = join_cols(lhs44,Z4_Z4A,Z4_Z5A);  
  //
  lhs55 = join_cols( Z5_Z1A, Z5_Z2A,Z5_Z3A);
  lhs5 =join_cols(lhs55,Z5_Z4A,Z5_Z5A);

  lhs_tmp = join_rows(lhs1,lhs2,lhs3);
  lhs=join_rows(lhs_tmp,lhs4,lhs5);

  cout<<"lhs done"<<endl;
  // very stupid process but i was rush 
  return lhs;
}


// [[Rcpp::export]]
// extremly slow but save lot of ram
MatrixXd solvepcg (SparseMatrix<double> lhs){
  //int core = 7;
  //#pragma omp parallel num_threads(core) {
   // int kk = omp_get_thread_num();
  
//  #pragma omp for schedule(static)
//Eigen::setNbThreads(n);
  cout<<"lhs inversion"<<endl;
   int size = lhs.cols();
   cout << "dim matrix: "<< size << endl;
   MatrixXd I(size,size)  ; 
   MatrixXd Caa(size,size);
   I =  MatrixXd::Identity(size, size);
 // Caa = lhs.llt().solve(I);
  //return Caa;
ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
cg.compute(lhs);
Caa = cg.solve(I);
//  }
std::cout << "#iterations:     " << cg.iterations() << std::endl;
std::cout << "estimated error: " << cg.error()      << std::endl;

return Caa;
cout<<"lhs done"<<endl;
} 



// [[Rcpp::export]]
MatrixXd solveLU (SparseMatrix<double> lhs){
    cout<<"lhs inversion"<<endl;
   int size = lhs.cols();
   cout << "dim matrix: "<< size << endl;
   MatrixXd I(size,size)  ; 
   MatrixXd Caa(size,size);
   I =  MatrixXd::Identity(size, size);
 // Caa = lhs.llt().solve(I);
  //return Caa;
// fill A and b
 SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >   solver;
// fill A and b;
// Compute the ordering permutation vector from the structural pattern of A
solver.analyzePattern(lhs); 
// Compute the numerical factorization 
solver.factorize(lhs); 
//Use the factors to solve the linear system 
Caa = solver.solve(I);
   return Caa;
} 

// [[Rcpp::export]]
mat invArm (mat lhs){
   cout<<"lhs inversion"<<endl;
    int size = std::sqrt(lhs.size());
    cout << "dim matrix: "<< size << endl;
    mat caa (size, size);
    caa = inv(lhs);
    return caa;
} 






/*
  //lhs_e = cast_eigen(lhs);
   Caa = lhs_e.llt().solve(I);
  //Caa = lhs_e;
   //Caa = PCG_solve(lhs_e);
  cout<<"lhs inverted"<<endl;
  return Caa;

}
*/



// [[Rcpp::export]]

Rcpp::List re_orderCaa(Rcpp::NumericMatrix TMP,int p, int n) {
    Rcpp::NumericMatrix A(p,p);
    int ANIMAL=0;
    Rcpp::List A_LIST(n);
    //cout << "number of animal:"<< n << endl;
    for(int ANIMAL=0;ANIMAL<n;ANIMAL++) {
        int riga = 0;
        for(int h=0;h<p;h++ ){
            int colonna=0;
            for (int i=0;i<p;i++) {
              double io =  TMP( ANIMAL+(n*h) , ANIMAL+(n*i) ) ;
              A(riga,colonna) = io;
              colonna= colonna +1 ;
            }
    riga = riga +1  ;
    }
    A_LIST[ANIMAL] = clone(A); // CHECK COLONE FUNZIONALITY IN RCCP
 //  cout << ANIMAL << endl; 
}//CAORE QUA
   return A_LIST;
    
}


// [[Rcpp::export]]
mat   SOLVEMME(mat pheno,mat e, mat Z1 , mat Z2,mat c22) {
   mat Z3 = Z1;
   mat Z4 = Z1;
  //sex limit traits
   mat Z5 = Z2;
   
   int size = std::sqrt(Z1.size());

   mat R1(size,1);
   mat R2(size,1);
   mat R3(size,1);
   mat R4(size,1);
   mat R5(size,1);
   mat rhs_tmp(size*3,1);
   mat rhs(size*5,1);
   mat sol(size*5,1);

   //cout << size << endl;
   R1 =Z1*(pheno.col(0)*e(0,0) + pheno.col(1)*e(1,0)+pheno.col(2)*e(2,0) +  pheno.col(3)*e(3,0) +  pheno.col(4)*e(4,0));
   R2 =Z2*(pheno.col(1)*e(1,1) + pheno.col(1)*e(1,1)+pheno.col(2)*e(2,1) +  pheno.col(3)*e(3,1) +  pheno.col(4)*e(4,1));
   R3 =Z3*(pheno.col(2)*e(2,2) + pheno.col(1)*e(1,2)+pheno.col(2)*e(2,2) +  pheno.col(3)*e(3,2) +  pheno.col(4)*e(4,2));
   R4 =Z4*(pheno.col(3)*e(3,3) + pheno.col(1)*e(1,3)+pheno.col(2)*e(2,3) +  pheno.col(3)*e(3,3) +  pheno.col(4)*e(4,3));
   R5 =Z5*(pheno.col(4)*e(4,4) + pheno.col(1)*e(1,4)+pheno.col(2)*e(2,4) +  pheno.col(3)*e(3,4) +  pheno.col(4)*e(4,4));
  
   rhs_tmp = join_cols( R1,R2,R3);
   rhs = join_cols(rhs_tmp,R4,R5);  
   sol = trans(rhs)*c22;
   return trans(sol);
}
