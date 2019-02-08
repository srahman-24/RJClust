
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>



// save it as .cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

 *//* Variable definitions:
     X      input data matrix with N rows and P columns
     XX     XX'/P, the R matrix from the article
     J      J matrix, derived from XX'
     N      number of items to be clustered
     P      number of features measured on each item
     Cov    Armadillo cube of dimension N+1,N+1,C that stores the
                estimated coveriace matrix for each cluster
     Mu     Armadillo matrix of dimension N+1,C storing the C N+1 mean
                vectors.
     I  Vector of labels of cluster membership
     Alpha  Vector of mixture weights, length C and summing to 1
     */

int Mstep(Mat<double> J,Mat<double> Mu, Cube<double> Cov, Col<double> Alpha, Col<int> I){
    // implements the M-step of the EM algorithm for given W and J
    
    int i, j, k, m;
    unsigned long long labeli, labelj, labelk;
    unsigned long long N = J.n_rows;
    unsigned long long C = Alpha.n_elem;
    Cube<double> CovCluster(C+1,C+1,C); CovCluster.zeros(); //
    Mat<double> MuCluster(C,C+1); MuCluster.zeros(); // work vectors, note that C index is first
                                             // unlike mu, Cov where cluster index is last
    
    // begin by estimating Alpha
    Alpha.zeros();
    for( i=0; i<N; i++ ){
        Alpha(I(i))++;
    }
    for( i=0; i<C; i++ ) Alpha(i) /= (double)N;
    if( sum(Alpha) != 1.0 ){
        cout << "Alpha sums to " << sum(Alpha) << "\n";
        exit(-9);
    }
    
    // Update Mu
    // Mu contains the vector means of the rows of J, but these means depend only on cluster
    // There are C choose 2 cominations, plus C diagonal elements to estimate.  Note that this
    // function assumes discrete membership labels are availabe in Labels
    
    Mat<int> tempN(C,C+1); tempN.zeros();
    
    for( i=0; i<N; i++){
        labeli = I(i);
        for( j=(i+1); j<(N+1); j++ ){
            if( j<N ){ // non-diagonal entry
                labelj = I(j);
                if( labeli <= labelj ){
                    MuCluster(labeli,labelj) += J(i,j);
                    tempN(labeli,labelj) += 1;
                } else {
                    MuCluster(labelj,labeli) += J(i,j);
                    tempN(labelj,labeli) += 1;
                }
            }
            else{  // this gets mean for final column of J
                MuCluster(labeli,C) += J(i,N);
                tempN(labeli,C) += 1;
            }
        }
    }
    // Divide by sum of weights to get cluster by cluster means
    for( k=0; k<C; k++){
        for( m= k; m<(C+1); m++){
            MuCluster(k,m) /= (double) tempN(k,m);
            if( m<C ) MuCluster(m,k)=MuCluster(k,m);
        }
    }
    
    // Update Cov
    // The first coordinate of Cov refers to the cluster label of the row of J, the next two coordinates
    // refer to the labels of the two components in that row
    
    // First estimate covariances.  Treat last column in each cluster as cluster C.  No diagonal elements
    // or squares of the same element are used.
    Cube<int> tempN2(C+1,C+1,C); tempN2.zeros();
    for( i=0; i<N; i++){
        labeli = I(i);
        for( j=0; j<N; j++){
            if( j==i ) continue;
            if( j<N ){
                labelj = I(j);
            } else {
                labelj = C;
            }
            for( k=(j+1); k<(N+1); k++){
                if( k==i )continue;
                if( k<N ){
                    labelk = I(k);
                } else {
                    labelk = C;
                }
                if( labelj <= labelk ){
                    CovCluster(labelj,labelk,labeli) += (J(i,j)-MuCluster(labeli,labelj))*(J(i,k)-MuCluster(labeli,labelk));
                    tempN2(labelj,labelk,labeli) += 1;
                } else{
                    CovCluster(labelk,labelj,labeli) += (J(i,j)-MuCluster(labeli,labelj))*(J(i,k)-MuCluster(labeli,labelk));
                    tempN2(labelk,labelj,labeli) += 1;
                    cout << "labeli, labelj, labelk, i, j, J(i,j) J(i,k) muC(labeli,labelj) muC(labeli,labelk) \n";
                    cout <<labeli << " " << labelj << " " << labelk << " " << i << " " << j << " " << J(i,j)<< " "  << J(i,k) << " " << MuCluster(labeli,labelj) << " " << MuCluster(labeli,labelk) << "\n";
                }
            }
        }
    }
    
    // Divide to get mean for covariances
    for( i=0; i<C; i++){
        for( j=0; j<C; j++ ){
            for( k=j; k<(C+1); k++){
                if( tempN2(j,k,i)<1 ){
                    cout << "Cant estimate a covariance because of 0 tempN2 in Mstep" << i << j << k << tempN2(i,j,k) << "\n";
                    cout << "MuCluster : " << MuCluster << "\n";
                    cout << "tempN2: " << tempN2 << "\n";
                    exit(9);
                }
                CovCluster(j,k,i) /= (double) tempN2(j,k,i);
                CovCluster(k,j,i) = CovCluster(j,k,i);
            }
        }
    }
    // Now estimate variance of elements of J(,N) and other
    mat varDiag(C+1,C+1); varDiag.zeros();
    Mat <int> tempVar(C+1,C+1); tempVar.zeros();
    for(i=0; i<N; i++ ){
        labeli = I(i);
        for(j=(i+1); j<(N+1); j++ ){
            if( j<N ){
                labelj = I(j);
            } else {
                labelj = C;
            }
            if( labeli<=labelj){
                varDiag(labeli,labelj) += (J(i,j)-MuCluster(labeli,labelj))*(J(i,j)-MuCluster(labeli,labelj));
                tempVar(labeli,labelj) += 1;
            }else{
                varDiag(labelj,labeli) += (J(i,j)-MuCluster(labelj,labeli))*(J(i,j)-MuCluster(labelj,labeli));
                tempVar(labelj,labeli) += 1;
            }
        }
    }
    
    // Divide to get mean for variances
    for( i=0; i<C; i++){
        for( j=i; j<(C+1); j++ ){
            if( tempVar(i,j)<1 ){
                cout << "0 count in calculation of variances in Mstep\n";
                exit(9);
            }
            varDiag(i,j) /= (double)tempVar(i,j);
            varDiag(j,i) = varDiag(i,j);
        }
    }
    
    
    // Transfer MuCluster to Mu
    for( i=0; i<C; i++ ){
        for( j=0; j<N; j++ ){
            labelj = I(j);
            Mu(j,i) = MuCluster(i,labelj);
        }
        Mu(N,i) = MuCluster(i,C);
    }
    
    
    // Transfer CovCluster to Cov
    for( i=0; i<C; i++ ){
        for( j=0; j<(N+1); j++){
            if( j<N ){
                labelj = I(j);
            } else {
                labelj = C;
            }
            if( j<N) Cov(j,j,i) = varDiag(labelj,i);
            else Cov(j,j,i) = varDiag(i,C);
            for( k=(j+1); k<(N+1); k++ ){
                if( k<N ){
                    labelk = I(k);
                } else {
                    labelk = C;
                }
                Cov(j,k,i) = Cov(k,j,i) = CovCluster(labelj,labelk,i);
            }
        }
    }
     
    return(0);
}
