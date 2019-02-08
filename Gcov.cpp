#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>



// save it as .cpp
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;
using namespace arma;
using namespace std;

// [[Rcpp::export]]
Cube<double> GcovCPP(Mat<double> GG, int C, Mat<double> z, int N)  
{
  
  Col<double>  ss_diag(C); ss_diag.ones(); //     = rep(1,C)    
  Col<double>  mean_diag(C); mean_diag.ones(); // = rep(1,C)
  Mat<double>  ss_off(C,C); ss_off.zeros();
  Mat<double>  ss_bound(C,C); ss_bound.zeros();
  Mat<double>  nn(C,C); nn.zeros(); //= ss.bound = nn = matrix(0, nrow = C, ncol = C)
  Mat<double>  mean_off(C,C); mean_off.zeros(); //  = matrix(0, nrow = C, ncol = C)
  Cube<double> sss(C,C,C); sss.zeros();
  Cube<double> nnn(C,C,C); nnn.zeros(); // = nnn = array(0,dim = c(C, C, C))
  Mat<double>  gamma(C,N); gamma.fill(-99); // was a list in R.  gamma(c,n) is the n+1'st element in cluster c, -99 indicates no element  = list()
  Cube<double> Cov(C,N+1,N+1); Cov.fill(0.0); // Cov(c,,) provides covariance of values in cluster c+1       = list()
  
  Col<int>     labels(N); labels.ones(); //  cluster number = rep(1,N)
  Mat<double>  GG_1(C,C);
  Mat<double>  GG_2(C,C);

  printf("at 0\n");
  
  Col<int> nCluster(C); nCluster.zeros(); // number of items in each cluster
  
  int ii, jj, kk, ll, mm, m1, n1, ictr, r1, cc;
  double M, S_homo, S_hetero, S; 
  double ss1, a, ss2;
  
  printf("at 1\n");
  
  for(jj=0; jj<C; jj++) // for(jj in 1:C)
  {
    ictr = 0;
    for( ii=0; ii<N; ii++){
      if( z(ii,jj)==1 ){
        gamma(jj,ictr) = ii;   //gamma[1] = vector of elements which belong to cluster 1. 
        ictr++;                 
        labels(ii) = jj;       //labels[ii] = indicator of which cluster the observation belongs to. 
      }
    }
    
        nCluster(jj) = ictr;
    // gamma([ii]]       =  which(z[,ii] == 1) 
    // labels[gamma[[ii]]] =  ii
  }
  
 
  
  
  // ###Diagonal Covariance Estimates:
  // ###Total number of Diagonal Parameters: (C+1)  
  for( ii=0; ii<C; ii++ ){ //for(ii in 1:C)
    
    for( jj=ii; jj<C; jj++ ){ // for(jj in ii:C)
      
      if(ii == jj){
        
        // homo-off diagonals
          GG_1.reshape(nCluster(ii),nCluster(jj));
          for( kk=0; kk<nCluster(ii); kk++ ){
              for( ll=0; ll<nCluster(jj); ll++ ){
                  GG_1(kk,ll) = GG(gamma(ii,kk),gamma(jj,ll));
              }
          } // GG_1                = GG[gamma[[ii]], gamma[[jj]]]
        
          nn(ii,jj)           = nCluster[ii]*(nCluster[jj]-1); //          = (length(GG_1) - nrow(GG_1))
          mean_diag(ii)       = mean(GG_1.diag()); //     = mean(diag(GG_1))   #diagonal mean
          ss_diag(ii)         = var(GG_1.diag());   //diagonal variance
          M = mean_off(ii,jj) = (accu(GG_1) - accu(GG_1.diag()))/nn(ii,jj) ; //#offdiagonals mean
          /////////  Okay to subtract constant from matrix ???????????????????????????????????????????????????
          //Mat<double> M1(); M1.ones();
          
          ss1 = ss2 = 0.0;
          for(m1=0; m1<nCluster[ii]; m1++ ){
            a = GG_1(m1,m1)-M;
            ss2 += a*a;
            for( n1=0; n1<nCluster[jj]; n1++ ){
              a = GG_1(m1,n1)-M;
              ss1+= a*a;
            }
          }
          ss_off(ii,jj)       = (ss1 - ss2)/nn(ii,jj); // #offdiagonals variance
        
        }
        else
        {
        //hetero-off diagonals
            GG_2.reshape(nCluster(ii),nCluster(jj));
            for( kk=0; kk<nCluster(ii); kk++ ){
                for( ll=0; ll<nCluster(jj); ll++ ){
                    GG_2(kk,ll) = GG(gamma(ii,kk),gamma(jj,ll));
                }
            } // GG_2                   = GG[gamma[[ii]], gamma[[jj]]]
            nn(ii,jj)              = nCluster(ii)*nCluster(jj);
            M = mean_off(ii,jj)    = accu(GG_2)/nn(ii,jj);
            
            ss1 = 0.0;
            for(m1=0; m1<nCluster[ii]; m1++ ){
              for( n1=0; n1<nCluster[jj]; n1++ ){
                 a  =  GG_2(m1,n1) - M;
                ss1+= a*a;
              }
            }
            ss_off(ii,jj)          = ss1/nn(ii,jj);
        }
        }
    }
  
  printf("at 2\n");
  ss_off =  ss_off + ss_off.t() - diagmat(ss_off.diag());    
  
  //####### ss.diag and ss.off are completely cluster specific and symmetric 
  //###Boundary Covariance Estimates:  non-symmetric
  //###Total number of Boundary parameters: C  
  
  for( ii=0; ii<C; ii++ ) // for(ii in 1:C)
  {
      GG_1.reshape(nCluster(ii),nCluster(ii));
      for( kk=0; kk<nCluster(ii); kk++ ){
          for( ll=0; ll<nCluster(ii); ll++ ){
              GG_1(kk,ll) = GG(gamma(ii,kk),gamma(ii,ll));
          }
      } //GG_1                = GG(gamma((ii)), gamma((ii)))
      for( jj=0; jj<C; jj++ ) // for(jj in 1:C)
      {
          if(ii == jj)
          {
              nn(ii,jj)           = nCluster(ii)*(nCluster(jj)-1); //      = (length(GG_1) - nrow(GG_1))
              S_homo = 0.0;
              for( kk=0; kk<nCluster(ii); kk++ ) // for(kk in 1:nrow(GG_1))
              {
                  for( ll=0; ll<nCluster(jj); ll++ ){
                      S_homo     +=   (GG_1(kk,ll) - mean_off(ii,jj))*(GG_1(kk,kk) - mean_diag(ii));
                  }
          S_homo -= (GG_1(kk,kk) - mean_off(ii,jj))*(GG_1(kk,kk) - mean_diag(ii));
              }
              ss_bound(ii,jj)    =  S_homo/nn(ii,jj)  ;
          }
          else
          {
              GG_2.reshape(nCluster(ii),nCluster(jj));
              for( kk=0; kk<nCluster(ii); kk++ ){
                  for( ll=0; ll<nCluster(jj); ll++ ){
                      GG_2(kk,ll) = GG(gamma(ii,kk),gamma(jj,ll)) - mean_off(ii,jj);
                  }
              } //GG_2               = GG(gamma((ii)), gamma((jj)))
              // GG_2               = GG_2 - mean.off(ii,jj)
              nn(ii,jj)          = nCluster(ii)*nCluster(jj); // length(GG_2)
              S_hetero = 0.0;
              for( kk=0; kk<nCluster(ii); kk++ ) // for(kk in 1:nrow(GG_2))
              {
                  for( ll=0; ll<nCluster(jj); ll++){
                      S_hetero   +=   GG_2(kk,ll)*(GG_1(kk,kk) - mean_diag(ii));
                  }
              }
              ss_bound(ii,jj)     =   S_hetero/nn(ii,jj);
          }
      }
    }
  // #Offdiagonals Estimates: symmetric 
  // ## 1<=ii <=jj <= kk <= C 
  // ## Total number of parameters = C choose 3 + 2*C choose 2 + C choose 1. 
  // ## Hardest part. 
  
  printf("at 3\n");
  for( ii=0; ii<C; ii++ ) //for(ii in 1:C)
  {
    for(jj=0; jj<C; jj++ ) //for(jj in ii:C)
    {
        GG_1.reshape(nCluster(ii),nCluster(jj));
        if(jj == ii)
            {
                for( kk=0; kk<nCluster(ii); kk++ ){
                    for( ll=0; ll<nCluster(jj); ll++ ){
                        GG_1(kk,ll) = GG(gamma(ii,kk),gamma(ii,ll)) - mean_off(ii,jj);
                    }
                } // GG_1      =    GG[gamma[[ii]], gamma[[jj]]]
                // GG_1      =    (GG_1 - mean.off[ii,jj]) + diag(diag(GG_1)-mean.diag[ii]) - diag(diag(GG_1)-mean.off[ii,jj])
                GG_1 = GG_1 + diagmat(GG_1.diag()-mean_diag(ii)) - diagmat(GG_1.diag()-mean_off(ii,jj)); // ????????? okay to subtract constants
            }
        else
        {
            for( kk=0; kk<nCluster(ii); kk++ ){
                for( ll=0; ll<nCluster(jj); ll++ ){
                    GG_1(kk,ll) = GG(gamma(ii,kk),gamma(jj,ll)) - mean_off(ii,jj);
                }
            }
            // GG_1      =   GG[gamma[[ii]], gamma[[jj]]]
            // GG_1      =   GG_1 - mean.off[ii,jj]
        }
        for( kk=0; kk<C; kk++) // for(kk in jj:C)
        {
            GG_2.reshape(nCluster(ii),nCluster(kk));
            for( mm=0; mm<nCluster(ii); mm++ ){
                for( ll=0; ll<nCluster(kk); ll++ ){
                    GG_2(mm,ll) = GG(gamma(ii,mm),gamma(kk,ll)) - mean_off(ii,jj);
                }
            }
            // GG_2      =   GG[gamma[[ii]], gamma[[kk]]]
            // GG_2      =   GG_2 - mean.off[ii,kk]
            if(ii == jj && jj == kk )
            {
                // GG_2      =    GG[gamma[[ii]], gamma[[kk]]]  IS THIS NECESSARY?
                GG_2 = GG_2 + diagmat(GG_2.diag() - mean_diag(ii)) - diagmat(GG_2.diag()-mean_off(ii,kk));  //  HUH ????
                //  GG_2      =    (GG_2 - mean.off[ii,kk]) + diag(diag(GG_2) - mean.diag[ii]) - diag(diag(GG_2) - mean.off[ii,kk])
                S=0.0;
                nnn(ii,jj,kk) = nCluster(ii)*nCluster(jj)*nCluster(kk) - nCluster(ii)*nCluster(jj) - nCluster(ii)*nCluster(kk) + nCluster(ii);
                //nnn(ii,jj,kk)   =   nrow(GG_1)*ncol(GG_1)*ncol(GG_2) - length(GG_1) - length(GG_2) + nrow(GG_1)
                for(r1=0; r1<nCluster(ii); r1++ )// for(r1 in 1:nrow(GG_1))
                {
                    //S += S + accu(kron(GG_1[r1,],GG_2[r1,])) - sum(GG_1*GG_2) - sum(GG_1[r1,r1]*GG_2[r1,]) +
                    //+ sum((GG_1[r1,r1])^2)
                    S += accu(kron(GG_1.row(r1),GG_2.row(r1))) -accu(GG_1 % GG_2 ) - GG_1(r1,r1)*accu(GG_2.row(r1));
            //                                           element-wise mult intended?       ?? does armadillo support this multiplication?
                }
                sss(ii,jj,kk) = S/nnn(ii,jj,kk);
            }
            if(ii == jj && jj < kk)
            {
                S = 0.0;
                // nnn[ii,jj,kk]  =  nrow(GG_1)*ncol(GG_1)*ncol(GG_2) - length(GG_2)
                nnn(ii,jj,kk) = nCluster(ii)*nCluster(jj)*nCluster(kk) - nCluster(ii)*nCluster(kk);
                for( r1=0; r1<nCluster(ii); r1++) // for(r1 in 1:nrow(GG_1))
                {
                    //S = S + sum(kronecker(GG_1[r1,],GG_2[r1,])) - sum(GG_1[r1,r1]*GG_2[r1,])
                    S += accu(kron(GG_1.row(r1),GG_2.row(r1))) - GG_1(r1,r1)*accu(GG_2.row(r1));
                }
                sss(ii,jj,kk) = S/nnn(ii,jj,kk);
            }
            if(ii < jj && jj == kk)
            {
                S=0.0;
                nnn(ii,jj,kk) = GG_1.n_rows*GG_1.n_cols*GG_2.n_cols - GG_1.n_elem ;
                for(int r1=0; r1<GG_1.n_rows; r1++)
                {
                    //S = S + sum(kronecker(GG_1[r1,],GG_2[r1,])) - sum(GG_1[r1,]*GG_2[r1,])
                    S += accu(kron(GG_1.row(r1),GG_2.row(r1))) - accu(GG_1.row(r1) % GG_2.row(r1) );  // element-wise multiplication?
                }
                sss(ii,jj,kk) = S/nnn(ii,jj,kk);
            }
            if(ii < jj && jj < kk)
            {
                S=0.0;
                nnn(ii,jj,kk) = GG_1.n_rows*GG_1.n_cols*GG_2.n_cols;
                for( r1=0; r1<GG_1.n_rows; r1++ ){ //for(r1 in 1:nrow(GG_1))
                    {
                        //S = S +  sum(kronecker(GG_1[r1,],GG_2[r1,]))
                        S += accu(kron(GG_1.row(r1),GG_2.row(r1)));
                    }
                    sss(ii,jj,kk) = S/nnn(ii,jj,kk);
                }
                sss(ii,kk,jj) = sss(jj,kk,ii) = sss(kk,jj,ii) = sss(kk,ii,jj) = sss(jj,ii,kk) = sss(ii,jj,kk);
            }
        }
        }
    }
  
  printf("at 4\n");
  
  // ### Estimating all the cluster Covariance Matrix: 
  for( cc=0; cc<C; cc++ ) //for(cc in 1: C)
  {
    //Cov[[cc]]               = matrix(0, nrow = (N+1), ncol = (N+1))
    //Cov[[cc]][(N+1),(N+1)]  = ss.diag[cc]                 #boundary
    printf("bb\n");
    Cov(cc,N,N) = ss_diag(cc);
    printf("aa\n");
    for( ii=0; ii<N; ii++) //for(ii in 1:N)
    {
      printf("12\n");
      if( ii != gamma(cc,0) )// if(ii != gamma[[cc]][1]) 
      {
        printf("13\n");
        for( jj=ii; jj<N; jj++ ){ // for(jj in ii:N) ?????????
          {
            if( jj != gamma(cc,0) ) //if(jj != gamma[[cc]][1])
            {
              if(ii == jj)                              //diagonals 
              {
                printf("a\n");
                Cov(cc,ii,jj) = ss_off(cc,labels(ii)); //Cov[[cc]][ii,jj] = ss.off[cc,labels[ii]] 
                printf("b\n");
              } 
              else{
                printf("20\n");
                Cov(cc,ii,jj) = sss(cc,labels(ii),labels(jj)); //Cov[[cc]][ii,jj] = sss[cc,labels[ii],labels[jj]]  #offdiagonals 
                printf("at 21\n");
              }
            }
            printf("c\n");
            Cov(cc,jj,ii) = Cov(cc,ii,jj); //Cov[[cc]][jj,ii] = Cov[[cc]][ii,jj]
          }
          printf("d\n");
          Cov(cc,N,ii) = Cov(cc,ii,N) = ss_bound(cc,labels(ii)); //Cov[[cc]][(N+1), ii] = Cov[[cc]][ii, (N+1)]  =  ss.bound[cc, labels[ii]]                #boundary
        }
      }
      printf("endif\n");
      //# Covariance on the Average Column 
      //Cov[[cc]][gamma[[cc]][1], 1:N]              = colSums(Cov[[cc]][1:N, 1:N])/(N-1)
      printf("at 5\n");
      //Cov.subcube(cc,gamma(cc,0),0, cc,gamma(cc,0),N-1) = sum(Cov.subcube(cc,0,0, cc,N-1,N-1))/(N-1);
      for(int i=0;i<N;i++){
        a=0;
        for(int j=0; j<N; j++){
          a+= Cov(cc,j,i);
        }
        Cov(cc,gamma(cc,0),i) = a/(N-1);
      }
      //Cov.subcube(cc,gamma(cc,0),0, cc,gamma(cc,0),N-1) = sum(Cov(cc))/(N-1);
      for( int i=0; i<N; i++){
        Cov(cc,i,gamma(cc,0)) = Cov(cc,gamma(cc,0),i);
      }
      // Cov[[cc]][1:N, gamma[[cc]][1]]              = t(Cov[[cc]][gamma[[cc]][1], 1:N])
      printf("at 6\n");
      //Cov.subcube(cc,0,gamma(cc,0), cc,N-1,gamma(cc,0)) = Cov.subcube(cc,gamma(cc,0),0, cc,gamma(cc,0), N-1);
      for(int i=0;i<N;i++){
       
        for(int j=0; j<N; j++){
       
        }
        Cov(cc,gamma(cc,0),i)=Cov(cc,i,gamma(cc,0)) = a/(N-1);
      }
      //Cov[[cc]][gamma[[cc]][1], gamma[[cc]][1]]   = sum(diag(Cov[[cc]])[1:N])/(N-1)
      printf("at 7\n");
      mat temp = Cov.subcube(cc,0,0, cc, N-1, N-1); 
      printf("at 8\n");
      Cov(cc,gamma(cc,0),gamma(cc,0))  = accu(temp.diag())/(N-1); // ???????
      printf("at 9\n");
      // Cov[[cc]][(N+1), gamma[[cc]][1]]            = Cov[[cc]][gamma[[cc]][1],(N+1)]  = sum(Cov[[cc]][1:N,(N+1)])/(N-1)
      Cov(cc,N,gamma(cc,0)) = Cov(cc,gamma(cc,0),N) = accu(Cov.subcube(cc,0,N, cc, N-1, N))/(N-1);
    }
    
   // return(Cov);
    
  }
  return(Cov);
}
