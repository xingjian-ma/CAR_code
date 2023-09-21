#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
rowvec rcpp_cumsum(rowvec ratio){

  int n = ratio.n_elem;
  rowvec v(n);
  v[0] = ratio[0];
  for(int i = 1; i <= n; i++){
    v[i] = v[i-1] + ratio[i];
  }
  return v;
}

// [[Rcpp::export]]
mat rcpp_multinom(int n,rowvec ratio){
  
  int l = ratio.n_elem;
  mat m(n,l);
  rowvec range(l+1);
  range.subvec(1,l) = rcpp_cumsum(ratio);
  for(int i = 0; i < n; i++){
    NumericVector random = runif(1);
    for(int j = 0; j < l; j++){
      if(random[0] < range[j+1] && random[0] >= range[j]){
        m(i,j) = 1;
      }
    }
  }
  return m;
}

// [[Rcpp::export]]
urowvec rcpp_find_index(rowvec a){
  int n = a.n_elem;
  urowvec res(n);
  rowvec s = sort(a,"descend");
  
  for(int i = 0; i < n; i++){
    uvec x = find(s == a[i]);
    res[i] = x[0];
    }
  return res;
}

// [[Rcpp::export]]
rowvec rcpp_new_ratio(rowvec imb, rowvec allocation){
  urowvec res = rcpp_find_index(imb);
  rowvec ratio = allocation(res).t();
  return ratio/sum(ratio);
}

// [[Rcpp::export]]
mat rcpp_subset(mat m, int col, int x){
  vec column = m.col(col);
  uvec v = find(column == x);
  return m.rows(v);
}

// [[Rcpp::export]]
mat rcpp_subset_all(mat m, rowvec x, urowvec select){
  vec row_index = zeros(m.n_rows);
  for(unsigned int i = 0; i < m.n_rows; i++){
    uvec row = {i};
    uvec v = find(x == m.submat(row,select));
    if (v.n_elem == x.n_elem){
      row_index(i) = 1;
    }
  }
  uvec index = find(row_index == 1);  
  return m.rows(index);
}

// [[Rcpp::export]]
mat rcpp_get_covariate_data(int size, rowvec prob){
  uvec i1 = {0,1,2,3,4,5,6,7};
  uvec i2 = {0,1,2,3,8,9,10,11};
  uvec i3 = {0,1,4,5,8,9,12,13};
  uvec i4 = {0,2,4,6,8,10,12,14};
  
  mat raw_data = rcpp_multinom(size, prob);
  mat data1 = raw_data.cols(i1);
  mat data2 = raw_data.cols(i2);
  mat data3 = raw_data.cols(i3);
  mat data4 = raw_data.cols(i4);
  
  vec x1 = sum(data1,1);
  vec x2 = sum(data2,1);
  vec x3 = sum(data3,1);
  vec x4 = sum(data4,1);
  
  mat data = join_rows(x1,x2,x3,x4);
  return data;
}


// [[Rcpp::export]]
mat rcpp_get_data(int size,rowvec ratio,int groups,rowvec prob,mat real_data){
  
  //mat covariate_data = rcpp_get_covariate_data(size,prob);
  mat covariate_data = real_data;
  mat treat_data(size,groups);
  treat_data.row(0) = rcpp_multinom(1,ratio);
  mat data = join_rows(covariate_data,treat_data);
  return data;
}

// [[Rcpp::export]]
double rcpp_cal_overall(mat data,int covariates,int groups,rowvec ratio){
  int n = data.n_rows;
  rowvec s = sum(data,0);
  s = s.subvec(covariates,covariates+groups-1);
  double res =  sum(pow(s-n*ratio,2));
  return res;
}

// [[Rcpp::export]]
rowvec rcpp_cal_margin(mat data,int covariates,int observed,int groups,rowvec ratio){
  rowvec imb_mar_list(observed);
  int rows = data.n_rows;
  for(int i = 0; i < observed; i++){
    mat data_mar = rcpp_subset(data, i, data(rows-1,i));
    double res = rcpp_cal_overall(data_mar,covariates,groups,ratio);
    imb_mar_list(i) = res;
  }
  return imb_mar_list;

}

// [[Rcpp::export]]
double rcpp_cal_stratum(mat data,int covariates,int observed,int groups,rowvec ratio){
  int rows = data.n_rows;
  rowvec patient = data(rows-1,span(0,observed-1));
  urowvec v = linspace<urowvec>(0,observed-1,observed);
  mat data_str = rcpp_subset_all(data,patient,v);
  double res = rcpp_cal_overall(data_str,covariates,groups,ratio);
  return res;
}

// [[Rcpp::export]]
double rcpp_get_imb(mat data, int covariates, int group, rowvec ratio){
  //group:1,2,3...
  int n = data.n_rows;
  rowvec s = sum(data,0);
  double imb = s(covariates+group-1)-n*ratio(group-1);
  return imb;
}



// [[Rcpp::export]]
mat rcpp_car(int size,int covariates,int observed,int groups,
                rowvec prob,rowvec ratio,rowvec allocation,
                double wo, vec wm, double ws, mat real_data){

  mat data = rcpp_get_data(size,ratio,groups,prob,real_data);

  for(int i = 1; i < size; i++){
    rowvec imb_list(groups);
    for(int j = 0; j < groups; j++){
      data(i,covariates+j) = 1;
      mat now_data = data.rows(0,i);
      double imb_overall = rcpp_cal_overall(now_data,covariates,groups,ratio)*wo;
      rowvec imb_margin_vec = rcpp_cal_margin(now_data,covariates,observed,groups,ratio)*wm;
      double imb_margin = imb_margin_vec[0];
      double imb_stratum = rcpp_cal_stratum(now_data,covariates,observed,groups,ratio)*ws;
      imb_list[j] = imb_overall + imb_margin + imb_stratum;
      data(i,covariates+j) = 0;
    }
    
    rowvec r = rcpp_new_ratio(imb_list,allocation);
    data(i,span(covariates,covariates+groups-1)) = rcpp_multinom(1,r);

  }
  return data;
}


// [[Rcpp::export]]
mat rcpp_loop(int n,int size,int covariates,int observed,int groups,int group,
              rowvec prob,rowvec ratio,rowvec allocation,
              double wo, vec wm, double ws, mat real_data){
  mat res(1,2);
  rowvec list1(n);
  rowvec list2(n);
  rowvec list3(n);
  rowvec list4(n);
  rowvec list5(n);
  rowvec list6(n);
  rowvec list7(n);
  

  for(int i = 0; i < n; i++){
    mat data = rcpp_car(size,covariates,observed,groups,prob,ratio,allocation,wo,wm,ws,real_data);
    // list1[i] = rcpp_get_imb(rcpp_subset_all(data,{0,0},{2,3}),covariates,group,ratio)/sqrt(size);
    // list2[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1},{2,3}),covariates,group,ratio)/sqrt(size);
    // list3[i] = rcpp_get_imb(rcpp_subset_all(data,{1,0},{2,3}),covariates,group,ratio)/sqrt(size);
    // list4[i] = rcpp_get_imb(rcpp_subset_all(data,{1,1},{2,3}),covariates,group,ratio)/sqrt(size);
    list5[i] = rcpp_get_imb(rcpp_subset(data,5,0),covariates,group,ratio)/sqrt(size);
    list6[i] = rcpp_get_imb(rcpp_subset(data,5,1),covariates,group,ratio)/sqrt(size);
    list7[i] = rcpp_get_imb(rcpp_subset(data,5,2),covariates,group,ratio)/sqrt(size);
  }
  
  
  //res(0,0) = var(list1,0)+var(list2,0)+var(list3,0)+var(list4,0);
  res(0,1) = var(list5,0)+var(list6,0)+var(list7,0);
  return res;
}



// [[Rcpp::export]]
int test(rowvec a, int b){
  Function print("print");
  uvec x = find(a == b);
  print(x);

  return 0;
}


