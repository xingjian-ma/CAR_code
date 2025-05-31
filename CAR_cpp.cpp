#include <RcppArmadillo.h>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double normalCDF(double value)
{
  return 0.5 * erfc(-value * M_SQRT1_2);
}


// [[Rcpp::export]]
double logit(double value)
{
  return 1/(1+exp(-value))  ;
}

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
double allocation_function(double value, double limit){
  return -2*limit*(1/(1+exp(-value)) - 0.5);
}


// [[Rcpp::export]]
rowvec rcpp_sym_ratio(rowvec imb, rowvec ratio){
  rowvec m1 = min(ratio,1-ratio);
  double min_value = min(m1);
  
  int k = imb.n_elem;
  double limit = min_value/(k-1);
  
  
  for(int i = 0; i < k; i++){
    for(int j = i+1; j < k; j++){
      double dif = allocation_function(imb[i]-imb[j], limit);
      ratio[i] = ratio[i] + dif;
      ratio[j] = ratio[j] - dif;
      
    }
    
  }
  
  return ratio;
}




// [[Rcpp::export]]
//col:第几列，从0开始
//x:取值
mat rcpp_subset(mat m, int col, int x){
  vec column = m.col(col);
  uvec v = find(column == x);
  return m.rows(v);
}

// [[Rcpp::export]]
mat rcpp_subset_all(mat m, urowvec select, rowvec x){
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

// // [[Rcpp::export]]
// mat rcpp_get_covariate_data(int size, rowvec prob){
//   uvec i1 = {0,1,2,3,4,5,6,7};
//   uvec i2 = {0,1,2,3,8,9,10,11};
//   uvec i3 = {0,1,4,5,8,9,12,13};
//   uvec i4 = {0,2,4,6,8,10,12,14};
//   
//   mat raw_data = rcpp_multinom(size, prob);
//   mat data1 = raw_data.cols(i1);
//   mat data2 = raw_data.cols(i2);
//   mat data3 = raw_data.cols(i3);
//   mat data4 = raw_data.cols(i4);
//   
//   vec x1 = sum(data1,1);
//   vec x2 = sum(data2,1);
//   vec x3 = sum(data3,1);
//   vec x4 = sum(data4,1);
//   
//   mat data = join_rows(x1,x2,x3,x4);
//   return data;
// }

// [[Rcpp::export]]
mat rcpp_get_covariate_data(int size, int covariates, int observed, rowvec beta){
  mat data(size, 12);
  for(int i = 0; i < size; i++ ){
    NumericVector random = runif(12);
    for(int j = 0; j < 10; j++){
      if(random[j] > 0.5){
        data(i,j) = 1;
      }
    }
   

      if(random[10] <  normalCDF(dot(data(i,span(0,9)),beta(span(0,9))))) {
        data(i,10) = 1;
      
    }

      if(random[11] <  logit(dot(data(i,span(0,9)),beta(span(0,9))))) {
        data(i,11) = 1;
        
      }
  }
  
  return data;
}




// [[Rcpp::export]]
mat rcpp_get_data(int size,rowvec ratio,int groups,rowvec beta, mat real_data, int covariates, int observed){
  
  mat covariate_data = rcpp_get_covariate_data(size, covariates, observed, beta);
  //mat covariate_data = real_data;
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
rowvec rcpp_cal_margin(mat data,int covariates,int observed,int groups, rowvec ratio){
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
  mat data_str = rcpp_subset_all(data,v,patient);
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
                rowvec ratio,rowvec allocation,
                double wo, vec wm, double ws, mat data){

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
mat rcpp_cr(int size,int covariates,int observed,int groups,
            rowvec ratio, mat data){
  data.submat(span(0,size-1),span(covariates,covariates+groups-1)) = rcpp_multinom(size,ratio);
  
  return data;
}


// [[Rcpp::export]]
mat rcpp_rblock(rowvec ratio, int groups){
  rowvec index(groups+1);
  index(span(1,groups)) = rcpp_cumsum(ratio)*10;
  mat block(10, groups);
  
  for(int i = 0; i < groups; i++){
    vec c(ratio[i]*10, fill::ones);
    block(span(index[i], index[i]+ratio[i]*10-1),i)  = c;
  }
  
  return shuffle(block,0);
}


// [[Rcpp::export]]
mat rcpp_str(int size,int covariates,int observed,int groups,
            rowvec ratio, mat data){

  int strata = pow(2, observed);
  
  int block_n = ceil(size/10.0);
  cube block(10*block_n, groups, strata);
  

  
  for(int i = 0; i < strata; i++){
    for(int j = 0; j < block_n; j++){
    block(span(10*j,10*j+9),span(0,groups-1),span(i,i)) = rcpp_rblock(ratio,groups);
    }
  }

  for(int i = 0; i < size; i++){
    rowvec patient = data(i,span(0,observed-1));
    rowvec v = {1,2,4,8,16,32,64,128,256,512};
    int slice = dot(v(span(0,observed-1)), patient);
    urowvec r = linspace<urowvec>(0,observed-1,observed);
    mat data_str = rcpp_subset_all(data.rows(0,i),r,patient);
    int index = data_str.n_rows-1;
    rowvec treat = block(span(index, index),span(0,groups-1),span(slice,slice));
    data(i,span(covariates,covariates+groups-1)) = treat;
  }
  return data;
}



// [[Rcpp::export]]
mat rcpp_loop(int n,int size,int covariates,int observed,int groups,int group,
              rowvec beta,rowvec ratio,rowvec allocation,
              double wo, vec wm, double ws, mat real_data){
  mat res(3,2);
  rowvec list1(n);
  rowvec list2(n);
  rowvec list3(n);
  rowvec list4(n);
  rowvec list5(n);
  rowvec list6(n);
  rowvec list7(n);
  rowvec list8(n);
  

  for(int i = 0; i < n; i++){
    mat data = rcpp_get_data(size,ratio,groups,beta,real_data,covariates,observed);
    data = rcpp_car(size,covariates,observed,groups,ratio,allocation,wo,wm,ws,data);
    
  
    
    // list1[i] = rcpp_get_imb(rcpp_subset_all(data,{0,0},{2,3}),covariates,group,ratio)/sqrt(size);
    // list2[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1},{2,3}),covariates,group,ratio)/sqrt(size);
    // list3[i] = rcpp_get_imb(rcpp_subset_all(data,{1,0},{2,3}),covariates,group,ratio)/sqrt(size);
    // list4[i] = rcpp_get_imb(rcpp_subset_all(data,{1,1},{2,3}),covariates,group,ratio)/sqrt(size);
    //list5[i] = rcpp_get_imb(rcpp_subset(data,5,0),covariates,group,ratio)/sqrt(size);
    //list6[i] = rcpp_get_imb(rcpp_subset(data,5,1),covariates,group,ratio)/sqrt(size);
    //list7[i] = rcpp_get_imb(rcpp_subset(data,5,2),covariates,group,ratio)/sqrt(size);

    // p1 - p3

    list1[i] = rcpp_get_imb(data,covariates,group,ratio)/sqrt(size);
    list2[i] = rcpp_get_imb(rcpp_subset_all(data,{0},{0}),covariates,group,ratio)/sqrt(size);
    list3[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1},{0,0}),covariates,group,ratio)/sqrt(size);
    
    
    //p4

    // list1[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1,10,11},{0,0,0,0}),covariates,group,ratio)/sqrt(size);
    // list2[i] = rcpp_get_imb(rcpp_subset_all(data,{12,13},{1,1}),covariates,group,ratio)/sqrt(size);

  }
  
  
  res(0,0) = mean(list1);
  res(0,1) = sqrt(var(list1,0));
  
  res(1,0) = mean(list2);
  res(1,1) = sqrt(var(list2,0));
  
  res(2,0) = mean(list3);
  res(2,1) = sqrt(var(list3,0));
  
  // res(3,0) = mean(list4);
  // res(3,1) = sqrt(var(list4,0));
  
  
  return res;
}



// [[Rcpp::export]]
mat rcpp_cr_loop(int n,int size,int covariates,int observed,int groups,int group,
              rowvec beta,rowvec ratio,rowvec allocation,
              double wo, vec wm, double ws, mat real_data){
  mat res(2,2);
  rowvec list1(n);
  rowvec list2(n);
  rowvec list3(n);
  rowvec list4(n);
  rowvec list5(n);
  rowvec list6(n);
  rowvec list7(n);
  rowvec list8(n);
  
  
  for(int i = 0; i < n; i++){
    mat data = rcpp_get_data(size,ratio,groups,beta,real_data,covariates,observed);
    data = rcpp_cr(size,covariates,observed,groups,ratio,data);
    
    
    
    //list4[i] = rcpp_get_imb(rcpp_subset_all(data,{1,1},{2,3}),covariates,group,ratio)/sqrt(size);
    //list5[i] = rcpp_get_imb(rcpp_subset(data,5,0),covariates,group,ratio)/sqrt(size);
    //list6[i] = rcpp_get_imb(rcpp_subset(data,5,1),covariates,group,ratio)/sqrt(size);
    //list7[i] = rcpp_get_imb(rcpp_subset(data,5,2),covariates,group,ratio)/sqrt(size);
    
    
    //p1 - p3
    
    // list1[i] = rcpp_get_imb(data,covariates,group,ratio)/sqrt(size);
    // list2[i] = rcpp_get_imb(rcpp_subset_all(data,{0},{0}),covariates,group,ratio)/sqrt(size);
    // list3[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1},{0,0}),covariates,group,ratio)/sqrt(size);
    
    
    //p4
    
    list1[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1,10,11},{0,0,0,0}),covariates,group,ratio)/sqrt(size);
    list2[i] = rcpp_get_imb(rcpp_subset_all(data,{12,13},{0,0}),covariates,group,ratio)/sqrt(size);
    
    
  }
  
  
  res(0,0) = mean(list1);
  res(0,1) = sqrt(var(list1,0));
  
  res(1,0) = mean(list2);
  res(1,1) = sqrt(var(list2,0));
  
  // res(2,0) = mean(list3);
  // res(2,1) = sqrt(var(list3,0));
  
  // res(3,0) = mean(list4);
  // res(3,1) = sqrt(var(list4,0));
  
  
  return res;
}


// [[Rcpp::export]]
mat rcpp_str_loop(int n,int size,int covariates,int observed,int groups,int group,
                 rowvec beta,rowvec ratio,rowvec allocation,
                 double wo, vec wm, double ws, mat real_data){
  mat res(2,2);
  rowvec list1(n);
  rowvec list2(n);
  rowvec list3(n);
  rowvec list4(n);
  rowvec list5(n);
  rowvec list6(n);
  rowvec list7(n);
  rowvec list8(n);


  for(int i = 0; i < n; i++){
    mat data = rcpp_get_data(size,ratio,groups,beta,real_data,covariates,observed);
    data = rcpp_str(size,covariates,observed,groups,ratio,data);




    // list4[i] = rcpp_get_imb(rcpp_subset_all(data,{1,1},{2,3}),covariates,group,ratio)/sqrt(size);
    //list5[i] = rcpp_get_imb(rcpp_subset(data,5,0),covariates,group,ratio)/sqrt(size);
    //list6[i] = rcpp_get_imb(rcpp_subset(data,5,1),covariates,group,ratio)/sqrt(size);
    //list7[i] = rcpp_get_imb(rcpp_subset(data,5,2),covariates,group,ratio)/sqrt(size);


    //p1 - p3

    // list1[i] = rcpp_get_imb(data,covariates,group,ratio)/sqrt(size);
    // list2[i] = rcpp_get_imb(rcpp_subset_all(data,{0},{0}),covariates,group,ratio)/sqrt(size);
    // list3[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1},{0,0}),covariates,group,ratio)/sqrt(size);


    //p4

    list1[i] = rcpp_get_imb(rcpp_subset_all(data,{0,1,12,13},{0,0,0,0}),covariates,group,ratio)/sqrt(size);
    list2[i] = rcpp_get_imb(rcpp_subset_all(data,{12,13},{0,0}),covariates,group,ratio)/sqrt(size);


  }


  res(0,0) = mean(list1);
  res(0,1) = sqrt(var(list1,0));

  res(1,0) = mean(list2);
  res(1,1) = sqrt(var(list2,0));

  // res(2,0) = mean(list3);
  // res(2,1) = sqrt(var(list3,0));

  // res(3,0) = mean(list4);
  // res(3,1) = sqrt(var(list4,0));


  return res;
}


// [[Rcpp::export]]
cube sd(cube input, int n){
  
  cube final = sqrt((sum(pow(input,2),2) - n*pow(mean(input,2),2))/(n-1));
  
  return final;
}

// [[Rcpp::export]]
cube rcpp_all_loop(int n,int size,int covariates,int observed,int groups,int group,
                  rowvec beta,rowvec ratio,rowvec allocation,
                  double wo, vec wm, double ws, vec wm_ps, mat real_data, List subset, List value){

  int length = subset.length();

  cube res(4,length,n);


  for(int i = 0; i < n; i++){
    mat data = rcpp_get_data(size,ratio,groups,beta,real_data,covariates,observed);
    mat data_cr = rcpp_cr(size,covariates,observed,groups,ratio,data);
    mat data_str = rcpp_str(size,covariates,observed,groups,ratio,data);
    mat data_ps = rcpp_car(size,covariates,observed,groups,ratio,allocation,0,wm_ps,0,data);
    mat data_car = rcpp_car(size,covariates,observed,groups,ratio,allocation,wo,wm,ws,data);


    for(int j = 0; j < length; j++){

      res(0,j,i) = rcpp_get_imb(rcpp_subset_all(data_cr,subset[j],value[j]),covariates,group,ratio)/sqrt(size);
      res(1,j,i) = rcpp_get_imb(rcpp_subset_all(data_str,subset[j],value[j]),covariates,group,ratio)/sqrt(size);
      res(2,j,i) = rcpp_get_imb(rcpp_subset_all(data_ps,subset[j],value[j]),covariates,group,ratio)/sqrt(size);
      res(3,j,i) = rcpp_get_imb(rcpp_subset_all(data_car,subset[j],value[j]),covariates,group,ratio)/sqrt(size);
    }
  }
  cube final(4,length,2);

  final(span(0,3),span(0,length-1),span(0,0)) = mean(res,2);
  final(span(0,3),span(0,length-1),span(1,1)) = sd(res,n);


  return final;
}






