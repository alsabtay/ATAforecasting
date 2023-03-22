#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
double meanIT(NumericVector x, int t){
  double total = 0.0;
  for(int i = 0; i < t; ++i)
    total += x[i];
  return 1.0 * total / t;
}

// [[Rcpp::export]]
double medianIT(NumericVector x, int t){
  NumericVector y = x[Rcpp::Range(0, t)];
  int size = y.size();
  std::sort(y.begin(), y.begin()+t);
  if (size % 2 != 0)
    return (double)y[size/2];
  return (double)(y[(size-1)/2] + y[size/2])/2.0;
}

// [[Rcpp::export]]
double calc_amse(NumericMatrix IAX){
  const int nrw = IAX.nrow();
  const int ncl = IAX.ncol();
  int i, j;
  double denom;
  double accmeasure = 0.0;
  NumericVector amse(nrw-1);

  for(i = 1; i < nrw; i++)
  {
    amse[i-1] = 0.0;
    denom = 0.0;
    for(j = 0; j < ncl; j++)
    {
      if (!NumericMatrix::is_na(IAX(i,j)))
      {
        denom += 1.0;
        amse[i-1] = ((amse[i-1] * (denom - 1.0)) + IAX(i,j)) / denom;
      }
    }
  }
  accmeasure = mean(amse);
  return accmeasure;
}

// [[Rcpp::export]]
double NaiveSD_Accry(NumericVector train_set, int frqx, int accry){
  int n = train_set.size();
  int j, fn;
  double accmeasure=0.0;
  NumericVector fitsn(n);
  NumericVector ITAcc(n);
  NumericVector ITErr(n);
  NumericVector ITsmape(n);
  NumericVector pe(n);
  fn = frqx - 1;
  for(j = 0; j < n; j++){
    if(fn < j)
      fitsn[j] = train_set[j - frqx];
  }
  ITErr = train_set - fitsn;
  pe = (ITErr / fitsn) * 100.0;
  ITsmape = (abs(ITErr) / (abs(train_set) + abs(fitsn))) * 200.0;
  if ((accry==1) | (accry==2) | (accry==12) | (accry==13))
    ITAcc = abs(ITErr);
  else if ((accry==3) | (accry==4) | (accry==11) | (accry==15))
    ITAcc = pow(ITErr, 2.0);
  else if ((accry==5) | (accry==6) )
    ITAcc = pe;
  else if ((accry==7) | (accry==8) )
    ITAcc = abs(pe);
  else if ((accry==9) | (accry==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(n-1,NumericVector::get_na());
  if ((accry==1) | (accry==3) | (accry==5) | (accry==7) | (accry==9) )
    accmeasure = mean(ITAcc);
  else if ((accry==2) | (accry==4) | (accry==6) | (accry==8) | (accry==10) )
    accmeasure = median(ITAcc);
  else if (accry==11)
    accmeasure = sqrt(mean(ITAcc));
  else if (accry==12)
    accmeasure = 1.0;
  else if (accry==13)
    accmeasure = 1.0 * (mean(ITsmape) + mean(ITAcc)) / 2;
  else if (accry==15)
    accmeasure = ITAcc.size() * log(sum(ITAcc));
  else if (accry==16)
    accmeasure = var(ITErr);
  else
    accmeasure = NA_REAL;
  return accmeasure;
}

// [[Rcpp::export]]
double NaiveSV_Accry(NumericVector train_set, NumericVector frqx, int accry ){
  int n = train_set.size();
  int p = frqx.length();
  int z, j, fn;
  double accmeasure=0.0;
  NumericVector fitsn(n);
  NumericVector deseas(n);
  NumericVector tfitsn(n);
  NumericVector ITAcc(n);
  NumericVector ITErr(n);
  NumericVector ITsmape(n);
  NumericVector pe(n);
  deseas = train_set;
  for(z = 0; z < p; z++){
    fn = frqx[z] - 1;
    for(j = 0; j < n; j++){
      if(fn < j)
        fitsn[j] = deseas[j - frqx[z]];
      tfitsn[j] = tfitsn[j] + fitsn[j];
    }
    deseas = deseas - fitsn;
  }
  ITErr = train_set - tfitsn;
  pe = (ITErr / tfitsn) * 100.0;
  ITsmape = (abs(ITErr) / (abs(train_set) + abs(tfitsn))) * 200.0;
  if ((accry==1) | (accry==2) | (accry==12) | (accry==13))
    ITAcc = abs(ITErr);
  else if ((accry==3) | (accry==4) | (accry==11) | (accry==15))
    ITAcc = pow(ITErr, 2.0);
  else if ((accry==5) | (accry==6) )
    ITAcc = pe;
  else if ((accry==7) | (accry==8) )
    ITAcc = abs(pe);
  else if ((accry==9) | (accry==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(n-1,NumericVector::get_na());
  if ((accry==1) | (accry==3) | (accry==5) | (accry==7) | (accry==9) )
    accmeasure = mean(ITAcc);
  else if ((accry==2) | (accry==4) | (accry==6) | (accry==8) | (accry==10) )
    accmeasure = median(ITAcc);
  else if (accry==11)
    accmeasure = sqrt(mean(ITAcc));
  else if (accry==12)
    accmeasure = 1.0;
  else if (accry==13)
    accmeasure = 1.0 * (mean(ITsmape) + mean(ITAcc)) / 2;
  else if (accry==15)
    accmeasure = ITAcc.size() * log(sum(ITAcc));
  else if (accry==16)
    accmeasure = var(ITErr);
  else
    accmeasure = NA_REAL;
  return accmeasure;
}

// [[Rcpp::export]]
double NaiveSD_Accry_hin(NumericVector train_set, double frqx, int accry, int h){
  int n = train_set.size();
  int j, fn;
  double accmeasure=0.0;
  NumericVector fitsn(n);
  NumericVector deseas(n);
  NumericVector ITAcc(n);
  NumericVector ITErr(n);
  NumericVector ITsmape(n);
  NumericVector pe(n);
  NumericVector hITAcc(h);
  NumericVector hITErr(h);
  NumericVector hITsmape(h);
  fn = frqx - 1;
  for(j = 0; j < n; j++){
    if(fn < j)
      fitsn[j] = train_set[j - frqx];
  }
  ITErr = train_set - fitsn;
  pe = (ITErr / fitsn) * 100.0;
  ITsmape = (abs(ITErr) / (abs(train_set) + abs(fitsn))) * 200.0;
  if ((accry==1) | (accry==2) | (accry==12) | (accry==13))
    ITAcc = abs(ITErr);
  else if ((accry==3) | (accry==4) | (accry==11) | (accry==15))
    ITAcc = pow(ITErr, 2.0);
  else if ((accry==5) | (accry==6) )
    ITAcc = pe;
  else if ((accry==7) | (accry==8) )
    ITAcc = abs(pe);
  else if ((accry==9) | (accry==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(n-1,NumericVector::get_na());
  hITAcc = head(ITAcc, (n - h));
  hITErr = head(ITErr, (n - h));
  hITsmape = head(ITsmape, (n - h));
  if ((accry==1) | (accry==3) | (accry==5) | (accry==7) | (accry==9) )
    accmeasure = mean(hITAcc);
  else if ((accry==2) | (accry==4) | (accry==6) | (accry==8) | (accry==10) )
    accmeasure = median(hITAcc);
  else if (accry==11)
    accmeasure = sqrt(mean(hITAcc));
  else if (accry==12)
    accmeasure = 1.0;
  else if (accry==13)
    accmeasure = 1.0 * (mean(hITsmape) + mean(hITAcc)) / 2;
  else if (accry==15)
    accmeasure = hITAcc.size() * log(sum(hITAcc));
  else if (accry==16)
    accmeasure = var(hITErr);
  else
    accmeasure = NA_REAL;
  return accmeasure;
}


// [[Rcpp::export]]
double NaiveSV_Accry_hin(NumericVector train_set, NumericVector frqx, int accry, int h){
  int n = train_set.size();
  int p = frqx.length();
  int z, j, fn;
  double accmeasure=0.0;
  NumericVector fitsn(n);
  NumericVector deseas(n);
  NumericVector tfitsn(n);
  NumericVector ITAcc(n);
  NumericVector ITErr(n);
  NumericVector ITsmape(n);
  NumericVector pe(n);
  NumericVector hITAcc(h);
  NumericVector hITErr(h);
  NumericVector hITsmape(h);
  deseas = train_set;
  for(z = 0; z < p; z++){
    fn = frqx[z] - 1;
    for(j = 0; j < n; j++){
      if(fn < j)
        fitsn[j] = deseas[j - frqx[z]];
      tfitsn[j] = tfitsn[j] + fitsn[j];
    }
    deseas = deseas - fitsn;
  }
  ITErr = train_set - tfitsn;
  pe = (ITErr / tfitsn) * 100.0;
  ITsmape = (abs(ITErr) / (abs(train_set) + abs(tfitsn))) * 200.0;
  if ((accry==1) | (accry==2) | (accry==12) | (accry==13))
    ITAcc = abs(ITErr);
  else if ((accry==3) | (accry==4) | (accry==11) | (accry==15))
    ITAcc = pow(ITErr, 2.0);
  else if ((accry==5) | (accry==6) )
    ITAcc = pe;
  else if ((accry==7) | (accry==8) )
    ITAcc = abs(pe);
  else if ((accry==9) | (accry==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(n-1,NumericVector::get_na());
  hITAcc = head(ITAcc, (n - h));
  hITErr = head(ITErr, (n - h));
  hITsmape = head(ITsmape, (n - h));
  if ((accry==1) | (accry==3) | (accry==5) | (accry==7) | (accry==9) )
    accmeasure = mean(hITAcc);
  else if ((accry==2) | (accry==4) | (accry==6) | (accry==8) | (accry==10) )
    accmeasure = median(hITAcc);
  else if (accry==11)
    accmeasure = sqrt(mean(hITAcc));
  else if (accry==12)
    accmeasure = 1.0;
  else if (accry==13)
    accmeasure = 1.0 * (mean(hITsmape) + mean(hITAcc)) / 2;
  else if (accry==15)
    accmeasure = hITAcc.size() * log(sum(hITAcc));
  else if (accry==16)
    accmeasure = var(hITErr);
  else
    accmeasure = NA_REAL;
  return accmeasure;
}

// [[Rcpp::export]]
double SubATACore(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, NumericVector IZFRQ, int IZNMSE) {
  int LENZ = IAZ.size();
  int f = IZFRQ.length();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, S, T, S_1, T_1, T_0;
  int i, indx, h;
  NumericVector ITErr(LENZ);
  NumericVector ITAcc(LENZ);
  double accmeasure=0.0;
  double phiTotal=0.0;
  NumericVector pe(LENZ);
  NumericVector ITsmape(LENZ);
  NumericMatrix FC(LENZ, IZNMSE);

  if (IZMO==1)
    IZT_0 = IZTA_0;
  else
    IZT_0 = IZTM_0;

  for(indx = 0; indx < LENZ-1; indx++) {
    i = indx + 1;

    if (indx==0)
    {
      Xlag = IAZ[indx];
      Xobs = IAZ[indx];
    }
    else
    {
      if ((IZIL==1) & (i<=IZP))
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else if ((IZIL==2) & (i<=IZP))
      {
        Xlag = medianIT(IAZ,indx-1);
        Xobs = medianIT(IAZ,indx);
      }
      else
      {
        if ((IZIL==1) & (indx<=IZP))
          Xlag = meanIT(IAZ,indx-1);
        else if ((IZIL==2) & (indx<=IZP))
          Xlag = medianIT(IAZ,indx-1);
        else
          Xlag = IAZ[indx-1];
        Xobs = IAZ[indx];
      }
    }

    if (IZMO==1)
      T_0 = 1.0 * Xobs - Xlag;
    else
      T_0 = 1.0 * Xobs / Xlag;

    if (i == 1) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        T = 0.0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        T = 1.0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i<=IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i<=IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i>IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else {
      IZFIT[i] = NA_REAL;
      S_1 = NA_REAL;
      T_1 = NA_REAL;
      for(h = 0; h < IZNMSE; h++)
        FC(i,h) = NA_REAL;
    }
  }

  ITErr = IAZ - IZFIT;
  pe = (ITErr / IAZ) * 100.0;
  ITsmape = (abs(ITErr) / (abs(IZFIT) + abs(IAZ))) * 200.0;
  if ((IZAC==1) | (IZAC==2) | (IZAC==12) | (IZAC==13))
    ITAcc = abs(ITErr);
  else if ((IZAC==3) | (IZAC==4) | (IZAC==11) | (IZAC==15))
    ITAcc = pow(ITErr, 2.0);
  else if ((IZAC==5) | (IZAC==6) )
    ITAcc = pe;
  else if ((IZAC==7) | (IZAC==8) )
    ITAcc = abs(pe);
  else if ((IZAC==9) | (IZAC==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(LENZ-1,NumericVector::get_na());
  if ((IZAC==1) | (IZAC==3) | (IZAC==5) | (IZAC==7) | (IZAC==9) )
    accmeasure = mean(ITAcc);
  else if ((IZAC==2) | (IZAC==4) | (IZAC==6) | (IZAC==8) | (IZAC==10) )
    accmeasure = median(ITAcc);
  else if (IZAC==11)
    accmeasure = sqrt(mean(ITAcc));
  else if (IZAC==12){
    if(f > 1)
      accmeasure = mean(ITAcc) / NaiveSV_Accry(IAZ, IZFRQ, 1);
    else
      accmeasure = mean(ITAcc) / NaiveSD_Accry(IAZ, IZFRQ[0], 1);
  }
  else if (IZAC==13){
    if(f > 1)
      accmeasure = 1.0 * ((mean(ITsmape) / NaiveSV_Accry(IAZ, IZFRQ, 9)) + (mean(ITAcc) / NaiveSV_Accry(IAZ, IZFRQ, 1))) / 2;
    else
      accmeasure = 1.0 * ((mean(ITsmape) / NaiveSD_Accry(IAZ, IZFRQ[0], 9)) + (mean(ITAcc) / NaiveSD_Accry(IAZ, IZFRQ[0], 1))) / 2;
  }
  else if ((IZAC==14) | (IZAC==17))
    accmeasure = calc_amse(FC);
  else if (IZAC==15)
    accmeasure = ITAcc.size() * log(sum(ITAcc));
  else if (IZAC==16)
    accmeasure = var(ITErr);
  else
    accmeasure = NA_REAL;
  return accmeasure;
}

// [[Rcpp::export]]
NumericVector SubATADamped(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0, NumericVector IXFRQ, int IXNMSE){
  int LENX = IAX.size();
  int  	d_opt_p;
  int  	d_opt_q;
  double  d_opt_phi;
  NumericVector out(4);
  long double optAccryStart;
  long double optAccryEnd;
  int m, i, j, mstart, mfinish ;
  double k;

  d_opt_phi = 1.0;
  optAccryStart = 9999999999999.9;
  optAccryEnd = 9999999999999.9;

  if (IXMO==0) {
    mstart = 1;
    mfinish = 2;
  }
  else if (IXMO==1) {
    mstart =1;
    mfinish = 1;
  }
  else {
    mstart =2;
    mfinish = 2;
  }

  if (IXLF==1) {
    d_opt_q = 0;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = SubATACore(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 0;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
        for(j = 0; j <= d_opt_p; j++) {
          optAccryEnd = SubATACore(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_q = j;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  else if (IXTF==1) {
    d_opt_q = 1;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = SubATACore(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 1;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  else if (IXTS==1) {
    d_opt_q = 1;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(j = 1; j <= LENX; j++) {
          for(i=j; ((j<=i) & (i<=LENX)); i++) {
            optAccryEnd = SubATACore(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_p = i;
              d_opt_q = j;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
  }
  else {
    if ( (IXP==-1) & (IXQ==-1) ) {
      d_opt_q = 0;
      d_opt_p = 1;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            for(j=0; j<=i; j++) {
              optAccryEnd = SubATACore(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
              if (optAccryEnd <= optAccryStart) {
                d_opt_phi = 1.0 * k;
                d_opt_p = i;
                d_opt_q = j;
                IXMO = m;
                optAccryStart = optAccryEnd;
              }
            }
          }
        }
      }
    }
    else if ( (IXP==-1) & (IXQ!=-1) ) {
      d_opt_q = IXQ;
      d_opt_p = 1;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            optAccryEnd = SubATACore(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_p = i;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
    else if ( (IXP!=-1) & (IXQ==-1) ) {
      d_opt_q = 0;
      d_opt_p = IXP;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(j = 0; j <= IXP; j++) {
            optAccryEnd = SubATACore(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_q = j;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
    else {
      d_opt_q = IXQ;
      d_opt_p = IXP;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          optAccryEnd = SubATACore(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  out[0] = d_opt_p;
  out[1] = d_opt_q;
  out[2] = d_opt_phi;
  out[3] = IXMO;
  return out;
}


// [[Rcpp::export]]
NumericVector SubATA(arma::mat IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, arma::mat IXTA_0, arma::mat IXTM_0, NumericVector IXSMO, NumericVector IXST, int max_smo, int max_st, int IXFRQ, int IXNMSE){
  int  d_opt_p;
  int  d_opt_q;
  double  d_opt_phi;
  int d_opt_mo;
  int d_opt_smo;
  int d_opt_st;
  int d_opt_clmn;
  int sm, st, LastIXSMO, LastIXST, mod_clmn;
  NumericVector out(7);
  NumericVector output;
  long double optAccryStart;
  long double optAccryEnd;

  d_opt_phi = 1.0;
  optAccryStart = 9999999999999.9;
  optAccryEnd = 9999999999999.9;

  for(sm = 1; sm <= max_smo; sm++) {
    for (st = 1; st <= max_st; st++) {
      LastIXSMO = IXSMO[sm-1];
      LastIXST = IXST[st-1];
      mod_clmn = (sm*max_st)-(st%max_st);
      NumericVector subIAX = wrap(IAX.col(mod_clmn-1));
      NumericVector subIXTA_0 = wrap(IXTA_0.col(mod_clmn-1));
      NumericVector subIXTM_0 = wrap(IXTM_0.col(mod_clmn-1));
      output = SubATADamped(subIAX, IXP, IXQ, IXMO, IXAC, IXLF, IXTF, IXTS, IXPHIS, IXPHIE, IXPHISS, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IXNMSE);
      optAccryEnd = SubATACore(subIAX, output[0], output[1], output[2], output[3], IXAC, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IXNMSE);
      if (optAccryEnd <= optAccryStart){
        d_opt_p = output[0];
        d_opt_q = output[1];
        d_opt_phi = output[2];
        d_opt_mo = output[3];
        d_opt_smo = LastIXSMO;
        d_opt_st = LastIXST;
        d_opt_clmn = mod_clmn;
        optAccryStart = optAccryEnd;
      }
    }
  }
  out[0] = d_opt_p;
  out[1] = d_opt_q;
  out[2] = d_opt_phi;
  out[3] = d_opt_mo;
  out[4] = d_opt_smo;
  out[5] = d_opt_st;
  out[6] = d_opt_clmn;
  return out;
}


// [[Rcpp::export]]
double SubATACoreHoldout(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, NumericVector IZFRQ, NumericVector IAZout, int onestep) {
  int LENZ = IAZ.size();
  int LENH = IAZout.size();
  int f = IZFRQ.length();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, phiTotal, S, T, S_1, T_1, T_0;
  int i, indx, h;
  NumericVector IZFRCST(LENH);
  NumericVector ITFrcstErr(LENH);
  NumericVector peOUT(LENH);
  NumericVector ITsmapeOUT(LENH);
  NumericVector ITAccOUT(LENH);
  NumericVector IAZonestep(LENZ+LENH);
  IAZonestep[Range(0,LENZ-1)] = IAZ;
  IAZonestep[Range(LENZ, LENZ+LENH-1)] = IAZout;

  double accmeasureOUT = 0;

  if (IZMO==1)
    IZT_0 = IZTA_0;
  else
    IZT_0 = IZTM_0;

  for(indx = 0; indx < LENZ-1; indx++) {
    i = indx + 1;

    if (indx==0)
    {
      Xlag = IAZ[indx];
      Xobs = IAZ[indx];
    }
    else
    {
      if ((IZIL==1) & (i<=IZP))
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else if ((IZIL==2) & (i<=IZP))
      {
        Xlag = medianIT(IAZ,indx-1);
        Xobs = medianIT(IAZ,indx);
      }
      else
      {
        if ((IZIL==1) & (indx<=IZP))
          Xlag = meanIT(IAZ,indx-1);
        else if ((IZIL==2) & (indx<=IZP))
          Xlag = medianIT(IAZ,indx-1);
        else
          Xlag = IAZ[indx-1];
        Xobs = IAZ[indx];
      }
    }

    if (IZMO==1)
      T_0 = 1.0 * Xobs - Xlag;
    else
      T_0 = 1.0 * Xobs / Xlag;

    if (i == 1) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        T = 0.0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        T = 1.0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i<=IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i<=IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else {
      IZFIT[i] = NA_REAL;
      S_1 = NA_REAL;
      T_1 = NA_REAL;
    }
  }

  if (onestep==0){
    Xobs = IAZ[LENZ-1];
    if (IZMO==1) {
      coefpk = 1.0 * IZP / LENZ;
      coefqk = 1.0 * IZQ / LENZ;
      S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
      T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
      IZFRCST[0] = S + (IZPHI * T);
      phiTotal = IZPHI;
      for(h = 1; h < LENH; h++) {
        phiTotal += pow(IZPHI, h);
        IZFRCST[h] = S + (phiTotal * T);
      }
    }
    if (IZMO==2) {
      coefpk = 1.0 * IZP / LENZ;
      coefqk = 1.0 * IZQ / LENZ;
      S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
      T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
      IZFRCST[0] = S * pow(T, IZPHI);
      phiTotal = IZPHI;
      for(h = 1; h < LENH; h++) {
        phiTotal += pow(IZPHI, h);
        IZFRCST[h] = S * pow(T, phiTotal);
      }
    }
  }
  else{
    h = 0;
    for(indx = LENZ-1; indx < LENZ+LENH-1; indx++) {
      Xobs = IAZonestep[indx];
      if (IZMO==1) {
        coefpk = 1.0 * IZP / (indx+1);
        coefqk = 1.0 * IZQ / (indx+1);
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFRCST[h] = S + (IZPHI * T);
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / (indx+1);
        coefqk = 1.0 * IZQ / (indx+1);
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFRCST[h] = S * pow(T, IZPHI);
      }
      h = h + 1;
    }
  }

  ITFrcstErr = IAZout - IZFRCST;
  peOUT = (ITFrcstErr / IAZout) * 100.0;
  ITsmapeOUT = (abs(ITFrcstErr) / (abs(IZFRCST) + abs(IAZout))) * 200.0;
  if ( (IZAC==1) | (IZAC==2) | (IZAC==12) | (IZAC==13) )
    ITAccOUT = abs(ITFrcstErr);
  else if ( (IZAC==3) | (IZAC==4) | (IZAC==11) | (IZAC==15))
    ITAccOUT = pow(ITFrcstErr, 2.0);
  else if ( (IZAC==5) | (IZAC==6) )
    ITAccOUT = peOUT;
  else if ( (IZAC==7) | (IZAC==8) )
    ITAccOUT = abs(peOUT);
  else if ( (IZAC==9) | (IZAC==10) )
    ITAccOUT = ITsmapeOUT;
  else
    NumericVector ITAccOUT(LENH-1,NumericVector::get_na());
  if ( (IZAC==1) | (IZAC==3) | (IZAC==5) | (IZAC==7) | (IZAC==9) )
    accmeasureOUT = mean(ITAccOUT);
  else if ( (IZAC==2) | (IZAC==4) | (IZAC==6) | (IZAC==8) | (IZAC==10) )
    accmeasureOUT = median(ITAccOUT);
  else if (IZAC==11)
    accmeasureOUT = sqrt(mean(ITAccOUT));
  else if (IZAC==12){
    if(f > 1)
      accmeasureOUT = mean(ITAccOUT) / NaiveSV_Accry(IAZ, IZFRQ, 1);
    else
      accmeasureOUT = mean(ITAccOUT) / NaiveSD_Accry(IAZ, IZFRQ[0], 1);
  }
  else if (IZAC==13){
    if(f > 1)
      accmeasureOUT = 1.0 * ((mean(ITsmapeOUT) / NaiveSV_Accry(IAZ, IZFRQ, 9)) + (mean(ITAccOUT) / NaiveSV_Accry(IAZ, IZFRQ, 1))) / 2;
    else
      accmeasureOUT = 1.0 * ((mean(ITsmapeOUT) / NaiveSD_Accry(IAZ, IZFRQ[0], 9)) + (mean(ITAccOUT) / NaiveSD_Accry(IAZ, IZFRQ[0], 1))) / 2;
  }
  else if (IZAC==15)
    accmeasureOUT = ITAccOUT.size() * log(sum(ITAccOUT));
  else if (IZAC==16)
    accmeasureOUT = var(ITFrcstErr);
  else
    accmeasureOUT = NA_REAL;
  return accmeasureOUT;
}

// [[Rcpp::export]]
NumericVector SubATADampedHoldout(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0, NumericVector IXFRQ, NumericVector IAXout, int onestep){
  int LENX = IAX.size();
  int  	d_opt_p;
  int  	d_opt_q;
  double  d_opt_phi;
  NumericVector out(5);
  long double optAccryStart;
  long double optAccryEnd;
  int m, i, j, mstart, mfinish ;
  double k;

  d_opt_phi = 1.0;
  optAccryStart = 9999999999999.9;
  optAccryEnd = 9999999999999.9;

  if (IXMO==0) {
    mstart = 1;
    mfinish = 2;
  }
  else if (IXMO==1) {
    mstart =1;
    mfinish = 1;
  }
  else {
    mstart =2;
    mfinish = 2;
  }

  if (IXLF==1) {
    d_opt_q = 0;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = SubATACoreHoldout(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 0;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
        for(j = 0; j <= d_opt_p; j++) {
          optAccryEnd = SubATACoreHoldout(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_q = j;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  else if (IXTF==1) {
    d_opt_q = 1;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = SubATACoreHoldout(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 1;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  else if (IXTS==1) {
    d_opt_q = 1;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(j = 1; j <= LENX; j++) {
          for(i=1; ((j<=i) & (i<=LENX)); i++) {
            optAccryEnd = SubATACoreHoldout(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_p = i;
              d_opt_q = j;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
  }else {
    if ( (IXP==-1) & (IXQ==-1) ) {
      d_opt_q = 0;
      d_opt_p = 1;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            for(j=0; j<=i; j++) {
              optAccryEnd = SubATACoreHoldout(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
              if (optAccryEnd <= optAccryStart) {
                d_opt_phi = 1.0 * k;
                d_opt_p = i;
                d_opt_q = j;
                IXMO = m;
                optAccryStart = optAccryEnd;
              }
            }
          }
        }
      }
    }
    else if ( (IXP==-1) & (IXQ!=-1) ) {
      d_opt_q = IXQ;
      d_opt_p = 1;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            optAccryEnd = SubATACoreHoldout(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_p = i;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
    else if ( (IXP!=-1) & (IXQ==-1) ) {
      d_opt_q = 0;
      d_opt_p = IXP;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(j = 0; j <= IXP; j++) {
            optAccryEnd = SubATACoreHoldout(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_q = j;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
    else {
      d_opt_q = IXQ;
      d_opt_p = IXP;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          optAccryEnd = SubATACoreHoldout(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout, onestep);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  out[0] = d_opt_p;
  out[1] = d_opt_q;
  out[2] = d_opt_phi;
  out[3] = IXMO;
  out[4] = optAccryStart;
  return out;
}


// [[Rcpp::export]]
NumericVector SubATAHoldout(arma::mat IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS,
                            int IXIL, int IXIT, arma::mat IXTA_0, arma::mat IXTM_0, NumericVector IXSMO, NumericVector IXST, int max_smo, int max_st,
                            NumericVector IXFRQ, arma::mat IAXout, int onestep){
  int  d_opt_p;
  int  d_opt_q;
  double  d_opt_phi;
  int d_opt_mo;
  int d_opt_smo;
  int d_opt_st;
  int d_opt_clmn;
  int sm, st, LastIXSMO, LastIXST, mod_clmn;
  NumericVector out(8);
  NumericVector output;
  long double optAccryStart;
  long double optAccryEnd;

  d_opt_phi = 1.0;
  optAccryStart = 9999999999999.9;
  optAccryEnd = 9999999999999.9;

  for(sm = 1; sm <= max_smo; sm++) {
    for (st = 1; st <= max_st; st++) {
      LastIXSMO = IXSMO[sm-1];
      LastIXST = IXST[st-1];
      mod_clmn = (sm*max_st)-(st%max_st);
      NumericVector subIAX = wrap(IAX.col(mod_clmn-1));
      NumericVector subIAXout = wrap(IAXout.col(mod_clmn-1));
      NumericVector subIXTA_0 = wrap(IXTA_0.col(mod_clmn-1));
      NumericVector subIXTM_0 = wrap(IXTM_0.col(mod_clmn-1));
      output = SubATADampedHoldout(subIAX, IXP, IXQ, IXMO, IXAC, IXLF, IXTF, IXTS, IXPHIS, IXPHIE, IXPHISS, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, subIAXout, onestep);
      optAccryEnd = SubATACoreHoldout(subIAX, output[0], output[1], output[2], output[3], IXAC, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, subIAXout, onestep);
      if (optAccryEnd <= optAccryStart){
        d_opt_p = output[0];
        d_opt_q = output[1];
        d_opt_phi = output[2];
        d_opt_mo = output[3];
        d_opt_smo = LastIXSMO;
        d_opt_st = LastIXST;
        d_opt_clmn = mod_clmn;
        optAccryStart = optAccryEnd;
      }
    }
  }
  out[0] = d_opt_p;
  out[1] = d_opt_q;
  out[2] = d_opt_phi;
  out[3] = d_opt_mo;
  out[4] = d_opt_smo;
  out[5] = d_opt_st;
  out[6] = d_opt_clmn;
  out[7] = optAccryStart;
  return out;
}


// [[Rcpp::export]]
NumericVector ATAHoldoutForecast(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, NumericVector IZFRQ, NumericVector IAZout, int onestep) {
  int LENZ = IAZ.size();
  int LENH = IAZout.size();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, phiTotal, S, T, S_1, T_1, T_0;
  int i, indx, h;
  NumericVector IZFRCST(LENH);
  NumericVector IAZonestep(LENZ+LENH);
  IAZonestep[Range(0,LENZ-1)] = IAZ;
  IAZonestep[Range(LENZ, LENZ+LENH-1)] = IAZout;

  if (IZMO==1)
    IZT_0 = IZTA_0;
  else
    IZT_0 = IZTM_0;

  for(indx = 0; indx < LENZ-1; indx++) {
    i = indx + 1;

    if (indx==0)
    {
      Xlag = IAZ[indx];
      Xobs = IAZ[indx];
    }
    else
    {
      if ((IZIL==1) & (i<=IZP))
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else if ((IZIL==2) & (i<=IZP))
      {
        Xlag = medianIT(IAZ,indx-1);
        Xobs = medianIT(IAZ,indx);
      }
      else
      {
        if ((IZIL==1) & (indx<=IZP))
          Xlag = meanIT(IAZ,indx-1);
        else if ((IZIL==2) & (indx<=IZP))
          Xlag = medianIT(IAZ,indx-1);
        else
          Xlag = IAZ[indx-1];
        Xobs = IAZ[indx];
      }
    }

    if (IZMO==1)
      T_0 = 1.0 * Xobs - Xlag;
    else
      T_0 = 1.0 * Xobs / Xlag;

    if (i == 1) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        T = 0.0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        T = 1.0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i<=IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i<=IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else {
      IZFIT[i] = NA_REAL;
      S_1 = NA_REAL;
      T_1 = NA_REAL;
    }
  }

  if (onestep==0){
    Xobs = IAZ[LENZ-1];
    if (IZMO==1) {
      coefpk = 1.0 * IZP / LENZ;
      coefqk = 1.0 * IZQ / LENZ;
      S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
      T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
      IZFRCST[0] = S + (IZPHI * T);
      phiTotal = IZPHI;
      for(h = 1; h < LENH; h++) {
        phiTotal += pow(IZPHI, h);
        IZFRCST[h] = S + (phiTotal * T);
      }
    }
    if (IZMO==2) {
      coefpk = 1.0 * IZP / LENZ;
      coefqk = 1.0 * IZQ / LENZ;
      S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
      T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
      IZFRCST[0] = S * pow(T, IZPHI);
      phiTotal = IZPHI;
      for(h = 1; h < LENH; h++) {
        phiTotal += pow(IZPHI, h);
        IZFRCST[h] = S * pow(T, phiTotal);
      }
    }
  }
  else{
    h = 0;
    for(indx = LENZ-1; indx < LENZ+LENH-1; indx++) {
      Xobs = IAZonestep[indx];
      if (IZMO==1) {
        coefpk = 1.0 * IZP / (indx+1);
        coefqk = 1.0 * IZQ / (indx+1);
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFRCST[h] = S + (IZPHI * T);
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / (indx+1);
        coefqk = 1.0 * IZQ / (indx+1);
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFRCST[h] = S * pow(T, IZPHI);
      }
      h = h + 1;
    }
  }
  return IZFRCST;
}


// [[Rcpp::export]]
double SubATACoreHoldhin(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, NumericVector IZFRQ, int IZH, int IZNMSE) {
  int LENZ = IAZ.size();
  int f = IZFRQ.length();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, S, T, S_1, T_1, T_0;
  int i, indx, h;
  NumericVector ITErr(LENZ);
  NumericVector ITAcc(LENZ);
  double accmeasure=0.0;
  double phiTotal=0.0;
  NumericVector pe(LENZ);
  NumericVector ITsmape(LENZ);
  NumericMatrix FC(LENZ, IZNMSE);
  NumericVector hITAcc(IZH);
  NumericVector hITErr(IZH);
  NumericVector hITsmape(IZH);
  NumericMatrix hFC(IZH, IZNMSE);

  if (IZMO==1)
    IZT_0 = IZTA_0;
  else
    IZT_0 = IZTM_0;

  for(indx = 0; indx < LENZ-1; indx++) {
    i = indx + 1;

    if (indx==0)
    {
      Xlag = IAZ[indx];
      Xobs = IAZ[indx];
    }
    else
    {
      if ((IZIL==1) & (i<=IZP))
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else if ((IZIL==2) & (i<=IZP))
      {
        Xlag = medianIT(IAZ,indx-1);
        Xobs = medianIT(IAZ,indx);
      }
      else
      {
        if ((IZIL==1) & (indx<=IZP))
          Xlag = meanIT(IAZ,indx-1);
        else if ((IZIL==2) & (indx<=IZP))
          Xlag = medianIT(IAZ,indx-1);
        else
          Xlag = IAZ[indx-1];
        Xobs = IAZ[indx];
      }
    }

    if (IZMO==1)
      T_0 = 1.0 * Xobs - Xlag;
    else
      T_0 = 1.0 * Xobs / Xlag;

    if (i == 1) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        T = 0.0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        T = 1.0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i<=IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        S = 1.0 * Xobs;
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i<=IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        coefqk = 1.0 * IZQ / i;
        S = 1.0 * Xobs;
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i>IZP) & (i<=IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
        else if (IZIT==2)
          T = medianIT(IZT_0,indx);
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
        IZFIT[i] = S + (IZPHI * T);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S + (phiTotal * T)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
      if (IZMO==2) {
        coefpk = 1.0 * IZP / i;
        coefqk = 1.0 * IZQ / i;
        S = coefpk * Xobs + (1-coefpk) * S_1 * pow(T_1, IZPHI);
        T = coefqk * (1.0 * S / S_1) + (1-coefqk) * pow(T_1, IZPHI);
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
        if (IZAC==14) {
          FC(i,0) = pow(IAZ[i] - IZFIT[i],2);
          phiTotal = IZPHI;
          for(h = 1; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              phiTotal += pow(IZPHI, (h+1));
              FC(i,h) = pow(IAZ[i+h] - (S * pow(T, phiTotal)),2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
        if (IZAC==17) {
          for(h = 0; h < IZNMSE; h++) {
            if((i+h) < LENZ){
              FC(i,h) = pow(IAZ[i+h] - IZFIT[i],2);
            }else{
              FC(i,h) = NA_REAL;
            }
          }
        }
      }
    }
    else {
      IZFIT[i] = NA_REAL;
      S_1 = NA_REAL;
      T_1 = NA_REAL;
      for(h = 0; h < IZNMSE; h++)
        FC(i,h) = NA_REAL;
    }
  }
  ITErr = IAZ - IZFIT;
  pe = (ITErr / IAZ) * 100.0;
  ITsmape = (abs(ITErr) / (abs(IZFIT) + abs(IAZ))) * 200.0;
  if ( (IZAC==1) | (IZAC==2) | (IZAC==12) | (IZAC==13) )
    ITAcc = abs(ITErr);
  else if ( (IZAC==3) | (IZAC==4) | (IZAC==11) | (IZAC==15))
    ITAcc = pow(ITErr, 2.0);
  else if ( (IZAC==5) | (IZAC==6) )
    ITAcc = pe;
  else if ( (IZAC==7) | (IZAC==8) )
    ITAcc = abs(pe);
  else if ( (IZAC==9) | (IZAC==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(LENZ-1,NumericVector::get_na());
  hITAcc = tail(ITAcc, IZH);
  hITErr = tail(ITErr, IZH);
  hITsmape = tail(ITsmape, IZH);
  if ( (IZAC==1) | (IZAC==3) | (IZAC==5) | (IZAC==7) | (IZAC==9) )
    accmeasure = mean(hITAcc);
  else if ( (IZAC==2) | (IZAC==4) | (IZAC==6) | (IZAC==8) | (IZAC==10) )
    accmeasure = median(hITAcc);
  else if (IZAC==11)
    accmeasure = sqrt(mean(hITAcc));
  else if (IZAC==12){
    if(f > 1)
      accmeasure = mean(hITAcc) / NaiveSV_Accry_hin(IAZ, IZFRQ, 1, IZH);
    else
      accmeasure = mean(hITAcc) / NaiveSD_Accry_hin(IAZ, IZFRQ[0], 1, IZH);
  }
  else if (IZAC==13){
    if(f > 1)
      accmeasure = 1.0 * ((mean(hITsmape) / NaiveSV_Accry_hin(IAZ, IZFRQ, 9, IZH)) + (mean(hITAcc) / NaiveSV_Accry_hin(IAZ, IZFRQ, 1, IZH))) / 2;
    else
      accmeasure = 1.0 * ((mean(hITsmape) / NaiveSD_Accry_hin(IAZ, IZFRQ[0], 9, IZH)) + (mean(hITAcc) / NaiveSD_Accry_hin(IAZ, IZFRQ[0], 1, IZH))) / 2;
  }
  else if ((IZAC==14) | (IZAC==17)){
    hFC = FC(Range(LENZ-IZH,LENZ-1), Range(0,IZNMSE-1));
    accmeasure = calc_amse(hFC);
  }
  else if (IZAC==15)
    accmeasure = hITAcc.size() * log(sum(hITAcc));
  else if (IZAC==16)
    accmeasure = var(hITErr);
  else
    accmeasure = 0.0;
  return accmeasure;
}

// [[Rcpp::export]]
NumericVector SubATADampedHoldhin(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0, NumericVector IXFRQ, int IXH, int IXNMSE){
  int LENX = IAX.size();
  int  	d_opt_p;
  int  	d_opt_q;
  double  d_opt_phi;
  NumericVector out(4);
  long double optAccryStart;
  long double optAccryEnd;
  int m, i, j, mstart, mfinish ;
  double k;

  d_opt_phi = 1.0;
  optAccryStart = 9999999999999.9;
  optAccryEnd = 9999999999999.9;

  if (IXMO==0) {
    mstart = 1;
    mfinish = 2;
  }
  else if (IXMO==1) {
    mstart =1;
    mfinish = 1;
  }
  else {
    mstart =2;
    mfinish = 2;
  }

  if (IXLF==1) {
    d_opt_q = 0;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = SubATACoreHoldhin(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 0;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
        for(j = 0; j <= d_opt_p; j++) {
          optAccryEnd = SubATACoreHoldhin(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_q = j;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  else if (IXTF==1) {
    d_opt_q = 1;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = SubATACoreHoldhin(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 1;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  else if (IXTS==1) {
    d_opt_q = 1;
    d_opt_p = 1;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(j = 1; j <= LENX; j++) {
          for(i=1; ((j<=i) & (i<=LENX)); i++) {
            optAccryEnd = SubATACoreHoldhin(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_p = i;
              d_opt_q = j;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
  }
  else {
    if ( (IXP==-1) & (IXQ==-1) ) {
      d_opt_q = 0;
      d_opt_p = 1;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            for(j=0; j<=i; j++) {
              optAccryEnd = SubATACoreHoldhin(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
              if (optAccryEnd <= optAccryStart) {
                d_opt_phi = 1.0 * k;
                d_opt_p = i;
                d_opt_q = j;
                IXMO = m;
                optAccryStart = optAccryEnd;
              }
            }
          }
        }
      }
    }
    else if ( (IXP==-1) & (IXQ!=-1) ) {
      d_opt_q = IXQ;
      d_opt_p = 1;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            optAccryEnd = SubATACoreHoldhin(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_p = i;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
    else if ( (IXP!=-1) & (IXQ==-1) ) {
      d_opt_q = 0;
      d_opt_p = IXP;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(j = 0; j <= IXP; j++) {
            optAccryEnd = SubATACoreHoldhin(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
            if (optAccryEnd <= optAccryStart) {
              d_opt_phi = 1.0 * k;
              d_opt_q = j;
              IXMO = m;
              optAccryStart = optAccryEnd;
            }
          }
        }
      }
    }
    else {
      d_opt_q = IXQ;
      d_opt_p = IXP;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          optAccryEnd = SubATACoreHoldhin(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH, IXNMSE);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
      }
    }
  }
  out[0] = d_opt_p;
  out[1] = d_opt_q;
  out[2] = d_opt_phi;
  out[3] = IXMO;
  return out;
}


// [[Rcpp::export]]
NumericVector SubATAHoldhin(arma::mat IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, arma::mat IXTA_0, arma::mat IXTM_0, NumericVector IXSMO, NumericVector IXST, int max_smo, int max_st, NumericVector IXFRQ, int IXH, int IXNMSE){
  int  d_opt_p;
  int  d_opt_q;
  double  d_opt_phi;
  int d_opt_mo;
  int d_opt_smo;
  int d_opt_st;
  int d_opt_clmn;
  int sm, st, LastIXSMO, LastIXST, mod_clmn;
  NumericVector out(7);
  NumericVector output;
  long double optAccryStart;
  long double optAccryEnd;

  d_opt_phi = 1.0;
  optAccryStart = 9999999999999.9;
  optAccryEnd = 9999999999999.9;

  for(sm = 1; sm <= max_smo; sm++) {
    for (st = 1; st <= max_st; st++) {
      LastIXSMO = IXSMO[sm-1];
      LastIXST = IXST[st-1];
      mod_clmn = (sm*max_st)-(st%max_st);
      NumericVector subIAX = wrap(IAX.col(mod_clmn-1));
      NumericVector subIXTA_0 = wrap(IXTA_0.col(mod_clmn-1));
      NumericVector subIXTM_0 = wrap(IXTM_0.col(mod_clmn-1));
      output = SubATADampedHoldhin(subIAX, IXP, IXQ, IXMO, IXAC, IXLF, IXTF, IXTS, IXPHIS, IXPHIE, IXPHISS, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IXH, IXNMSE);
      optAccryEnd = SubATACoreHoldhin(subIAX, output[0], output[1], output[2], output[3], IXAC, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IXH, IXNMSE);
      if (optAccryEnd <= optAccryStart){
        d_opt_p = output[0];
        d_opt_q = output[1];
        d_opt_phi = output[2];
        d_opt_mo = output[3];
        d_opt_smo = LastIXSMO;
        d_opt_st = LastIXST;
        d_opt_clmn = mod_clmn;
        optAccryStart = optAccryEnd;
      }
    }
  }
  out[0] = d_opt_p;
  out[1] = d_opt_q;
  out[2] = d_opt_phi;
  out[3] = d_opt_mo;
  out[4] = d_opt_smo;
  out[5] = d_opt_st;
  out[6] = d_opt_clmn;
  return out;
}
