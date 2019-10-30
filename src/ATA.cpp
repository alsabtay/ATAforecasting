#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double NaiveSD(NumericVector insmp, int frqx){
  int n = insmp.size();
  NumericVector fNaiveSD(n);
  double sum = 0.0;
  double naiveAcc = 0.0;
  int i,j,jx,fn;

  fn = frqx - 1;
  for(j = 0; j < n; j++){
    if (j <= fn)
      fNaiveSD[j] = 0;
    else
    {
      jx = j - frqx;
      fNaiveSD[j] = insmp[jx];
    }
  }
  for(i = frqx; i < n; i++)
    sum += fabs(insmp[i] - fNaiveSD[i]);
  naiveAcc = 1.0 * sum / (n - frqx);
  return naiveAcc;
}

// [[Rcpp::export]]
double NaiveSDholdin(NumericVector insmp, int frqx, int frqh){
  int n = insmp.size();
  NumericVector fNaiveSD(n-frqh);
  double sum = 0.0;
  double naiveAcc = 0.0;
  int i,j,jx,fn;

  fn = frqx - frqh - 1;
  for(j = 0; j < (n-frqh); j++){
    if (j <= fn)
      fNaiveSD[j] = 0;
    else
    {
      jx = j - frqx;
      fNaiveSD[j] = insmp[jx];
    }
  }
  for(i = frqx; i < (n-frqh); i++)
    sum += fabs(insmp[i] - fNaiveSD[i]);
  naiveAcc = 1.0 * sum / (n - frqx - frqh);
  return naiveAcc;
}

// [[Rcpp::export]]
double outMASE(NumericVector insmp, NumericVector outsmp, NumericVector frcst, int frqx){
  int n = insmp.size();
  int m = outsmp.size();
  NumericVector fNaiveSD(n);
  double sum = 0.0;
  double pMASE, sumMASE, oMASE = 0.0;
  int i,j,jx,fn;

  fn = frqx - 1;
  for(j = 0; j < n; j++){
    if (j <= fn)
      fNaiveSD[j] = 0;
    else
    {
      jx = j - frqx;
      fNaiveSD[j] = insmp[jx];
    }
  }
  for(i = frqx; i < n; i++)
    sum += fabs(insmp[i] - fNaiveSD[i]);
  pMASE = 1.0 * sum / (n - frqx);
  for(i = 1; i < m; i++)
    sumMASE += fabs(outsmp[i]-frcst[i]);
  oMASE = 1.0 * sumMASE / (m-1);
  return 1.0 * oMASE / pMASE;
}

// [[Rcpp::export]]
double inMASE(NumericVector insmp, NumericVector fitsmp, int frqx){
  int n = insmp.size();
  int m = fitsmp.size();
  NumericVector fNaiveSD(n);
  double sum = 0.0;
  double pMASE, sumMASE, nMASE = 0.0;
  int i,j,jx,fn;

  fn = frqx - 1;
  for(j = 0; j < n; j++){
    if (j <= fn)
      fNaiveSD[j] = 0;
    else
    {
      jx = j - frqx;
      fNaiveSD[j] = insmp[jx];
    }
  }
  for(i = frqx; i < n; i++)
    sum += fabs(insmp[i]-fNaiveSD[i]);
  pMASE = 1.0 * sum / (n - frqx);
  for(i = 1; i < m; i++)
    sumMASE += fabs(insmp[i]-fitsmp[i]);
  nMASE = 1.0 * sumMASE / (m-1);
  return 1.0 * nMASE / pMASE;
}


// [[Rcpp::export]]
double inMASEholdin(NumericVector insmp, NumericVector fitsmp, int frqx, int frqh){
  int n = insmp.size();
  int m = fitsmp.size();
  NumericVector fNaiveSD(n-frqh);
  double sum = 0.0;
  double pMASE, sumMASE, nMASE = 0.0;
  int i,j,jx,fn;

  fn = frqx - frqh - 1;
  for(j = 0; j < (n-frqh); j++){
    if (j <= fn)
      fNaiveSD[j] = 0;
    else
    {
      jx = j - frqx;
      fNaiveSD[j] = insmp[jx];
    }
  }
  for(i = frqx; i < (n-frqh); i++)
    sum += fabs(insmp[i]-fNaiveSD[i]);
  pMASE = 1.0 * sum / (n - frqx - frqh);
  for(i = 1; i < (m-frqh); i++)
    sumMASE += fabs(insmp[i]-fitsmp[i]);
  nMASE = 1.0 * sumMASE / (m-1-frqh);
  return 1.0 * nMASE / pMASE;
}


// [[Rcpp::export]]
double meanIT(NumericVector x, int t){
  double total = 0.0;
  for(int i = 0; i < t; ++i)
    total += x[i];
  return 1.0 * total / t;
}

// [[Rcpp::export]]
double AutoATACore(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, int IZFRQ) {
  int LENZ = IAZ.size();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, S, T, S_1, T_1, T_0;
  int i, indx;
  NumericVector ITErr(LENZ);
  NumericVector ITAcc(LENZ);
  double accmeasure=0.0;
  NumericVector pe(LENZ);
  NumericVector ITsmape(LENZ);

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
      if (IZIL==1)
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else
      {
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
      Xobs = IAZ[indx];
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
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
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      Xobs = IAZ[indx];
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
  ITErr = IAZ - IZFIT;
  pe = (ITErr / IAZ) * 100.0;
  ITsmape = (abs(ITErr) / (abs(IZFIT) + abs(IAZ))) * 200.0;
  if ( (IZAC==1) | (IZAC==2) | (IZAC==12) | (IZAC==13) )
    ITAcc = abs(ITErr);
  else if ( (IZAC==3) | (IZAC==4) | (IZAC==11) )
    ITAcc = pow(ITErr, 2.0);
  else if ( (IZAC==5) | (IZAC==6) )
    ITAcc = pe;
  else if ( (IZAC==7) | (IZAC==8) )
    ITAcc = abs(pe);
  else if ( (IZAC==9) | (IZAC==10) )
    ITAcc = ITsmape;
  else
    NumericVector ITAcc(LENZ-1,NumericVector::get_na());
  if ( (IZAC==1) | (IZAC==3) | (IZAC==5) | (IZAC==7) | (IZAC==9) )
    accmeasure = mean(ITAcc);
  else if ( (IZAC==2) | (IZAC==4) | (IZAC==6) | (IZAC==8) | (IZAC==10) )
    accmeasure = median(ITAcc);
  else if (IZAC==11)
    accmeasure = sqrt(mean(ITAcc));
  else if (IZAC==12)
    accmeasure = inMASE(IAZ, IZFIT, IZFRQ);
  else if (IZAC==13)
    accmeasure = 1.0 * ((mean(ITsmape)/NaiveSD(IAZ, IZFRQ)) + (inMASE(IAZ, IZFIT, IZFRQ)/NaiveSD(IAZ, IZFRQ))) / 2;
  else
    accmeasure = 0;
  return accmeasure;
}

// [[Rcpp::export]]
NumericVector AutoATADamped(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0, int IXFRQ){
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = AutoATACore(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 0;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
        for(j = 0; j <= d_opt_p; j++) {
          optAccryEnd = AutoATACore(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = AutoATACore(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
            optAccryEnd = AutoATACore(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
      d_opt_p = 0;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            for(j=0; j<=i; j++) {
              optAccryEnd = AutoATACore(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
      d_opt_p = 0;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            optAccryEnd = AutoATACore(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
            optAccryEnd = AutoATACore(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
          optAccryEnd = AutoATACore(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ);
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
NumericVector AutoATA(arma::mat IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, arma::mat IXTA_0, arma::mat IXTM_0, NumericVector IXSMO, NumericVector IXST, int max_smo, int max_st, int IXFRQ){
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
      output = AutoATADamped(subIAX, IXP, IXQ, IXMO, IXAC, IXLF, IXTF, IXTS, IXPHIS, IXPHIE, IXPHISS, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ);
      optAccryEnd = AutoATACore(subIAX, output[0], output[1], output[2], output[3], IXAC, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ);
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
double AutoATACoreHoldout(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, int IZFRQ, NumericVector IAZout) {
  int LENZ = IAZ.size();
  int LENH = IAZout.size();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, phiTotal, S, T, S_1, T_1, T_0;
  int i, indx, h;
  NumericVector IZFRCST(LENH);
  NumericVector ITFrcstErr(LENH);
  NumericVector peOUT(LENH);
  NumericVector ITsmapeOUT(LENH);
  NumericVector ITAccOUT(LENH);
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
      if (IZIL==1)
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else
      {
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
      Xobs = IAZ[indx];
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
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
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      Xobs = IAZ[indx];
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

  if (IZIL==1)
    Xobs = meanIT(IAZ,LENZ-1);
  else
    Xobs = IAZ[LENZ-1];

  if (IZMO==1) {
    coefpk = 1.0 * IZP / LENZ;
    coefqk = 1.0 * IZQ / LENZ;
    S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
    T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
    IZFRCST[0] = S + (IZPHI * T);
    phiTotal = IZPHI;
    for(h = 1; h < LENH; h++) {
      phiTotal = phiTotal + pow(IZPHI, h);
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
      phiTotal = phiTotal + pow(IZPHI, h);
      IZFRCST[h] = S * pow(T, phiTotal);
    }
  }
  ITFrcstErr = IAZout - IZFRCST;
  peOUT = (ITFrcstErr / IAZout) * 100.0;
  ITsmapeOUT = (abs(ITFrcstErr) / (abs(IZFRCST) + abs(IAZout))) * 200.0;
  if ( (IZAC==1) | (IZAC==2) | (IZAC==12) | (IZAC==13) )
    ITAccOUT = abs(ITFrcstErr);
  else if ( (IZAC==3) | (IZAC==4) | (IZAC==11) )
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
  else if (IZAC==12)
    accmeasureOUT = outMASE(IAZ, IAZout, IZFRCST, IZFRQ);
  else if (IZAC==13)
    accmeasureOUT = 1.0 * ((mean(ITsmapeOUT)/NaiveSD(IAZ, IZFRQ)) + (outMASE(IAZ, IAZout, IZFRCST, IZFRQ)/NaiveSD(IAZ, IZFRQ))) / 2;
  else
    accmeasureOUT = 0;
  return accmeasureOUT;
}

// [[Rcpp::export]]
NumericVector AutoATADampedHoldout(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0, int IXFRQ, NumericVector IAXout){
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = AutoATACoreHoldout(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 0;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
        for(j = 0; j <= d_opt_p; j++) {
          optAccryEnd = AutoATACoreHoldout(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = AutoATACoreHoldout(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(j = 1; j <= LENX; j++) {
          for(i=1; ((j<=i) & (i<=LENX)); i++) {
            optAccryEnd = AutoATACoreHoldout(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
      d_opt_p = 0;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            for(j=0; j<=i; j++) {
              optAccryEnd = AutoATACoreHoldout(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
      d_opt_p = 0;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            optAccryEnd = AutoATACoreHoldout(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
            optAccryEnd = AutoATACoreHoldout(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
          optAccryEnd = AutoATACoreHoldout(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IAXout);
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
NumericVector AutoATAHoldout(arma::mat IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, arma::mat IXTA_0, arma::mat IXTM_0, NumericVector IXSMO, NumericVector IXST, int max_smo, int max_st, int IXFRQ, NumericVector IAXout){
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
      NumericVector subIXTA_0 = wrap(IXTA_0.col(mod_clmn-1));
      NumericVector subIXTM_0 = wrap(IXTM_0.col(mod_clmn-1));
      output = AutoATADampedHoldout(subIAX, IXP, IXQ, IXMO, IXAC, IXLF, IXTF, IXTS, IXPHIS, IXPHIE, IXPHISS, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IAXout);
      optAccryEnd = AutoATACoreHoldout(subIAX, output[0], output[1], output[2], output[3], IXAC, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IAXout);
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
NumericVector ATAHoldoutForecast(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, int IZFRQ, int LENH) {
  int LENZ = IAZ.size();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, phiTotal, S, T, S_1, T_1, T_0;
  int i, indx, h;
  NumericVector IZFRCST(LENH);

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
      if (IZIL==1)
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else
      {
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
      Xobs = IAZ[indx];
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
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
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      Xobs = IAZ[indx];
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
  if (IZIL==1)
    Xobs = meanIT(IAZ,LENZ-1);
  else
    Xobs = IAZ[LENZ-1];
  if (IZMO==1) {
    coefpk = 1.0 * IZP / LENZ;
    coefqk = 1.0 * IZQ / LENZ;
    S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
    T = coefqk * (S - S_1) + (1-coefqk) * (IZPHI * T_1);
    IZFRCST[0] = S + (IZPHI * T);
    phiTotal = IZPHI;
    for(h = 1; h < LENH; h++) {
      phiTotal = phiTotal + pow(IZPHI, h);
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
      phiTotal = phiTotal + pow(IZPHI, h);
      IZFRCST[h] = S * pow(T, phiTotal);
    }
  }
  return IZFRCST;
}


// [[Rcpp::export]]
double AutoATACoreHoldhin(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0, int IZFRQ, int IZH) {
  int LENZ = IAZ.size();
  NumericVector IZT_0(LENZ);
  NumericVector IZFIT(LENZ);
  double coefpk, coefqk, Xobs, Xlag, S, T, S_1, T_1, T_0;
  int i, indx;
  NumericVector ITErr(LENZ);
  NumericVector ITAcc(LENZ);
  double accmeasure=0.0;
  NumericVector pe(LENZ);
  NumericVector ITsmape(LENZ);
  NumericVector hITAcc(IZH);

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
      if (IZIL==1)
      {
        Xlag = meanIT(IAZ,indx-1);
        Xobs = meanIT(IAZ,indx);
      }
      else
      {
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
      Xobs = IAZ[indx];
      if (IZMO==1) {
        coefpk = 1.0 * IZP / i;
        S = coefpk * Xobs + (1-coefpk) * (S_1 + (IZPHI * T_1));
        if (IZIT==1)
          T = meanIT(IZT_0,indx);
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
        else
          T = T_0;
        IZFIT[i] = S * pow(T, IZPHI);
        S_1 = S;
        T_1 = T;
      }
    }
    else if ( (i>IZP) & (i>IZQ) & (IZP>=IZQ) ) {
      Xobs = IAZ[indx];
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
  ITErr = IAZ - IZFIT;
  pe = (ITErr / IAZ) * 100.0;
  ITsmape = (abs(ITErr) / (abs(IZFIT) + abs(IAZ))) * 200.0;
  if ( (IZAC==1) | (IZAC==2) | (IZAC==12) | (IZAC==13) )
    ITAcc = abs(ITErr);
  else if ( (IZAC==3) | (IZAC==4) | (IZAC==11) )
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
  if ( (IZAC==1) | (IZAC==3) | (IZAC==5) | (IZAC==7) | (IZAC==9) )
    accmeasure = mean(hITAcc);
  else if ( (IZAC==2) | (IZAC==4) | (IZAC==6) | (IZAC==8) | (IZAC==10) )
    accmeasure = median(hITAcc);
  else if (IZAC==11)
    accmeasure = sqrt(mean(hITAcc));
  else if (IZAC==12)
    accmeasure = inMASEholdin(IAZ, IZFIT, IZFRQ, IZH);
  else if (IZAC==13)
    accmeasure = 1.0 * ((mean(tail(ITsmape, IZH))/NaiveSDholdin(IAZ, IZFRQ, IZH)) + (inMASEholdin(IAZ, IZFIT, IZFRQ, IZH)/NaiveSDholdin(IAZ, IZFRQ, IZH))) / 2;
  else
    accmeasure = 0;
  return accmeasure;
}

// [[Rcpp::export]]
NumericVector AutoATADampedHoldhin(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0, int IXFRQ, int IXH){
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = AutoATACoreHoldhin(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
          if (optAccryEnd <= optAccryStart) {
            d_opt_phi = 1.0 * k;
            d_opt_p = i;
            d_opt_q = 0;
            IXMO = m;
            optAccryStart = optAccryEnd;
          }
        }
        for(j = 0; j <= d_opt_p; j++) {
          optAccryEnd = AutoATACoreHoldhin(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(i = 1; i <= LENX; i++) {
          optAccryEnd = AutoATACoreHoldhin(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
    d_opt_p = 0;
    for(m = mstart; m <= mfinish; m++) {
      for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
        for(j = 1; j <= LENX; j++) {
          for(i=1; ((j<=i) & (i<=LENX)); i++) {
            optAccryEnd = AutoATACoreHoldhin(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
      d_opt_p = 0;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            for(j=0; j<=i; j++) {
              optAccryEnd = AutoATACoreHoldhin(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
      d_opt_p = 0;
      for(m = mstart; m <= mfinish; m++) {
        for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
          for(i = 1; i <= LENX; i++) {
            optAccryEnd = AutoATACoreHoldhin(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
            optAccryEnd = AutoATACoreHoldhin(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
          optAccryEnd = AutoATACoreHoldhin(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0, IXFRQ, IXH);
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
NumericVector AutoATAHoldhin(arma::mat IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, int IXTS, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, arma::mat IXTA_0, arma::mat IXTM_0, NumericVector IXSMO, NumericVector IXST, int max_smo, int max_st, int IXFRQ, int IXH){
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
      output = AutoATADampedHoldhin(subIAX, IXP, IXQ, IXMO, IXAC, IXLF, IXTF, IXTS, IXPHIS, IXPHIE, IXPHISS, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IXH);
      optAccryEnd = AutoATACoreHoldhin(subIAX, output[0], output[1], output[2], output[3], IXAC, IXIL, IXIT, subIXTA_0, subIXTM_0, IXFRQ, IXH);
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
