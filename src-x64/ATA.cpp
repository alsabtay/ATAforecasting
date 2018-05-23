#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
 double meanMASE(NumericVector x, double meanNB) {
    double sum = 0.0;
	double mMASE = 0.0;
	int i;
	int n = x.size();
	for(i = 0; i < n; i++)
        sum += fabs(x[i]-meanNB);
    mMASE = 1.0 * sum / n;
	return mMASE;
}

// [[Rcpp::export]]
double meanIT(NumericVector x, int t){
  double total = 0.0;
  for(int i = 0; i < t; ++i)
    total += x[i];
  return 1.0 * total / t;
}


// [[Rcpp::export]]
double AutoATACore(NumericVector IAZ, int IZP, int IZQ, double IZPHI, int IZMO, int IZAC, int IZIL, int IZIT, NumericVector IZTA_0, NumericVector IZTM_0) {	
	int LENZ = IAZ.size();
	NumericVector IZT_0(LENZ);
	NumericVector IZFIT(LENZ);
	double coefpk, coefqk, Xobs, Xlag, S, T, S_1, T_1, T_0;
	int i, indx;
	NumericVector ITErr(LENZ);
	NumericVector ITAcc(LENZ);
	double accmeasure=0.0;
	double naive_benchmark=0.0;
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
	if ( (IZAC==1) | (IZAC==2) | (IZAC==12) )
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
	{
		naive_benchmark = meanMASE(IAZ,mean(IAZ));
		accmeasure = meanMASE(ITErr,naive_benchmark);
	}
	else
		accmeasure = 0;
	return accmeasure;
}

// [[Rcpp::export]]
NumericVector AutoATADamped(NumericVector IAX, int IXP, int IXQ, int IXMO, int IXAC, int IXLF, int IXTF, double IXPHIS, double IXPHIE, double IXPHISS, int IXIL, int IXIT, NumericVector IXTA_0, NumericVector IXTM_0){	
	int LENX = IAX.size();
	int  	d_opt_p;
	int  	d_opt_q;
	double  d_opt_phi;
	NumericVector out(4);
	double optAccryStart;
	double optAccryEnd;
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
					optAccryEnd = AutoATACore(IAX, i, 0, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
					if (optAccryEnd <= optAccryStart) {
						d_opt_phi = 1.0 * k;
						d_opt_p = i;
						d_opt_q = 0;
						IXMO = m;
						optAccryStart = optAccryEnd;
					}
				}
				for(j = 0; j <= d_opt_p; j++) {
					optAccryEnd = AutoATACore(IAX, d_opt_p, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
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
					optAccryEnd = AutoATACore(IAX, i, 1, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
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
	else {
		if ( (IXP==-1) & (IXQ==-1) ) {
			d_opt_q = 0;
			d_opt_p = 0;
			for(m = mstart; m <= mfinish; m++) {
				for(k = IXPHIS; k < IXPHIE+IXPHISS; k = k+IXPHISS) {
					for(i = 1; i <= LENX; i++) {
						for(j=0; j<=i; j++) {
							optAccryEnd = AutoATACore(IAX,i, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
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
						optAccryEnd = AutoATACore(IAX, i, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
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
						optAccryEnd = AutoATACore(IAX, IXP, j, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
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
					optAccryEnd = AutoATACore(IAX, IXP, IXQ, k, m, IXAC, IXIL, IXIT, IXTA_0, IXTM_0);
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