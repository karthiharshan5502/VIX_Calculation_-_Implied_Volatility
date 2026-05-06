proc import out = option
datafile= "/home/u63799767/MATH802/option(Jan2015-Feb2015).csv"
dbms="csv";
delimiter=',';
getnames='yes';
run;
proc print data=option(obs = 10);
run;

/* Convert Date variable to SAS date format */
data options;
set option;
format Date date9.;
run;
proc print data = options(obs=10);
run;

/*Filtering data for a particular day */
data jan2;
set options;
where date = '02JAN2015'd;
run;
proc print data = jan2(obs=10);
run;

/* VIX */
proc iml;
start main;
use jan2;
read all var {date tau k call_option_price s r put_option_price};
/* Calculate T1 and T2 */
T1 = min(tau)/365;/*Near term*/
T2 = max(tau)/365;/*Next term */
/* Initialize Tnew */
Tnew = j(nrow(tau), 1, .);
/* Assign T1 or T2 based on the condition */
do i = 1 to nrow(tau);
if tau[i] = min(tau) then Tnew[i] = T1;
else Tnew[i] = T2;
end;
/* Initialize matrix to store ΔKi */
Delta_K = j(nrow(k), 1, .);
/* Locate indices where tau = min(tau) and tau = max(tau) */
idx_tau_min = loc(tau=min(tau));/*Near term*/
idx_tau_max = loc(tau=max(tau));/*Next term*/
/* Calculate ΔKi for tau = min(tau) */
do i = 2 to ncol(idx_tau_min) - 1;
Delta_K[idx_tau_min[i]] = (k[idx_tau_min[i]+1] - k[idx_tau_min[i]-1]) / 2;
end;
Delta_K[idx_tau_min[1]] = k[idx_tau_min[2]] - k[idx_tau_min[1]];
Delta_K[idx_tau_min[ncol(idx_tau_min)]] = k[idx_tau_min[ncol(idx_tau_min)]] - k[idx_tau_min[ncol(idx_tau_min)]-1];
/* Calculate ΔKi for tau = max(tau) */
do i = 2 to ncol(idx_tau_max) - 1;
Delta_K[idx_tau_max[i]] = (k[idx_tau_max[i]+1] - k[idx_tau_max[i]-1]) / 2;
end;
Delta_K[idx_tau_max[1]] = k[idx_tau_max[2]] - k[idx_tau_max[1]];
Delta_K[idx_tau_max[ncol(idx_tau_max)]] = k[idx_tau_max[ncol(idx_tau_max)]] - k[idx_tau_max[ncol(idx_tau_max)]-1];

/* Calculate the absolute difference */
Difference = abs(call_option_price - put_option_price);

/* Calculate the midquote */
midquote = (call_option_price + put_option_price)/2;

/* Initialize matrix to store F */
F = j(nrow(k), 1, .);
/* Calculate F for every value of k */
e = exp(1);
do i = 1 to nrow(k);
F[i] = k[i] + (e##(r[i] * Tnew[i])) * (call_option_price[i] - put_option_price[i]);
 end;
/* Find the index where Difference is minimum and tau is minimum */
min_diff1 = min(Difference[idx_tau_min]);
idx_min_diff1 = loc(Difference = min_diff1 & tau = min(tau));
min_diff2 = min(difference[idx_tau_max]);
idx_min_diff2 = loc(difference = min_diff2 & tau = max(tau));
/* Extract the corresponding F1 value */
F1 = F[idx_min_diff1];/*Near term*/
F2 = F[idx_min_diff2];/*Next term*/
/* Determine K0 for near-term (K01) and next-term (K02) options */
K01 = .;
K02 = .;
do i = 1 to nrow(k);
if k[i] <= F1 then K01 = k[i];
if k[i] <= F2 then K02 = k[i];
end;
conbystrike = (Delta_K / (k##2)) # (e##(r # Tnew) # midquote);;
p1near = (2/T1)*sum(conbystrike[loc(tau=min(tau))]);
p1next = (2/T2)*sum(conbystrike[loc(tau=max(tau))]);
p2near = (1/T1) * ((F1/K01 - 1)**2);
p2next = (1/T2) * ((F2/K02 - 1)**2);

sigma1 = p1near - p2near; /* VOLATILITY FOR NEAR TERM*/

sigma2 = p1next - p2next; /*VOLATILITY FOR NEXT TERM*/
/* Initialize sigma */
sigma = j(nrow(tau), 1, .);
/* Assign T1 or T2 based on the condition */
do i = 1 to nrow(tau);
if tau[i] = min(tau) then sigma[i] = sigma1;
else sigma[i] = sigma2;
end;

N1 = min(tau);
N2 = max(tau);
N30 = 30;
N365 = 365;

wa = ((T1*sigma1*((N2-N30)/(N2-N1)))+((T2*sigma2*((N30-N1)/(N2-N1)))))*(N365/N30);
vix = 100*sqrt(wa);
print T1 T2 F1 F2 K01 K02 p1near p1next p2near p2next sigma1 sigma2 wa vix;
/* Create jan2new dataset with the calculated values */
create jan2new var {date tau k call_option_price s r put_option_price Difference midquote Tnew Delta_K conbystrike sigma};
append;
close jan2new;
finish main;
run;
quit;



data jan2new; /* Convert Date variable to SAS date format */
set jan2new;
format Date date9.;
run;
proc print data = jan2new (obs=10);
run;

/* IMPLIED VOLATILITY NEAR TERM */
PROC IML;
start main;
use jan2new;
read all var{date tau k call_option_price s r put_option_price Difference midquote Tnew Delta_K conbystrike sigma};
/* Find K closest to S0 */
idx = loc((abs(k - S) = min(abs(k - S))) & (tau = min(tau)));

/* Extracting values corresponding to the closest strike price and minimum tau */
K = k[idx];
T = Tnew[idx];
r = r[idx];
sigma = sigma[idx];
cop = call_option_price[idx];
S0 = s[idx]; 
pi = constant("pi");
total_run=3; /* 3 steps */
do i = 1 to total_run;
d1=(log(S0/K)+(r+sigma**2/2)*T)/(sigma*sqrt(T));
d2=d1-sigma*sqrt(T);
c=S0*cdf("Normal",d1)-K*(exp(-r*T)*cdf("Normal",d2));
vega=S0*1/sqrt(2*pi)*exp(-d1**2/2)*sqrt(T);
g=c-cop;
sigma=sigma-g/vega; /* Implied volatility*/
print k T r cop s0 d1 d2 c vega g sigma;
end;
finish main;
run;
quit;

/* IMPLIED VOLATILITY NEXT TERM */
PROC IML;
start main;
use jan2new;
read all var{date tau k call_option_price s r put_option_price Difference midquote Tnew Delta_K conbystrike sigma};
/* Find K closest to S0 */
idx = loc((abs(k - S) = min(abs(k - S))) & (tau = max(tau)));

/* Extracting values corresponding to the closest strike price and maximum tau */
K = k[idx];
T = Tnew[idx];
r = r[idx];
sigma = sigma[idx];
cop = call_option_price[idx];
S0 = s[idx]; 
pi = constant("pi");
total_run=3; /* 3 steps */
do i = 1 to total_run;
d1=(log(S0/K)+(r+sigma**2/2)*T)/(sigma*sqrt(T));
d2=d1-sigma*sqrt(T);
c=S0*cdf("Normal",d1)-K*(exp(-r*T)*cdf("Normal",d2));
vega=S0*1/sqrt(2*pi)*exp(-d1**2/2)*sqrt(T);
g=c-cop;
sigma=sigma-g/vega; /* Implied volatility*/
print k T r cop s0 d1 d2 c vega g sigma;
end;
finish main;
run;
quit;


/* PART 3 */
DATA part3;
    /* Define the variables in the dataset */
    INPUT DATE : $9.  VIX  NEARTERMIMPLIEDVOLATILITY  NEXTTERMIMPLIEDVOLATILITY;
DATALINES;
02JAN2015  20.089742   0.1301737  0.1368587
15JAN2015  24.994294   0.2157911  0.1774996
30JAN2015  25.189093   0.1882955  0.1703417
02FEB2015  22.507105   0.1446644  0.1487932
17FEB2015  19.715015   0.1078536  0.120316
26FEB2015  18.530757   0.0929367  0.1100007
;
RUN;

/* Sort the dataset by VIX */
PROC SORT DATA=part3;
    BY VIX;
RUN;

proc sgplot data=part3;
title 'VIX AND IMPLIED VOLATILITY';
series x=vix y=neartermimpliedvolatility / legendlabel='Near term IV'lineattrs=(color=blue);
series x=vix y=nexttermimpliedvolatility / legendlabel='Next term IV' lineattrs=(color=red);
xaxis label='VIX';
yaxis label='IMPLIED VOLATILITY';
run;

