## S functions for calculating power and sample size requirements
## when testing the general linear hypothesis.
##   Daniel F. Heitjan, 25 November 1991
##
pfnc_function(q,df1,df2,lm,iprec=c(6)) {
##   Daniel F. Heitjan, 15 May 1991
 iprec_iprec[1]
 nrw_max(c(length(q),length(df1),length(df2),length(lm)))
 if (length(q)<nrw) q_rep(q[1],nrw)
 if (length(df1)<nrw) df1_rep(df1[1],nrw)
 if (length(df2)<nrw) df2_rep(df2[1],nrw)
 if (length(lm)<nrw) lm_rep(lm[1],nrw)
 mlm_max(lm)
 xprec_floor(mlm/2+iprec*(mlm/2)^0.5)+1
 qa_rep(q,xprec+1)
 df1a_rep(df1,xprec+1)
 df2a_rep(df2,xprec+1)
 lma_rep(lm,xprec+1)
 ilma_floor(lma/2+iprec*(lma/2)^0.5)+1
 j_rep((0:xprec),rep(nrw,xprec+1))
 wt_0+1*(j==0)
 wt[lma!=0]_exp(-lma[lma!=0]/2+j[lma!=0]*log(lma[lma!=0]/2)-
  lgamma(j[lma!=0]+1))
 wt_wt*pbeta(df1a*qa/(df1a*qa+df2a),j+df1a/2,df2a/2)*(j<=ilma)
 wt_matrix(wt,nrw)
 c(apply(wt,1,sum)) }
##
pchisqnc_function(q,df,lm,iprec=c(6)) {
##   Daniel F. Heitjan, 15 May 1991
 iprec_iprec[1]
 nrw_max(c(length(q),length(df),length(lm)))
 if (length(q)<nrw) q_rep(q[1],nrw)
 if (length(df)<nrw) df_rep(df[1],nrw)
 if (length(lm)<nrw) lm_rep(lm[1],nrw)
 mlm_max(lm)
 xprec_floor(mlm/2+iprec*(mlm/2)^0.5)+1
 qa_rep(q,xprec+1)
 dfa_rep(df,xprec+1)
 lma_rep(lm,xprec+1)
 ilma_floor(lma/2+iprec*(lma/2)^0.5)+1
 j_rep((0:xprec),rep(nrw,xprec+1))
 wt_0+1*(j==0)
 wt[lma!=0]_exp(-lma[lma!=0]/2+j[lma!=0]*log(lma[lma!=0]/2)-
  lgamma(j[lma!=0]+1))
 wt_wt*pchisq(qa,2*j+dfa)*(j<=ilma)
 wt_matrix(wt,nrw)
 c(apply(wt,1,sum)) }
##
glhpwr_function(alpha,x,bt1,sigma,cc,u,th0,test=1,tol=1.e-8) {
##   Daniel F. Heitjan, 25 November 1991
 theta_cc%*%bt1%*%u-th0
 hpop_t(theta)%*%solve(cc%*%solve(t(x)%*%x)%*%t(cc))%*%theta
 n_nrow(x)-qr(x,tol)$rank
 epop_n*t(u)%*%sigma%*%u
 a_nrow(cc)
 b_ncol(u)
 s_min(a,b)
 if (test==1) {
## Calculations for Wilks's lambda test (W).
  W_prod(eigen(epop)$values)/prod(eigen(hpop+epop)$values)
  gnum_a^2*b^2-4
  gden_a^2+b^2-5
  if ((gnum==0)&(gden==0)) {g_1} else {g_(gnum/gden)^0.5}
  df1_a*b
  df2_g*(n-(b-a+1)/2)-(a*b-2)/2
  lm_(1-W^(1/g))*df2/W^(1/g)
 } else if (test==2) {
## Calculations for the Hotelling-Lawley trace test (HLT).
  HLT_sum(diag(hpop%*%solve(epop)))
  df1_a*b
  df2_s*(n-b-1)+2
  lm_df2*HLT/s
 } else if (test==3) {
## Calculations for the Pillai-Bartlett trace test (PB).
  PB_hpop%*%solve(hpop+epop)
  PB_sum(diag(PB))
  df1_a*b
  df2_s*(n+s-b)
  lm_df2*PB/(s-PB)
 }
 if ((df1<=0)|(df2<=0)) { 
  0.0
 } else { 
  fcrit_qf(1-alpha,df1,df2)
  1-pfnc(fcrit,df1,df2,lm)
 } }
##
sscomp_function(alpha,power,x,bt1,sigma,cc,u,th0,test=1,ninit=1) {
##   Daniel F. Heitjan, 14 May 1991
 xx_x
 nrx_nrow(x)
 if (ninit>1) {
  for (i in 2:ninit) xx_rbind(xx,x) }
 n_ninit
 pwr_glhpwr(alpha,xx,bt1,sigma,cc,u,th0,test)
 if (pwr>power) {
  while ((pwr>power)&(n>0)) {
   xx_matrix(xx[(-1:(-nrx)),],ncol=ncol(xx))
   n_n-1
   opwr_pwr
   pwr_glhpwr(alpha,xx,bt1,sigma,cc,u,th0,test)
  }
  pwr_opwr
  n_n+1 }
 else if (pwr<power) {
  while (pwr<power) {
   xx_rbind(xx,x)
   n_n+1
   pwr_glhpwr(alpha,xx,bt1,sigma,cc,u,th0,test)
  }
 }
 list(n=n,power=pwr) }
##
