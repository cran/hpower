## Test the power and sample size routines in power.fcn.  This
## takes about five minutes on a SPARCstation 1.
##   Daniel F. Heitjan, 25 November 1991
##
## First source the file containing the functions.
library('hpower')
## Try to recompute a portion of Table 4.2 of Kraemer and Thiemann.
## Results differ slightly because K&T use one-sided 5% tests.
alpha_c(.05)
betavec_c(.90,.99)
x_matrix(1,1,1)
bt0vec_c(6:10)/10
sigma_matrix(1,1,1)
cc_matrix(1,1,1)
u_matrix(1,1,1)
th0_matrix(0,1,1)
tbl_matrix(0,5,2)
dimnames(tbl)_list(c('.6','.7','.8','.9','1.0'),c('.90','.99'))
init_matrix(c(31,25,18,17,1,53,41,32,24,21),5,2)
## Compute the table.
for (i in 1:5) {
 for (j in 1:2) {
  tbl[i,j]_sscomp(alpha,betavec[j],x,matrix(bt0vec[i],1,1),
   sigma,cc,u,th0,1,init[i,j])$n }}
cat('\nPart of Table 4.2 of Kraemer and Thiemann\n')
cat('Results differ slightly because the tests differ\n')
print(tbl)
##
## Now try to recompute the ANOVA example on p. 454 of Neter and
## Wasserman (1974).
cat('\n\nANOVA power calculations from Neter&Wasserman (1974)\n')
tbl_matrix(c(.72,0,.56,0,.36,0),2,3)
x_matrix(c(1,1,rep(0,8),
         0,0,1,1,1,0,0,0,0,0,
         0,0,0,0,0,1,1,1,0,0,
         0,0,0,0,0,0,0,0,1,1),10,4)
dimnames(tbl)_list(c('N&W','Ours'),c('Ex. 1','Ex. 2','Ex. 3'))
bt0_matrix(c(-3.5,-3,2,5),4,1)
sigma_matrix(2.5^2,1,1)
cc_t(matrix(c(1,-1,0,0,0,1,-1,0,0,0,1,-1),4,3))
u_diag(1)
th0_matrix(0,3,1)
tbl[2,1]_glhpwr(alpha,x,bt0,sigma,cc,u,th0)
sigma_matrix(3.0^2,1,1)
tbl[2,2]_glhpwr(alpha,x,bt0,sigma,cc,u,th0)
alpha_.01
sigma_matrix(2.5^2,1,1)
tbl[2,3]_glhpwr(alpha,x,bt0,sigma,cc,u,th0)
print(round(tbl,3))
##
## Now let's try a repeated-measures example from p.161 of
## Morrison (1976).
cat('\n\nA repeated-measures example from Morrison (1976)\n')
cat('Morrison gives the power as 93%; the difference is due to\n')
cat('approximation of the noncentrality parameter.\n')
tbl_matrix(c(0),1,3)
dimnames(tbl)_list(c('Power'),c('Test W','Test HLT','Test PB'))
alpha_.05
x_matrix(1,8,1)
bt0_matrix(c(6,2,0),1,3)
sigma_4*matrix(1,3,3)+6*diag(3)
cc_diag(1)
u_matrix(c(1,-1,0,0,1,-1),3,2)
th0_0*cc%*%bt0%*%u
tbl[1,]_c(glhpwr(alpha,x,bt0,sigma,cc,u,th0,1),
          glhpwr(alpha,x,bt0,sigma,cc,u,th0,2),
          glhpwr(alpha,x,bt0,sigma,cc,u,th0,3))
print(round(tbl,3))
##
## Here's example 5B.1, p. 217 of Johnson and Wichern (1982), who
## consider the power of Hotelling's T^2 test.
cat(
 '\n\nA repeated-measures example from Johnson&Wichern (1982)\n')
cat('Again differences are due to the approximation\n')
alpha_.05
x_matrix(1,26,1)
bt0.big_matrix(1,1,4)
bt0.sml_matrix(1,1,2)
sig.11_matrix(c(1,.2,.2,1),2,2)
sig.22_matrix(c(1,.1,.1,1),2,2)
sig.12_diag(2)
cc_diag(1)
u.big_diag(4)
u.sml_diag(2)
th0.big_0*cc%*%bt0.big%*%u.big
th0.sml_0*cc%*%bt0.sml%*%u.sml
mu_c(.2,.3,.4,.5)
rho_c(.7,.8,.9)
tbl_matrix(0,4,4)
dimnames(tbl)_list(c('.2','.3','.4','.5'),
                   c('2v','4v,r=.7','4v,r=.8','4v,r=.9'))
## Make their table, using the W test.  The results differ because
## Muller and Peterson compute the noncentrality parameter slightly
## differently -- it's a little smaller than J&W's.  Note that
## J&W use an exact power, whereas M&P use an approximation.
for (i in 1:4) {
 tbl[i,1]_glhpwr(alpha,x,mu[i]*bt0.sml,sig.11,cc,u.sml,th0.sml)
 for (j in 1:3) {
  sigma_rbind(cbind(sig.11,rho[j]*sig.12),
              cbind(rho[j]*sig.12,sig.22))
  tbl[i,j+1]_glhpwr(alpha,x,mu[i]*bt0.big,
   sigma,cc,u.big,th0.big) } }
cat('Table 5B.1, p. 218 of J&W, using the W test\n')
print(round(tbl,3))
##
