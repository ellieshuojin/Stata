use gss_nkids
label drop f91f
drop if childs >= 9
label def clab 0 "0" 1 "1" 2 "2" 3 "3" 4 "4" 5 "5" 6 "6" 7 "7" 8 "8 or more"
label val childs clab
twoway histogram childs, xlabel(,value) discrete xtitle(Number of Children (per Family)) percent xlabel(0(1)10)
// Poisson Regression
// 
// NULL model
glm childs, f(p)
glm, eform
poisson childs
poisson, irr
gen y = childs 
scalar lampois = exp(_b[_cons])
di lampois
scalar lampois = 2
di lampois
// generate the predicted distribution for this value of lambda and plot against empirical distribution
gen PrY_pois = exp(-lampois)*lampois^y/exp(lngamma(y+1))
//
twoway (histogram y, discrete fraction fcolor(gs14) ytitle(Pr(Y=y)) ) ///
       (line PrY_pois y, lpattern(solid) lwidth(.5) xtitle("Number of Children") sort legend(off)) 
// tables 
tab y 
tab y, sum(PrY_pois) mean
//
//fit full model
//
gen married = 0
 replace married = 1 if marital==1
gen nonwht  = 0
 replace nonwht = 1 if race > 1
glm childs coh10 married i.deg nonwht income, family(p) eform
//or another way
poisson childs coh10 married i.deg nonwht income, irr 
est store m1
estat ic
fitstat
//
gen boomer = 0
 replace boomer = 1 if coh10 > 7
poisson childs i.boomer i.married ib2.deg i.nonwht c.income, irr 
est store m2, title(Poisson)
estat ic
fitstat
margins married#boomer#nonwht, atmeans

// do Poisson Assumptions hold?
//  zero inflation NULL model accounts for excess 0 counts --> compare to Poisson NULL model
//
zip childs, inflate(_cons) irr
scal lamzip = exp(_b[childs:_cons])
scal p1 = invlogit(_b[inflate:_cons])
scal p0 = 1 - p1
gen PrY_zip = p1 + p0*( exp(-lamzip)*lamzip^y/exp(lngamma(y+1)) ) if y==0
replace PrY_zip = p0*( exp(-lamzip)*lamzip^y/exp(lngamma(y+1)) ) if (y ~=0)
//
 twoway (histogram y, discrete fraction fcolor(gs14) xtitle("Number of Children") ytitle(Pr(Y=y)) ) ///
        (line PrY_pois y, lwidth(.75) lpattern(solid) sort)  ///
		(line PrY_zip y,  lwidth(.75) lpattern(dash) sort legend(label(1 "Observed") ///
		                                                         label(2 "Poisson" ) ///
																 label(3 "Zip") rows(1))  )
//
//
// fit zip (which variables affect inflation?)
zip childs i.boomer i.married ib2.deg i.nonwht c.income, inflate(i.boomer i.married ib2.deg i.nonwht c.income)
//
poisson childs i.boomer i.married ib2.deg i.nonwht c.income
scalar ll = e(ll)
scalar npar = e(k)
scalar nobs = e(N)
scalar AIC = -2*ll + 2*npar
scalar BIC = -2*ll + log(nobs)*npar 
scalar list AIC
scalar list BIC
estat ic
fitstat
// trimmed zip
// fix issue with this model
zip childs i.boomer i.married ib2.deg i.nonwht c.income, inflate(i.boomer ib2.deg i.nonwht c.income)
est store m3, title(ZIP)
estat ic
fitstat
scalar ll = e(ll)
scalar npar = e(k)
scalar nobs = e(N)
scalar AIC = -2*ll + 2*npar
scalar BIC = -2*ll + log(nobs)*npar 
scalar list AIC
scalar list BIC
//
// Negative Binomial Regression
//
nbreg childs
nbvargr childs
* or
gen idnum = _n
xtpois childs, i(idnum)
xtpois childs i.boomer i.married ib2.deg i.nonwht c.income, i(idnum)
est store m4, title(NB)
xtpois, irr
predict xb, xb
gen lambda1 = exp(xb)
scalar alpha = 1/(exp([lnalpha]_cons))
gen vhat1 = (y + alpha)/(lambda1 + alpha)
sum vhat1, detail
histogram vhat1
gen off1 = log(vhat1)
glm childs i.boomer i.married ib2.deg i.nonwht c.income, f(p)
glm childs i.boomer i.married ib2.deg i.nonwht c.income, f(p) offset(vhat1)

estout m2 m3 m4 , label cells(b(star fmt(4))) stats(ll chi2 df_m ,  fmt(%9.2f %9.2f %9.0g) labels(logL "Chi-Sqr" df))  drop(0b.* 2b.*) ///
  varlabels(childs:_cons Constant inflate:_cons Constant 1.married "Married" 1.boomer "Baby Boomer" 1.nonwht "Nonwhite" 1.deg "< 12" 3.deg ">12"  ///
   c.income "Income" lnalpha:_cons "Var(u)")  ///
  varwidth(40)  title(Count Models) transform(lns1_1_1:_cons exp(@)^2 exp(@)^2) ///
  eqlabels("Fixed Effects" "Inflation" "Variance Components")
