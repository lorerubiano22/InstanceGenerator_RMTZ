$inlinecom   /* */
$offlisting
$offsymxref offsymlist
option optcr = 0;
option Solver=gurobi;
*option ResLim=10;
option OptCR=0.7;
option lp=gurobi;
option limrow=0;
option limcol=0;
option SolPrint=off;
option SolveLink=5;

Scalar
   start 'date + time at start of GAMS job'
   now   'current date + time'
   elapsed ' run + time'
   ;

*set           bn, be(bn,bn), zn, ze(zn,zn), p, t, sp(zn,zn,zn,zn), zebe(zn,zn,bn,bn), g(zn,zn);
*
set           bn, be, zn, ze, p, t, sp, zebe, g, revzones, citytaz, sftaz, sjtaz, oataz, zz(zn,zn), speb, check, list, jL, befix, beout, snn, enn;

*parameter     rev(zn,zn,p,t);
*parameter zoneborder(i,j), lambdax;
* singleton set cbd, pp(p);
singleton set cbd, pp(p);
parameter     rev(zn,zn,p,t), PTGweight, pi, TranTrips, trips, bnX, bnY, ntz, izr, lambdax, zone, zoneborder,  marketshare, avgprice, maxprice, bnd;
*scalar        r, xdim, ydim;

scalar        r, xdim, ydim, totalrev, counter, limW, degress;
alias         (zn,h,i,j,k);
alias         (bn,u,v,n,m);
set           iLj(i,j), uLv(u,v), zzpt(i,j,p,t);

Set           iter  /1*1000/;
parameter     uij(i,j), uij1(i,j), uij_iter(i,j,iter), uij_iter1(i,j,iter);

parameter ntz(i,j),zone;
sets check,jL,list;

*$GDXIN  "instance_zn_121_bn_144_nc_0.5_2055.gdx"
*$LOADDC bn be zn ze p t rev sp zebe xdim ydim g
*$GDXIN
$GDXIN  "SFBA_data_07_2026.gdx"
$LOADDC bn be zn t sp zebe g bnX bnY  p pi rev izr sftaz sjtaz oataz trantrips trips befix degress
$GDXIN



*Setting price systems
pp(p) = ord(p)$(ord(p)=1);


iLj(i,j)      = (ord(i) < ord(j) );
uLv(u,v)      = (ord(u) < ord(v) and be(u,v));
zzpt(iLj,p,t) = (rev(iLj,p,t));
revzones(j)   = (sum((i,pp,t)$zzpt(i,j,pp,t), 1) > 0);
speb(i,j,iLj(h,k),uLv(u,v))
              = ((sp(i,j,h,k) or sp(i,j,k,h)) and zebe(h,k,u,v));
citytaz(i,j)  = (iLj(i,j) and
                                   (sftaz(i) and sftaz(j))
                                or (sjtaz(i) and sjtaz(j))
                                or (oataz(i) and oataz(j))
                                    );

*set zz(zn,zn), speb(i,i,i,i,u,u);
zz(iLj) = (sum(t, rev(iLj,pp,t)));
speb(i,j,iLj(h,k),u,v) = ((sp(i,j,h,k) or sp(i,j,k,h)) and zebe(h,k,u,v));

bnd(u) = sum(v$be(u,v), 1);                 /* border node degree */ 
beout(be) =(not sum((h,k)$zebe(h,k,be),1)); /* border edges without corresponding PT edge */
limW  = 1;             /* constraining (outer rim ) flowvariable W upper bound 0 or 1 */








Parameter
init_u(i,j)
Uy(u,v)
Uy1(u,v);

init_u(i,j)=0;

uij_iter(i,j,iter) =init_u(i,j);
uij_iter1(i,j,iter) =init_u(i,j);
uij(i,j)=smin(iter,uij_iter(i,j,iter));
uij1(i,j)=smin(iter,uij_iter(i,j,iter));
Uy(u,v)=0;
Uy1(u,v)=0;
**********+*******
* Problem Lower *
**********+*******
variables        F;
binary variables Y;
Positive Variable W;

equations obj_lb, border, flow, degree,nntz,coupling1,coupling2;
obj_lb            .. F=e=sum((u,v), uy(u,v)*(Y(u,v)));
border(uLv(u,v))  .. Y(v,u) + Y(u,v)                     =l= 1;
flow(u)           .. sum(be(u,v), W(be)) - sum(be(v,u), W(be))    =e= 0;
*degree(u)         .. sum(be(u,v), W(be)) + sum(be(v,u), W(be))    =l= 2;
nntz(zz) .. sum(speb(zz,h,k,u,v), Y(u,v))+1              =l= card(t);
coupling1(u,v) ..card(be)*Y(u,v)=g=W(u,v);
coupling2(u,v) .. Y(u,v)=l=W(u,v);

model netFlow_lb /obj_lb, border, flow, nntz,coupling1,coupling2/;

**********+*******

******** problem upper bound ****
variables        F_up;
binary variables X;



equations obj_up,tariffstairs;
obj_up            .. F_up=e= sum(zzpt(zz(ilj),pp,t),rev(ilj,pp,t)*X(zzpt))-sum((ilj,pp,t),uij(ilj)*ord(t)*X(ilj,pp,t))+sum((ilj),uij(ilj));
tariffstairs(zz)  .. sum(zzpt(zz,pp,t), X(zzpt))                  =e= 1;
model rev_ub /obj_up,tariffstairs/;


*****************


***********************
* Lagrangian relaxation

*Relaxed constraints:
*tariffborders(zz) .. sum(zzpt(zz,pp,t), ord(t)*X(zzpt))- sum(speb(zz,h,k,u,v), Y(u,v))              =e= 1;
***********************

****************************************

alias(iter,iter2);

set
next_Iter(iter,iter2) 'Tuple of the current iteration and the subsequent iteration'
remaining_Iter(iter,iter2) ;

Scalar
abs_Viol
sum_Check
status_netFlow_lb
status_rev_ub
epsilon
factorGradient
steps
gradientLR
percen
a
la
;


Parameter
pc_lb(p)
pc_ub(p)
LB_bestP(p)
UB_bestP(p)
UB_best(iter)
LB_best(iter)
delta(iter)
lambda(iter)
UB(iter)
LB(iter)
maxR(i,j)
tij_lbSol(i,j)
tij(i,j)
rev2(i,j,p,t)
check_viol(i,j)
check_viol_control(i,j,iter)
check_viol1(i,j)
x_ub(i,j,p,t)
x_lb(i,j,p,t)
x_lbSol(i,j,p,t)
Y_lbSol(n,m)
W_lbSol(n,m)
Y_upSol(n,m)
x_ubIter(i,j,p,t,iter)
check_penaltyReward(iter)
check_violIter(iter)
LR_NEW_OF(iter)
LR_sum_U(iter)
controlVar(iter)
Uy(u,v)
gap(iter)
gapP(p)
rateViol(iter)
LR_penalty1(iter)
LR_penalty2(iter)
MaxIter(p)
stopParam
totalRunTime(p) 'total run + time'
totalCPU;


la=2;
steps=50;
stopParam=steps*3;
gradientLR=0.5;
LR_sum_U(iter)=inf;
UB_best(iter)=inf;
LB_best(iter)=-inf;
UB(iter)=inf;
LB(iter)=-inf;
delta(iter)=la;
lambda(iter)=la;
epsilon=0.1;
pc_lb(p)=0;
pc_ub(p)=0;
**********+*******
MaxIter(p)=max(stopParam,card(iter));
totalRunTime(p)=0;
totalCPU=0;
next_Iter(iter,iter2)$(ord(iter2)=ord(iter)+1)=yes;
remaining_Iter(iter,iter2)$(ord(iter2)>=ord(iter))=yes;
elapsed=0;




**********+*******
MaxIter(p)=max(stopParam,card(iter));
totalRunTime(p)=0;

next_Iter(iter,iter2)$(ord(iter2)=ord(iter)+1)=yes;
remaining_Iter(iter,iter2)$(ord(iter2)>=ord(iter))=yes;

*********************************
file rev_result /LR_result_07_noPattern_2026.csv/;
put  rev_result;
put "price degression parameter: "degress /;
put "FP.it, basePrice, avgPrice, maxPrice, PTshare, Revenue, Num.Zones, CPU, Opt, Gap" /;


*********************************


loop(p,
pp(p) = ord(p);
start=jnow;
*********************************
iLj(i,j)      = NO;
uLv(u,v)      = NO;
zzpt(iLj,pp,t) = NO;
zz(iLj) =  NO;
speb(i,j,iLj(h,k),u,v) = NO;
init_u(i,j)=0;

****

iLj(i,j)      = (ord(i) < ord(j) );
uLv(u,v)      = (ord(u) < ord(v) and be(u,v));
zzpt(iLj,p,t) = (rev(iLj,p,t));
revzones(j)   = (sum((i,pp,t)$zzpt(i,j,pp,t), 1) > 0);
speb(i,j,iLj(h,k),uLv(u,v))
              = ((sp(i,j,h,k) or sp(i,j,k,h)) and zebe(h,k,u,v));
citytaz(i,j)  = (iLj(i,j) and
                                   (sftaz(i) and sftaz(j))
                                or (sjtaz(i) and sjtaz(j))
                                or (oataz(i) and oataz(j))
                                    );
*$offText

*init_u(i,j)=0; tariffstairs
init_u(i,j)=0;
*********************************
*********************************
*********************************

uij_iter(i,j,'1') =init_u(i,j);
uij_iter1(i,j,'1') =init_u(i,j);
Uy(u,v)=0;
Uy1(u,v)=0;
LR_sum_U(iter)=inf;
UB_best(iter)=inf;
LB_best(iter)=-inf;
UB(iter)=inf;
LB(iter)=-inf;
delta(iter)=la;
lambda(iter)=la;
tij(ilj)=0;
tij_lbSol(ilj)=0;
x_lb(ilj,pp,t)=0;
x_lbSol(ilj,pp,t)=0;
Y_lbSol(n,m)=0;
Y_upSol(n,m)=0;
W_lbSol(n,m)=0;gap(iter)=0;
check_viol(ilj)=0;
check_viol1(ilj)=0;
check_violIter(iter)=0;
check_viol_control(ilj,iter)=0;
controlVar(iter)=0;
abs_Viol=0;
factorGradient=0;
rateViol(iter)=0;
LR_penalty1(iter)=0;
LR_penalty2(iter)=0;
x_ub(ilj,pp,t)=0;
x.l(i,j,pp,t)=0;
y.l(n,m)=0;
w.l(n,m)=0;

*********************************


loop(iter$(ord(iter)<=MaxIter(p) and totalRunTime(p)<5000),
   
    uij(ilj)=0;
    Uy(u,v)=0;
    x_lb(ilj,pp,t)=0;
    tij(ilj)=0;
    
    
    uij(ilj)=uij_iter(ilj,iter);    
    Uy(u,v)=sum(ilj$(sum(speb(ilj,h,k,u,v),1)>0),uij(ilj));
    
**upper bound problem
solve rev_ub max F_up using mip;

$onecho > gurobi.opt
threads 1
optfile=1; 
$offEcho

*lower bound problem
solve netFlow_lb max F using mip;
        
* Obtaining the values of x from the values of the variables y of the netFlow_lb Model
    tij(ilj)=sum(speb(ilj,h,k,u,v), Y.l(u,v))+1;
    x_lb(ilj,pp,t)$(tij(ilj)>0 and ord(t)= tij(ilj))=1;
    
********
LB(iter)=sum(zzpt(zz(ilj),pp,t),rev(ilj,pp,t)*x_lb(zzpt));


***************
** Computing violation
abs_Viol=0;
check_viol(ilj)=sum((t,pp),ord(t)*x.l(ilj,pp,t))-sum(speb(ilj,h,k,u,v),Y.l(u,v))-1;
check_viol_control(ilj,iter)=check_viol(ilj);
abs_Viol=sum((ilj),sqr(check_viol(ilj)));
check_violIter(iter)=abs_Viol;


if(ord(iter)>1,
rateViol(iter)=(check_violIter(iter-1)-check_violIter(iter))/check_violIter(iter-1);
);

if(LB(iter)>LB_best(iter),
    LB_best(iter2)$(remaining_Iter(iter,iter2))=LB(iter);
    pc_lb(p)=(jnow - start)*24*3600;
* clean
    tij_lbSol(ilj)=0;
    x_lbSol(ilj,pp,t)=0;
    Y_lbSol(n,m)=0;
    Y_upSol(n,m)=0;
    W_lbSol(n,m)=0;

*storing
    tij_lbSol(ilj)=tij(ilj);
    x_lbSol(ilj,pp,t)=x_lb(ilj,pp,t);
    Y_lbSol(n,m)=y.l(n,m);
    Y_upSol(n,m)=y.up(n,m);
    W_lbSol(n,m)=w.l(n,m);
    check_viol1(ilj)=sum((t,pp),ord(t)*x_lbSol(ilj,pp,t))-sum(speb(ilj,h,k,u,v),Y_lbSol(u,v))-1;

);

************
  UB(iter)=F_up.L+f.l;
    

    if(UB(iter)<UB_best(iter),
        pc_ub(p)=(jnow - start)*24*3600;
        UB_best(iter2)$(remaining_Iter(iter,iter2))=UB(iter);
        x_ub(ilj,pp,t)=0;
        x_ub(ilj,pp,t)=x.l(ilj,pp,t);
    );


gap(iter)=(UB_best(iter) - LB_best(iter))/LB_best(iter)

if(gap(iter)<=epsilon,
MaxIter(p)=ord(iter);
);

********************


    if(mod(ord(iter), steps) = 0,
    lambda(iter2)$(remaining_Iter(iter,iter2))=lambda(iter)*gradientLR;
    );
    

    
    if(ord(iter)>(steps/2) and mod(ord(iter)-(steps/2), steps) = 0 and gap(iter)=gap(iter-(steps/2)),
    if(lambda(iter)<0.01,
    lambda(iter2)$(remaining_Iter(iter,iter2))=la/2;
    else
    lambda(iter2)$(remaining_Iter(iter,iter2))=lambda(iter)*gradientLR;);
    );
    



   
   
    if(abs_Viol<>0,
        controlVar(iter)=UB_best(iter) - LB_best(iter);
        delta(iter2)$(remaining_Iter(iter,iter2))=lambda(iter)*((controlVar(iter))/abs_Viol);
    );
    
*** Multiplier


  uij_iter(ilj,iter2)$(next_Iter(iter,iter2))=(uij_iter(ilj,iter)+delta(iter)*(check_viol(ilj)));

*********************
gapP(pp)=gap(iter);
******
*ctime
 elapsed = (jnow - start)*24*3600;
);
    UB_bestP(pp)=smin(iter,UB_best(iter));
    LB_bestP(pp)=smax(iter,LB_best(iter));
    elapsed = (jnow - start)*24*3600;
    totalRunTime(p)=elapsed+totalRunTime(p);
    totalCPU=totalCPU+elapsed;



************* OUTPUT **********************
***** revenue
 ntz(iLj) = 1 + sum(speb(iLj,h,k,u,v), Y_lbSol(u,v));
 totalrev = sum((iLj(i,j),t)$(ord(t) = ntz(iLj)), rev(i,j,pp,t))+ sum(i, izr(i,pp));
 
****** compute price
 avgprice = sum((iLj(i,j),t)$(ord(t) = ntz(iLj)), pi(pp,t) ) / card(iLj);
 maxprice = smax((iLj,t)$(ord(t) = ntz(iLj)), pi(pp,t) );


****** compute market share
 marketshare = (sum((iLj,t)$(ord(t) = ntz(iLj)), TranTrips(iLj,pp,t)) + sum(i,TranTrips(i,i,pp,"1"))
                    ) / sum((iLj), trips(iLj))
                    ;



****** Compute Zones
zoneborder(i,j)$(g(i,j) and sum(be(n,m)$zebe(i,j,n,m), W_lbSol(n,m) + W_lbSol(m,n) ) >0) = 1;

lambdax = 1; check(i) = no;
jL(i)=yes;
loop(i$jL(i), 
  jL(i)    = no; 
  list(j)  = no;
  check(j) = no;
  list(i)  = yes; 
  while((card(list)>card(check)),
   loop(k$(list(k) and not check(k) ),
    check(k) = yes; 
    loop(j$(g(k,j) and
            jL(j)   and not
            list(j) and 
            (zoneborder(k,j) + zoneborder(j,k) = 0)
           ),
    list(j) = yes; 
    jL(j)   = no; 
    );
   );
  );
  loop(k$list(k),
   zone(k) = lambdax;
  ); 
  lambdax = lambdax + 1;
);
zoneborder(i,j) =0;
                    
 put ord(p) ","
     pi(pp,"1") ","
     avgprice ","
     maxprice ","
     marketshare:6:6 ","
    totalrev ","
    (lambdax -1) ","
    totalRunTime(p) ',' netFlow_lb.modelstat ',' gapP(pp):6:6
    /;
);





parameter grenzen; grenzen(n,m)$(Y_lbSol(n,m) or Y_lbSol(m,n)) = 1;
execute_unload "LR_result_07_noPattern_2026_all.gdx" ;

