$inlinecom   /* */
$offlisting
$offsymxref offsymlist
option optcr = 0;
option solver=gurobi;
option ResLim =3600;
option lp=gurobi;
option limrow=0;
option limcol=0;
option SolPrint=off;
option    optcr    = 0.0;

* Saving datos computational timeclose
Scalar
   start 'date + time at start of GAMS job'
   now   'current date + time'
   elapsed ' run + time'
   totalRunTime 'total run + time';


* Problem information
* bn borders nodes
* zn bus stations
* be(bn,bn) district borders


set           bn, be(bn,bn), zn, ze(zn,zn), p, t, sp(zn,zn,zn,zn), zebe(zn,zn,bn,bn), g(zn,zn);
parameter     rev(zn,zn,p,t), maxRevP(p);
singleton set cbd, pp(p);
alias         (zn,h,i,j,k);
alias         (bn,u,v,n,m);
set           iLj(i,j), uLv(u,v), zzpt(i,j,p,t),zz(zn,zn), speb(i,i,i,i,u,u);
scalar        r, xdim, ydim;

set           bn, be(bn,bn), zn, ze(zn,zn), p, t, sp, zebe, g, revzones, citytaz, sftaz, sjtaz, oataz, check, list, jL, befix, beout, snn, enn;
singleton set cbd, pp(p);
parameter     rev, PTGweight, pi, TranTrips, trips, bnX, bnY, ntz, izr, lambdax, zone, zoneborder,  marketshare, avgprice, maxprice, bnd;
scalar        r, xdim, ydim, totalrev, counter, limW, degress;
alias         (zn,h,i,j,k);
alias         (bn,u,v,n,m);
set           iLj(i,j), uLv(u,v), zzpt(i,j,p,t),zz(zn,zn), speb(i,i,i,i,u,u);



$GDXIN  "SFBA_data_09_2026.gdx"
$LOADDC bn be zn t sp zebe g bnX bnY  p pi rev izr sftaz sjtaz oataz trantrips trips befix degress
$GDXIN




r     = card(zn) * 1;  /* select subset size for heuristic */
pp(p) = (ord(p) = 1);  /* select price scenario */

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



* Result variables and sets
variables        handle,bestW,handleTotalrev;
parameter ntz(i,j),zone;
sets check,jL,list;
scalar maxRevenue, accumualtive,totalrev;
scalar bestP;
parameter zoneborder, lambda;
parameter Ycopy(v,u),Wcopy(v,u);
* Problem variables
variables        F;
binary variables X,Y;
Positive Variable W;

equations        obj, selectzone, sumselect, eachzone;
obj                    ..  F =e= sum(zzpt(i,j,pp,t), rev(i,j,pp,t) * X(zzpt));
selectzone(iLj(i,j))   ..  sum(zzpt(i,j,pp,t), X(zzpt)) =l= 1;
sumselect              ..  sum(zzpt(i,j,pp,t), X(zzpt)) =l= r;
eachzone(j)$(ord(j)>1
      and revzones(j)) ..  sum(zzpt(i,j,pp,t), X(zzpt)) =g= 1;

model selectrelations /obj, selectzone, sumselect, eachzone/;

equations maxrev, tariffborders, tariffstairs, border, flow,coupling1,coupling2;

maxrev                 .. F - sum(zzpt(zz(i,j),pp,t),
                          (rev(i,j,pp,t))*X(zzpt)) =e= 0;
tariffborders(zz) .. sum(zzpt(zz,pp,t), ord(t)*X(zzpt))
                     - sum(speb(zz,h,k,u,v), Y(u,v))              =e= 1;
tariffstairs(zz)  .. sum(zzpt(zz,pp,t), X(zzpt))                  =e= 1;
border(uLv(u,v))  .. Y(v,u) + Y(u,v)                     =l= 1;
flow(u)           .. sum(be(u,v), W(be)) - sum(be(v,u), W(be))    =e= 0;
*degree(u)         .. sum(be(u,v), W(be)) + sum(be(v,u), W(be))    =l= 2;
coupling1(u,v) ..card(be)*Y(u,v)=g=W(u,v);
coupling2(u,v) .. Y(u,v)=l=W(u,v);



model tariffzones /maxrev, tariffborders, tariffstairs, border, flow,coupling1,coupling2/;


maxRevenue=0;
accumualtive=0;
totalRunTime=0;


file rev_result /Heu_result_09_noPattern_2026.csv/;
put  rev_result;
put "price degression parameter: "degress /;
put "FP.it, basePrice, avgPrice, maxPrice, PTshare, Revenue, Num.Zones, CPU, Opt, Gap" /;



loop(p,
pp(p) = ord(p);

*Saving computatinal time
$onecho > gurobi.opt
threads 0
parallelmode 0
$offecho

start=jnow;

solve selectrelations maximizing F using rmip;
zz(iLj) = 0;
zz(iLj) = (sum(t, X.l(iLj,pp,t)));
speb(i,j,iLj(h,k),u,v) = ((sp(i,j,h,k) or sp(i,j,k,h)) and zebe(h,k,u,v));

solve tariffzones max F using mip;

handle.l(p)= (F.l)
if(F.l>maxRevenue,
Ycopy(v,u)=0;
Wcopy(u,v)=0;
bestP=ord(p);
maxRevP(p)=F.l;
maxRevenue=F.l;
Wcopy(u,v)=W.l(u,v);
Ycopy(v,u)=Y.l(u,v);
ntz(iLj) = 0;
ntz(iLj) = 1 + sum(speb(iLj,h,k,u,v), Y.l(u,v));
totalrev = sum((iLj(i,j),t)$(ord(t) = ntz(iLj)), rev(i,j,pp,t));
);
elapsed = (jnow - start)*24*3600;
totalRunTime=elapsed+totalRunTime;

************* OUTPUT **********************
***** revenue
 ntz(iLj) = 1 + sum(speb(iLj,h,k,u,v), Y.l(u,v));
 totalrev = sum((iLj(i,j),t)$(ord(t) = ntz(iLj)), rev(i,j,pp,t))+ sum(i, izr(i,pp));
****** compute price
 avgprice = sum((iLj(i,j),t)$(ord(t) = ntz(iLj)), pi(pp,t) ) / card(iLj);
 
 maxprice = smax((iLj,t)$(ord(t) = ntz(iLj)), pi(pp,t) );
****** compute market share
 marketshare = (sum((iLj,t)$(ord(t) = ntz(iLj)), TranTrips(iLj,pp,t)) + sum(i,TranTrips(i,i,pp,"1"))
                    ) / sum((iLj), trips(iLj))
                    ;
****** Compute Zones
zoneborder(i,j)$(g(i,j) and sum(be(n,m)$zebe(i,j,n,m), Wcopy(n,m) + Wcopy(m,n) ) >0) = 1;

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
    totalRunTime ',' tariffzones.modelstat ',' (abs(tariffzones.objest-tariffzones.objval)/(1e-10+abs(tariffzones.objval))):6:6
    /;
************ END OUTPUT ******************* 
 X.lo(i,j,pp,"1") =0;

);


parameter grenzen; grenzen(n,m)$(Ycopy(n,m) or Ycopy(m,n)) = 1;
execute_unload "Heu_result_09_noPattern_2026_all.gdx" ;

