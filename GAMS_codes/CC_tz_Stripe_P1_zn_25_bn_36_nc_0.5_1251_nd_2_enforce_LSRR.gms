$inlinecom   /* */
$offlisting
$offsymxref offsymlist
option optcr = 0;
option solver=gurobi;
option lp=gurobi;
option ResLim=3600;
* Saving datos computational time
Scalar
   start 'date + time at start of GAMS job'
   now   'current date + time'
   elapsed ' run + time'
   totalRunTime 'total run + time';


* Problem information
* bn borders nodes
* zn bus stations
* be(bn,bn) district borders
* file name

set           bn, be(bn,bn), zn, ze(zn,zn), p, t, sp(zn,zn,zn,zn), zebe(zn,zn,bn,bn), g(zn,zn), ijij, SetTrip, beout, snn, enn;
parameter     rev(zn,zn,p,t), dij(zn,zn), price(p,t), s2, bnd, maxRevP(p);
singleton set cbd, pp(p);
alias         (zn,h,i,j,k);
alias         (bn,u,v,n,m);
set           iLj(i,j), uLv(u,v), zzpt(i,j,p,t);
scalar        r, xdim, ydim,counter;

$GDXIN  "instance_zn_25_bn_36_nc_0.5_1251.gdx"
$LOADDC bn be zn ze p t rev sp zebe xdim ydim g
$GDXIN

* Setting price system
*pp(p) = ord(p)$(ord(p)=10);

snn("bn19") = yes;  
snn("bn13") = yes;

enn("bn24") = yes;   
enn("bn18") = yes;


bnd(u) = sum(v$be(u,v), 1);                 /* border node degree */



iLj(i,j)      = (ord(i) < ord(j));
uLv(u,v)      = (ord(u) < ord(v) and be(u,v));
zzpt(iLj,p,t) = (rev(iLj,p,t));

* Result variables and sets
variables        handle,bestW,handleTotalrev;
parameter ntz(i,j),zone;
sets Ycopy(v,u),Wcopy(v,u),check,jL,list;
scalar maxRevenue, accumualtive,totalrev;
scalar bestP;
parameter zoneborder, lambda;
* Problem variables
variables        F;
binary variables X,Y,W;


set zz(zn,zn), speb(i,i,i,i,u,u);
zz(iLj) =  yes;
speb(i,j,iLj(h,k),u,v) = ((sp(i,j,h,k) or sp(i,j,k,h)) and zebe(h,k,u,v));


equations maxrev, tariffborders, tariffstairs, border, flow, degree;

maxrev            .. F - sum(zzpt(zz(i,j),pp,t),
                                         rev(i,j,pp,t) * X(zzpt)) =e= 0;
tariffstairs(zz)  .. sum(zzpt(zz,pp,t), X(zzpt))                  =e= 1;
tariffborders(zz) .. sum(zzpt(zz,pp,t), ord(t)*X(zzpt))
                     - sum(speb(zz,h,k,u,v), Y(u,v))              =e= 1;
border(uLv(u,v))  .. W(u,v) + W(v,u) - Y(u,v)                     =e= 0;
flow(u)           .. sum(be(u,v), W(be)) - sum(be(v,u), W(be))    =e= 0
                                                                     +2$snn(u)
                                                                     -2$enn(u)
                                                                       ;
degree(u)         .. sum(be(u,v), W(be)) + sum(be(v,u), W(be))    =l= 2
                                                                    + bnd(u)$(snn(u) or enn(u))
                                                                    ;

model tariffzones /all/;
*model tariffzones /maxrev, tariffborders, tariffstairs, border, flow, degree/;
Parameter status(*);

*execute_unload "wstripe_P1.gdx" W, be;

maxRevenue=0;
accumualtive=0;
totalRunTime=0;
*for(counter = 1 to 1 by 1,
*pp(p) = (ord(p)=counter);
loop(p,
pp(p) = ord(p);
option MIP = gurobi;
$onecho > gurobi.opt
threads 0
parallelmode 0
$offecho
*Saving computatinal time
start=jnow;

* Problem

tariffzones.OptFile = 1;
solve tariffzones max F using mip;
handle.l(p)= (F.l);
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
totalrev = sum((iLj(i,j),t)$(ord(t) = ntz(iLj)), rev(i,j,pp,t)+ rev(j,i,pp,t));
status('objval') = tariffzones.objval;status('objest') = tariffzones.objest;status('gapsolver') = ((tariffzones.objest-tariffzones.objval)/tariffzones.objest);

);
elapsed = (jnow - start)*24*3600;
totalRunTime=elapsed+totalRunTime;

);

execute_unload 'tmp_Stripe_p1_nd_2_zn_25_bn_36_nc_0.5_1251.gdx';  

** selecting arcs
zoneborder(i,j)$(g(i,j) and sum(be(n,m)$(zebe(i,j,n,m) and (Ycopy(n,m) or Ycopy(m,n))),1 ) >0) = 1;

lambda = 1; check(i) = no;
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
   zone(k) = lambda;
  ); 
  lambda = lambda + 1;
);
zoneborder(i,j) =0;

                



********

scalar scale,xc,yc;
scale = 12.25/(xdim-1);

scalar noz, r,gg,b;
parameter first(zn), second(zn),third(zn);

noz = smax(zn, zone(zn) );

while (noz > 0,
 r = uniform(0,1);
 gg = uniform(0,1);
 b = uniform(0,1);
 loop (zn$(zone(zn) = noz),
 first(zn)=r;
 second(zn)=gg;
 third(zn)=b;
 );
 noz = noz - 1;
);

********
file tzres /tzresult_Stripe_nd_2_zn_25_bn_36_nc_0.5_1251_p1.tex/; 
put  tzres;

put '\documentclass{article}' /;
put '\usepackage{tikz}' /;
put '\usetikzlibrary{arrows.meta}' /;
put '\begin{document}' /;
put '\begin{tikzpicture}[scale = ' scale:7:5 ']' /;
yc = 0;
loop(bn,
 xc = mod(ord(bn),xdim +1);
 if (xc = 0, xc = xdim +1);
 if (xc = 1, yc = yc + 1);
 put '\coordinate (' bn.tl:5:0 ') at ('xc:5:0',' yc:5:0 ');' /;
);
yc = 0;

loop(zn,
 xc = mod(ord(zn),xdim);
 if (xc = 0, xc = xdim);
 if (xc = 1, yc = yc + 1);
** this is new here

******
 put '\definecolor{currentcolor}{rgb}{'(first(zn))', '(second(zn))', '(third(zn))'};' /;  
 put '\node (' zn.tl:4:0 ') at ('(xc+0.5):4:1',' (yc+0.5):4:1 ') [fill=currentcolor, inner sep=1.5pt, draw, circle, text=blue, minimum height = 1cm]{{\tiny ' zn.tl:4:0'}};'  /;
);
loop(be(u,v)$(W.up(be)>0),
  put '\draw[arrows = -{Stealth[inset=0pt]}, blue] (' u.tl:5:0')--('v.tl:5:0 ');'/;
);
loop(g(i,j),
  put '\draw[dashed,line width = 0.1cm,gray] (' i.tl:5:0')--('j.tl:5:0 ');'/;
);
loop(be(u,v)$(Wcopy(u,v) or Wcopy(v,u)),
put '\draw[inner sep=0.5pt,line width = 0.1cm,red] (' u.tl:5:0')--('v.tl:5:0 ');'/;

*  put '\draw[-,red] (' u.tl:5:0')--('v.tl:5:0 ');'/;
);
loop(bn,
  put '  \node at (' bn.tl:5:0 ') [inner sep=1.5pt,blue, text=blue, fill=white,draw]{{\tiny ' bn.tl:5:0'}};'/;
);
put '\end{tikzpicture}' /;
put '\end{document}' /;
putclose;    
$exit
