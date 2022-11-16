
$onExternalInput

set    i    zones         /z01  *   z36 /;
set    b    border nodes  /bn01  *   bn49 /;
set    p    price stages  /p1   *   p10/;
set    t    tariff zones  /t1   *   t7/;
$offExternalInput

$onExternalInput
Parameter
netConnect    fraction rho  /1/
xdim /6/


*execute_unload 'parameters.gdx' a,xdim,ydim, boderNodes;
$offExternalInput
ydim
Scalar boderNodes;
ydim=xdim;
boderNodes  =  card(i)+xdim+ydim+1;

option limrow = 0, limcol = 0, solprint = off, sysout = off,
seed=5793;


alias (b,n,m), (u,v,i,j,kk), (t,t2);

sets
       e(b,b)      border edge (two border nodes) - Dual graph
       a(u,u)      potential public transport network arcs - Original gragh
       ae(u,u,b,b) incident border edge of transport arc -  simplify the model??
    
;


*  geographic settings
*  -------------------
*  xy coordinates of border edges and public transport arcs
parameter xb(b), yb(b), xu(u), yu(u);
  xb(b) = mod(ord(b)-1,xdim+1) + 1;
  yb(b) = floor(div(ord(b)-1,xdim+1)) + 1;
  xu(u) = mod(ord(u)-1,xdim) + 1.5;
  yu(u) = floor(div(ord(u)-1,xdim)) + 1.5;

* borders and potential transport arcs
  e(m,n) = (xb(n) - xb(m) = 1 and yb(n) - yb(m) = 0 ) or
           (xb(n) - xb(m) = 0 and yb(n) - yb(m) = 1 );
  a(i,j) = (xu(i) - xu(j) = 1 and yu(i) - yu(j) = 0 ) or
           (xu(i) - xu(j) = 0 and yu(i) - yu(j) = 1 );
  a(i,j) = (a(i,j) or a(j,i));



set iLj(i,j) auxilliary set i less than j;
    iLj(i,j) = (ord(i) < ord(j));
$onExternalOutput
parameter d(i,j) length of transport network arc;
$offExternalOutput
  d(a(iLj))        = (10 + uniform(-1,1))/xdim;
  d(a(i,j))$d(j,i) = d(j,i);
  


* central business district located in the centre
$onExternalOutput
singleton set cbd(i);
$offExternalOutput
cbd(i) = (ord(i) = round((card(i)+1)/2));

* definition of public network graph
* shortest paths from central business district
parameter  s(i) supply of node;
s(i) = - 1 + card(u)$cbd(i);



variable          F;

positive variable X;

equations  obj, flow;
obj     .. F - sum(a, d(a)*X(a))                 =e=    0;
flow(i) .. sum(a(i,j), X(a)) - sum(a(j,i), X(a)) =e= s(i);

model spath /obj, flow /;

display a;

solve spath using lp minimizing F;


$onExternalOutput
set g(i,j) transport network graph;
$offExternalOutput
g(a(i,j)) = (X.l(i,j) or X.l(j,i));

* selecting additional arcs
equation add_arcs;
  add_arcs .. sum(a(iLj)$(not g(a)), X(a)) =e= round((card(a) - card(g))*netConnect/2);

model extension /obj, add_arcs/;

X.up(a)   = 1;

solve extension using lp minimizing F;


g(a(i,j)) = (g(a) or X.l(i,j) or X.l(j,i));
x.up(a)   = inf;
a(i,j)    = (g(i,j));



ae(a(i,j),e(n,m)) =
           abs( xu(i)- xb(n)) = 0.5 and abs( xu(i)- xb(m)) = 0.5 and
           abs( xu(j)- xb(n)) = 0.5 and abs( xu(j)- xb(m)) = 0.5 and
           abs( yu(i)- yb(n)) = 0.5 and abs( yu(i)- yb(m)) = 0.5 and
           abs( yu(j)- yb(n)) = 0.5 and abs( yu(j)- yb(m)) = 0.5;

* distance matrix transport network
parameter  s2(u,i) supply of node;
s2(u,i) = - 1;
s2(u,u) = - 1 + card(u);

positive variables X2;

equations  obj2, flow2;
obj2        .. F - sum((u,a), d(a)*X2(u,a))                =e=    0;
flow2(u,i)  .. sum(a(i,j), X2(u,a)) - sum(a(j,i), X2(u,a)) =e= s2(u,i);


model spath2 /obj2, flow2 /;

solve spath2 using lp minimizing F;


parameter d2(i,j);
d2(i,j) = d(i,j);

positive variables X3;
equation dist, obj3;
dist(u,a(i,j))$X2.l(u,a) .. X3(u,j) =g= X3(u,i) + d2(i,j);
obj3                     .. F       =e= sum((u,i), X3(u,i));
model spath3 /flow2, dist, obj3/;
solve spath3 using lp minimizing F;



$onExternalOutput
parameters dij;
$offExternalOutput
Parameter d2ij;
dij(i,j) = X3.l(i,j);

* Max number of tariff zones between i and j
d2(i,j) = 1;
solve spath3 using lp minimizing F;

d2ij(i,j) = X3.l(i,j)+1;

set uv(u,v), ua(u,i,j);

uv(u,v) = (ord(u)<ord(v));
ua(u,a) = (X2.l(u,a));

* arcs of each shortest path relation from u to v
positive variable X4;
equation obj4, flow4;
obj4             .. F - sum((uv(u,v),ua(u,a)), X4(uv,a))  =e= 0;
flow4(uv(u,v),i) .. sum(ua(u,i,j), X4(u,v,i,j)) -
                    sum(ua(u,j,i), X4(u,v,j,i))
                    =e= 0 + 1$(ord(i) = ord(u)) - 1$(ord(i) = ord(v));
model spath4 /obj4, flow4/;
solve spath4 using lp minimizing F;


variable           F;
positive  variable y, q;


equations objsrt, flowsrt;
objsrt       ..  F =e= sum((i,a), d(a) * y(i,a));
flowsrt(i,j)..
              sum(a(j,u), y(i,j,u)) - sum(a(u,j), y(i,u,j))
              =e= - 1 + card(i)$(ord(i) = ord(j));
model spat /objsrt, flowsrt/;
spat.optfile=1;
solve spat minimizing F using LP;

set ia;
* /* shortest path tree for TAZ i */
ia(i,a)   = y.l(i,a);   

y.lo(i,a) = y.l(i,a);

* q <- total areas touched by the shortest route

equation border;
border(i,j,u)$(ia(i,j,u)) .. q(i,u) =e= q(i,j) + 1;
model max_number_borders /obj2, border/;
solve max_number_borders minimizing F using LP;

parameter maxzones;
maxzones(i,j) = 1 + q.l(i,j);
maxzones(i,j) = max(maxzones(i,j), maxzones(j,i));


$onExternalOutput
parameter  rev(i,j,p,t);
set  sp;
$offExternalOutput
set ia, ijij, rel, sp;

set ijpt;
set taz2;
*/* only taz that are connected by PT */

taz2(j)                 = (sum(i, a(i,j) + a(j,i)));
set SetTrip; SetTrip(i,j)$(taz2(i) and taz2(j) and ord(i)< ord(j) )=yes;


ijpt(i,j,p,t) = (ord(t) <= maxzones(i,j));






rev(i,j,p,t)$( ijpt(i,j,p,t)) = 0 ;

rel(i,j)=(sum(u$ia(i,u,j),1)>0);



ijij(i,j,u,v)$( ia(i,u,v) and rel(i,j) ) = yes;



$onExternalOutput
positive variable W;
$offExternalOutput
W.fx(i,j,u,v)$(ijij(i,j,u,v) and not ia(i,u,v)) = 0;




variable F;
equation objsp, flowsp;
objsp          .. F - sum(ijij, W(ijij))  =e= 0;

flowsp(i,j,u)$(SetTrip(i,j) and ord(i) < ord(j) and taz2(i) and taz2(j) and taz2(u) and rel(i,j))
              .. sum(ijij(i,j,u,v), W(i,j,u,v))
               - sum(ijij(i,j,v,u), W(i,j,v,u))
             =e= 0
                 + 1$(ord(i) = ord(u))
                 - 1$(ord(j) = ord(u));

model shortestpaths /objsp, flowsp/;
shortestpaths.optfile=1;
solve shortestpaths using lp minimizing F;


sp(ijij)$W.l(ijij) = (W.l(ijij));



*execute_unload 'distancia.gdx' dij;

* Gravity model to determine number of trips
scalar total_att, tinh;

parameters inh(u)     inhabitants
           trips(i,j) trips
;
inh(i)   = round(uniform(10000,20000));
*execute_unload 'addInf.gdx' inh;

tinh     = sum(i,inh(i));
inh(i)   = round(inh(i)/tinh*500000);
*execute_unload 'addInf1.gdx' inh, tinh;

loop(i,
   total_att  = sum(j, inh(j)/ (1/xdim + dij(i,j) - 0.75*dij(i,j)$cbd(j)));
   trips(i,j) = inh(i)*inh(j)/ (1/xdim + dij(i,j) - 0.75*dij(i,j)$cbd(j)) / total_att;
);

execute_unload 'trips.gdx' trips;
*$exit
scalar ttt;
ttt = sum((i,j), trips(i,j));
display ttt;


* Utility to use public transport
* Minimum distance mind(i,j)
* Utility util(p,i,j,t)
* Fraction of public transport users frac(p,i,j,t)

parameter price;

price(p,t)  = (1 + 0.2*(ord(p)-1)) * ord(t);


$onExternalOutput
parameter mind(i,j), util(p,i,j,t), frac(p,i,j,t);

$offExternalOutput
mind(i,j)     = (abs(xu(j) - xu(i)) + abs(yu(j) - yu(i)))*10 / xdim;

util(p,i,j,t) = - price(p,t)/(1/xdim + mind(i,j)) - 0.15*(dij(i,j) - 0.9*mind(i,j)$cbd(j));



frac(p,i,j,t) = exp(util(p,i,j,t))/( 1 +  exp(util(p,i,j,t)));

* Number of public transport trips
parameter  ptt(p,i,j,t);

ptt(p,i,j,t)$(ord(t) <= d2ij(i,j)) = frac(p,i,j,t) * trips(i,j);
rev(i,j,p,t) = ptt( p,i,j,t) * price(p,t);


parameter rconst(p), tconst(p);

tconst(p)    = sum(i, frac(p,i,i,'t1') * trips(i,i));
rconst(p)    = tconst(p) * price(p,'t1');


rev(i,j,p,t) = round(100*rev(i,j,p,t));
rconst(p)    = round(100*rconst(p));




scalar scale;

scale = 12/(xdim-1);



set tb(b,b) potential tariff border;

loop(ae(i,j,n,m)$a(i,j), tb(n,m) = yes);
e(n,m) = (e(n,m) or e(m,n));

* parametrization tz_ring
*$LOADDC bn be zn ze p t rev sp zebe xdim ydim g

$onExternalOutput
alias (b,bn), (u,v,i,j,zn), (t,t2);
alias (e,be), (iLj,ze), (ae,zebe)
$offExternalOutput
execute_unload 'informationSofar.gdx' bn, be, d, dij, zn,ze, p, rev, t, zebe,sp,cbd, xdim ydim,W,  g;
execute_unload 'instance_zn_36_bn_49_NC_1_5793.gdx' bn, be, zn,ze, p, rev, t, zebe,sp,cbd, xdim ydim,  g;


file log /Illusttratrion.tex/;
put  log;
*put '\documentclass{article}' /;
put '\usepackage{tikz}' /;
put '\usetikzlibrary{arrows.meta}' /;
put '\begin{document}' /;
put '\begin{tikzpicture}[scale = ' scale:7:5 ']' /;
put '  \draw[step=1cm,gray,thick] (1,1) grid (' (xdim)'+1 , ' (ydim) '+1 );' /;
*put '  \draw[step=1cm,gray,thick] (1,1) grid (' (xdim):6:0 ',' (ydim):6:0 ');' /;
loop(a(iLj(i,j)),
 put '  \draw[dashdotted,blue, line width = 0.05cm] (';
 put xu(i):6:1 ',' yu(i):6:1  ') -- (' xu(j):6:1 ',' yu(j):6:1 ');' /;
);
loop(b,
 put '  \draw (' xb(b):3:0 ',' yb(b):3:0 ') node [rounded corners=0.2pt, inner sep=1pt, text=white, fill=gray]{{\tiny ' b.tl:<6:0'}};'/;
);
loop(u,
 put '  \draw (' xu(u):4:1 ',' yu(u):4:1 ') node [circle, inner sep=0.5pt, text=white, fill=blue]{{\tiny \bf ' u.tl:<4:0'}};'/;
);

put '\end{tikzpicture}' /;
*put '\end{document}'/;
putclose;
execute_unload 'taz.gdx';
$exit

