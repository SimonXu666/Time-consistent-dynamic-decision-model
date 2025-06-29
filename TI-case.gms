* 预处理算法
$set T 3

* 计算集合元素数量
$eval N 2**%T%-1
$eval S 2**(%T%-1)

Set
    x       /x1*x2/
    n       /n1*n%N%/
    s(n)    /n%S%*n%N%/
    t       /t1*t%T%/;

ALIAS   (n,parent,child);
SET     tn(t,n)
        anc(parent,child) ancestor mapping
        nroot(t,n,s);

tn(t,n)$(2**(ord(t)-1)<=ord(n) and ord(n) < 2**ord(t)) = yes;
anc(parent,child)$(2*ord(parent)<=ord(child) and ord(child)<2*ord(parent)+2) = yes;
nroot(t,n,s)$(tn(t,n) and 2**(%T%-ord(t))*(ord(n)-2**(ord(t)-1))<ord(s)
and ord(s)<2**(%T%-ord(t))*(ord(n)-2**(ord(t)-1)+1)+1) = yes;

parameter return(n);
return(n)$(mod(ord(n),2) eq 1) = -0.5;
return(n)$(mod(ord(n),2) eq 0) = 1;

display anc,return,nroot;

Parameter
    lamda /0.5/
    theta /0.95/;

Variable OF(t,n),w(n),gama;
NonNegative Variables I(n,x),js(s);

Equation const1,const2,const3,constOF;

w.fx('n1') = 1;

const1(n)$(ord(n)<2**(%T%-1)).. sum(x,I(n,x)) =e= w(n);

const2(t,parent,child)$(anc(parent,child) and tn(t,child))..    
    w(child) =e= I(parent,'x1')+(1+return(child))*I(parent,'x2');
    
const3(s)..
    js(s) =g= gama - w(s);

Equation constOF;
constOF(t,n)$(tn(t,n) and ord(t)<%T%)..
    OF(t,n) =e= ((1-lamda)*sum(s$(nroot(t,n,s)), w(s))/2**(%T%-ord(t))
    + lamda*(gama-sum(s$(nroot(t,n,s)),js(s))/2**(%T%-ord(t))/(1-theta)));

parameter W_T_plan(n,s);

$set i 1
$label loop
Variable OF%i%;
Equation constOF%i%;
constOF%i%..  OF%i% =e= sum((t,n)$(ord(n)=%i% and tn(t,n)),OF(t,n));
Model node%i% / const1,const2,const3,constOF,constOF%i% /;
solve node%i% maximizing OF%i% using LP;
W_T_plan('n%i%',s) = w.l(s);
I.fx('n%i%','x1')=I.l('n%i%','x1');
$eval i %i%+1
$if not %i%==%S%  $goto loop

execute_unload "output%T%.gdx" W_T_plan;
execute 'gdxxrw.exe output%T%.gdx par=W_T_plan rng=Sheet1!A1'

