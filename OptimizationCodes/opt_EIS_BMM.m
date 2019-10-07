
rand('twister', sum(100*clock)); %New random seed every time
 opts=optimset('MaxFunEvals',100000,'TolCon',1.e-14,'TolFun',1.e-12,'TolX',1.e-12,'GradObj','on','MaxIter',10000,'Diagnostics','on','Display','on'...
        ,'UseParallel','never','Algorithm','sqp');

if restart==0
    clear X x
    restart=0;      
end



if type ==1
    n=2*k^2-(k-1)^2 ;
elseif type==2
   n= (5*k^2)/2 - k/2 -(k-1)^2
elseif type==3
    n= 2*k^2 + k-(k-1)^2
elseif type==4       
n= 2*k^2 + k -(k-1)^2
elseif type==5
    (5*k^2)/2 + k/2-(k-1)^2;
end
Aeq=[]; beq=[];
bd=2*k;
lb=[zeros(1,n)-bd]; lb(end)=-10^12*k;
ub=zeros(1,n)+bd; ub(end)=0;
%==============================================
count=0
info=-2;
Failed=0;
while (info==-2 || r<minr||  info==0)
    if restart==0
        if count==100
            Failed=1;
            ('exceed count')
            break
        end
    end
    
    if restart==1
        X1=X;
        x=X;
        if count==100;
            Failed=1;
            ('exceed count')
            x=X1
            break
        end
    elseif restart==2
        if count==0
            X1=X;
        elseif count ~=0
            X=X1;
        end
        x=X+.05*rand(1,length(X));
        if count==100
            Failed=1;
            ('exceed count')
            x=X1
            break
        end
        
    else
        x(1:n)=(2*rand(1,n));
        x(n)=-.01;
    end
    
    %==============================================
    %The optimization call:
    [X,FVAL,info]=fmincon(@Scheme_obj,x,[],[],Aeq,beq,lb,ub,@(x) nlc_EIS_BMM(x,k,p,type),opts);
    
    r=-FVAL;
    count=count+1;
    
end %while loop

[D,A,R,c] = unpackScheme(X,k,type);
[con,coneq,tau] = nlc_EIS_BMM(X,k,p,type);

