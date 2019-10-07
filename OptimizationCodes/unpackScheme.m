function[D,A,R,c] = unpackScheme(X,k,type);    
 %This case uses the fact that weights on previous function values are zero
 %do to some restrictions it may cause with SSP coeffients
 
 %c - vector of abscissa
 %D-Weights on previous steps
 %A-Weights on function evaluations of previous steps
 %R-Weights on stages 
 
;
        %c=[3/2;1/2]; %uncomment and comment below for predetermined abscissa 
    %c=[-k+2:1];     %Spread out Nodes . **Makes Filtering Harder
    c=[1/k:1/k:1];  %Equidistant Nodes [0,1)
    count=k-1;
    %c=X(1:count)'; c(k)=1; 
    c=c(:);  
    X=X(count+1:end);        
     
   if type==4
    count=k*(k-1);    
    D=reshape(X(1:count),k,k-1); 
    D(:,k)=1-sum(D,2)
   else
    count=k-1;    
    D=repmat(X(1:count),k,1); 
    D(:,k)=1-sum(D,2);    
    
    X=X(count+1:end);

    count=k^2;
    A=reshape(X(1:count),k,k);
    X=X(count+1:end);
       

if type==1  %Explicit Parallel methods No using information from other blocks
      R = zeros(k); 

elseif type==2              %Explicit Method , which require serial implementation     
    count=.5*(k^2-k);
    ind  = tril(true(k),-1);
    R = zeros(k);
    R(ind) = X(1:count);
    X=X(count+1:end);
    
elseif type==3              % Implicit Parallel methods 
    count=k;
    R=diag(X(1:count));
    X=X(count+1:end);           
    
elseif type==4               % One Implicit Solve with Reuse Information          
    count=k;   
    R = zeros(k);
    R(:,1)=X(1:count)';
    X=X(count+1:end);  
    
elseif type==5              %K implicit solves with Reuse Information
    count=.5*(k^2+k);
    ind  = tril(true(k),0);
    R = zeros(k);
    R(ind) = X(1:count);
    X=X(count+1:end);  
    
end 

end