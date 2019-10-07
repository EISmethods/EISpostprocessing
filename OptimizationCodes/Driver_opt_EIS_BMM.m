
%Different parameters to try when runing the optimization routine
%k-Number of steps we are using
%p-Desired order of accuracy

%type =  1- Explicit Parallel i.e classic BMM
%        2- Explicit Method , which uses function evaluations of other stages (requires serial implementation)
%        3- Implicit Parallel methods i.e Implicit BMM
%        4- One Implicit Solve with Reuse Information
%        5- Implicit Method , which uses function evaluations of other stages (requires serial implementation)

%restart = 0 - Random starting vector of unknowns
%          1 - Start from a currently defined method 
%          2 - Start from a perturbation of a currently defined method. Hels get out of local minima

%minr - Smallest acceptable value for objective function
  %if restart= 0 start with minr=0 or some small value and grow using restart =1/2
  %if restart = 1/2 start with minr=r (this avoids getting a worse method  than you currently have.
  
  
  type=1
  k=2;
  p=3;
  minr=0
  restart=0;
  opt_EIS_BMM
  tau,coneq,A, D , R ,c,   %the outputs I tend to look at
  
  