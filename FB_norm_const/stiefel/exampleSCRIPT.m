clc;

d=3;
p=3; 
A = diag([10 5 1]); 
B = magic(d*p)/30; 

display(' ');
display('Calculating SP approx');
logC_SP = logNormConstSP(A,B,4);
display(' ');
display('Calculating MC approx (this may take a few minutes)');
logC_MC = logNormConstMC(A,B,1e6);

display(' ');
display('Using');
display(A);
display(B);
display(' ');
display('The SP and MC approximations to the log normalising const are as follows:');
display(' ');
display(['First order SP approx: ' num2str(logC_SP(1))]);
display(['Second order SP ("1+T" variant) approx: ' num2str(logC_SP(2))]);
display(['Second order SP ("exp(T)" variant) approx: ' num2str(logC_SP(3))]);
display(['By Monte Carlo integration: ' num2str(logC_MC)]);

display(' ');
display(' ');
display('Note: in the Fisher matrix case, i.e. when B=0, use');
display('logNormConstSP_matrixFisher for much faster calculations.');

display(' ');
display('E.g., for A as above and B=0, the first order SP and the two');
display('second-order SP approxes are: ');
solSP_matrixFisher = logNormConstSP_matrixFisher(A,4); 
display(num2str(solSP_matrixFisher));
