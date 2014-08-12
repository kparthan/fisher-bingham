clc;

d = [2 2]';
a = [10 0 5 0]';
B = zeros(4); B(2,4) = 30;

display(' ');
display('Calculating SP approx');
logC_SP = logNormConstSP(d,a,B,4);
display(' ');
display('Calculating MC approx');
logC_MC = logNormConstMC(d,a,B,1e6);

display(' ');
display('Using');
display(d);
display(a);
display(B);
display(' ');
display('The SP and MC approximations to the log normalising const are as follows:');
display(' ');
display(['First order SP approx: ' num2str(logC_SP(1))]);
display(['Second order SP ("1+T" variant) approx: ' num2str(logC_SP(2))]);
display(['Second order SP ("exp(T)" variant) approx: ' num2str(logC_SP(3))]);
display(['By Monte Carlo integration: ' num2str(logC_MC)]);
display(' ');

