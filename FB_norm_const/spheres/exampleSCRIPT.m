clc;

%d = [2 2]';
%a = [10 0 5 0]';
%B = zeros(4); B(2,4) = 30;

d = [3]';
a = [100 0 0]';                                % kmu
B = zeros(3); B(2,2) = 14.5; B(3,3) = -14.5;    % beta * (y2 y2' - y3 y3')

%a = [2.383908e+01 4.496280e+00 9.701279e+01]';
%B = [-3.417925e+01 -2.940266e+01 9.761645e+00; 
%     -2.940266e+01 3.683373e+01 5.518010e+00;
%     9.761645e+00 5.518010e+00 -2.654487e+00]

%a = [-7.923010e+01 -3.885625e+01 4.704024e+01]';
%B = [1.689740e+01 -2.055747e+01 1.147946e+01; 
%     -2.055747e+01 -4.536178e+00 -3.837203e+01;
%     1.147946e+01 -3.837203e+01 -1.236122e+01]
 
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

