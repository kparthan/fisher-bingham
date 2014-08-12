function out = sumOverPerms(fh,CHatInv,p,dcum,inds,lambda)

switch fh
   case {'lambdaCClambda', 'lambdaCCClambda', 'lambdaCCCClambda', 'lambdaCCCCClambda'}
      % for these, only provide the middle indices, which are then
      % permuted; the first and last indices are summed over within this
      % function
      permedInds = perms(inds);
      numPerms = size(permedInds,1);
   case {'trCC', 'trCCC', 'trCCCC', 'trCCCCC'}
      % permute all but the first index
      permedInds = [repmat(inds(1),[factorial(length(inds)-1) 1]) perms(inds(2:end))];
      numPerms = size(permedInds,1);
end

runSum = 0;

funcString = [fh '(CHatInv,p,dcum,permedInds(ind_,:),lambda)'];
for ind_ = 1:numPerms
   runSum = runSum + eval(funcString);
end

out = runSum;


% ---------------------------------------------------------------------

function out = trCC(C,p,dcum,inds,lambda)

out = trace(...
   subm(C,dcum,inds(1),inds(2)) * ...
   subm(C,dcum,inds(2),inds(1)));

% ---------------------------------------------------------------------

function out = trCCC(C,p,dcum,inds,lambda)

out = trace(...
   subm(C,dcum,inds(1),inds(2)) * ...
   subm(C,dcum,inds(2),inds(3)) * ...
   subm(C,dcum,inds(3),inds(1)));

% ---------------------------------------------------------------------

function out = trCCCC(C,p,dcum,inds,lambda)

out = trace(...
   subm(C,dcum,inds(1),inds(2)) * ...
   subm(C,dcum,inds(2),inds(3)) * ...
   subm(C,dcum,inds(3),inds(4)) * ...
   subm(C,dcum,inds(4),inds(1)));

% ---------------------------------------------------------------------

function out = lambdaCClambda(C,p,dcum,inds,lambda)

runSum = 0;

for u = 1:p
   for v = 1:p
      runSum = runSum + subv(lambda,dcum,u)' * ...
         subm(C,dcum,u,inds(1))* ...
         subm(C,dcum,inds(1),v) * ...
         subv(lambda,dcum,v);
   end
end

out = runSum;

% ---------------------------------------------------------------------

function out = lambdaCCClambda(C,p,dcum,inds,lambda)

runSum = 0;

for u = 1:p
   for v = 1:p
      runSum = runSum + subv(lambda,dcum,u)' * ...
         subm(C,dcum,u,inds(1))* ...
         subm(C,dcum,inds(1),inds(2))* ...
         subm(C,dcum,inds(2),v) * ...
         subv(lambda,dcum,v);
   end
end

out = runSum;

% ---------------------------------------------------------------------

function out = lambdaCCCClambda(C,p,dcum,inds,lambda)

runSum = 0;

for u = 1:p
   for v = 1:p
      runSum = runSum + subv(lambda,dcum,u)' * ...
         subm(C,dcum,u,inds(1))* ...
         subm(C,dcum,inds(1),inds(2))* ...
         subm(C,dcum,inds(2),inds(3)) * ...
         subm(C,dcum,inds(3),v) * ...
         subv(lambda,dcum,v);
   end
end

out = runSum;

% ---------------------------------------------------------------------

function out = lambdaCCCCClambda(C,p,dcum,inds,lambda)

runSum = 0;

for u = 1:p
   for v = 1:p
      runSum = runSum + subv(lambda,dcum,u)' * ...
         subm(C,dcum,u,inds(1))* ...
         subm(C,dcum,inds(1),inds(2))* ...
         subm(C,dcum,inds(2),inds(3)) * ...
         subm(C,dcum,inds(3),inds(4)) * ...
         subm(C,dcum,inds(4),v) * ...
         subv(lambda,dcum,v);
   end
end

out = runSum;

