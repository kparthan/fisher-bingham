function out = calc_KHat_n(n,CHatInv,p,dcum,inds,lambda)

if length(inds) ~= n, error('Check number of elements of inds'); end

switch n
   case 1
      out = trace(subm(CHatInv,dcum,inds(1),inds(1))) + ...
         sumOverPerms('lambdaCClambda',CHatInv,p,dcum,inds,lambda);
   case 2
      out = 2* ...
         (sumOverPerms('trCC',CHatInv,p,dcum,inds,lambda) + ...
         sumOverPerms('lambdaCCClambda',CHatInv,p,dcum,inds,lambda));
   case 3
      out = 4* ...
         (sumOverPerms('trCCC',CHatInv,p,dcum,inds,lambda) + ...
         sumOverPerms('lambdaCCCClambda',CHatInv,p,dcum,inds,lambda));
   case 4
      out = 8* ...
         (sumOverPerms('trCCCC',CHatInv,p,dcum,inds,lambda) + ...
         sumOverPerms('lambdaCCCCClambda',CHatInv,p,dcum,inds,lambda));
end


