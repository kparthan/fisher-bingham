function out = calc_KHat_n(n,CHatInv,p,q,ind,lambda)

% NOTE: ind specifies the index required [need numel(ind) = n].
%    If isempty(ind) then the whole n-dimensional array is returned.
%
% For n=3 and n=4, this only calculates the "upper-diagonal" elements -
% other elements can be deduced by symmetry (i.e. invariance w.r.t. 
% subscript ordering)


if isempty(ind)
   m = p*(p+1)/2;
   switch n
      case 1
         out = NaN(m,1);
         for k1 = 1:m
            out(k1) = calc_KHat_n(n,CHatInv,p,q,k1,lambda);
         end
      case 2
         out = NaN(m,m);
         for k1 = 1:m
            for k2 = 1:m
               out(k1,k2) = calc_KHat_n(n,CHatInv,p,q,[k1 k2],lambda);
            end
         end
      case 3  % calculates only the "upper-diagonal" elements
         out = NaN(m,m,m);
         for k1 = 1:m
            for k2 = k1:m
               for k3 = k2:m
                  out(k1,k2,k3) = calc_KHat_n(n,CHatInv,p,q,[k1 k2 k3],lambda);
               end
            end
         end
      case 4  % calculates only the "upper-diagonal" elements
         out = NaN(m,m,m,m);
         for k1 = 1:m
            for k2 = k1:m
               for k3 = k2:m
                  for k4 = k3:m
                     out(k1,k2,k3,k4) = calc_KHat_n(n,CHatInv,p,q,[k1 k2 k3 k4],lambda);
                  end
               end
            end
         end
   end
elseif numel(ind) ~= n
   error('ind must be a vector of mth n');
else
   
   s = zeros(n,1);
   r = zeros(n,1);
   
   for k = 1:n
      [r(k), s(k)] = get_r_s(ind(k));
   end
   
   switch n
      
      case 1   %%%%%%%%%%%%%% CASE 1 %%%%%%%%%%%%%%
         t1_list = [r(1) s(1); s(1) r(1)];
         term1 = 0;
         for i1 = [1 2]
            t1 = t1_list(i1,1); t1c = t1_list(i1,2);
            term1 = term1 + trace(subm(CHatInv,q,t1,t1c));
         end
         term2 = 0;
         for u = 1:p
            for v = 1:p
               CCsum = 0;
               for i1 = [1 2]
                  t1 = t1_list(i1,1);
                  t1c = t1_list(i1,2);
                  CCsum = CCsum + subm(CHatInv,q,u,t1)*subm(CHatInv,q,t1c,v);
               end
               term2 = term2 + subv(lambda,q,u)'*CCsum*subv(lambda,q,v);
            end
         end
         out = 1/2*(term1 + term2);
         
      case 2   %%%%%%%%%%%%%% CASE 2 %%%%%%%%%%%%%%
         t1_list = [r(1) s(1); s(1) r(1)];
         t2_list = [r(2) s(2); s(2) r(2)];
         term1 = 0;
         for i1 = [1 2]
            for i2 = [1 2]
               t1 = t1_list(i1,1); t1c = t1_list(i1,2);
               t2 = t2_list(i2,1); t2c = t2_list(i2,2);
               term1 = term1 + trace(...
                  subm(CHatInv,q,t1,t2)*subm(CHatInv,q,t2c,t1c));
            end
         end
         term2 = 0;
         for u = 1:p
            for v = 1:p
               CCCperm_sum_over_t = 0;
               for i1 = [1 2]
                  for i2 = [1 2]
                     t = zeros(2,2);
                     t(1,:) = t1_list(i1,:);
                     t(2,:) = t2_list(i2,:);
                     CCCperm_sum_over_t = CCCperm_sum_over_t + ...
                        sumOverPermsTerm2(3,CHatInv,p,q,t,u,v);
                  end
               end
               term2 = term2 + ...
                  subv(lambda,q,u)'*CCCperm_sum_over_t*subv(lambda,q,v);
            end
         end
         out = 1/2*(term1 + term2);
         
      case 3   %%%%%%%%%%%%%% CASE 3 %%%%%%%%%%%%%%
         t1_list = [r(1) s(1); s(1) r(1)];
         t2_list = [r(2) s(2); s(2) r(2)];
         t3_list = [r(3) s(3); s(3) r(3)];
         term1 = 0;
         for i1 = [1 2]
            for i2 = [1 2]
               for i3 = [1 2]
                  t = zeros(3,2);
                  t(1,:) = t1_list(i1,:);
                  t(2,:) = t2_list(i2,:);
                  t(3,:) = t3_list(i3,:);
                  term1 = term1 + sumOverPermsTerm1(3,CHatInv,p,q,t);
               end
            end
         end
         term2 = 0;
         for u = 1:p
            for v = 1:p
               CCCCperm_sum_over_t = 0;
               for i1 = [1 2]
                  for i2 = [1 2]
                     for i3 = [1 2]
                        t = zeros(3,2);
                        t(1,:) = t1_list(i1,:);
                        t(2,:) = t2_list(i2,:);
                        t(3,:) = t3_list(i3,:);
                        CCCCperm_sum_over_t = CCCCperm_sum_over_t + ...
                           sumOverPermsTerm2(4,CHatInv,p,q,t,u,v);
                     end
                  end
               end
               term2 = term2 + ...
                  subv(lambda,q,u)'*CCCCperm_sum_over_t*subv(lambda,q,v);
            end
         end
         out = 1/2*(term1 + term2);
         
      case 4   %%%%%%%%%%%%%% CASE 4 %%%%%%%%%%%%%%
         t1_list = [r(1) s(1); s(1) r(1)];
         t2_list = [r(2) s(2); s(2) r(2)];
         t3_list = [r(3) s(3); s(3) r(3)];
         t4_list = [r(4) s(4); s(4) r(4)];
         term1 = 0;
         for i1 = [1 2]
            for i2 = [1 2]
               for i3 = [1 2]
                  for i4 = [1 2]
                     t = zeros(4,2);
                     t(1,:) = t1_list(i1,:);
                     t(2,:) = t2_list(i2,:);
                     t(3,:) = t3_list(i3,:);
                     t(4,:) = t4_list(i4,:);
                     term1 = term1 + sumOverPermsTerm1(4,CHatInv,p,q,t);
                  end
               end
            end
         end
         term2 = 0;
         for u = 1:p
            for v = 1:p
               CCCCCperm_sum_over_t = 0;
               for i1 = [1 2]
                  for i2 = [1 2]
                     for i3 = [1 2]
                        for i4 = [1 2]
                           t = zeros(4,2);
                           t(1,:) = t1_list(i1,:);
                           t(2,:) = t2_list(i2,:);
                           t(3,:) = t3_list(i3,:);
                           t(4,:) = t4_list(i4,:);
                           CCCCCperm_sum_over_t = CCCCCperm_sum_over_t + ...
                              sumOverPermsTerm2(5,CHatInv,p,q,t,u,v);
                        end
                     end
                  end
               end
               term2 = term2 + ...
                  subv(lambda,q,u)'*CCCCCperm_sum_over_t*subv(lambda,q,v);
            end
         end
         out = 1/2*(term1 + term2);
         
   end

%    t1 = term1/2
%    t2 = term2/2
   
end


% ---------------------------------------------------------------------

function out = sumOverPermsTerm1(termType,CHatInv,p,q,t)

% This is for 'trCCC', 'trCCCC', specified by termType:
%   'trCCC'   - termType = 3
%   'trCCCC'  - termType = 4

n = size(t,1);
inds = 1:n;

% permute all but first index
permedInds = [...
   repmat(inds(1),[factorial(length(inds)-1) 1]) ...
   perms(inds(2:end))];
numPerms = size(permedInds,1);

runSum = 0;

for ind_ = 1:numPerms
   t_tmp = t(permedInds(ind_,:),:);
   switch termType
      case 3 % 'trCCC'
         runSum = runSum + trCCC(CHatInv,p,q,t_tmp);
      case 4 % 'trCCCC'
         runSum = runSum + trCCCC(CHatInv,p,q,t_tmp);
   end
end

out = runSum;


% ---------------------------------------------------------------------

function out = sumOverPermsTerm2(termType,CHatInv,p,q,t,u,v)

% This is for 'CCC', 'CCCC', and 'CCCCC' terms, specified by termType:
%   'CCC'   - termType = 3
%   'CCCC'  - termType = 4
%   'CCCCC' - termType = 5

n = size(t,1);
inds = 1:n;

permedInds = flipud(perms(inds));
numPerms = size(permedInds,1);

runSum = 0;

for ind_ = 1:numPerms
   t_tmp = t(permedInds(ind_,:),:);
   switch termType
      case 3 % 'CCC'
         runSum = runSum + CCC(CHatInv,p,q,t_tmp,u,v);
      case 4 % 'CCCC'
         runSum = runSum + CCCC(CHatInv,p,q,t_tmp,u,v);
      case 5 % 'CCCCC'
         runSum = runSum + CCCCC(CHatInv,p,q,t_tmp,u,v);
   end
end

out = runSum;


% ---------------------------------------------------------------------

function out = trCCC(C,p,q,t) %#ok<*INUSD>

% display([t1 t2; t2c t3; t3c t4; t4c t1c]);

t1 =  t(1,1);
t1c = t(1,2);
t2 =  t(2,1);
t2c = t(2,2);
t3 =  t(3,1);
t3c = t(3,2);

out = trace(...
   subm(C,q,t1,t2) * ...
   subm(C,q,t2c,t3) * ...
   subm(C,q,t3c,t1c));


% ---------------------------------------------------------------------

function out = trCCCC(C,p,q,t) %#ok<*INUSL,*INUSD>

% display([t1 t2; t2c t3; t3c t4; t4c t1c]);

t1 =  t(1,1);
t1c = t(1,2);
t2 =  t(2,1);
t2c = t(2,2);
t3 =  t(3,1);
t3c = t(3,2);
t4 =  t(4,1);
t4c = t(4,2);

out = trace(...
   subm(C,q,t1,t2) * ...
   subm(C,q,t2c,t3) * ...
   subm(C,q,t3c,t4) * ...
   subm(C,q,t4c,t1c));

a= 1;


% ---------------------------------------------------------------------

function out = CCC(C,p,q,t,u,v)

t1 =  t(1,1);
t1c = t(1,2);
t2 =  t(2,1);
t2c = t(2,2);

out = ...
   subm(C,q,u,t1) * ...
   subm(C,q,t1c,t2) * ...
   subm(C,q,t2c,v);

% ---------------------------------------------------------------------

function out = CCCC(C,p,q,t,u,v)

t1 =  t(1,1);
t1c = t(1,2);
t2 =  t(2,1);
t2c = t(2,2);
t3 =  t(3,1);
t3c = t(3,2);

out = ...
   subm(C,q,u,t1) * ...
   subm(C,q,t1c,t2) * ...
   subm(C,q,t2c,t3) * ...
   subm(C,q,t3c,v);

% ---------------------------------------------------------------------

function out = CCCCC(C,p,q,t,u,v)

t1 =  t(1,1);
t1c = t(1,2);
t2 =  t(2,1);
t2c = t(2,2);
t3 =  t(3,1);
t3c = t(3,2);
t4 =  t(4,1);
t4c = t(4,2);

out = ...
   subm(C,q,u,t1) * ...
   subm(C,q,t1c,t2) * ...
   subm(C,q,t2c,t3) * ...
   subm(C,q,t3c,t4) * ...
   subm(C,q,t4c,v);


