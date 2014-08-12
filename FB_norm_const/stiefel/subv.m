function out = subv(v,d,i)

% returns the i'th sub-vector of v (each sub-vector has length d)

out = v(1+(i-1)*d:i*d);