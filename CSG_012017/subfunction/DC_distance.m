function d = DC_distance(pos)

% DC_distance returns a matrix d(N,N) with the element d(i,j) corresponding to the distance between the ith and jth point in pos 
% check the dimensionality: pos is a (N,M) matrix composed of N points and M
% coordinates
% distance are given in cm

d = zeros(size(pos,1));
M = ones(size(pos,1),1);
for k = 1 : size(pos,1)
    d(k,:) = sqrt(sum((M*pos(k,:) - pos).^2,2));
end    
