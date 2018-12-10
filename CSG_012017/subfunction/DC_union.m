function xy = DC_union(x,y)

% union function for cells

% check dimension 
if size(x) == size(y)
    xy = cell(size(x));
    for i1 = 1 : size(x,1)
        for i2 = 1 : size(x,2)
            xy{i1,i2} = union(x{i1,i2},y{i1,i2});
        end
    end
else 
    error('Dimension error: the two cells must have the same dimension')
end