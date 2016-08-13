%reads a vector
function x = getVector(fileName) 
    x = dlmread(fileName, ' ', 2, 0);
end