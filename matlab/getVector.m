%reads a vector
function x = getVector(fileName) 
    x = dlmread(fileName, ' ', 1, 0);
end