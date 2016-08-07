%reads a vector
function x = getVector(fileName) 
    x = dlmread(fileName, ' ', 3, 0);
end