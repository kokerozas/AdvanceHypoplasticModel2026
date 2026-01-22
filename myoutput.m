function stop = myoutput(~, optimValues, ~)
if optimValues.fval < 0.03
    stop = true;
else
    stop = false;
end
end