function mae = maerr(x, x_err)
%MAERR Calculate mean absolute error
%   Detailed explanation goes here

    err = x-x_err;
    a_err = abs(err);
    mae = mean(a_err);
end

