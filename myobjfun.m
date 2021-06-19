function [residual, SSE] = myobjfun(rates, te, ye, tin, tfin , Xiniz)

    % [t,y] = ode45(@(t,X) enzReaction(t,X, rates), [tin tfin], Xiniz);
    [t,y]  = ode15s(@(t,X) enzReaction(t,X, rates), [tin tfin], Xiniz);
    ysimul = interp1(t,y,te);

    %unweighted residual
    residual = (ysimul-ye);

    %weighted residual
    % residual= (ysimul-ye)./ye;
    % Caution: there are zeros!

    %residual as vector
    residual = residual(:);

    % sum of squares 
    SSE = sum(residual.^2);

end