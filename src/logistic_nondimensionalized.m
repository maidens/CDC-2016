function [X, C, I, signals] = logistic_nondimensionalized(inp, par)

    r = par(1);   % growth rate
    K = par(2);   % carrying capacity
    
    % state update 
    X{1} = inp.X{1}.*(1 - inp.U{1}) + r*(1 - inp.U{1}).*inp.X{1}.*(1 - inp.X{1}.*(1-inp.U{1})); 
    
    % sensitivity w.r.t. r 
    Df_r = (1 - inp.U{1}).*inp.X{1}.*(1 - inp.X{1}.*(1-inp.U{1})); 
    Df_xt = (1+r)*(1 - inp.U{1}) - 2*r*(1 - inp.U{1}).^2.*inp.X{1}; 
    X{2} = Df_r + Df_xt.*inp.X{2}; 
    
    
    % objective function 
    C{1} = - X{1}.*(1 - inp.U{1}) - (K./X{1}).^2.*(1 - inp.U{1}).*X{2}.^2; 
    % infeasibility
    I=0;

    % input signals 
    signals.U{1} = inp.U{1};

end