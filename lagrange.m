x = 2;
x0 = [2; 3; 4; 5];
y0 = [2; 6; 24; 120];

    % x0 - vector containing inputs (x values)
    % y0 - vector containing outputs (results for these x values
    % x - value you want to compute, for interpolation
    % y - computed value

    n = size(x0, 1); 
    y = 0;

    for i=1:n
        p = 1;
        for j=1:n
            if j == i   % avoiding fancy division by 0
                continue;
            endif;
            p *= (x-x0(j)) / (x0(i)-x0(j));
        endfor;
        y += y0(i) * p;   
    endfor;
    disp("--------Resultado--------")
    disp(y);
