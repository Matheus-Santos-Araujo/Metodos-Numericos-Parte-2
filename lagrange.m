x = 2;
x0 = [2; 3; 4; 5];
y0 = [2; 6; 24; 120];

    n = size(x0, 1); 
    y = 0;

    for i=1:n
        p = 1;
        for j=1:n
            if j == i   
                continue;
            endif;
            p *= (x-x0(j)) / (x0(i)-x0(j));
        endfor;
        y += y0(i) * p;   
    endfor;
    disp("--------Resultado--------")
    disp(y);
