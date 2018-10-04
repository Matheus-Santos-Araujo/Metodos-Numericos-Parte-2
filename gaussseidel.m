A = [6 0 -1 2 0; 0 8 -4 0 2; 0 2 -8 4 0; -2 0 0 -5 1; 0 1 0 9 14];
b = [ 265; 100; 234; -400; 721];
Niter = 100;
tol = 1e-4;

  n=length(A);
  D=diag(diag(A));
  L=D-tril(A);
  U=D-triu(A);
  invDL=inv(D-L);
  x=zeros(n,1);
  xout=x;
  oldx=ones(n,1);
  for i=1:Niter
    if (max(abs(x-oldx))<tol)
      disp("")
      disp ("Metodo completado com sucesso!")
      disp("--------Numero de Iteracoes--------")
      disp(i)
      disp("--------Resultado--------")
      disp(x)
      disp("")
      return;
    else
      oldx=x;
      x=invDL*U*x+invDL*b;
      xout=[xout x];
    endif
  endfor
  disp ("Maximo de iteracoes atingido!!")
