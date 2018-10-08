  opcao = input("\n Digite o numero de acordo com a opcao desejada: \n 1) LU.\n 2) Gauss-Seidel.\n 3) Lagrange.\n");

if( opcao == 1 )

  A = load('S1Matriz.txt');
  b = load('S1Vetor.txt');

 %A = [6 0 -1 2 0; 0 8 -4 0 2; 0 2 -8 4 0; -2 0 0 -5 1; 0 1 0 9 14];
 %b = [ 265; 100; 234; -400; 721];
 %pega o tamanho
  [numLin numCol] = size(A);
  % gera uma matriz identidade
  L = eye(numLin, numLin);

  % Calcula matriz L
  for k=1:numLin
     if(A(k,k) ==0)
       L(k+1:numLin,k) = 0;
     elseif(A(k,k)!=0)
       L(k+1:numLin,k)=A(k+1:numLin,k)/A(k,k);
     endif
     
      % Calcula matriz U
      for j=k+1:numLin
        A(j,:)=A(j,:)-(L(j,k)*A(k,:));
      end
  end
  % iguala U a A
  U = A;

  %calcula a conscistencia
  for i=1:numCol
    if (U(i,i)==0)
      inc=1;
    else
      inc=0;
    endif
  endfor
  if(inc==0)
  % Calcula resultado de x
  y = inv(L)*b;
  x = inv(U)*y;
  endif

  %calcula o determinante
  detA=1;
  for i=1:numLin
    detA = detA*U(i,i);
  endfor

  %calcula a conscistencia
  for i=1:numCol
    if (U(i,:)==0)
      inc=1;
    else
      inc=0;
    endif
  endfor

  %Imprime matriz L
  disp("\nMatriz L:")
  disp(L);

  %Imprime matriz U
  disp("\nMatriz U:")
  disp(U);

  if (inc==0)
    %imprime vetor x
    disp("Vetor x");
    disp(x);
  elseif(inc==1)
    disp("\nImpossivel Calcular o X");
  endif

elseif( opcao == 2 )

  A = load("S2Matriz.txt");
  b = load("S2Vetor.txt");

function beta=sassenfeld(A)
  [m n]=size(A);
  beta=zeros(m,1);
  for i=1:m
   for j=1:i-1
     beta(i)=beta(i)+abs(A(i,j))/abs(A(i,i))*beta(j);
   endfor;
   for j=i+1:n
    beta(i)=beta(i)+abs(A(i,j))/abs(A(i,i));
   endfor;
  endfor;
  endfunction


%A = [6 0 -1 2 0; 0 8 -4 0 2; 0 2 -8 4 0; -2 0 0 -5 1; 0 1 0 9 14];
%b = [ 265; 100; 234; -400; 721];
Niter = 100;
tol = 1e-4;
  beta=sassenfeld(A);
  if(max(beta)>1) 
    fprintf("Nao vai convergir!\n");
    break;
  endif
  disp(beta);
  n=length(A);
  D=diag(diag(A));
  L=D-tril(A);
  U=D-triu(A);
  invDL=inv(D-L);
  x2=zeros(n,1);
  xout=x2;
  oldx=ones(n,1);
  for i=1:Niter
    if (max(abs(x2-oldx))<tol)
      disp("")
      disp ("Metodo completado com sucesso!")
      disp("--------Numero de Iteracoes--------")
      disp(i)
      disp("--------Resultado--------")
      disp(x2)
      disp("")
      return;
    else
      oldx=x2;
      x2=invDL*U*x2+invDL*b;
      xout=[xout x2];
    endif
  endfor
  disp ("Maximo de iteracoes atingido!!")

elseif( opcao == 3 )

x = 2;
x0 = [23.0000; 7.0000; 30.0000; -50.0000; -12.0000];
y0 = [2.0000; 24.0000; 72.0000; 9.0000; 57.0000];

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

else 
disp("Opcao invalida\n");
endif
