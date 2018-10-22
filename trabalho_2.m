  Aname = input(" Digite o nome do arquivo da matriz do S1: ");
  bname = input(" Digite o nome do arquivo do vetor do S1: ");
  
  A2name = input(" Digite o nome do arquivo da matriz do S2: ");
  b2name = input(" Digite o nome do arquivo do vetor do S2: ");
  
  tol = input("Digite a presisao desejada: ");
  
  A = load(Aname);
  b = load(bname);
  
  A2 = load(A2name);
  b2 = load(b2name);
  
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

  if (inc==0)
    %imprime vetor x
    disp("--------Pontos Amostrais--------");
   disp(sort(x));
  elseif(inc==1)
    disp("\nImpossivel Calcular o X");
  endif

%------------ Gauss Seidel--------------

function beta=sassenfeld(A2)
  [m n]=size(A2);
  beta=zeros(m,1);
  for i=1:m
   for j=1:i-1
     beta(i)=beta(i)+abs(A2(i,j))/abs(A2(i,i))*beta(j);
   endfor;
   for j=i+1:n
    beta(i)=beta(i)+abs(A2(i,j))/abs(A2(i,i));
   endfor;
  endfor;
  endfunction


%A = [6 0 -1 2 0; 0 8 -4 0 2; 0 2 -8 4 0; -2 0 0 -5 1; 0 1 0 9 14];
%b = [ 265; 100; 234; -400; 721];
Niter = 100;

  beta=sassenfeld(A2);
  if(max(beta)>1) 
    fprintf("Nao vai convergir!\n");
    break;
  endif
  n=length(A2);
  D=diag(diag(A2));
  L=D-tril(A2);
  U=D-triu(A2);
  invDL=inv(D-L);
  x2=zeros(n,1);
  xout=x2;
  oldx=ones(n,1);
  for i=1:Niter
    if (max(abs(x2-oldx))<tol)
      disp("--------Poluiçao--------")
      disp(x2)
      disp("")
      % LAGRANGE --------------------------------
      % 
      x0 = x;
      y0 = x2;
       xaux = median(x0);
          % Grau n do polinomio
          n = size(x0, 1); 
          y = 0;

          for i=1:n
              p = 1;
              for j=1:n
                  if j == i   
                      continue;
                  endif;
                  p *= (xaux-x0(j)) / (x0(i)-x0(j));
              endfor;
              y += y0(i) * p;   
          endfor;
          disp("--------Local da Cidade--------")
          disp(median(x0));
          disp("--------Poluiçao--------")
          disp(y);
       % ---------------------------------------
      return;
    else
      oldx=x2;
      x2=invDL*U*x2+invDL*b;
      xout=[xout x2];
    endif
  endfor
  disp ("Maximo de iteracoes atingido!!")
