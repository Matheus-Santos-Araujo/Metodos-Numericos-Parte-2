% PEGA E DEFINE VALORES -----------------------------------------
  Aname = input(" Digite o nome do arquivo da matriz do S1: ");
  bname = input(" Digite o nome do arquivo do vetor do S1: ");
  
  A2name = input(" Digite o nome do arquivo da matriz do S2: ");
  b2name = input(" Digite o nome do arquivo do vetor do S2: ");
  
  tol = input("Digite a presisao desejada: ");
 
  A1 = load(Aname);
  b1 = load(bname);
  
  A = load(A2name);
  b = load(b2name);
  
  Niter = 100;
 
  % --------------------------------------------------------
  
  % DECOMPOSIÇAO LU ----------------------------------------
  
 %pega o tamanho
  [numLin numCol] = size(A1);
  % gera uma matriz identidade
  L = eye(numLin, numLin);

  % Calcula matriz L
  for k=1:numLin
     if(A1(k,k) ==0)
       L(k+1:numLin,k) = 0;
     elseif(A1(k,k)!=0)
       L(k+1:numLin,k)=A1(k+1:numLin,k)/A1(k,k);
     endif
     
      % Calcula matriz U
      for j=k+1:numLin
        A1(j,:)=A1(j,:)-(L(j,k)*A1(k,:));
      end
  end
  % iguala U a A
  U = A1;

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
  y = inv(L)*b1;
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

  % -------------------------------------------------
  
  % VERIFICAÇAO DE SASSENFELD -----------------------
  
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
  
   beta=sassenfeld(A);
  if(max(beta)>1) 
    fprintf("Nao vai convergir!\n");
    break;
  endif
  
  % -------------------------------------------------
  
  % METODO DE GAUSS SEIDEL --------------------------
  
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
      disp("--------Poluiçao--------")
      disp(sort(x2))
      disp("")
      
   % METODO DE LAGRANGE ----------------------------
     
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
