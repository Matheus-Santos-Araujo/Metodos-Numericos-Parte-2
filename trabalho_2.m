% SUPRIME AVISOS -----------------------------------------
warning("off","all"); 

% PEGA E DEFINE VALORES -----------------------------------------
  Aname = input(" Digite o nome do arquivo da matriz do S1: ");
  bname = input(" Digite o nome do arquivo do vetor do S1: ");
  
  A2name = input(" Digite o nome do arquivo da matriz do S2: ");
  b2name = input(" Digite o nome do arquivo do vetor do S2: ");
  
  tol = input("Digite a precisao desejada: ");
  disp("")
  
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
    fprintf("O sistema 2 nao vai convergir, sassenfeld negado!\n");
    break;
  endif
  
  % -------------------------------------------------
  
  % METODO DE GAUSS SEIDEL --------------------------
  
  % Pega o tamanho de A
  n=length(A);
  % Pega a diagonal de A
  D=diag(diag(A));
  % Pega o L
  L=D-tril(A);
  % Pega o U
  U=D-triu(A);
  % Pega a inversa da diagonal menos L
  invDL=inv(D-L);
  % Inicia o X2
  x2=zeros(n,1);
  % Atualiza X2
  xout=x2;
  % Atualiza oldx
  oldx=ones(n,1);
  
  % Inicia o laço
  for i=1:Niter
    % Verifica se convergiu
    if (max(abs(x2-oldx))<tol || max(abs(x2-oldx))/max(abs(x2))<tol)
         % Verifica a conscistencia
         if (inc==0)
            disp("--Ponto---- Poluicao ---")
    % Associa os pontos amostrais ordenados com as poluicoes
    [sorted, indices] = sort(x);
    % Imprime o resultado
    disp([x(indices) x2(indices)]);
    disp("")
        % Se n~ao for consistente retorna erro
        elseif(inc==1)
    disp("\nImpossivel Calcular os pontos amostrais");
  endif
      
   % METODO DE LAGRANGE ----------------------------
      % Pega os pontos amostrais
      x0 = x;
      % Pega as poluicoes
      y0 = x2;
      % Pega o ponto da cidade
      xaux = median(x0);
      
          % Grau n do polinomio
          n = size(x0, 1);
          % Inicia y 
          y = 0;

          % Inicia o primeiro laço com grau do polinomio
          for i=1:n
              % Inicializa p
              p = 1;
              % Inicia o segundo laço com grau do polinomio
              for j=1:n
                  % Se j e i percorridos forem iguais continua
                  if j == i   
                      continue;
                  endif;
                  % Executa formula de Lagrange
                  p *= (xaux-x0(j)) / (x0(i)-x0(j));
              endfor;
              % Atualiza y
              y += y0(i) * p;   
          endfor;
          % Imprime o resultado 
          disp("---Cidade--- Poluicao ----")
         disp([ median(x0)  y])
       % ---------------------------------------
       
      return;
    else
      % Atualiza oldx
      oldx=x2;
      % Executa formula de Lagrange
      x2=invDL*U*x2+invDL*b;
      % Atualiza xout
      xout=[xout x2];
    endif
  endfor
  disp ("Maximo de iteracoes atingido!!")
