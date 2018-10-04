 A = [6 0 -1 2 0; 0 8 -4 0 2; 0 2 -8 4 0; -2 0 0 -5 1; 0 1 0 9 14];
 
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
 