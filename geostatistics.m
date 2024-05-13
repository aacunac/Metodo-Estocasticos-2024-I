x = [60.6,	36.2;
60.9,	17.7;];

t = [60.4, 37.2;
60.9,   36.1;
56.3,	34.6;
76.4,	23.3;
63.9,	21.5;
77.4,   25.3;
73.9,   34.8;
69.4,	35.6;];

vert_cuadro = [57.5, 22.5;
57.5,32.5;
72.5,22.5;
72.5,32.5];

%M
b = kriging_1(x,t,vert_cuadro,4);
%%
r = kriging_2(x,t,vert_cuadro,4);

%%
function [var,mejor_comb,menor_var] = kriging_1(datos_1, datos_variables,vertices,n)
       comb_posibles = generar_indices(datos_variables,n);
       comb_n = datos_1;
       for i_1 = 1:length(comb_posibles)
           for j_1 = 1:size(comb_posibles,2)
               a = comb_posibles(i_1,j_1);
               b = datos_variables(a,1);
               c = datos_variables(a,2);
               comb_n(length(datos_1)+j_1,1) = b;
               comb_n(length(datos_1)+j_1,2) = c;
           end
       m_d =calcular_matriz_distancias(comb_n);
       m_d_v = calcular_distancia_vertices(comb_n,vertices);
       m_s_v = calcula_semivarianza(m_d);
       m_s_v_v = calcula_semivarianza(m_d_v);
       v_s_v_p_v = zeros(size(m_s_v_v,2),1);
       for i_2 = 1:size(m_s_v_v,2)
           v_s_v_p_v(i_2,1) = mean(m_s_v_v(:,i_2));
       end
       m_kriging = matriz_final(m_s_v); 
       v_kriging = vertcat(v_s_v_p_v,1);
       pesos = m_kriging\v_kriging;
       multiplicacion = zeros(size(v_s_v_p_v,1),1);
       for i = 1:size(v_s_v_p_v,1)
           multiplicacion(i,1) = v_s_v_p_v(i,1)* pesos(i,1);
       end
       
       varianza_kringing = 0.18-sum(multiplicacion,1)-pesos(size(pesos,1),1);
       var(i_1) = varianza_kringing;
       comb_n = datos_1;
       end 
     [min_valor, indice] = min(var);
     menor_var = min_valor;
     mejor_comb = [comb_posibles(indice,:)];
     disp('Menor varianza posible con 4 puntos extras que se puede obtener es:')
     disp(menor_var);
     disp('Empleando la convinación de puntos:')
     disp(mejor_comb)
end


%%
function [varianzas,mejor_comb,menor_varianza] = kriging_2(datos_1, datos_variables,vertices,n)
      datos_adicionales = 1:length(datos_variables);
      mejor_comb = [];
      for i = 1:n
          for j = 1:length(datos_adicionales)
             combo = [datos_1;datos_variables(datos_adicionales(j),:)];
             m_d = calcular_matriz_distancias(combo);
             m_d_v = calcular_distancia_vertices(combo,vertices);
             m_s_v = calcula_semivarianza(m_d);
             m_s_v_v = calcula_semivarianza(m_d_v);
             v_s_v_p_v = zeros(size(m_s_v_v,2),1);
             for i_2 = 1:size(m_s_v_v,2)
                v_s_v_p_v(i_2,1) = mean(m_s_v_v(:,i_2));
             end
             m_kriging = matriz_final(m_s_v); 
             v_kriging = vertcat(v_s_v_p_v,1);
             pesos = m_kriging\v_kriging;
             multiplicacion = zeros(size(v_s_v_p_v,1),1);
             for i_1 = 1:size(v_s_v_p_v,1)
                multiplicacion(i_1,1) = v_s_v_p_v(i_1,1)* pesos(i_1,1);
             end
             varianza_kringing = 0.18-sum(multiplicacion,1)-pesos(size(pesos,1),1);
             varianzas(j) = varianza_kringing;
          end
          indice_minimo = find(varianzas == min(varianzas));
          if i == n
              menor_varianza = min(varianzas);
          end
          varianzas = [];
          datos_1 = [datos_1;datos_variables(datos_adicionales(indice_minimo),:)];
          mejor_comb(i) = datos_adicionales(indice_minimo);
          datos_adicionales(indice_minimo) = [];
      end
     disp('Menor varianza posible con 4 puntos extras que se puede obtener es:')
     disp(menor_varianza);
     disp('Empleando la convinación de puntos:')
     disp(mejor_comb)
end

%%
function gamma = variograma(h)
    % Definir los intervalos y las expresiones correspondientes
    if h == 0
        gamma = 0; 
    elseif h > 0 && h < 40
        gamma = 0.08+0.1*((3/2)*(h/40)-(1/2)*(h/40)^3); % Expresión para 0 <= x < 1
    elseif h >= 40
        gamma =0.18 ; % Expresión para x >= 1
    end
end

%%
function combinaciones = generar_indices(puntos,n)
        N = length(puntos);
        if n <= 0 || n > N
        error('El valor de n debe ser un entero positivo menor o igual que la longitud del vector de datos.');
        end

        indice = combnk(1:N, n);
        combinaciones = flipud(indice);
        %disp(combinaciones)
end   

%%
function matriz_distancias = calcular_matriz_distancias(datos)
    % Número de puntos
    num_puntos = size(datos, 1);

    % Matriz de distancias
    matriz_distancias = zeros(num_puntos);

    % Calcular distancias entre cada par de puntos
    for i = 1:num_puntos
        for j = 1:num_puntos
            % Calcular la distancia euclidiana entre los puntos i y j
            matriz_distancias(i, j) = sqrt((datos(i, 1) - datos(j, 1))^2 + (datos(i, 2) - datos(j, 2))^2);
        end
    end

    % Mostrar la matriz de distancias
    %disp('Matriz de Distancias:')
    %disp(matriz_distancias);
end

%%
function matriz_distancia_vertices = calcular_distancia_vertices(puntos,vertices)
         num_vertices = size(vertices,1);
         num_puntos = size(puntos,1);
         matriz_distancia_vertices = zeros(num_vertices,num_puntos);

         for i = 1:num_vertices
             for j = 1:num_puntos
                 matriz_distancia_vertices(i,j) = sqrt((vertices(i,1)-puntos(j,1))^2+(vertices(i,2)-puntos(j,2))^2);
             end
         end
%disp('Matriz de Distancias a vertices:')         
%disp(matriz_distancia_vertices);         
end

%%
function matriz_semivarianza = calcula_semivarianza(x)
 for i = 1:size(x,1)
    for j = 1:size(x,2)
        matriz_semivarianza(i,j) = variograma(x(i,j));
    end
 end
end

%%
function matriz_kriging = matriz_final(matriz_semivarianza)

for i = 1:size(matriz_semivarianza,1)+1
    for j = 1:size(matriz_semivarianza,2)+1
        if i <= size(matriz_semivarianza,1) && j<= size(matriz_semivarianza,2)
            matriz_kriging(i,j) = matriz_semivarianza(i,j);
        elseif i <= size(matriz_semivarianza,1) && j > size(matriz_semivarianza,2)
            matriz_kriging(i,j) = 1;
        elseif i > size(matriz_semivarianza,1) && j <= size(matriz_semivarianza,2)
            matriz_kriging(i,j) = 1;
        elseif i > size(matriz_semivarianza,1) && j > size(matriz_semivarianza,2)
            matriz_kriging(i,j) = 0;
        end
    end
end
end

      