%creo vettore che nel nostro caso sara il vettore m composto dalle varie
%celle ocupate dai veicoli
a=[ 3 5 6 8 8 3 2 4 5 3 1 1 ];

% calcolo la lunghezza del vettore
n=length(a);

%primo ciclo for per ogni componenete del vettore
for i=1:n
    %secondo ciclo for fare la verifica 
    for k=1:n
        %creo un vetto x composta da 1 quando è soddisfatta l'uguaglianza e
        %0 quando non è soddisfatta
        x(:,k) = (a(i)==a(k));
    end
    %creo la matrice composta dai vettori x
    X(i,:)=x;
end

%elimino le righe di X che si ripetono
X1=unique(X,'rows');
disp(X)
disp(X1)
