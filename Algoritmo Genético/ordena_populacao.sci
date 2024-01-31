function populacao = ordena_populacao(populacao, n)
    x=0;
    while(x==0)
        x=1;
        for i=1:n-1
            if (populacao(i+1,4)<populacao(i,4)) then
                x=0;
                aux1(i,1) = populacao(i,1);
                aux1(i,2) = populacao(i,2);
                aux1(i,3) = populacao(i,3);
                aux1(i,4) = populacao(i,4);
                populacao(i,1) = populacao(i+1,1);
                populacao(i,2) = populacao(i+1,2);
                populacao(i,3) = populacao(i+1,3);
                populacao(i,4) = populacao(i+1,4);
                populacao(i+1,1) = aux1(i,1);
                populacao(i+1,2) = aux1(i,2);
                populacao(i+1,3) = aux1(i,3);
                populacao(i+1,4) = aux1(i,4);
            end
        end
    end
endfunction
