function populacao = mutacao(n,populacao,taxa_mutacao,ranking)
    for i=ranking+1:n
        z = int(grand(1,1,"unf",0,n+1));
        if z<taxa_mutacao then
            individuo = int(grand(1,1,"unf",ranking+1,n+1));
            gene = int(grand(1,1,"unf",1,4));
            if gene==1 then
                populacao(individuo,gene)=(rand()*0.05);
            elseif gene==2 then
                populacao(individuo,gene)=((rand()/10)+0.9);
            else
                populacao(individuo,gene)=(rand()*0.15);
                while populacao(individuo,gene)<0.001
                    populacao(individuo,gene)=(rand()*0.15);
                end
            end
        end
    end
endfunction
