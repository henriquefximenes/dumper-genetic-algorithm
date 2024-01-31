function populacao_cruzamento = cruzamento(n,populacao_torneio)
    for i=1:n
        x=grand(1,1,"unf",0,101);
        if x<taxa_cruzamento then
            pai = int(grand(1,1,"unf",1,n+1));
            mae = int(grand(1,1,"unf",1,n+1));
            while pai==mae
                pai = int(grand(1,1,"unf",1,n+1));
                mae = int(grand(1,1,"unf",1,n+1));
            end
            //sorteando as porcentagens que definem a dominância dos genes
            pct1 = rand();
            pct2 = rand();
            pct3 = rand();
            
            //cromossomos novos do filho
            populacao_cruzamento(i,1) = (1-pct1)*populacao_torneio(pai,1)+pct1*populacao_torneio(mae,1);
            populacao_cruzamento(i,2) = (1-pct2)*populacao_torneio(pai,2)+pct2*populacao_torneio(mae,2);
            populacao_cruzamento(i,3) = (1-pct3)*populacao_torneio(pai,3)+pct2*populacao_torneio(mae,3);
                 
            //se a probabilidade calculada por "x" for menor que taxa_cruzamento, cria-se um novo indivíduo.
        else
            populacao_cruzamento(i,1) = (rand()*0.05);
            populacao_cruzamento(i,2) = ((rand()/10)+0.9);
            populacao_cruzamento(i,3) = (rand()*0.15);
            while populacao_cruzamento(i,3)<0.001
                populacao_cruzamento(i,3)=(rand()*0.15);
            end
        end
    end
endfunction
