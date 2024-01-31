function [populacao,populacao_torneio] = selecao_torneio(populacao,n,ranking)
    i=1;
    while i<=n
        num1 = int(grand(1,1,"unf",1,n+1)); //função grand desta maneira não retorna o valor máximo
        num2 = int(grand(1,1,"unf",1,n+1));
        if num2~=num1 then
            if populacao(num1,4)<=populacao(num2,4) then
                for j=1:4
                    populacao_torneio(i,j) = populacao(num1,j);
                end
            else
                for j=1:4
                    populacao_torneio(i,j) = populacao(num2,j);
                end
            end
            i=i+1;
        else
            i=i;
        end
    end
endfunction
