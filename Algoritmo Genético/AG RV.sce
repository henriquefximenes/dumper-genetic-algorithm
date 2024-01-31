//Início do Programa


//Comandos para limpar a tela do console:

clc;
clear;


format(8); //mostrar as 5 casas depois da virgula


//----------Parâmetros dos operadores genéticos----------

taxa_cruzamento = 60;
taxa_mutacao = 30;         //Taxas para configurar as probabilidades de cruzamento e mutação.
max_geracao = 3000;       //numero de gerações desejadas.
ranking = 8;               //quantidade de indivíduos selecionados como melhores por geração
n = 60;                    //quantidade de individuos por população
tol = 1E-6;
//-------------------------------------------------------


//-----------Parâmetros de calibração da função objetivo

q1 = 10;
q2 = 7; 
q3 = 40;
q4 = 1.07;
//-------------------------------------------------------


//---------Entrada de dados da Laje-----------

L = 9;                          //comprimento da laje
mp = 13576.6;                   //massa da laje
Ecdin = 28560*10^6;             //módulo de elasticidade dinâmico
Ieq = 0.0051321;                //Inercia equivalente
amort_estr = 0.04;              //amortecimento da estrutura
alim = 0.049;                   //aceleração limite
FS = 1.4;                       //Fator de segurança para conforto à vibração
//--------------------------------------------


//---------Dados de excitação------------

w = 33.3291881;                //frequencia de excitação
m=80;                          //Massa da pessoa
//---------------------------------------


//--------Aceleração da gravidade--------

acel_g = 9.80665;              //gravidade
//---------------------------------------


//-------Cálculos Iniciais---------------

//carregamento
Pp = m*acel_g;                //peso da pessoa
p3 = 0.12*Pp;                  //peso no terceiro harmonico

//propriedades da laje
kp = 2*((48*Ecdin*Ieq)/(L^3)); //rigidez da estrutura
wnp = sqrt(kp/mp);             //frequencia natural da estrutura
cp = 2*amort_estr*sqrt(kp*mp); //coeficiente de amortecimento da estrutura

//limites de vibração
uest = p3/kp;
ulim = alim/(wnp^2);
Glim = (ulim/uest)/q4;      //amplificação limite da estrutura
//---------------------------------------


//---------Declarando funções utilizadas pelo AG---------

exec("ordena_populacao.sci");
exec("selecao_torneio.sci");
exec("mutacao.sci");
exec("cruzamento.sci");
//-------------------------------------------------------


//--------------------------------------------------------------------
//========================Início da programação=======================
//--------------------------------------------------------------------


//----------Definindo a quantidade de indivíduos na população------------



//Fase de inicialização do programa: 
//Onde são criados os indivíduos que serão utilizados como a primeira população
//a qual será submetida aos vários processos sequenciais de aprimoramento com a
//função objetivo, elitismo, seleção, cruzamento e mutação.


//Criando os indivíduos para a população inicial:


for i=1:n
    for j=1:3
        if j==1 then
            populacao(i,j)=(rand()*0.05);
        elseif j==2 then
            populacao(i,j)=((rand()/10)+0.9);
        else
            populacao(i,j)=(rand()*0.15);
            while populacao(i,j)<0.001
                populacao(i,j)=(rand()*0.15);
            end
        end
    end
end


//Laço de repetição while baseado no fato em que o programa repetirá os processos seguintes
//de geração (número 1 inicialmente) até o número máximo de gerações informadas pelo usuário
//(max_geracao).


//Necessário declaração da variável Gr(1,4) para o primeiro laço.

Gr(1,4) = 1;
geracao = 1;               //Para o laço de repetição.

while geracao<=max_geracao
    
    
    if abs(Gr(1,4)-Glim)>tol then
        
        
        
        //-----------------******----------------------
        //Calculo de r1, r2 e alguns parametros do AMS:
        //-----------------******----------------------
    
    
        for i=1:n
            ma(i) = mp*populacao(i,1);              //massa do AMS
            wa(i) = wnp*populacao(i,2);             //frequencia do AMS
            ka(i) = ma(i)*(wa(i)^2);                //rigidez do AMS
            m = [mp,0; 0,ma(i)];                    //matriz de massa do sistema
            k = [kp+ka(i),-ka(i);-ka(i),ka(i)];     //matriz de rigidez do sistema
            A = (m^-1)*k;                           //matriz caracteristica
            [Ve,Va] = spec(A);                      //determina os autovetores (Ve) e autovalores (Va) da matriz A
            r(i,1) = sqrt(real(Va(1,1)))/wnp;       //matriz que armazena os valores de r1, r2 e r3
            r(i,2) = sqrt(real(Va(2,2)))/wnp;
            r(i,3) = 1;
        end
    
    
        //Salvando os valores de massa, frequência e rigidez do melhor dimensionamento
        //encontrado pelo AG na última geração.
    
        if geracao==max_geracao then
            MA = ma(1);
            WA = wa(1);
            KA = ka(1);
        end
    
    
        //--------------------********----------------------
        //Calculo dos Gr's para utilizar na função objetivo:
        //--------------------********----------------------
    
    
        for i=1:n
            for j=1:3
                r2 = r(i,j);
                B = wa(i)/wnp;   //beta da função de amplificação do artigo.
                A1 = (1-(r2^2/B^2))^2;
                A2 = 4*((populacao(i,3)*(r2/B))^2);
                A3 = ((r2^4/B^2)-(((4*populacao(i,3)*amort_estr)/B)+(1/B^2)+(populacao(i,1)+1))*(r2^2)+1)^2;
                A4 = 4*(r2*(amort_estr+(populacao(i,3)/B))-((r2^3/B)*(populacao(i,3)+(amort_estr/B)))-((r2^3/B)*populacao(i,3)*populacao(i,1)))^2;
                Gr(i,j) = sqrt((A1+A2)/(A3+A4)); //matriz com calculo dos r's para colocar na função objetivo
            end
        end
    
    
        //--------------------------********-----------------------------
        //Definindo Grmax para a função objetivo e aplicando a penalidade
        //--------------------------********-----------------------------
    
    
        for i=1:n
            if Gr(i,1)>Gr(i,2) & Gr(i,1)>Gr(i,3) then
                Gr(i,4) = Gr(i,1);
            elseif Gr(i,2)>Gr(i,1) & Gr(i,2)>Gr(i,3) then
                Gr(i,4) = Gr(i,2);
            else                   //penalidade caso G3 maior que G2 e G1
                Gr(i,3)=1000;
            end
        end
    
        //----------------------********----------------------
        //Corrigindo penalidade aplicada se G3 fosse a máxima.
        //----------------------********----------------------
    
    
        //Obs: Penalidade criada pois foi observado que ouveram muitos indivíduos que
        //sofriam a penalidade de que se a amplificação calculada com a variável r3 
        //(sempre igual a 1) fosse dada como a amplificação máxima, Gr(i,3) receberia
        //o valor de 1000 para resultar em valores altos e serem descartadas durante o
        //processo de otimização.
    
        for i=1:n
            if Gr(i,3)==1000 then
                while Gr(i,3)==1000 then
                    for y=1:3
                        if y==1 then
                            populacao(i,y) = (rand()*0.05);
                        elseif y==2 then
                            populacao(i,y) = ((rand()/10)+0.9);
                        else
                            populacao(i,y) = (rand()*0.15);
                        end
                    end
                    ma(i) = mp*populacao(i,1);
                    wa(i) = wnp*populacao(i,2);
                    ka(i) = ma(i)*(wa(i)^2);
                    m = [mp,0; 0,ma(i)];           //matriz de massa do sistema
                    k = [kp+ka(i),-ka(i);-ka(i),ka(i)];     //matriz de rigidez do sistema
                    A = (m^-1)*k;               //matriz caracteristica
                    [Ve,Va] = spec(A);          //determina os autovetores (Ve) e autovalores (Va) da matriz A
                    r(i,1) = sqrt(real(Va(1,1)))/wnp;
                    r(i,2) = sqrt(real(Va(2,2)))/wnp;
                    r(i,3) = 1
                    for x=1:3
                        r2 = r(i,x);
                        B = wa(i)/wnp;
                        A1 = (1-(r2^2/B^2))^2;
                        A2 = 4*(populacao(i,3)*(r2/B))^2;
                        A3 = ((r2^4/B^2)-(((4*populacao(i,3)*amort_estr)/B)+(1/B^2)+(populacao(i,1)+1))*(r2^2)+1)^2;
                        A4 = 4*(r2*(amort_estr+(populacao(i,3)/B)-(r2^3/B)*(populacao(i,3)+(amort_estr/B))-(r2^3/B)*populacao(i,3)*populacao(i,1)))^2;
                        Gr(i,x) = sqrt((A1+A2)/(A3+A4));
                    end
                    if Gr(i,1)>Gr(i,2) & Gr(i,1)>Gr(i,3) then
                        Gr(i,4) = Gr(i,1);
                    elseif Gr(i,2)>Gr(i,1) & Gr(i,2)>Gr(i,3) then
                        Gr(i,4) = Gr(i,2);
                    else                   //penalidade caso G3 maior que G2 e G1
                        Gr(i,3)=1000;
                    end
                end
            end
        end
    
    
        //----------------------------******-------------------------------------
        //Aplicando a FUNÇÃO OBJETIVO com as repectivas amplificações encontradas
        //----------------------------******-------------------------------------
    
    
        //A função objetivo será adicionada na quarta coluna da população.
    
    
        for i=1:n
            populacao(i,4)=(q1*abs(Gr(i,1)- Gr(i,2)) + q2*abs(Gr(i,4)-Glim) + q3*populacao(i,1) + q2*abs(Gr(i,3)-Glim/FS)); //Função objetivo
        end
    
    
        //----------------------------******----------------------------------
        //Ordenando população por melhores respostas perante a função objetivo
        //----------------------------******----------------------------------
    
    
        populacao = ordena_populacao(populacao, n); //utilizando a função ordena_populacao.sci
    
    
        //-----------------------------------******-------------------------------------
        //Definindo o elitismo da populacao para que não se perca os melhores resultados
        //-----------------------------------******-------------------------------------
    
    
        for i=1:ranking
            for j=1:4
                ranking1(i,j) = populacao(i,j);
            end
        end
    
    
    
        //Algumas impressões na tela.
        //===============================================================//
    
    
        if geracao==1 then
            mprintf("----------Populacao Inicial----------\n\n      mi        qopt     amort_AMS     fobj");
            disp(populacao());
        end
        if geracao==max_geracao then
            mprintf("\n\n----------Populacao Final----------\n\n      mi        qopt     amort_AMS     fobj");
            disp (populacao());
        end
        for j=1:4
            melhores_individuos(geracao,j)=populacao(1,j);
            melhores_individuos(geracao,5)=geracao;
        end
        if geracao==max_geracao   then
            mprintf("\n--------------------------------------------------------");
            mprintf("\n\nOs melhores indivíduos e suas respectivas gerações foram:\n\n     mi        qopt     amort_AMS     fobj       Geração");
            mprintf("\n  %f   %f   %f   %f   %d",melhores_individuos(1,1),melhores_individuos(1,2),melhores_individuos(1,3),melhores_individuos(1,4),melhores_individuos(1,5));
            for j=2:max_geracao-1
                if melhores_individuos(j+1,4)~=melhores_individuos(j,4);
                    mprintf("\n  %f   %f   %f   %f   %d",melhores_individuos(j+1,1),melhores_individuos(j+1,2),melhores_individuos(j+1,3),melhores_individuos(j+1,4),melhores_individuos(j+1,5));
                end
            end
            mprintf("\n------------------------------------------------------\n");
        end
        
    
    //===============================================================//
    
    
        //----------------------------******---------------------------------
        //If utilizado para que o ultimo laço não seja realizado erroneamente.
        //----------------------------******---------------------------------
    
    
        if geracao<=max_geracao then
        
        
            //----------------------------******---------------------------------
            //                Aplicando a seleção por torneio:
            //----------------------------******---------------------------------
        
            //Obs:
            //Seleção realizada da seguinte maneira: é selecionado dois indivíduos aleatoriamente
            //dentre a população atual e são comparadas seus valores encontrados para a FOBJ. O
            //selecionado sempre será o que estiver com o melhor resultado, que no caso é mais 
            //próximo de zero.
        
        
            [populacao,populacao_torneio] = selecao_torneio(populacao,n,ranking) //cria-se uma nova populacao, definida por populacao_torneio
        
        
            //----------------------------******---------------------------------
            //                  Aplicando os cruzamentos:
            //----------------------------******---------------------------------
        
        
            //A populacao cruzamento é criada da seguinte maneira: se o numero sorteado for abaixo da probabilidade
            //da taxa de cruzamento, então acontece o cruzamento entre 2 indivíduos da populacao_torneio 
            //diferentes sorteados ao acaso.
        
        
            populacao_cruzamento = cruzamento(n,populacao_torneio);


        //----------------------------------------******------------------------------------------------------
        //    Adicionando os indivíduos do elitismo e da populacao_cruzamento para a populacao principal.
        //----------------------------------------******------------------------------------------------------
        
        
            for i=1:n
                for j=1:3
                    if i<=ranking then
                        populacao(i,j) = ranking1(i,j);
                    else
                        populacao(i,j) = populacao_cruzamento(i-ranking,j);
                    end
                end
            end
        
        //----------------------******-----------------------------
        //                      Mutação
        //----------------------******-----------------------------
        
        
            populacao = mutacao(n,populacao,taxa_mutacao,ranking);
        
        end
        geracao=geracao+1; //Partindo para a próxima geração.
    else
        break
    end
end //fim do laço

//-------------------------****------------------------------
//Gerando os gráficos de G x r da solução encontrada pelo AG:
//-------------------------****------------------------------


//Adicionando o melhor individuo na variável ams_ag

for i=1:4
    ams_ag(1,i) = ranking1(1,i);
end

//Calculando as amplificações num dado intervalo.

r = linspace(0.8,1.2,1000);
for i=1:1000
    B = WA/wnp;   //beta da função de amplificação do artigo.
    A1 = (1-(r(i)^2/B^2))^2;
    A2 = 4*((ams_ag(1,3)*(r(i)/B))^2);
    A3 = ((r(i)^4/B^2)-(((4*ams_ag(1,3)*amort_estr)/B)+(1/B^2)+(ams_ag(1,1)+1))*(r(i)^2)+1)^2;
    A4 = 4*(r(i)*(amort_estr+(ams_ag(1,3)/B))-((r(i)^3/B)*(ams_ag(1,3)+(amort_estr/B)))-((r(i)^3/B)*ams_ag(1,3)*ams_ag(1,1)))^2;
    G(i) = sqrt((A1+A2)/(A3+A4));
end

//Calculando e plotando o gráfico com as respostas do ábaco

ams_abaco(1,1) = 0.017;
ams_abaco(1,2) = 0.978;
ams_abaco(1,3) = 0.058;
WA_abaco = wnp*ams_abaco(1,2);

for i=1:1000
    B = WA_abaco/wnp;   //beta da função de amplificação do artigo.
    A1 = (1-(r(i)^2/B^2))^2;
    A2 = 4*((ams_abaco(1,3)*(r(i)/B))^2);
    A3 = ((r(i)^4/B^2)-(((4*ams_abaco(1,3)*amort_estr)/B)+(1/B^2)+(ams_abaco(1,1)+1))*(r(i)^2)+1)^2;
    A4 = 4*(r(i)*(amort_estr+(ams_abaco(1,3)/B))-((r(i)^3/B)*(ams_abaco(1,3)+(amort_estr/B)))-((r(i)^3/B)*ams_abaco(1,3)*ams_abaco(1,1)))^2;

    Ga(i) = sqrt((A1+A2)/(A3+A4));
end

//Plotando os gráficos no mesmo esboço.


scf(0);
plot(r',G,"blue","thickness",1.5);
plot(r',Ga,"red","thickness",1.5);
xgrid(5);
xlabel("r","fontsize",3,"color","black");
ylabel("G","fontsize",3,"color","black");
title("Amplificação Dinâmica","fontsize",4);
legend(["Algoritmo Genético";"Ábacos"]);
    
//Gráfico para evolução das respostas/geração


evolucao(1) = melhores_individuos(1,4); //iniciando os valores da evolucao
Geracao(1) = 1; //numero da geracao dos melhores individuos diferentes
x(1) = 1;
cont=1;
for i=2:max_geracao
    if melhores_individuos(i,4)~=melhores_individuos(i-1,4) then //adicionando os melhores (se diferentes) para evolucao
        cont = cont+1;
        evolucao(cont) = melhores_individuos(i,4);
        Geracao(cont) = melhores_individuos(i,5);
    end
end

scf(1);
plot2d(Geracao,evolucao,-4);
xlabel("Número de Gerações otimizadas","fontsize",3,"color","black");
ylabel("Resultados da Função Objetivo","fontsize",3,"color","black");
title("Otimização por Geração","fontsize",4);


//--------------------------------Impressões Finais - Relatório------------------------------------------

mprintf("\n\n======================Relatório da Simulação===========================\n\n");
mprintf("************Parâmetros da Simulação************\n\n");
mprintf("Tamanho da população = %d \nNúmero de gerações = %d \nNúmeros de gerações completas = %d \nValor da função objetivo = %f \nTaxa de cruzamento = %d \nTaxa de mutação = %d \nNúmeros de indivíduos de elite = %d \nTolerância = %f \nPonderação q1 = %d \nPonderação q2 = %d \nPonderação q3 = %d \nPonderação q4 = %.3f \nFator de segurança (FS) = %.2f",n,max_geracao,geracao-1,ams_ag(1,4),taxa_cruzamento,taxa_mutacao,ranking,tol,q1,q2,q3,q4,FS);
mprintf("\n***********************************************");
mprintf("\n\n***************Dados da laje***************\n\n");
mprintf("Massa(mp) = %.2f kg \nRigidez equivalente (kp) = %.2f N/m \nAmortecimento adimensional = %.2f \nFrequência natural(wnp) = %.2f rad/s \nMáxima amplificação admissível = %.2f", mp,kp,amort_estr,wnp,Glim*q4);
mprintf("\n*********************************************");
mprintf("\n\n**********Dados obtidos para o AMS*********\n\n");
mprintf("Massa(ma) = %.2f kg \nRigidez equivalente (ka) = %.2f N/m \nAmortecimento adimensional = %f \nFrequência natural(wa) = %.2f rad/s",MA,KA,ams_ag(1,3),WA);
mprintf("\n*********************************************");
