
id.outliers = function (x, quant, group=0, id="box", NUMBER=10, visual="box", res="LOW", csv=FALSE)
{
  if (group!=0) #Para quando houver subgrupos para analisar
  {
    fatores=as.factor(unique(x[,group])) #extrair os valores únicos da coluna que sera usada como subgrupo e transforma-los em fatores 
    for (fator in fatores)#fazer as analises para cada fator, em ciclo 
    {
      df= x[x[,group]==fator,] #extrair todas as linhas do data frame original que contem o fator de interesse e salvar cada uma delas como 'df'
      resultados = etapa2 (df, quant, fator, id, NUMBER, csv) #usar a funcao auxiliar 'etapa2' para a identificacao dos outliers e salvar a lista resultante em 'resultados'
      etapa3 (df, quant, fator, visual, res) # usar a funcao auxiliar 'etapa3' para a elaboracao dos graficos
    }
  }
  else #Quando nao houver subgrupos para analisar
  {
    df=x #considerar a tabela toda como df
    fator="tabela_completa" #considerar o fator com esse nome
    resultados = etapa2 (df, quant, fator, id, NUMBER, csv) #funcao auxiliar para identificacao dos outliers
    etapa3 (df, quant, fator, visual, res)#funcao auxiliar para a identificacao dos outliers
  }
  return(resultados) #retorna uma lista dos os possiveis outliers identificados
}

################################################################
######### FUNCOES AUXILIARES ###################################

mudarSubDir = function (fator, subfolder)
  {#funcao auxiliar para mudar o diretorio antes de salvar graficos ou csv
      
      if (dir.exists(file.path(getwd(), fator)) == FALSE) #se a pasta com o nome do fator nao existir
        {
        dir.create(fator, showWarnings = FALSE)# crar a pasta com o nome do fator
        }
      setwd(file.path(getwd(), fator)) #mudar para a pasta criada
  
      if (dir.exists(file.path(getwd(), subfolder)) == FALSE) #se o subfolder nao existir
        {
        dir.create(subfolder, showWarnings = FALSE)# criar subfolder
        }
      setwd(file.path(getwd(), subfolder)) #mudar o diretorio para o subfolder
  }

#######################################################################################

rval = function(y) #funcao auxiliar para calcular o ESD
  {
  ares = abs(y - mean(y, na.rm = TRUE))/sd(y, na.rm = TRUE)
  tab = data.frame(y, ares)
  r = max(tab$ares)
  list(r, tab)
  }

#######################################################################################

etapa2= function (df, quant, fator, id="box", NUMBER=10, csv=FALSE)
{
  resultado.box=list()#criar listas para salvar os resultados dos tres metodos
  resultado.z=list()
  resultado.ESD=list()
 
  if (id=="box" | id=="ALL")#se o metodo for 'box'
    {
      fatorSemObs=TRUE #variavel que vai indicar a existencia ou nao de outliers 
      loop=0 #variavel que vai indicar a qual e o ciclo do loop 
      
      for (i in quant)#ciclo para cada coluna quantitativa indicada no inicio da funcao
        {
          loop=loop+1#incremento para o numero de loop
          
          if (length(df[,i])>=5) #se o numero de amostras para tal coluna for menor do que 5 nao e possivel indicar os outliers, vai para o proximo 'i'
              {
                graph=boxplot(df[i], data=df, main= colnames(df[i]), na.omit(df[i])) #cria o boxplot em um objeto, uma lista dentro dele chamado 'out' tem os valores dos possiveis outliers
                posi.out= which (df[,i] %in% graph$out) #descobrir em qual posicao da tabela esta esses valores de outliers
                nobs=length(graph$out) #numero de outliers encontrados
              
                if (nobs != 0) #se o numero de outliers for diferente de 0
                  {
                  fatorSemObs=FALSE#muda o valor da variavel indicativa
                  out.lines= as.data.frame(matrix(NA, nrow= nobs, ncol=ncol(df))) #cria uma tabela para salvar os resultados
                  colnames(out.lines)= colnames(df) #dar os nomes para as colunas da tabela
                  
                      for (j in 1:nobs) #para cada observacao de outliers
                        {
                        out.lines[j,]= df[posi.out[j], ] #salvar a linha inteira da tabela original "df" que contenha o valor de outlier na tabela de resultados
                        }
                  resultado.box[[loop]]= out.lines#salvar a tabela de resultados no objeto lista final na posicao indicada pelo loop. 1ºciclo posicao 1, 2ºciclo posicao 2, etc.
                  names(resultado.box)[loop]= colnames(df[i]) #dar o nome para cada tabela salva na lista
                  
                  if (csv==TRUE) #se foi pedido para salvar csv
                        {
                        mudarSubDir(fator, "box") #usar funcao auxiliar para mudar o diretorio
                        write.table(out.lines, file=paste(colnames(df[i]),"csv", sep="."), sep=",") #salvar resultado como csv 
                        setwd('../..')#voltar ao diretorio original
                        }
                  }
              }
            
          else #se nao tiver mais 4 amostras
              {
                print(paste(colnames(df[i])," do agrupamento ",fator," nao possui numero de amostras suficientes para realizar as analises (n>=5)", sep="")) #print a mensagem do erro
                next #ir para o proximo ciclo
              }
        }
    
      if (fatorSemObs==TRUE) print (paste(fator, "nao tem outliers para metodo 'box'", sep=" "))#se nao teve outliers identificados imprimir essa mensagem
    }
  
  if (id=="z" | id=="ALL")#se for o metodo 'z'
    {
      fatorSemObs=TRUE#variavel que vai indicar a existencia ou nao de outliers
      loop=0#variavel que vai indicar a qual e o ciclo do loop
      for (i in quant)#ciclo para cada coluna quantitativa indicada no inicio da funcao
          { 
            loop=loop+1#incremento para o ciclo do loop
            obser=sapply(df[,i], as.numeric)#salvar valores da coluna do data frame original como 'numeric'
            mod.z= abs((0.6745*(obser-median(obser, na.rm=TRUE)))/ mad(df[i], na.rm=TRUE))#calculo da estatistica
            
            if (any(mod.z > 3.5, na.rm = TRUE)==TRUE)#se o valor da estatistica for maior que 3.5, ha outliers 
            {
              fatorSemObs=FALSE#mudar variavel indicativa de outliers
              positions=which(mod.z> 3.5)#descobrir as posicoes dos outliers no vetor numerico. e a mesma na tabela original
              out.lines=df[positions,]#salvar as linhas correspondentes as posicoes em um novo arquivo
              resultado.z[[loop]]= out.lines#salvar arquivo em uma lista
              names(resultado.z)[loop]= colnames(df[i])#nomear elemento da lista com o nome da coluna do arquivo original
              
              if (csv==TRUE) #salvar em csv se isso foi pedido
                {
                  mudarSubDir(fator, "z")
                  write.table(out.lines, file=paste(colnames(df[i]),"csv", sep="."), sep=",") 
                  setwd('../..')
                }
            }
          }
    if (fatorSemObs==TRUE) print (paste(fator, "nao tem outliers para metodo 'z'", sep=" "))#se nao teve outliers identificados imprimir essa mensagem
    }
  
  if (id=="ESD"| id=="ALL") #se for o metodo 'ESD'
    {    
      fatorSemObs=TRUE
      loop=0
      for (j in quant)
        {
          loop=loop+1
          #y=as.vector(df[j][-1,])
          y=sapply(df[,j], as.numeric)
          y=y[!is.na(y)] #retirar os NA da tabela
          n = length(y) #comprimento do vetor
          alpha = 0.05
          NUMBER2= min(NUMBER,ceiling(n/2)) #numero maximo de outliers testado nunca podera ser maior do que 50% dos dados. Por exemplo é pedido que se teste 20 outleirs, para um subgrupo com n=50 não tem problema na fórmula do teste. Mas se um dos vários subgrupos tiver só 6 amostras (acontece com os meus dados), e pede para retirar 20 outleirs, a fórmula dará erro (NaN's produced). O teste para ser eficaz precisa ter n maior que 15 segundo os autores.
          lam = 1:NUMBER2 #cria objeto para salvar o valor de lam
          R = 1:NUMBER2 #cria objeto para salvaro valor de R
          
              if(n>=5) #se tiver numero de amostra sufuciente
                {
                    for (i in 1:NUMBER2)#calcular o teste estatistico, retirei da referência que pus na proposta, só fiz algumas mudancas. Essa parte do codigo nao e de minha autoria  
                        {
                            if(i==1)
                                {
                                  rt = rval(y) #funcao auxiliar de autoria dos autores do algoritmo
                                  R[i] = unlist(rt[1])
                                  tab = data.frame(rt[2])
                                  newtab = tab[tab$ares!=max(tab$ares),]
                                }
                            else if(i!=1)
                                {
                                  rt = rval(newtab$y)
                                  R[i] = unlist(rt[1])
                                  tab = data.frame(rt[2])
                                  newtab = tab[tab$ares!=max(tab$ares),]
                                }
                            #Computar o valor critico.
                            p = 1 - alpha/(2*(n-i+1))
                            t = qt(p,(n-i-1))
                            lam[i] = t*(n-i) / sqrt((n-i-1+t**2)*(n-i+1))
                        }
                }
              
              else #se o numero de amostra for pequeno 
                {
                  print(paste(colnames(df[j])," do agrupamento ",fator," nao possui numero de amostras suficientes para realizar as analises (n>4)", sep=""))#imprimir o aviso
                  next#ir para o proximo 'i'
                }
              
          newtab = data.frame(c(1:NUMBER2),R,lam) #criar um tabela com o resultado
          names(newtab)=c("No. Outliers","Test Stat.", "Critical Val.")#acrescentar os nomes das colunas
          resultado.ESD[[loop]]= newtab#salvar tabela na lista
          names(resultado.ESD)[loop]= colnames(df[j]) #salvar o nome da coluna da tabela original usada para o calculo
      
          if (any(newtab[2] > newtab[3], na.rm = TRUE)==TRUE)#se houve algum valor que ultrapassou o ponto critico
              {
                fatorSemObs=FALSE #mudar a variavel de indicacao de outliers
                numb.out=which(newtab[2] > newtab[3]) #descobrir quais numeros ultrapassaram
                print(paste(colnames(df[j])," do agrupamento ",fator," possui ", max(numb.out), " possivel(is) outlier(s)", sep="")) #imprimir a mensagem informando o numero de possiveis ouliers
                
                if (csv==TRUE) #salvar em csv se foi requisitado
                    {
                      mudarSubDir(fator, "ESD")
                      write.table(newtab, file=paste(colnames(df[j]),"csv", sep="."), sep=",") 
                      setwd('../..')
                    }
              }
    }
    
    if (fatorSemObs==TRUE) print (paste(fator, "nao tem outliers para metodo 'ESD'", sep=" ")) #imprimir essa mensagem se nao ha outliers  
    }
  
  return (list(resultado.box, resultado.z, resultado.ESD)) #retornar os resultados dos tres metodos, se um deles nao for escolhido a lista retornara como NULL
}

#######################################################################################

etapa3 = function (df, quant, fator, visual= "boxplot", res="LOW")
  {
      #funcao para elaboracao dos graficos 
      #parametros definidos como padrao, para o 'LOW'
      width = 1000
      height = 480*ceiling(length(quant)/2)
      resol = 150
      
      if (res=="MED")#parametros para 'MED'
          {
          width = 2000
          height = 960*ceiling(length(quant)/2)
          resol = 300
          }
  
      if (res=="HIGH")#parametros para 'HIGH'
          {
            width = 4000
            height = 1920*ceiling(length(quant)/2)
            resol = 300
          }
          
      if (visual=="boxplot" | visual=="ALL") #para graficos boxplot
          {
            mudarSubDir(fator, "boxplot")#funcao auxiliar para mudar o diretorio antes de salvar o grafico
            jpeg(filename = paste("graph_boxplot_", res, ".jpg", sep=""), width = width, height = height, units = "px", pointsize = 12, quality = 100, bg = "white", res = resol, family = "") #abrir o dispositivo e padronizar o formato jpeg para salvar o grafico
            par(mfrow= c(ceiling(length(quant)/2), 2)) #dividir a janela grafica   
            
            for (k in quant) #Ciclo para todas as variaveis quantitativas, identificadas pela posicao, um grafico por variavel 
              {
                  if (length(df[,k])>=5)#se houver mais de cinco amostras, contruir o grafico
                  {
                    graph=boxplot(df[k], data=df, main= colnames(df[k]), na.omit(df[k]))#desenhar o grafico
                    posi.out= which (df[,k] %in% graph$out)#saber a posicao dos outliers
                        
                      if(length(posi.out)>0) #se tiver outliers
                          {text(graph$group, graph$out,posi.out, pos = 4, col= "red")} #identifica-los no grafico com o numero da posicao da linha na tabela original 'x'
                          
                  }
                  
                  else #se o numero de amostra for pequeno, nao realizar a analise 
                  {print(paste(colnames(df[k])," do agrupamento ",fator," nao possui numero de amostras suficientes para realizar as analises (n>4)", sep=""))#imprimir essa mensagem
                    next#ir para o proximo k
                  }
              }
            
            par(mfrow= c(1,1))#voltar ao padrão original da janela grafica
            dev.off() #desligar o dispositivo
            setwd('../..')#voltar para o diretorio original
          }
          
  
      if (visual=="pontos" | visual=="ALL") #para graficos dotchart
          {
            mudarSubDir(fator, "pontos")
            jpeg(filename = paste("graph_pontos_", res, ".jpg", sep=""), width = width, height = height, units = "px", pointsize = 12, quality = 100, bg = "white", res = resol, family = "") #padronizar o formato jpeg para salver o gráfico
            par(mfrow= c(ceiling(length(quant)/2), 2)) #dividir a janela grafica   
    
                for (k in quant) #Ciclo para todas as variaveis quantitativas, identificadas pela posição 
                  {
                      if (length(df[,k])>=5)
                      {
                      dotchart(as.numeric(df[,k]), main= colnames(df[k]), pch=16)
                      text(df[,k], 1:length(as.vector(df[,k])),rownames(df[k]), pos=4, cex=1, col="red") #nome das linhas nos graficos
                      }
                      
                      else 
                      {print(paste(colnames(df[k])," do agrupamento ",fator," nao possui numero de amostras suficientes para realizar as analises (n>4)", sep=""))
                      next
                      }
                  }
              
            par(mfrow= c(1,1))#voltar ao padrão original da janela grafica
            dev.off() #desligar o dispositivo
            setwd('../..') #voltar para diretorio original
          }

  
      if (visual=="biplot" | visual=="ALL") #para graficos com duas variaveis
          {
              mudarSubDir(fator, "biplot") #usar funcao auxiliar para mudar o diretorio antes de salvar os graficos
              jpeg(filename = paste("graph_biplot_", res, ".jpg", sep=""), width = width, height = height, units = "px", pointsize = 12, quality = 100, bg = "white", res = resol, family = "") #padronizar o gráfico
              par(mfrow= c(ceiling(length(quant)/2), 2)) #dividir a janela grafica   
              biplot=quant[-length(quant)] #retirar a posicao do ultimo 'quant'
            
                  for(k in biplot)
                      {
                         if (length(df[,k])>=5) #se tiver amostras suficientes
                            {
                              plot(df[,k], df[,k+1], xlab= colnames(df[k]), ylab= colnames(df[k+1])) #desenhar os graficos
                              text(df[,k], df[,k+1],rownames(df[k]), pos=4, cex=1, col="red")#colocar a posicao das linhas nos pontos
                            }
                         else #se não tiver amostras suficientes
                            {print(paste(colnames(df[k])," do agrupamento ",fator," nao possui numero de amostras suficientes para realizar as analises (n>=5)", sep=""))
                            next
                            }
                      }
            
            par(mfrow= c(1,1))#voltar ao padrão original da janela grafica
            dev.off() #desligar o dispositivo
            setwd('../..') #voltar ao diretorio original
          }
        }

########################################################################
