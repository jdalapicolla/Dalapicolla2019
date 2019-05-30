contar.NA = function (x, var, group=0)
{
      
      if (group == 0)
      {
          tabela=apply(x[var], 2, function(x)sum(is.na(x)))
          
          if (dir.exists(file.path(getwd(), "NA")) == FALSE)
              {
              dir.create("NA", showWarnings = FALSE)
              }
          setwd(file.path(getwd(), "NA"))
          write.table(tabela, file=paste(paste0("NA_", "colunas"),"csv", sep="."), sep=",", col.names = TRUE)
          setwd('..')
  
  
          amostras=x[1]
          names(amostras)= "Observacao"
          linhas=apply(x[var], 1, function(x)sum(is.na(x)))
          tabela.linhas=data.frame(amostras, linhas)
          setwd(file.path(getwd(), "NA"))
          write.table(tabela.linhas, file=paste(paste0("NA_", "observacoes"),"csv", sep="."), sep=",", col.names = TRUE)
          setwd('..')
      }
          
      if (group != 0)
        {
          fatores=sort(as.factor(unique(x[,group])))
          nrow=length(fatores)
          ncol=length(var)
          tabela= as.data.frame(matrix(NA, nrow, ncol))
          names(tabela)= colnames(x[var])
          loop=0
                
          for (i in var)
              {
              loop=loop+1
              linha= aggregate.data.frame(x[i], x[group], function(x)sum(is.na(x)))
              linha=linha[,-1]
              tabela[,loop]= linha
              }
                
          newrow= c(1:ncol)
          tabela=rbind(tabela,newrow)
          nomes= c(as.character(fatores), "TOTAL")
          row.names(tabela)= nomes
          tabela[nrow(tabela),]=apply(x[var], 2, function(x)sum(is.na(x)))
                
          if (dir.exists(file.path(getwd(), "NA")) == FALSE) #se a pasta com o nome NA nao existir
              {
              dir.create("NA", showWarnings = FALSE)# criar a pasta com o nome NA
              }
          
          setwd(file.path(getwd(), "NA")) 
          write.table(tabela, file=paste(paste0("NA_", colnames(x[group])),"csv", sep="."), sep=",", col.names = TRUE)
          setwd('..')
          
          amostras=x[1]
          names(amostras)= "Observacao"
          linhas=apply(x[var], 1, function(x)sum(is.na(x)))
          tabela.linhas=data.frame(amostras, linhas)
          setwd(file.path(getwd(), "NA")) 
          write.table(tabela.linhas, file=paste(paste0("NA_", "observacoes"),"csv", sep="."), sep=",", col.names = TRUE)
          setwd('..')
      }
      
      resultados=list(tabela, tabela.linhas)
      return (resultados)
   }