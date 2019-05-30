repetibilidade = function (x, data)
{
  
  
  if("ICC" %in% rownames(installed.packages()) == FALSE){install.packages("ICC")
  } else {print (paste0("'ICC' ja instalado na biblioteca"))}
  
  library(ICC)
  
  if (dir.exists(file.path(getwd(), "Repetibilidade")) == FALSE) 
  {
    dir.create("Repetibilidade", showWarnings = FALSE)
  }
  
  setwd(file.path(getwd(), "Repetibilidade"))
  
  
  variavel=colnames(data)
  variavel=variavel[-1]
  n=length(variavel)
  
  resultados=matrix(data = NA, nrow = 7, ncol= n, byrow = FALSE, dimnames = NULL)
  row.names(resultados) = c("ICC", "LowerCI", "UpperCI", "N", "k", "varw", "vara")
  colnames(resultados) = c(variavel)
  
  for (i in variavel){
    resul=ICCest(x, i, data, alpha = 0.05, CI.type = "Smith")
    resultados[1,i] = resul$ICC
    resultados[2,i] = resul$LowerCI
    resultados[3,i] = resul$UpperCI
    resultados[4,i] = resul$N
    resultados[5,i] = resul$k
    resultados[6,i] = resul$varw
    resultados[7,i] = resul$vara
  }
  write.table(resultados, file=paste("Repetibilidade_Tabela", "csv", sep="."), sep=",", col.names = NA, row.names = T)
  setwd('..')
  
  return(resultados)
}
