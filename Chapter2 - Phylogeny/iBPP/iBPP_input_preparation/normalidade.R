normalidade = function (x, quant, group=0){

  if (dir.exists(file.path(getwd(), "Normalidade")) == FALSE) 
      {
      dir.create("Normalidade", showWarnings = FALSE)
      }
      setwd(file.path(getwd(), "Normalidade"))
  
  
  if (group==0)
      {
      loop=0
      w= 1:length(quant)
      pvalue= 1:length(quant)
  
      for (i in quant)
          {
          loop=loop+1
          teste= shapiro.test(x[,i])
          w[loop]=as.numeric(teste$statistic)
          pvalue[loop]=as.numeric(teste$p.value)
     
          jpeg(filename = paste("Normalidade_", colnames(x[i]), ".jpg", sep=""), width = 2000, height = 2000, units = "px", pointsize = 12, quality = 100, bg = "white", res = 300, family = "") #padronizar o gráfico
          par(mar= c(8, 4, 4, 2))
          qqnorm(x[,i])
          qqline(x[,i])
          mtext(paste("Shapiro-Wilk:", as.numeric(teste$statistic)), 1, line=5)
          mtext(paste("p-valor:", as.numeric(teste$p.value)), 1, line=6)
          dev.off()
          }
      
      tabela=data.frame(w, pvalue, row.names = colnames(x[quant]))
      write.table(tabela, file=paste("Normalidade_Tabela", "csv", sep="."), sep=",")
      setwd('..')
      resultados=tabela
      }
  
  if (group!=0){
      fatores=as.factor(unique(x[,group]))
      resultados=list()
      fac=0
    
      for (fator in fatores)
          {
          df= x[x[,group]==fator,]
          fac=fac+1
          loop=0
          w= 1:length(quant)
          pvalue= 1:length(quant)
     
          if (dir.exists(file.path(getwd(), fator)) == FALSE) 
              {
              dir.create(fator, showWarnings = FALSE)
              }
      
          
          setwd(file.path(getwd(), fator))
    
          for (i in quant)
              {
              loop=loop+1
              teste= shapiro.test(df[,i])
              w[loop]=as.numeric(teste$statistic)
              pvalue[loop]=as.numeric(teste$p.value)
      
              jpeg(filename = paste("Normalidade_", colnames(df[i]), ".jpg", sep=""), width = 2000, height = 2000, units = "px", pointsize = 12, quality = 100, bg = "white", res = 300, family = "") #padronizar o gráfico
              par(mar= c(8, 4, 4, 2))
              qqnorm(df[,i])
              qqline(df[,i])
              mtext(paste("Shapiro-Wilk:", as.numeric(teste$statistic)), 1, line=5)
              mtext(paste("p-valor:", as.numeric(teste$p.value)), 1, line=6)
              dev.off()
              }
    
          tabela = data.frame(w, pvalue, row.names = colnames(df[quant]))
          resultados[[fac]]= tabela
          names(resultados)[fac]= fator
          write.table (tabela, file=paste(paste0("Normalidade_", fator),"csv", sep="."), sep=",")
          setwd('..')
      }
      setwd('..')
      }
  
  return (resultados)
  
}