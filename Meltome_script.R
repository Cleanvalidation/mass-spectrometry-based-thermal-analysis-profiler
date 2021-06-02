MEK1<-read_xlsx("Q02750.xlsx")
MEK2<-read_xlsx("P36507.xlsx")
MEK1<-MEK1 %>% dplyr::rename("temperature"="...1")
MEK2<-MEK2 %>% dplyr::rename("temperature"="...1")
MEK1<-MEK1 %>% pivot_longer(names(MEK1)[!names(MEK1)=="temperature"]) %>%
  dplyr::rename("Replicate"="name") %>%
  dplyr::mutate(Replicate=as.factor(Replicate))

MEK2<-MEK2 %>% pivot_longer(names(MEK2)[!names(MEK2)=="temperature"]) %>%
  dplyr::rename("Replicate"="name")%>%
  dplyr::mutate(Replicate=as.factor(Replicate))

MEK1<-MEK1 %>% 
  dplyr::mutate(fit=list(cetsa_fit(MEK1,norm=FALSE))) %>% 
  dplyr::mutate(fitted_values=predict(fit[[1]],se.fit=TRUE))
MEK2<-MEK2 %>% 
  dplyr::mutate(fit=list(cetsa_fit(MEK2,norm=FALSE))) %>%
  dplyr::mutate(fitted_values=predict(fit[[1]],se.fit=TRUE))
#generate data for confidence intervals
MEK1<-MEK1 %>% dplyr::ungroup(.) %>%
  tidyr::separate(Replicate,c("dataset","replicate")) %>% 
  dplyr::group_by(dataset,temperature) %>%  
  dplyr::mutate(mvalue=mean(value,na.rm=TRUE),
                sdvalue=sd(value,na.rm=TRUE),
                se=qnorm(.975)*(sdvalue/sqrt(as.numeric(max(replicate,na.rm=TRUE)))),
                ymin=mvalue-se,
                ymax=mvalue+se)
MEK2<-MEK2 %>% dplyr::ungroup(.) %>%
  tidyr::separate(Replicate,c("dataset","replicate")) %>% 
  dplyr::group_by(dataset,temperature) %>%  
  dplyr::mutate(mvalue=mean(value,na.rm=TRUE),
                sdvalue=sd(value,na.rm=TRUE),
                se=qnorm(.975)*(sdvalue/sqrt(as.numeric(max(replicate,na.rm=TRUE)))),
                ymin=mvalue-se,
                ymax=mvalue+se)
                                                        
                                                        

check<-ggplot2::ggplot(MEK1,mapping=aes(x=temperature,y=value)) +
  geom_line(mapping=aes(x=temperature,y=fitted_values),color="blue",size=2)+
  geom_errorbar(mapping=aes(ymin=ymin,ymax=ymax),size=2,color="black",width=1)+
  xlab('Temperature (\u00B0C)')+
  ylab("normalized intensity")+
  ggtitle("MEK1 (Jurkat Cell, Meltome)")+ 
  theme_classic()+
  theme(text = element_text(size = 30),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
        axis.text.x = element_text(face="bold", color="black", 
                                   size=30),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=30),
        axis.line = element_line(size=2),
        axis.ticks= element_line(size=2),
        axis.ticks.length=unit(.25, "cm"))+
  xlim(30,70)+ylim(0,1.5)
pdf("Meltome_MEK1.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
check
dev.off()
check<-ggplot2::ggplot(MEK2,mapping=aes(x=temperature,y=value)) +
  geom_line(mapping=aes(x=temperature,y=fitted_values),color="blue",size=2)+
  geom_errorbar(mapping=aes(ymin=ymin,ymax=ymax),size=2,color="black",width=1)+
  xlab('Temperature (\u00B0C)')+
  ylab("normalized intensity")+
  ggtitle("MEK2 (Jurkat Cell, Meltome)")+ 
  theme_classic()+
  theme(text = element_text(size = 30),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),
        axis.text.x = element_text(face="bold", color="black", 
                                   size=30),
        axis.text.y = element_text(face="bold", color="black", 
                                   size=30),
        axis.line = element_line(size=2),
        axis.ticks= element_line(size=2),
        axis.ticks.length=unit(.25, "cm"))+
  xlim(30,70)+ylim(0,1.5)
pdf("Meltome_MEK2.pdf",encoding="CP1253.enc",compress=FALSE,width=12.13,height=7.93)
check
dev.off()
