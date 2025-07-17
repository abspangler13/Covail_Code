#Note: this portion of the pipeline is an optional comparison so we can look at how the flow data
#stacks up to CITESeq

library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)
library(stringr)

#read in excel file of flow data
df <- read_xlsx(here::here("01_raw-data", "FlowData", "AllCOVAILMetadata_240314.xlsx"))

#filter down to only stage 1 samples present at 4 timepoints
infected <- unique(df$`Subject ID`[df$infect_flag == "Y"])
df$infect_flag <- ifelse(df$`Subject ID` %in% infected, "Y", 0)

df <- filter(df, df$Stage == 1, !is.na(df$infect_flag)) %>% 
      mutate(BoostInfect = paste0(Treatment, "_", infect_flag),
             Treatment = str_replace(Treatment, "  ", " ")) %>%
      filter(`Subject ID` %in% count(df, `Subject ID`)[count(df, `Subject ID`)$n == 4, 1][[1]])
  
summary <- df %>%
            group_by(df$BoostInfect) %>%
            summarize(TotalPropSum = mean(sum(`Live/IgG/Beta+/Beta | Freq. of IgG`,
                                              `Live/IgG/Beta+/Beta-BA1 | Freq. of IgG`,
                                              `Live/IgG/Beta+/Beta-BA1-XBB | Freq. of IgG`,
                                              `Live/IgG/Beta+/Beta-XBB | Freq. of IgG`,
                                              `Live/IgG/Proto neg Beta neg/BA1 | Freq. of IgG`,
                                              `Live/IgG/Proto neg Beta neg/BA1-XBB | Freq. of IgG`,
                                              `Live/IgG/Proto neg Beta neg/XBB | Freq. of IgG`,
                                              `Live/IgG/Proto+/Proto | Freq. of IgG`,
                                              `Live/IgG/Proto+/Proto-BA1 | Freq. of IgG`,
                                              `Live/IgG/Proto+/Proto-BA1-XBB | Freq. of IgG`,
                                              `Live/IgG/Proto+/Proto-XBB | Freq. of IgG`,
                                              `Live/IgG/Proto+Beta+/Proto-Beta | Freq. of IgG`,
                                              `Live/IgG/Proto+Beta+/Proto-Beta-BA1 | Freq. of IgG`,
                                              `Live/IgG/Proto+Beta+/Proto-Beta-BA1-XBB | Freq. of IgG`,
                                              `Live/IgG/Proto+Beta+/Proto-Beta-XBB | Freq. of IgG`)),
                      ProtoPOmiP = mean(sum(df$`Live/IgG/Proto+/Proto-BA1 | Freq. of IgG`,
                                       `Live/IgG/Proto+Beta+/Proto-Beta-XBB | Freq. of IgG`,
                                       `Live/IgG/Proto+Beta+/Proto-Beta-BA1-XBB | Freq. of IgG`,
                                       `Live/IgG/Proto+/Proto-BA1-XBB | Freq. of IgG`,
                                       `Live/IgG/Proto+/Proto-XBB | Freq. of IgG`,
                                       `Live/IgG/Proto+Beta+/Proto-Beta-BA1 | Freq. of IgG`
                                       )),
                      ProtoPOmiN = mean(sum(`Live/IgG/Proto+/Proto | Freq. of IgG`,
                                            `Live/IgG/Proto+Beta+/Proto-Beta | Freq. of IgG`
                                        )),
                      ProtoNOmiP = mean(sum(`Live/IgG/Beta+/Beta-BA1 | Freq. of IgG`,
                                            `Live/IgG/Beta+/Beta-BA1-XBB | Freq. of IgG`,
                                            `Live/IgG/Beta+/Beta-XBB | Freq. of IgG`,
                                            `Live/IgG/Proto neg Beta neg/BA1 | Freq. of IgG`,
                                            `Live/IgG/Proto neg Beta neg/BA1-XBB | Freq. of IgG`,
                                            `Live/IgG/Proto neg Beta neg/XBB | Freq. of IgG`
                                        )),
                      ProtoNOmiN = mean(sum(`Live/IgG/Beta+/Beta | Freq. of IgG`)),
                      )

summary$sumTotal <- summary$ProtoPOmiP + summary$ProtoPOmiN + summary$ProtoNOmiP

summary <- summary %>%
            mutate(ProtoPOmiP = ProtoPOmiP / sumTotal,
                   ProtoPOmiN = ProtoPOmiN / sumTotal,
                   ProtoNOmiP = ProtoNOmiP / sumTotal)

#plot
colnames(summary)[1] <- "Group"
summary <- pivot_longer(summary[,!colnames(summary) %in% c("TotalPropSum", "ProtoNOmiN", "sumTotal")], cols = starts_with("Proto"), names_to="Type", values_to="Percentage")

pdf(here("04_Analysis", "plots", "04_probe", "FlowDataReference_AveragePerGroup.pdf"))
ggplot(summary, aes(fill=Type, y=Percentage, x=Group))+
  geom_bar(position="stack",stat="identity")+
  ylab("Percentage of RBD+ Populations")+
  xlab("Panel")+
  ggtitle("Average Flow Data Populations - Stage 1")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12, angle = 90), axis.text.y=element_text(size=12))
dev.off() 
