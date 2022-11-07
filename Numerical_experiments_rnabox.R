
# R code with numerical experiments from paper (in preparation):
# Wesolowski J., Wieczorkowski R., Wojciak W. (2023?), 
# Adjusting the recursive Neyman algorithm for two-sided bounds on sample strata sizes.


library(stratallo) # implementation of 'rnabox' algorithm
library(dplyr)
library(microbenchmark)
library(bench)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(patchwork)

# R codes for used algorithms

# integer allocation algorithms
source("CapacityScaling.R")
source("SimpleGreedy.R")

# fixed point iteration algorithm
source("fpia.R")


# random rounding
ran_round<-function(x)
{
  return( floor(x)+(runif(length(x))<(x-floor(x))) )
}


# rounding based on article:
# Rama Count, Massoud Heidari, Optimal rounding under integer constraints
# December 2014, arxiv
round_oric <- function(x)
{
  n <- round(sum(x))
  m <- floor(x)
  y <- x - m
  Ix <- sum(y)
  
  if (Ix==0) return(x)
  else {
    iy <- order(-y)
    u <- unique(y[iy])
    z <- integer(length(x))
    for (i in 1:length(u)) z[iy] <- z[iy] + (y[iy]==u[i])*i
    iy2 <- order(-y,z,-m)
    #m[iy][iy2][1:Ix] <- ceiling(x[iy][iy2][1:Ix])
    m[iy2][1:Ix] <- (m[iy2][1:Ix])+1
    return(m)
  }
  }




# generation of artificial populations and box constraints

gen_population_parameters <- function(pop_n=1) {
# pop_n - numbering of populations, used values 1, 2 or 3 

if (pop_n==1) { Nrep <- 100; set.seed(2234)}  # parameter used in generation
if (pop_n==2) { Nrep <- 200; set.seed(876)}

source("gen_population.R")

if (pop_n %in% c(1,2)) {
  
  pop <- gen_population(Nrep=Nrep)
  Nh<-pop$Nh
  Sh<-pop$Sh
  NROW(Nh)
  # plot(Nh*Sh)
  # hist(Sh*Nh)
  
# Generation of lower and upper bounds - method used in JSSaM article
# improving population to have allocation with values greater than 0
# using integer allocation algorithm
mh<-rep(0,length(Nh)) # lower bounds
Mh <- Nh # upper bounds
(n <- round(0.1*sum(Nh)))
(alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh))
sum(alc)

length(alc[alc>2])/length(alc)
ix <- which(alc>2)
Nh <- Nh[ix]
Sh <- Sh[ix]
dh <- Nh*Sh
mh <- mh[ix]
Mh <- Mh[ix]
mh <- pmin(100,Mh) # additional lower constraints
(s1<-sum(mh))
(s2<-sum(Mh))
}

if (pop_n==3) {

# Generation simple population used in article published in JSSaM (2021, fig.2)
# added lower constraints
set.seed(4321)
Nh <- rep(1000,20)
Sh <- 10^(1:20)
Sh <- sample(Sh,20)
mh <- rep(100,20)
Mh <- Nh
dh <- Sh*Nh
}

return(list(Nh=Nh,Sh=Sh,mh=mh,Mh=Mh,dh=dh))
}



pop1 <- gen_population_parameters(pop_n = 1)
pop2 <- gen_population_parameters(pop_n = 2)
pop3 <- gen_population_parameters(pop_n = 3)




#N <- sum(Nh) # population size
#NROW(Nh)




# Time comparison
# creating data with times for selected algorithms and different fractions

tab_with_times <- function(pop_n=1) {
# pop_n - number of population (1,2 or 3)
pop <- gen_population_parameters(pop_n)
Nh <- pop$Nh
Sh <- pop$Sh
mh <- pop$mh
Mh <- pop$Mh
dh <- pop$dh


tab<-NULL

(s1<-sum(mh))
(s2<-sum(Mh))
(N <- sum(Nh))
dh <- Nh*Sh

#for (f in seq(0.01,0.8,0.05)) {
for (f in seq(s1/N,s2/N,0.1)) {
print(f)
N<-sum(Nh)
n<-round(f*N)

if (s1<n && n<s2) {

alc<-CapacityScaling(n, Nh, Sh, mh = mh, Mh = Mh)
h_under <- sum(alc<=mh)  # number of take-all strata  
h_over <- sum(alc>=Mh)  # number of take-all strata  


al_rnabox <- round(stratallo::dopt(n, dh, mh, Mh))
if (max(abs((al_rnabox-alc)))>1) { stop("Bad allocation rnabox!") }


al_fpi <- round(fpia(n,Nh,Sh, mh, Mh)$nh)
if (max(abs((al_fpi-alc)))>1) { stop("Bad allocation fpia!") }

options(digits=10)

ex<-microbenchmark(times=10,unit="ms",
                   fpia=fpia(n,Nh,Sh, mh, Mh)$nh,
                   rnabox=dopt(n, dh, mh, Mh)
)


summary(ex)
autoplot(ex)

                                            

exi<-group_by(ex,expr) %>% 
  summarise(Median_time=median(time)/1e6,
            Mean_time=mean(time)/1e6) # from nanoseconds to miliseconds
exi<-mutate(exi,s1=s1,s2=s2,N=N,f=f,H=length(Nh), h_under=h_under,h_over=h_over)

tab<-bind_rows(tab,exi)

}
}
 return(tab)
}


tab1 <- tab_with_times(pop_n = 1)
saveRDS(tab1,"tab1.rds") # results for population generated with Nrep=100

tab2 <- tab_with_times(pop_n = 2)
saveRDS(tab2,"tab2.rds") # results for population generated with Nrep=200

tab3 <- tab_with_times(pop_n = 3)
saveRDS(tab3,"tab3.rds") # results for small population



### Creation of plots 
options(digits=6)

# fig_n=1 uses populations 1 and 2
# fig_n=2 uses small population 3

fig_n <- 2

if (fig_n==1) {
  tab1 <- readRDS("tab1.rds")
  tab2 <- readRDS("tab2.rds")
  tab <- bind_rows(tab1,tab2)
} else tab <- readRDS("tab3.rds")


tab<-mutate(tab,H=as.factor(H),algorithm=expr , 
            flab=paste0("[",as.character(h_under),",",as.character(h_over),"]")
            )

as.vector(table(tab$N))

xN<-count(tab,N)$N

levels(tab$H)<-paste(levels(tab$H),"strata, N =", 
                     #sprintf("%f6.1",xN))
                     xN)
round(table(tab$N))
count(tab,H)

levels(tab$algorithm)

tab <- mutate(tab,f=round(f,3),algorithm=relevel(algorithm,c("fpia")))

group_by(tab,algorithm) %>% summarise(mean(Median_time),mean(Mean_time))

s1 <- tab$s1[1]
s2 <- tab$s2[1]
N <- tab$N[1]

p1<-
  ggplot(data=tab,aes(x=f,y=Median_time, shape=algorithm)) +
  #scale_y_log10() + 
  geom_point(size=2) +
  geom_line(data=tab,aes(x=f,y=Median_time,linetype=algorithm)) +
  facet_wrap(~H,scale="free") +
  #coord_cartesian(xlim=c(0.15,0.7)) +
  coord_cartesian(xlim=c(s1/N,s2/N)) +
  #scale_x_continuous(breaks = seq(0.0,0.6,0.1)) +
  #scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  labs(
       #x=" ",
       x="Sample fraction",
       #y = expression(paste("Time [miliseconds] (", log[10],  " scale)")),
       y="Time [miliseconds]" ,
       color="Algorithms: ", 
       title="Time comparison of selected algorithms"
      #, subtitle = "using microbenchmark package from R"
      )  +
  theme_bw(base_size=12) + theme(axis.title.x=element_blank(),
                                 axis.text.x=element_blank(),
                                 axis.ticks.x=element_blank(),
                                 legend.position = "right")

  #theme_bw(base_size=12) + theme(legend.position = "right") 
  #theme(#panel.background = element_rect(fill = NA, colour = "black"),
        #panel.grid = element_line(colour = "grey")
        #,strip.background = element_blank(), strip.text.x = element_blank())
  #theme(legend.position = "right",legend.text=element_text(size=rel(1.2)))
  #coord_flip()

##ggsave("fig_times.png",p,device="png", dpi=600, width = 8, height = 8/1.618)


### for additional parts of graphs

if (fig_n==1) {
  tab1 <- readRDS("tab1.rds")
  tab2 <- readRDS("tab2.rds")
  tab <- bind_rows(tab1,tab2)
} else tab <- readRDS("tab3.rds")


tab<-mutate(tab,algorithm=expr)
count(tab,algorithm)

df <- as.data.frame(tab)
df <- tab

# preapre population column, and algorithms facet
population <- paste0(df$H, " strata, N = ", df$N) # population info (no of strata, pop size)
population <- factor(population, levels = unique(population)[order(df$H)]) # levels order influences plotting order

df <- data.frame(facet = "algorithms",
                 plyr::rename(df[, c("expr", "Median_time", "f", "h_under","h_over")], c('expr' = 'series', 'Median_time' = 'value')),
                 population = population,
                 takeall_under_pct = round((df$h_under/df$H )*100,1),
                 takeall_over_pct = round((df$h_over/df$H )*100,1)) # Take-strata in %

# prepare takeall facet
df_takeall_under <- df[, c("f", "population", "takeall_under_pct", "h_under")]
df_takeall_under <- df_takeall_under[!duplicated(df_takeall_under), ]
df_takeall_under <- data.frame(facet = "takeall_under",
                         series = "takeall_under_pct",
                         plyr::rename(df_takeall_under, c('takeall_under_pct' = 'value')))

df_takeall_over <- df[, c("f", "population", "takeall_over_pct", "h_over")]
df_takeall_over <- df_takeall_over[!duplicated(df_takeall_over), ]
df_takeall_over <- data.frame(facet = "takeall_over",
                               series = "takeall_over_pct",
                               plyr::rename(df_takeall_over, c('takeall_over_pct' = 'value')))

df_takeall <- bind_rows(df_takeall_under,df_takeall_over)

df_takeall <- mutate(df_takeall,nover=ifelse(!is.na(h_under),h_under,h_over))

# rbind algorithms and takeall facets
#df_plot <- plyr::rbind.fill(df[, colnames(df_takeall)], df_takeall)

df_takeallG <- filter(df_takeall,facet=="takeall_over")
df_takeallD <- filter(df_takeall,facet=="takeall_under")


df_takeallD <- mutate(df_takeallD, f=round(f,1))
df_takeallG <- mutate(df_takeallG, f=round(f,1))


p2 <- 
  ggplot(data = df_takeallG, mapping = aes(x = f, y = value)) +
  #geom_text(data = subset(df_plot, facet == "algorithms" & series == "rNa"), aes(label = niter), vjust = -.6) +
  geom_bar(data = df_takeallG, mapping = aes(y = value), stat = "identity") +
  geom_text(data = df_takeallG, mapping = aes(label = nover), size = 3.0, vjust = 1.1, color="white") +
  #geom_blank(data = subset(df_plot, facet == "takeall"), mapping = aes(y = 100)) +
  facet_wrap(~population) + 
  #theme_bw(base_size=12) + 
  theme(legend.position = "right",
        panel.spacing.x=unit(2.5, "lines"),
        axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    #panel.background = element_rect(fill = NA, colour = "black"),
    #panel.grid = element_line(colour = "grey"),
    strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x=" ", y= 'Take-max \nstrata [%]') +
  #coord_cartesian(xlim=c(0.1,0.7))
  coord_cartesian(xlim=c(s1/N,s2/N)) 
  #scale_x_continuous(breaks = seq(0.0,0.6,0.1)) +
  #scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  #scale_y_continuous(breaks = scales::pretty_breaks(5))


p3 <- 
  ggplot(data = df_takeallD, mapping = aes(x = f, y = value)) +
  #geom_text(data = subset(df_plot, facet == "algorithms" & series == "rNa"), aes(label = niter), vjust = -.6) +
  geom_bar(data = df_takeallD, mapping = aes(y = value), stat = "identity") +
  geom_text(data = df_takeallD, mapping = aes(label = nover), size = 3.0, vjust = 1.1, color="white") +
  #geom_blank(data = subset(df_plot, facet == "takeall"), mapping = aes(y = 100)) +
  facet_wrap(~population) + 
  #theme_bw(base_size=12) + 
  theme(legend.position = "right",
        panel.spacing.x=unit(2.5, "lines"),
        #panel.background = element_rect(fill = NA, colour = "black"),
        #panel.grid = element_line(colour = "grey"),
        strip.background = element_blank(), strip.text.x = element_blank()) +
  labs(x="Sample fraction", y= 'Take-min \nstrata [%]') +
  #coord_cartesian(xlim=c(0.1,0.7)) +
  coord_cartesian(xlim=c(s1/N,s2/N)) +
  #scale_x_continuous(breaks = seq(s1/N,s2/N,0.1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(5)) +
  scale_y_continuous(breaks = scales::pretty_breaks(7))


pp <- 
  p1 + p2 + p3  + plot_layout(ncol = 1, heights = c(3, 1, 1))

ggsave(paste0("fig_",fig_n,"_rnabox.pdf"),pp,device="pdf", dpi=600, 
           width = 8, height = 9/1.618)


group_by(tab,algorithm) %>% summarise(mean(Median_time),mean(Mean_time))

