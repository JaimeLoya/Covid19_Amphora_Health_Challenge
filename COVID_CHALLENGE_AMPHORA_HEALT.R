# LOAD NECESSARY LIBRARIES TO PROCESS THE DATA 
library(dplyr)
library(magrittr)
library(lubridate)
library(ggplot2)
library(stringr)
library(cluster)
library(factoextra)

# SET DE WORKING DIRECTORY
setwd("D:\\AMPHORA")

################ 
################ 
##########     PROCESS DATA FROM MEXICO
###############
################ 
# READS THE MEXICO DATA OBTAINED FROM:
# https://www.gob.mx/salud/documentos/datos-abiertos-bases-historicas-direccion-general-de-epidemiologia
data_mexico<-read.csv("COVID_MEXICO.csv",encoding = "UTF-8")

# SELECT SPECIFIC COLUMNS WITH THE AIM TO REDUCE THE SIZE OF DATABASE AND CREATE A NEW DATA FRAME
# ALSO FILTER BY SYMPTOM DATE, AGE AND POSITIVE TEST
data_mexico<-select(data_mexico,ID_REGISTRO,SEXO,EDAD,
                    FECHA_SINTOMAS,FECHA_DEF,CLASIFICACION_FINAL,FECHA_ACTUALIZACION)%>%filter(FECHA_SINTOMAS>="2021-07-01"&
                                                                                                 FECHA_SINTOMAS<="2021-07-31"&
                                                                                                 EDAD>=0 & 
                                                                                                 EDAD<=100 &
                                                                                                 FECHA_DEF!="9999-99-99" &
                                                                                                 CLASIFICACION_FINAL==c(1,2,3))
# ESTIMATES THE DAYS BEFORE THE PATIENT DIE. IT WILL BE USED IN CLUSTER ANALYSIS
data_mexico$days_before_die<-yday(data_mexico$FECHA_DEF)-yday(data_mexico$FECHA_SINTOMAS)
# WRITE THE NEW DATA FRAME LIKE A CSV FORMAT 
write.csv(data_mexico,file="DATA_MEXICO_JULY2021.csv")


################ 
################ 
##########     PROCESS DATA FROM ARGENTINA
###############
################ 

# READS THE ARGENTINA DATA OBTAINED FROM:
# http://datos.salud.gob.ar/dataset/covid-19-casos-registrados-en-la-republica-argentina/archivo/fd657d02-a33a-498b-a91b-2ef1a68b8d16
data_arg<-read.csv("COVID_ARGENTINA.csv",encoding = "UTF-8")

# SELECT SPECIFIC COLUMNS AND FILTER BY DATE OF SYMPTOMS, DATE OF DEATH, AGE AND IF THE PERSON DIED
data_arg<-select(data_arg,id_evento_caso, sexo, edad, edad_años_meses,fecha_inicio_sintomas,fecha_fallecimiento,
                   clasificacion_resumen, ultima_actualizacion,fallecido)%>%filter(fecha_inicio_sintomas!="" &
                                                                         fecha_inicio_sintomas>="2021-07-01"&
                                                                         fecha_inicio_sintomas<="2021-07-31"&
                                                                         clasificacion_resumen=="Confirmado" &
                                                                         fallecido=="SI"&
                                                                         edad>=0 & edad <=100)
# ESTIMATES THE DAYS BEFORE THE PATIENT DIE. IT WILL BE USED IN CLUSTER ANALYSIS
data_arg$days_before_die<-yday(data_arg$fecha_fallecimiento)-yday(data_arg$fecha_inicio_sintomas)
# WRITE THE NEW DATA FRAME LIKE A CSV FORMAT 
write.csv(data_arg,file="DATA_ARGENTINA_JULY2021.csv")


################ 
################ 
##########     PROCESS DATA FROM COLOMBIA
###############
################ 

# READS THE COLOMBIA DATA OBTAINED FROM:
# https://www.datos.gov.co/Salud-y-Protecci-n-Social/Casos-positivos-de-COVID-19-en-Colombia/gt2j-8ykr
data_colombia<-read.csv("COVID_COLOMBIA.csv",encoding = "UTF-8")
# SELECT SPECIFIC COLUMNS WITH THE AIM TO REDUCE THE SIZE OF DATABASE AND CREATE A NEW DATA FRAME
# THE COLOMBIAN DATABASE REGISTER THE DATE WITH TIME, I NEED PROCESS IT IDEPENDIENTLY
data_colombia<-select(data_colombia,ID.de.caso,Sexo,Edad,Unidad.de.medida.de.edad,Fecha.de.inicio.de.síntomas,Fecha.de.muerte,
                      Recuperado,fecha.reporte.web)%>%filter(Recuperado=="Fallecido"|
                                                               Recuperado=="fallecido",Edad>=0 & Edad<=100)

#  CREATE A NEW COLUMN WITH THE DATE OF SYMTOMPS
data_colombia$fecha_sintomas<-as.Date(str_split_fixed(data_colombia$Fecha.de.inicio.de.síntomas," ",2)[,1],format ="%d/%m/%Y" )
# CREATE A NEW COLUMN WITH THE DATE OF DIE
data_colombia$fecha_muerte<-as.Date(str_split_fixed(data_colombia$Fecha.de.muerte," ",2)[,1],format="%d/%m/%Y")
# SELECT ONLY THE DATA FROM JULY 2021
data_colombia<-filter(data_colombia,fecha_sintomas=="2021-07-01" & fecha_sintomas<="2021-07-31")
# ESTIMATE THE DAYS BEFORE DIE
data_colombia$days_before_die<-yday(data_colombia$fecha_muerte)-yday(data_colombia$fecha_sintomas)
# WRITE THE NEW DATA FRAME LIKE A CSV FORMAT
write.csv(data_colombia,file="DATA_COLOMBIA_JULY2021.csv")

####READ THE SUBSET DATA FROM MEXICO, ARGENTINA AND COLOMBIA

data_mexico<-read.csv("DATA_MEXICO_JULY2021.csv")%>%select(EDAD,days_before_die)
colnames(data_mexico)<-c("age","days_before_die")
data_argentina<-read.csv("DATA_ARGENTINA_JULY2021.csv")%>%select(edad,days_before_die)
colnames(data_argentina)<-c("age","days_before_die")
data_colombia<-read.csv("DATA_COLOMBIA_JULY2021.csv")%>%select(Edad,days_before_die)
colnames(data_colombia)<-c("age","days_before_die")

data_merge_non_standarized<-rbind(data_mexico,data_argentina,data_colombia)

data_merge_standarized<-data.frame(scale(data_merge))


# Set the seed of R's random number generator, which is useful for creating simulations or random 
# objects that can be reproduced.
# Con el mismo número de seed el resultado siempre es el mismo
estimate_wcss<-function(data){
  set.seed(1234)
  wcss <- vector()
  for(i in 1:20){
    wcss[i] <- sum(kmeans(data, i)$withinss)
  }
  return(wcss)
}
# ESTIMATE THE WITHIN-CLUSTER SUM OF SQUARES (WCSS) IN STANDARIZED AND NON-STANDARIZED DATA
wcss_non_standarized<-estimate_wcss(data_merge_non_standarized)
wcss_standarized<-estimate_wcss(data_merge_standarized)
# SPECIFY HOW MANY GRAPHS TO PUT IN THE SAME IMAGE
par(mfrow=c(1,2)) 

#PLOT WCSS AND CUSTOMIZE THE GRAPHS
plot(wcss_non_standarized,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares",
     col = "dodgerblue",
     main="Elbow method for optimal K (non standardized values)",
     xaxt="n",
     yaxt="n")
     axis(1, at=seq(1, 20, 1),cex.axis=1)
     axis(2, at=seq(0, 1400000, 200000),cex.axis=1)


plot(wcss_standarized,
         type="b", pch = 19, frame = FALSE, 
         xlab="Number of clusters K",
         xlim=c(0,20),
         ylab="Total within-clusters sum of squares",
         col = "dodgerblue",
         main="Elbow method for optimal K (standardized values)",
     xaxt="n",
     yaxt="n")
axis(1, at=seq(1, 20, 1),cex.axis=1)
axis(2, at=seq(0, 10000, 2000),cex.axis=1)


set.seed(1234)
kmeans_non_standarized <- kmeans(data_merge_non_standarized, 4, iter.max = 10, nstart = 10)
data_merge_non_standarized$cluster <- kmeans_non_standarized$cluster


set.seed(1234)
kmeans_standarized <- kmeans(data_merge_standarized, 4, iter.max = 10, nstart = 10)
data_merge_standarized$cluster <- kmeans_standarized$cluster

# SAVE PLOT WITH NON STANDARDIZED VALUES
cluster_non_standarized<-ggplot() + geom_point(aes(x = age, y = days_before_die, color = cluster), data = data_merge_non_standarized, size = 2) +
  scale_colour_gradientn(colours=rainbow(4)) +
  geom_point(aes(alpha="Centroid",x = kmeans_non_standarized$centers[, 1], y = kmeans_non_standarized$centers[, 2])) +
  scale_alpha_manual(name=NULL, values=c(1,1),breaks = c("Centroid"))+
  ggtitle('Data clusters with Kmeans non-standarized (K=4)') + 
  scale_x_continuous("Age years", breaks = seq(0,100, by = 5))+
  scale_y_continuous("Days elapsed from symptoms to death", breaks = seq(0,32,2))+ #Draw y left axis 
  theme_bw()

# SAVE PLOT IN A VARIABLE WITH STANDARDIZED VALUES
cluster_standarized<-ggplot() + geom_point(aes(x = age, y = days_before_die, color = cluster), data = data_merge_standarized, size = 2) +
  scale_colour_gradientn(colours=rainbow(4)) +
  geom_point(aes(alpha="Centroid",x = kmeans_standarized$centers[, 1], y = kmeans_standarized$centers[, 2])) +
  scale_alpha_manual(name=NULL, values=c(1,1),breaks = c("Centroid"))+
  ggtitle('Data clusters with Kmeans standarized (K=4)') + 
  scale_x_continuous("Age years standarized", breaks = seq(-4,4, by = 1))+
  scale_y_continuous("Days standarized elapsed from symptoms to death", breaks = seq(-2,4,1))+ #Draw y left axis 
  theme_bw()

# PLOT BOTH GRAPHS IN THE SAME IMAGE
cluster_non_standarized+cluster_standarized

#CLUSTERING VALIDATION USING SILHOUETTE COEFFICIENT
#SILHOUETTE > 0 MEANS THAT THE OBSERVATIOS IS WELL CLUSTERED. THE CLOSEST IT IS TO 1, THE BEST IT IS CLUSTERED.
#SILHOUETTE < 0 EANS THAT THE OBSERVATIOS WAS PLACED IN THE WRONG CLUSTER.
#SILHOUETTE = 0 EANS THAT THE OBSERVATIOS IS BETWEEN TWO CLUSTERS.

# SILHOUETTE COEFFICIENTE TO CLUSTER VALIDATION (NON-STANDARIZED VALUES) 
sil_non_standarized <- silhouette(data_merge_non_standarized$cluster, dist(data_merge_non_standarized))
fviz_silhouette(sil_non_standarized)

# SILHOUETTE COEFFICIENTE TO CLUSTER VALIDATION (STANDARIZED VALUES) 
sil_standarized <- silhouette(data_merge_standarized$cluster, dist(data_merge_standarized))
fviz_silhouette(sil_standarized)

######################


# READS THE MEXICO DATA OBTAINED FROM:
# https://www.gob.mx/salud/documentos/datos-abiertos-bases-historicas-direccion-general-de-epidemiologia
data_mexico<-read.csv("COVID_MEXICO.csv",encoding = "UTF-8")

# SELECT SPECIFIC COLUMNS WITH THE AIM TO REDUCE THE SIZE OF DATABASE AND CREATE A NEW DATA FRAME
# ALSO FILTER BY SYMPTOM DATE, AGE AND POSITIVE TEST
data_mexico<-read.csv("COVID_MEXICO.csv",encoding = "UTF-8")
data_mexico<-select(data_mexico,ID_REGISTRO,SEXO,EDAD,NEUMONIA,DIABETES,
                    EPOC,ASMA,INMUSUPR,HIPERTENSION, CARDIOVASCULAR, OBESIDAD, RENAL_CRONICA,TABAQUISMO,OTRO_CASO,
                    FECHA_SINTOMAS,FECHA_DEF,CLASIFICACION_FINAL,FECHA_ACTUALIZACION)%>%filter(FECHA_SINTOMAS>="2021-07-01"&
                                                                                                 FECHA_SINTOMAS<="2021-07-31"&
                                                                                                 EDAD>=0 & 
                                                                                                 EDAD<=100 &
                                                                                                 FECHA_DEF!="9999-99-99" &
                                                                                                 CLASIFICACION_FINAL==c(1,2,3))
# REPLACE DATA VALUES 97, 98 AND 99 BY NA, IT IS SIMILAR TO  NULL
data_mexico[data_mexico==97|data_mexico==98|data_mexico==99]<-NA
# OMMIT NA VALUES
data_mexico<-na.omit(data_mexico)
data_mexico<-select(data_mexico,ID_REGISTRO,SEXO,NEUMONIA,DIABETES,EPOC,ASMA,INMUSUPR,
                           HIPERTENSION,CARDIOVASCULAR,OBESIDAD,RENAL_CRONICA,TABAQUISMO,OTRO_CASO)


# CORRESPONDENCE ANALYSIS ACCORDIG TO Nenadi´c & Greenacre (2007)
# https://www.jstatsoft.org/index.php/jss/article/view/v020i03/v20i03.pdf

# SELECT VARIABLES FROM COVID DATABASE
var_interes<-colnames(data_mexico)[c(2,3,5)]
data_interes<-data_mexico%>%select(var_interes)

# LOAD THE LIBRARY fastDummies TO DO THE DUMMIES VARIABLES
library(fastDummies)
N<-dummy_cols(data_interes,var_interes)
N<-N[,-c(1:3)]
# VERIFY THAT THE SUM OF THE COLUMN IS NEVER ZERO, OTHERWISE WE MUST DELETE SUCH COLUMN
# IN THIS CASE THE SUM OF COLUMN ALWAYS IS GREATER THAN ZERO
apply(N,2,sum)

#STEP 1
n<-sum(N) # SUM THE ELEMENTS IN THE MATRIX
P<-N/n  # ESTIMATES PROPORTIONS
P<-as.matrix(P) # CONVERT THE PROPORTION DATAFRAME AS MATRIX
addmargins(P,c(1,2)) # ADD THE MARGINALS 
rr<-margin.table(P,1) # CREATE A TABLE WITHE  THE MARGINALS IN ROWS
cc<-margin.table(P,2) # CREATE A TABLE WITHE  THE MARGINALS IN COLUMNS
S<-diag(rr^(-0.5))%*%(P-rr %*%t (cc)) %*% diag(cc^(-0.5)) # ESTIMATE THE RESIUDAL MATRIX (IT STANDARIZE THE DATA)

# STEP 2
u<-svd(S)$u # VECTOR SINGULAR OF ROWS
v<-svd(S)$v # VECTOR SINGULAR OF COLS
Da<-diag(svd(S)$d) # INERTIA = VARIANCE

#STEP 3
FF<-diag(rr^(-0.5)) %*% u %*% Da #ESTIMATE COORDINATES OF ROWS

# STEP 4 
GG<-diag(cc^(-0.5)) %*% v %*% Da #ESTIMATE COORDINATES OF COLUMNS

# PLOT
labels<-c("Mujer","Hombre","Con neumonía","Sin neumonía",
          "Con EPOC","Sin EPOC")

plot(GG[,1],GG[,2],type="n",xlim = c(-5,5),ylim=c(-3,3))
text(GG[,1],GG[,2],labels=labels,cex=0.7)
points(FF[,1],FF[,2],col="black")


######### ALTERNATIVE LIBRARY TO DO THE CORRESPONDECE ANALYSIS MORE DIRECT
# ONLY NEEDED CREATE THE N MATRIX
library(ca)
ca(N)
plot(ca(N))

######### ALTERNATIVE LIBRARY TO DO THE CORRESPONDENCE ANALYSIS, IN A DIRECT AND INTERACTIVE WAY
# ONLY NEEDED CREATE THE N MATRIX

library(FactoMineR)
library(explor)

## SELECT VARIABLES
var_interes<-colnames(data_mexico)[c(2:13)]     # HERE CAN SELECT SPECIFIC COLUMN POSITION
data_interes<-data_mexico%>%select(var_interes)


#library(fastDummies)
N<-dummy_cols(data_interes,var_interes)
N<-N[,-c(1:12)]

colnames(N)<-c("Mujer","Hombre","Con.neumonía","Sin.neumonía","Con.diabetes","Sin.diabetes",
               "Con.EPOC","Sin.EPOC","Con.asma","Sin.asma","Con.inmunosupresión","Sin.inmunosupresión",
               "Con.hipertension","Sin.hipertension","Enfermedad.cardiovascular","Sin.enfermedad.cardiovascular",
               "Con.obesidad","Sin.obesidad","Insuficiencia.renal","Sin.insuficiencia.renal",
               "Tabaquismo","Sin.tabaquismo","Contacto.COVID", "Sin.contacto.COVID")
N<-lapply(N, factor)
N<-data.frame(N)
row.names(N)<-data_mexico$ID_REGISTRO
mca<-MCA(N[,c(1:16)])
explor(mca)
