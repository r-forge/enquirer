\name{ENellipse}

\alias{ENellipse}

\title{Confidence ellipses}

\description{
This function allows to draw ellipses around the categories of a chosen variable.
}

\usage{
ENellipse(matP, cluster, resampling="population", iter=500, alpha = 0.05, coord = c(1, 2), eig, cex = 1, color = NULL)
}

\arguments{
  \item{matP}{factorial coordinates from MCA}
  \item{cluster}{chosen categorical variable}
  \item{resampling}{If "population", the resampling is done within the whole dataset. If "cluster", it is done within the clusters}
  \item{iter}{number of iterations for the Bootstrap method, by default iter=500}
  \item{alpha}{the confidence level of the ellipses. By default alpha=0.05}
  \item{coord}{selected axes}
  \item{eig}{eigenvalues of the MCA}
  \item{cex}{graphical option: character size}
  \item{color}{graphical option: ellipses' colours}
}

\details{
This function can be used to draw ellipses around the centers of gravity of each cluster,
when the chosen categorical variable is the variable obtained after the cluster analysis.
}


\author{
Jerome Bouche \email{jerome.bouche@yahoo.fr} \cr
Marine Cadoret \email{marine.cadoret@agrocampus-ouest.fr} \cr
Gwenaelle Fournier \email{gwenaelle.fournier@gmail.com} \cr
Olivier Fournier \email{olivier.fournier76@gmail.com} \cr
Francois Le Poder \email{francois.lepoder@gmail.com}
}

\examples{
\dontrun{
data(tea)
res.enmca=ENMCA(tea[,1:18])
x<-res.enmca$MCA$ind$coord[,1:2]
b<-res.enmca$data[,ncol(res.enmca$data)]
ENellipse(matP=x, cluster=b, resampling="population", iter=500, alpha = 0.05, coord = c(1, 2), eig=res.enmca$MCA$eig, cex = 1, color = NULL)
}
}