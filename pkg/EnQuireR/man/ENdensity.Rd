\name{ENdensity}

\alias{ENdensity}

\title{Density levels}

\description{
  Display density levels around clusters on the cloud of the individuals or around a cloud of individuals. 
}

\usage{
ENdensity(res,classe=NULL,nlevels=20)
}

\arguments{
  \item{res}{the result of a factorial analysis}
  \item{classe}{the variable of interest. Density levels will be drawn around its categories}
  \item{nlevels}{number of contour levels desired}
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
cluster<-res.enmca$data[,ncol(res.enmca$data)]
X<-res.enmca$MCA
ENdensity(X,cluster,10)
}
}