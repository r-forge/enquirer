\name{ENlisib}

\alias{ENlisib}

\title{Graph readability}

\description{
This function allows to improve the MCA factor map readability.
}

\usage{
ENlisib(res.mca,nbvar,nbind,axes=c(1,2))
}

\arguments{
  \item{res.mca}{results from an MCA}
  \item{nbvar}{chosen percentage or number for categories}
  \item{nbind}{chosen percentage or number for individuals}
  \item{axes}{selected axes, by default the two first axes are selected}
}

\details{
Categories with high v-test and individuals with high coordinates are kept. The chosen number and percentage concern individuals and categories beyond a threshold.
}


\author{
Sebastien Le \email{sebastien.le@agrocampus-ouest.fr} \cr
Jerome Bouche \email{jerome.bouche@yahoo.fr} \cr
Marine Cadoret \email{marine.cadoret@agrocampus-ouest.fr} \cr
Gwenaelle Fournier \email{gwenaelle.fournier@gmail.com} \cr
Olivier Fournier \email{olivier.fournier76@gmail.com} \cr
Francois Le Poder \email{francois.lepoder@gmail.com}
}

\examples{
\dontrun{
data(tea)
res.enmca<-ENMCA(tea[,1:18])
ENlisib(res.enmca$MCA,0.05,50,c(1,2))

res<-MCA(tea[,1:18])
ENlisib(res,0.05,50,c(1,2))
}
}