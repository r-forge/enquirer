\name{chisq.desc}

\alias{chisq.desc}

\title{Display the results of a chi-square test}

\description{
Color the cells of the table containing the results from the Chi-square test, 
crossing all the selected categorical variables, according to a selected threshold.
}

\usage{
chisq.desc(dataset,Y,X,method="proba",print=TRUE,report=FALSE,language="english")
}

\arguments{
  \item{dataset}{a dataframe}
  \item{Y}{first group of variable}
  \item{X}{second group of variable (to be crossed with the first one)}
  \item{method}{if "proba" : the colour depends on the p-values given by the hypergeometric law. If "mean" : the colour depends on the value of the contribution to the chi-square distance }
  \item{print}{a logical value. If TRUE, tables of the chi-square distance are printed with cells coloured according to the chosen threshold}
  \item{report}{a logical value.  If TRUE, a report is automatically generated}
  \item{language}{english or french}
}

\details{
This function is useful when comparing categorical variables.
If report=TRUE, a directory named "EnQuireR" is created in your working directory (which can be known thanks to the getwd() function). This directory contains the final pdf report and also the teX file used to produce it.
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
chisq.desc(tea,c(1,2),c(3,4))

chisq.desc(tea,c(1,2),c(3,4),report=TRUE)
}
}
