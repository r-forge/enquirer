\name{ENbarplot}

\alias{ENbarplot}

\title{Plots of categorical variables}

\description{
 This function allows to obtain barplots either sorted by frequency or by alphabetical order. 
The percentage of missing values is indicated at the top of the graphical device.
}

\usage{
ENbarplot(dataset,X,spl=FALSE,numr=NULL,numc=NULL,report=FALSE)
}

\arguments{
  \item{dataset}{a data frame}
  \item{X}{the numbers of the variables to be plotted}
  \item{spl}{if TRUE, the barplots will be sorted per length}
  \item{numr}{the number of plots per row (by default 2)}
  \item{numc}{the number of plots per column (by default 2)}
  \item{report}{a logical value. If TRUE, a report is automatically generated}
}

\details{
The colour changes depending on the number of individuals of the category. 
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

\seealso{ 
\code{\link{XvsYbarplot}} which gives barplots of a categorical variable depending on the categoricals of another categorical variable.
}

\examples{
\dontrun{
data(tea)
ENbarplot(tea,c(18,20,21,23))

ENbarplot(tea,c(18,20,21,23),report=TRUE)
}
}

\keyword{univar}