\name{p_inertia}

\alias{p_inertia}

\title{Test of MCA's percentages of inertia}

\description{
 This function creates new data sets and performs MCA on them. It then keeps the percentages of inertia and calculates the means and the associated p-values on the first two dimensions.
}

\usage{
p_inertia(dataset)
}

\arguments{
  \item{dataset}{a data frame}
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
p_inertia(tea[,1:5])
}
}
