\name{GRmarking}

\alias{GRmarking}

\title{Graphical outputs for a semantic marking}

\description{
 This function gives graphical outputs for the variables, pairs or triplets of categories obtained with a semantic marking.
}

\usage{
GRmarking(res,lev,col.neg="lightblue",col.pos="pink",colour="black")
}

\arguments{
  \item{res}{the result of a semantic marking}
  \item{lev}{level of the semantic marking (1=variables, 2=pairs and 3=triplets)}
  \item{col.neg}{colour for the negative v-tests}
  \item{col.pos}{colour for the positive v-tests}
  \item{colour}{colour and font to be used for the characters}
}


\author{
Jerome Bouche \email{jerome.bouche@yahoo.fr} \cr
Gwenaelle Fournier \email{gwenaelle.fournier@gmail.com} \cr
Olivier Fournier \email{olivier.fournier76@gmail.com} \cr
Francois Le Poder \email{francois.lepoder@gmail.com}
}

\examples{
\dontrun{
data(tea)
res=ENmarking(tea[,1:21],21)
GRmarking(res,lev=2)
}
}
