\name{ENmarking}

\alias{ENmarking}

\title{Semantic marking on categorical variables}

\description{
  Returns a semantic marking of categorical variables with 3 different levels (variables, pairs and triplets of categories). 
}

\usage{
ENmarking(dataset,var.int,proba=0.05)
}

\arguments{
  \item{dataset}{a data frame}
  \item{var.int}{the variable to mark}
  \item{proba}{Threshold to select variables, pairs and triplets. By default 0.05}
}


\author{
Jerome Bouche \email{jerome.bouche@yahoo.fr} \cr
Gwenaelle Fournier \email{gwenaelle.fournier@gmail.com} \cr
Olivier Fournier \email{olivier.fournier76@gmail.com} \cr
Francois Le Poder \email{francois.lepoder@gmail.com}
}

\seealso{ 
\code{\link{GRmarking}} which is a graphical representation of the semantic marking.
}

\examples{
\dontrun{
data(tea)
ENmarking(tea[,1:21],21)

res.enmca=ENMCA(tea[,1:18])
res.enmark=ENmarking(res.enmca$data,ncol(res.enmca$data))
GRmarking(res.enmark,1)
}
}