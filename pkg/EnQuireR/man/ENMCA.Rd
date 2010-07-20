\name{ENMCA}

\alias{ENMCA}

\title{MCA and cluster analysis}

\description{
This function allows to perform a cluster analysis following a MCA.
}

\usage{
ENMCA(dataset,ind.sup = NULL, quanti.sup = NULL, quali.sup = NULL, axes = c(1, 2), row.w = NULL,level.ventil = 0,signif=FALSE,proba=0.05,report=FALSE,language="english")
}

\arguments{
  \item{dataset}{current dataset}
  \item{ind.sup}{a vector indicating the indexes of the supplementary individuals}
  \item{quanti.sup}{a vector indicating the indexes of the supplementary quantitatives variables}
  \item{quali.sup}{a vector indicating the indexes of the supplementary qualitatives variables} 
  \item{axes}{a length 2 vector specifying the components to plot} 
  \item{row.w}{an optional row weights (by default, uniform row weights)}
  \item{level.ventil}{a number corresponding to the level under which the category is ventilated; by default 0, and no ventilation is done}  
  \item{signif}{a logical value. If TRUE, a p-value is associated to the 2 first dimensions}
  \item{proba}{Threshold used for the automatic description of the clusters. By default, proba=0.05}
  \item{report}{a logical value. If TRUE, a pdf report is produced using the Sweave package}
  \item{language}{a value of c("english","french") which determines the language to use for the report}
}


\value{
A list containing the following elements:
  \item{data}{dataset withe the new categorical variable, obtained with the cluster analysis}
  \item{catdes}{result of the catdes of the clusters}
  \item{MCA}{normal results of an MCA. On the MCA factor map, each cluster is displayed in a different colour}
}

\details{
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
res.enmca=ENMCA(tea[,1:18])
}
}