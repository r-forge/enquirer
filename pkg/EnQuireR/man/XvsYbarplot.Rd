\name{XvsYbarplot}

\alias{XvsYbarplot}

\title{Barplot per category with respect to a categorical variable}

\description{
  This function allows to obtain coloured barplots of a categorical variable, depending on the categories of another categorical variable.
}

\usage{
XvsYbarplot(var1, var2, dataset, width = 1, space = NULL, names.arg = NULL, 
    legend.text = NULL, horiz = FALSE, density = NULL, angle = 45, 
    col = NULL, border = par("fg"), main = NULL, sub = NULL, 
    xlab = NULL, ylab = NULL, xlim = NULL, ylim = NULL, xpd = TRUE, 
    log = "", axes = TRUE, axisnames = TRUE, cex.axis = par("cex.axis"), 
    cex.names = par("cex.axis"), inside = TRUE, plot = TRUE, 
    axis.lty = 0, offset = 0, add = FALSE, ...) 

}

\arguments{
  \item{var1}{the name of the first variable.}
  \item{var2}{the name of the second variable.}
  \item{dataset}{the dataset to be studied.}
  \item{width}{optional vector of bar widths. Re-cycled to length the number of bars drawn. Specifying a single value will no visible effect unless xlim is specified.}
  \item{space}{the amount of space (as a fraction of the average bar width) left before each bar. May be given as a single number or one number per bar. If height is a matrix and beside is TRUE, space may be specified by two numbers, where the first is the space between bars in the same group, and the second the space between the groups. If not given explicitly, it defaults to c(0,1) if height is a matrix and beside is TRUE, and to 0.2 otherwise.}
  \item{names.arg}{a vector of names to be plotted below each bar or group of bars. If this argument is omitted, then the names are taken from the names attribute of height if this is a vector, or the column names if it is a matrix.}
  \item{legend.text}{a vector of text used to construct a legend for the plot, or a logical indicating whether a legend should be included. This is only useful when height is a matrix. In that case given legend labels should correspond to the rows of height; if legend.text is true, the row names of height will be used as labels if they are non-null.}
  \item{horiz}{a logical value. If FALSE, the bars are drawn vertically with the first bar to the left. If TRUE, the bars are drawn horizontally with the first at the bottom.}
  \item{density}{a vector giving the density of shading lines, in lines per inch, for the bars or bar components. The default value of NULL means that no shading lines are drawn. Non-positive values of density also inhibit the drawing of shading lines..}
  \item{angle}{the slope of shading lines, given as an angle in degrees (counter-clockwise), for the bars or bar components.}
  \item{col}{a vector of colors for the bars or bar components. By default, grey is used if height is a vector, and a gamma-corrected grey palette if height is a matrix.}
  \item{border}{the color to be used for the border of the bars. Use border = NA to omit borders. If there are shading lines, border = TRUE means use the same colour for the border as for the shading lines.}
  \item{main,sub}{overall and sub title for the plot.}
  \item{xlab}{a label for the x axis.}
  \item{ylab}{a label for the y axis.}
  \item{xlim}{limits for the x axis.}
  \item{ylim}{limits for the y axis.}
  \item{xpd}{logical. Should bars be allowed to go outside region?}
  \item{log}{string specifying if axis scales should be logarithmic; see \code{\link{plot.default}}.}
  \item{axes}{logical. If TRUE, a vertical (or horizontal, if horiz is true) axis is drawn.}
  \item{axisnames}{logical. If TRUE, and if there are names.arg (see above), the other axis is drawn (with lty=0) and labeled.}
  \item{cex.axis}{expansion factor for numeric axis labels.}
  \item{cex.names}{expansion factor for axis names (bar labels).}
  \item{inside}{logical. If TRUE, the lines which divide adjacent (non-stacked!) bars will be drawn. Only applies when space = 0 (which it partly is when beside = TRUE).}
  \item{plot}{logical. If FALSE, nothing is plotted.}
  \item{axis.lty}{the graphics parameter lty applied to the axis and tick marks of the categorical (default horizontal) axis. Note that by default the axis is suppressed.}
  \item{offset}{a vector indicating how much the bars should be shifted relative to the x axis.}
  \item{add}{logical specifying if bars should be added to an already existing plot; defaults to FALSE.}
  \item{...}{arguments to be passed to/from other methods. For the default method these can include further arguments (such as axes, asp and main) and graphical parameters (see par) which are passed to plot.window(), title() and axis.}
}

\details{
Each category has its own colour which is shaded inside the category 
according to the values of the categories of the other variable.
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
\code{\link{ENbarplot}} which does barplots for categorical variables.
}

\examples{
\dontrun{
data(tea)
XvsYbarplot(tea[,15],tea[,20],tea, legend.text=TRUE)
XvsYbarplot("socio.professional.category","sex",tea, legend.text=TRUE)
}
}
