## Quick start

```r
library(ccTools)

dat <- loadData()

ccBoxplots(dat)

alr <- analearn(dat)

ccPlotMI(alr)
