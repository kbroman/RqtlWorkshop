## R/qtl Workshop

This repository contains materials for a one-day workshop on
[R/qtl](http://rqtl.org) and [R/qtl2](http://kbroman.org/qtl2), to be
held at the
[Complex Trait Community meeting](http://complextrait.org/ctc2017/) in
Memphis, TN, on 13 June 2017.


### Preparations

Participants will need to bring their own laptop computer. Mac users
will need to be running a very recent version of Mac OS X: El Capitan
(10.11) or Sierra (10.12). To determine your installed version, click
the  in the upper-left and select "About this Mac".

- Install the latest version of [R](https://cran.r-project.org), 3.4.0

- Install the latest version of
  [RStudio](https://www.rstudio.com/products/rstudio/download/),
  1.0.143

- Open RStudio and install the latest version of
  [R/qtl](http://rqtl.org) and
  [R/qtlcharts](http://kbroman.org/qtlcharts) by typing, at the
  command prompt (`>`),

  ```{r}
  install.packages(c("qtl", "qtlcharts"))
  ```

- Install [R/qtl2](http://kbroman.org/qtl2) by typing

  ```{r}
  install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
  ```

### Schedule

The workshop is 9-3 with a lunch break 12-1. In the morning, we'll
cover [R/qtl](http://rqtl.org) and the foundations of QTL mapping. In
the afternoon, we'll cover [R/qtl2](http://kbroman.org/qtl2)

---

### License

To the extent possible under law,
[Karl Broman](http://github.com/kbroman) has waived all copyright and
related or neighboring rights to
&ldquo;[R/qtl workshop](https://github.com/kbroman/RqtlWorkshop)&rdquo;.
This work is published from the United States.
<br/>
[![CC0](http://i.creativecommons.org/p/zero/1.0/88x31.png)](http://creativecommons.org/publicdomain/zero/1.0/)

Code is licensed under the
[MIT license](https://cran.r-project.org/web/licenses/MIT).
([More information here](https://en.wikipedia.org/wiki/MIT_License).)
