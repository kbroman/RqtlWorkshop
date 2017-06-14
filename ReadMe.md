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

The workshop will take place on Tuesday, 13 June, 9am-3pm, in Room 227
of the FedEx Institute. In the morning, we'll cover
[R/qtl](http://rqtl.org) and the foundations of QTL mapping. In the
afternoon, we'll cover [R/qtl2](http://kbroman.org/qtl2).

- 8:00-9:00  Light breakfast, 2nd floor FedEx Institute
- 9:00-noon  R/qtl (_break 10:20 - 10:40_)
- noon-1:00  Light lunch, 2nd floor FedEx Institute
- 1:00-3:00  R/qtl2 (_break 1:50 - 2:10_)

### Course materials

We'll switch back-and-forth between slides and follow-along software
demonstration.

- Source for the slides for the R/qtl part in [`Slides_Rqtl/`](Slides_Rqtl/)
  (PDF: [`Slides_Rqtl/rqtl_workshop.pdf`](Slides_Rqtl/rqtl_workshop.pdf))
- In demonstrating R/qtl, we'll follow the
  [Shorter tour of R/qtl](http://rqtl.org/tutorials/rqtltour2.pdf).
  My basic script is [`Scripts/rqtl_workshop.R`](Scripts/rqtl_workshop.R).
  Here's the version I actually typed at the workshop: [`Scripts/rqtl_workshop_live.R`](Scripts/rqtl_workshop_live.R).
- Source for the slides for the R/qtl2 part will be in [`Slides_Rqtl2/`](Slides_Rqtl2/)
  (PDF: [`Slides_Rqtl2/rqtl2_workshop.pdf`](Slides_Rqtl2/rqtl2_workshop.pdf))
- In demonstrating R/qtl2, we'll follow the
  [User guide](http://kbroman.org/qtl2/assets/vignettes/user_guide.html).
  My basic script is [`Scripts/rqtl2_workshop.R`](Scripts/rqtl2_workshop.R).
  Here's the version I actually typed at the workshop: [`Scripts/rqtl2_workshop_live.R`](Scripts/rqtl2_workshop_live.R).

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
