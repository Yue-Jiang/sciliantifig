# sciliantifig

This is an R package that keeps track of analysis scripts, functions, `Rmarkdown`s and immediate input data for the final figure generation in the sci-LIANTI manuscript ([bioRxiv](https://www.biorxiv.org/content/early/2018/06/12/338053), [published version](https://www.cell.com/molecular-cell/pdfExtended/S1097-2765(19)30618-5)).

# Install

1. Clone the remote repo to current the working directory. `git clone https://github.com/Yue-Jiang/sciliantifig.git`.

2. Install.

- Option 1, with local R installation. `install.packages("sciliantifig/", repos=NULL, type="source")`.

- Option 2, with Docker. This removes the need for installing `R`, `Rstudio` and dependent libraries, but adds the need to have a working docker installation.

```
cd sciliantifig
pull yuejiang/sciliantifig:3.5.1 # pull image from dockerhub
docker run --rm -p 8787:8787 --name sciliantifig -v `pwd`:/home/rstudio/sciliantifig yuejiang/sciliantifig:3.5.1 # run the image and make the current working directory available to docker
```

Go to `http://localhost:8787/` in a web browser and it should provide a running instance of `Rstudio`. The username and password are both `rstudio`. Click on `sciliantifig.Rproj` in the files panel to open this project. Then press Ctrl/Cmd + Shift + B. This builds and installs the package, then restarts R and (re)loads the package.

# Usage

The figures in the manuscript can be reproduced by the `Rmarkdown` file `inst/sci-lianti-figures.Rmd`. A sample output is included as `inst/sci-lianti-figures.html`. See also https://yue-jiang.github.io/sciliantifig/. Its input data are in `inst/extdata`.
