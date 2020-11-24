# ATAforecasting
Forecasting Time Series by ATA Method.

The ATA Method is a new alternative forecasting method. This method is alternative to two major forecasting approaches: Exponential Smoothing and ARIMA.

# Installation

Before installing ATAforecasting package, below packages are required for stR package:

Download lastest version of pandoc from https://pandoc.org/

Firstly, open Command Prompt (CMD.exe) and write down
```
pacman -S mingw-w64-{i686,x86_64}-freetype
```

Secondly, install below R packages.
```
install.packages("rgl", type = "source", INSTALL_opts = "--force-biarch")
```
or
```
install.packages("rgl", type = "source", INSTALL_opts = "--merge-multiarch")
```

Development version with latest features:
```
devtools::install_github("alsabtay/ATAforecasting")
devtools::install_github("alsabtay/ATAforecasting", upgrade_dependencies=FALSE)
```

# License
This package is free and open source software, licensed under GPL-3.
