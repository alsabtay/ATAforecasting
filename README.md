# ATAforecasting <img src="man/figures/logo.png" align="right" alt="ATAforecasting logo" />

[![CRAN](https://www.r-pkg.org/badges/version/ATAforecasting)](https://cran.r-project.org/package=ATAforecasting)
[![Downloads](https://cranlogs.r-pkg.org/badges/ATAforecasting)](https://cran.r-project.org/package=ATAforecasting)
[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://opensource.org/licenses/GPL-3.0)
## Synopsis

Automatic Time Series Analysis and Forecasting using Ata Method with Box-Cox Power Transformations Family and Seasonal Decomposition Techniques.

## Description

The Ata Method is a new alternative forecasting method. This method is alternative to two major forecasting approaches: Exponential Smoothing and ARIMA.
The Ata method based on the modified simple exponential smoothing as described in Yapar, G. (2016) [doi:10.15672/HJMS.201614320580](https://doi.org/10.15672/HJMS.201614320580), 
Yapar G., Capar, S., Selamlar, H. T., Yavuz, I. (2017) [doi:10.15672/HJMS.2017.493](https://doi.org/10.15672/HJMS.2017.493) and Yapar G., Selamlar, H. T., Capar, S., Yavuz, I. (2019) 
[doi:10.15672/hujms.461032](https://doi.org/10.15672/hujms.461032) is a new univariate time series forecasting method which provides innovative solutions to issues faced during 
the initialization and optimization stages of existing methods.

Forecasting performance of the Ata method is superior to existing methods both in terms of easy implementation and accurate forecasting. 
It can be applied to non-seasonal or seasonal time series which can be decomposed into four components (remainder, level, trend and seasonal).
This methodology performed well on the M3 and M4-competition dat

## Installation

You can install the **stable** version from
[CRAN](https://cran.r-project.org/package=ATAforecasting):

``` 
install.packages("ATAforecasting")
```

Development version with latest features:
```
devtools::install_github("alsabtay/ATAforecasting")
```

Fable Modelling Wrappers for ATAforecasting Package
```
devtools::install_github("alsabtay/fable.ata")
```

## Example

USAccDeaths: Accidental Deaths in the US 1973--1978

``` r
library(ATAforecasting)
ATA(USAccDeaths, h = 18, model.type = "A", seasonal.type = "A", seasonal.model = "stl")
``` 

## Links

[Github page](https://github.com/alsabtay/ATAforecasting)

[Github.io page](https://alsabtay.github.io/ATAforecasting/index.html)

[Project team website](https://atamethod.wordpress.com/)

[Github - Fable Modelling Wrappers for ATAforecasting Package](https://github.com/alsabtay/fable.ata)

[Github.io - Fable Modelling Wrappers for ATAforecasting Package](https://alsabtay.github.io/fable.ata/index.html)

[Github - Intermittent Ata Method Package](https://github.com/alsabtay/intermittentATA)

[Github.io Intermittent Ata Method Package](https://alsabtay.github.io/intermittentATA/index.html)


## Optional Dependencies

ATAforecasting provides multiple seasonal decomposition methods:

**Always available:**
- `decomp` - Classical decomposition (stats package)
- `stl` - STL decomposition (stats/forecast packages)
- `stlplus` - Enhanced STL (stlplus package)
- `tbats` - TBATS model (forecast package)
- `stR` - Regression-based STR (stR package)

**Optional (requires additional installation):**
- `x13` - X-13ARIMA-SEATS (requires: seasonal, x13binary)
- `x11` - X11 decomposition (requires: seasonal, x13binary)

To use X13/X11 methods:
```r
install.packages(c("seasonal", "x13binary"))
```

```
**Note:** x13binary contains software from the U.S. Census Bureau with 
[special licensing terms](https://github.com/x13org/x13binary). 
Users outside the United States should review the license before installation.
```

## License Information and Disclaimer
The ATAforecasting R Package was created by Ali Sabri Taylan, Hanife Taylan Selamlar, Guckan Yapar and is licensed under GPL-3.

Optional Dependencies
=====================
The x13 and x11 seasonal decomposition features optionally use x13binary,
which contains X-13ARIMA-SEATS software from the U.S. Census Bureau.

The X-13ARIMA-SEATS software is not subject to copyright in the United States
but has special licensing terms for use outside the U.S.

Users of x13/x11 features should review the licensing terms at:
https://github.com/x13org/x13binary

All other features of ATAforecasting do not depend on x13binary.

As stated in the manual of
[X-13ARIMA-SEATS](https://web.archive.org/web/20250412173420/https://www2.census.gov/software/x-13arima-seats/x13as/windows/documentation/docx13as.pdf)
(June 12, 2023):

> This Software was created by U.S. Government employees and therefore is not
> subject to copyright in the United States (17 U.S.C. §105). The United
> States/U.S. Department of Commerce (“Commerce”) reserve all rights to seek and
> obtain copyright protection in countries other than the United States. The
> United States/Commerce hereby grant to User a royalty-free, nonexclusive
> license to use, copy, and create derivative works of the Software outside of
> the United States.

> The Software is provided to the User and those who may take by, through or
> under it, “as is,” without any warranty (whether express or implied) or
> representation whatsoever, including but not limited to any warranty of
> merchantability. The Software is taken hereunder without any right to support
> or to any improvements, extensions, or modifications, except as may be agreed
> to separately, in writing, by Commerce.

> User, on behalf of itself and all others who take by, through or under it,
> hereby and forever waives, releases, and discharges the United States/Commerce
> and all its instrumentalities from any and all liabilities and obligations in
> connection with the use, application, sale or conveyance of the Software. User
> shall indemnify and hold harmless the United States/Commerce and its
> instrumentalities from all claims, liabilities, demands, damages, expenses,
> and losses arising from or in connection with User's use, application, sale or
> conveyance of the Software, including those who take by, through or under User
> whether or not User was directly involved. This provision will survive
> termination of this Agreement and will include any and all claims or
> liabilities arising under intellectual property rights, such as patents,
> copyrights, trademarks, and trade secrets. If User of software is an Executive
> Agency of the United States, this clause is not applicable.

> The construction, validity, performance, and effect of this Agreement for all
> purposes will be governed by Federal law of the United States.

> User agrees to make a good faith effort to use the Software in a way that does
> not cause damage, harm, or embarrassment to the United States/Commerce. The
> United States/Commerce expressly reserve all rights and remedies.

