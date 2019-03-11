
<!-- README.md is generated from README.Rmd. Please edit that file -->
Statistical Hypothesis Testing in R
===================================

`SHT` aims at providing a casket of tools for hypothesis testing procedures ranging from classical to modern techniques. We hope it not be used as a primary means of [***p*-hacking**](https://en.wikipedia.org/wiki/Data_dredging).

Installation
------------

`SHT` released version can be obtained from [CRAN](https://CRAN.R-project.org/package=SHT) with:

``` r
install.packages("SHT")
```

or the up-to-date development version from github:

``` r
## install.packages("devtools")
## library(devtools)
devtools::install_github("kisungyou/SHT")
```

List of Available Tests
-----------------------

We categorized available functions by their object of interest for better navigation.

-   Notations ***x*** and ***y*** refer to samples.
-   Authors are referred by last names. See the help page of each function for complete references.
-   ***k*-sample** means that the test is checking the *homogeneity* across multiple samples.
-   Function naming convention is {`type of test`.`test name`}, or {`type of test`.`year` `authors`}, where there are two or three authors, we took their initials as abbreviation or simply the last name of the first author otherwise.
-   **`usek1d`** and **`useknd`** lets a user to apply any *k*-sample tests for two-sample testings.
-   When ℝ<sup>*p*</sup> notation is used, it denotes ***multivariate*** procedures.

### 0. utilities

| function name | description                                               |
|---------------|-----------------------------------------------------------|
| `usek1d`      | apply *k*-sample test method for two univariate samples   |
| `useknd`      | apply *k*-sample test method for two multivariate samples |

### 1. tests for univariate mean *μ* ∈ ℝ

<table style="width:71%;">
<colgroup>
<col width="11%" />
<col width="34%" />
<col width="25%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>mean1.ttest</code></td>
<td><a href="https://www.jstor.org/stable/2331554">Student (1908)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> {≤, =, ≥} <em>μ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="even">
<td><code>mean2.ttest</code></td>
<td><a href="https://www.jstor.org/stable/2331554">Student (1908)</a> &amp; <a href="https://doi.org/10.1093/biomet/34.1-2.28">Welch (1947)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> {≤, =, ≥} <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>meank.anova</code></td>
<td>-</td>
<td align="left"><span class="math inline"><em>μ</em><sub>1</sub> = ⋯ = <em>μ</em><sub><em>k</em></sub></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
</tbody>
</table>

### 2. tests for multivariate mean *μ* ∈ ℝ<sup>*p*</sup>

<table style="width:82%;">
<colgroup>
<col width="22%" />
<col width="34%" />
<col width="25%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>mean1.1931Hotelling</code></td>
<td><a href="https://projecteuclid.org/euclid.aoms/1177732979">Hotelling (1931)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="even">
<td><code>mean1.1958Dempster</code></td>
<td><a href="https://projecteuclid.org/euclid.aoms/1177706437">Dempster (1958, 1960)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="odd">
<td><code>mean1.1996BS</code></td>
<td><a href="http://www3.stat.sinica.edu.tw/statistica/j6n2/j6n21/j6n21.htm">Bai and Saranadasa (1996)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="even">
<td><code>mean1.2008SD</code></td>
<td><a href="https://doi.org/10.1016/j.jmva.2006.11.002">Srivastava and Du (2008)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="odd">
<td><code>mean2.1931Hotelling</code></td>
<td><a href="https://projecteuclid.org/euclid.aoms/1177732979">Hotelling (1931)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>mean2.1958Dempster</code></td>
<td><a href="https://projecteuclid.org/euclid.aoms/1177706437">Dempster (1958, 1960)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>mean2.1965Yao</code></td>
<td><a href="https://www.jstor.org/stable/2333819">Yao (1965)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>mean2.1980Johansen</code></td>
<td><a href="https://doi.org/10.1093/biomet/67.1.85">Johansen (1980)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>mean2.1986NVM</code></td>
<td><a href="https://doi.org/10.1080/03610928608829342">Nel and Van der Merwe (1986)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>mean2.1996BS</code></td>
<td><a href="http://www3.stat.sinica.edu.tw/statistica/j6n2/j6n21/j6n21.htm">Bai and Saranadasa (1996)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>mean2.2004KY</code></td>
<td><a href="https://doi.org/10.1016/j.spl.2003.10.012">Krishnamoorthy and Yu (2004)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>mean2.2008SD</code></td>
<td><a href="https://doi.org/10.1016/j.jmva.2006.11.002">Srivastava and Du (2008)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>mean2.2011LJW</code></td>
<td><a href="https://papers.nips.cc/paper/4260-a-more-powerful-two-sample-test-in-high-dimensions-using-random-projection">Lopes, Jacob, and Wainwright (2011)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>mean2.2014CLX</code></td>
<td><a href="https://doi.org/10.1111/rssb.12034">Cai, Liu, and Xia (2014)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>mean2.2014Thulin</code></td>
<td><a href="https://doi.org/10.1016/j.csda.2013.12.003">Thulin (2014)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>meank.2007Schott</code></td>
<td><a href="https://doi.org/10.1016/j.jmva.2006.11.007">Schott (2007)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub>1</sub> = ⋯ = <em>μ</em><sub><em>k</em></sub></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
<tr class="odd">
<td><code>meank.2009ZX</code></td>
<td><a href="https://doi.org/10.1007/s11425-009-0091-x">Zhang and Xu (2009)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub>1</sub> = ⋯ = <em>μ</em><sub><em>k</em></sub></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
<tr class="even">
<td><code>meank.2019CPH</code></td>
<td><a href="https://doi.org/10.1016/j.jspi.2018.12.002">Cao, Park, and He (2019)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub>1</sub> = ⋯ = <em>μ</em><sub><em>k</em></sub></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
</tbody>
</table>

### 3. tests for variance *σ*<sup>2</sup>

<table style="width:78%;">
<colgroup>
<col width="22%" />
<col width="26%" />
<col width="29%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>var1.chisq</code></td>
<td>-</td>
<td align="left"><span class="math inline"><em>σ</em><sub><em>x</em></sub><sup>2</sup> {≤, =, ≥} <em>σ</em><sub>0</sub><sup>2</sup></span> (1-sample)</td>
</tr>
<tr class="even">
<td><code>var2.F</code></td>
<td>-</td>
<td align="left"><span class="math inline"><em>σ</em><sub><em>x</em></sub><sup>2</sup> {≤, =, ≥} <em>σ</em><sub><em>y</em></sub><sup>2</sup></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>vark.1937Bartlett</code></td>
<td><a href="https://doi.org/10.1098/rspa.1937.0109">Bartlett (1937)</a></td>
<td align="left"><span class="math inline"><em>σ</em><sub>1</sub><sup>2</sup> = ⋯ = <em>σ</em><sub><em>k</em></sub><sup>2</sup></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
<tr class="even">
<td><code>vark.1960Levene</code></td>
<td>Levene (1960)</td>
<td align="left"><span class="math inline"><em>σ</em><sub>1</sub><sup>2</sup> = ⋯ = <em>σ</em><sub><em>k</em></sub><sup>2</sup></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
<tr class="odd">
<td><code>vark.1974BF</code></td>
<td><a href="https://www.jstor.org/stable/2285659">Brown and Forsythe (1974)</a></td>
<td align="left"><span class="math inline"><em>σ</em><sub>1</sub><sup>2</sup> = ⋯ = <em>σ</em><sub><em>k</em></sub><sup>2</sup></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
</tbody>
</table>

### 4. tests for covariance *Σ*

<table style="width:72%;">
<colgroup>
<col width="22%" />
<col width="25%" />
<col width="25%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>cov1.2012Fisher</code></td>
<td><a href="https://doi.org/10.1016/j.jspi.2011.07.019">Fisher (2012)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="even">
<td><code>cov1.2015WL</code></td>
<td><a href="https://arxiv.org/abs/1511.01611">Wu and Li (2015)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub>0</sub></span> (1-sample)</td>
</tr>
<tr class="odd">
<td><code>cov2.2012LC</code></td>
<td><a href="https://projecteuclid.org/euclid.aos/1338515142">Li and Chen (2012)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>cov2.2013CLX</code></td>
<td><a href="https://doi.org/10.1080/01621459.2012.758041">Cai, Liu, and Xia (2013)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="odd">
<td><code>cov2.2015WL</code></td>
<td><a href="https://arxiv.org/abs/1511.01611">Wu and Li (2015)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
<tr class="even">
<td><code>covk.2001Schott</code></td>
<td><a href="https://doi.org/10.1016/j.csda.2007.03.004">Schott (2001)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub>1</sub> = ⋯ = <em>Σ</em><sub><em>k</em></sub></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
<tr class="odd">
<td><code>covk.2007Schott</code></td>
<td><a href="https://doi.org/10.1016/j.csda.2007.03.004">Schott (2007)</a></td>
<td align="left"><span class="math inline"><em>Σ</em><sub>1</sub> = ⋯ = <em>Σ</em><sub><em>k</em></sub></span> (<span class="math inline"><em>k</em></span>-sample)</td>
</tr>
</tbody>
</table>

### 5. simultaneous tests for mean *μ* and covariance *Σ*

<table style="width:82%;">
<colgroup>
<col width="22%" />
<col width="34%" />
<col width="25%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>sim1.2017Liu</code></td>
<td><a href="https://doi.org/10.1016/j.jspi.2017.03.009">Liu et al. (2017)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub>,  <em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub><em>y</em></sub></span> (1-sample)</td>
</tr>
<tr class="even">
<td><code>sim2.2018HN</code></td>
<td><a href="https://doi.org/10.1007/s11749-017-0567-x">Hyodo and Nishiyama (2018)</a></td>
<td align="left"><span class="math inline"><em>μ</em><sub><em>x</em></sub> = <em>μ</em><sub><em>y</em></sub>,  <em>Σ</em><sub><em>x</em></sub> = <em>Σ</em><sub><em>y</em></sub></span> (2-sample)</td>
</tr>
</tbody>
</table>

### 6. tests for equality of distributions

<table style="width:81%;">
<colgroup>
<col width="22%" />
<col width="30%" />
<col width="27%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>eqdist.2014BG</code></td>
<td><a href="https://doi.org/10.1016/j.jmva.2013.09.004">Biswas and Ghosh (2014)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = <em>F</em><sub><em>Y</em></sub> ∈ ℝ<sup>1</sup> &amp; ℝ<sup><em>p</em></sup></span> (2-sample)</td>
</tr>
</tbody>
</table>

### 7. goodness-of-fit tests of normality

<table style="width:82%;">
<colgroup>
<col width="22%" />
<col width="34%" />
<col width="25%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>norm.1965SW</code></td>
<td><a href="https://doi.org/10.1093/biomet/52.3-4.591">Shapiro and Wilk (1965)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Normal ∈ ℝ<sup>1</sup></span></td>
</tr>
<tr class="even">
<td><code>norm.1972SF</code></td>
<td><a href="https://doi.org/10.1080/01621459.1972.10481232">Shapiro and Francia (1972)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Normal ∈ ℝ<sup>1</sup></span></td>
</tr>
<tr class="odd">
<td><code>norm.1980JB</code></td>
<td><a href="https://doi.org/10.1016/0165-1765(80)90024-5">Jarque and Bera (1980)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Normal ∈ ℝ<sup>1</sup></span></td>
</tr>
<tr class="even">
<td><code>norm.1996AJB</code></td>
<td><a href="https://doi.org/10.1016/S0165-1765(96)00923-8">Urzua (1996)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Normal ∈ ℝ<sup>1</sup></span></td>
</tr>
<tr class="odd">
<td><code>norm.2008RJB</code></td>
<td><a href="https://doi.org/10.1016/j.econlet.2007.05.022">Gel and Gastwirth (2008)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Normal ∈ ℝ<sup>1</sup></span></td>
</tr>
</tbody>
</table>

### 8. goodness-of-fit tests of uniformity

<table style="width:82%;">
<colgroup>
<col width="22%" />
<col width="34%" />
<col width="25%" />
</colgroup>
<thead>
<tr class="header">
<th>function name</th>
<th>authors</th>
<th align="left">description of <span class="math inline"><em>H</em><sub>0</sub></span></th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td><code>unif.2017YMi</code></td>
<td><a href="https://doi.org/10.1007/s00362-015-0715-x">Yang and Modarres (2017)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Uniform ∈ ℝ<sup><em>p</em></sup></span></td>
</tr>
<tr class="even">
<td><code>unif.2017YMq</code></td>
<td><a href="https://doi.org/10.1007/s00362-015-0715-x">Yang and Modarres (2017)</a></td>
<td align="left"><span class="math inline"><em>F</em><sub><em>X</em></sub> = Uniform ∈ ℝ<sup><em>p</em></sup></span></td>
</tr>
</tbody>
</table>

<!---
your comment goes here
and here
| `cov1.mxPBF`  | - | $\Sigma_x = \Sigma_0$ (1-sample) |
| `cov2.mxPBF`    | - | $\Sigma_x = \Sigma_y$ (2-sample) |
| `mean2.mxPBF` | - | $\mu_x = \mu_y$ (2-sample) |
-->
