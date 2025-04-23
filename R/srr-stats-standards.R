#' srr_stats
#'
#' @srrstatsVerbose TRUE
#'
#' @srrstatsTODO {G2.4} *Provide appropriate mechanisms to convert between
#' different data types, potentially including:*
#' @srrstatsTODO {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstatsTODO {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstatsTODO {G2.4c} *explicit conversion to character via `as.character()`
#' (and not `paste` or `paste0`)*
#' @srrstatsTODO {G2.4d} *explicit conversion to factor via `as.factor()`*

#' @srrstatsTODO {G2.6} Software which accepts one-dimensional input should
#' ensure values are appropriately pre-processed regardless of class structures.
#' @srrstatsTODO {G2.7} *Software should accept as input as many of the above
#' standard tabular forms as possible, including extension to domain-specific
#' forms.*
#' @srrstatsTODO {G2.9} *Software should issue diagnostic messages for type
#' conversion in which information is lost (such as conversion of variables from
#' factor to character; standardisation of variable names; or removal of
#' meta-data such as those associated with
#' [`sf`-format](https://r-spatial.github.io/sf/) data) or added (such as
#' insertion of variable or column names where none were provided).*
#' @srrstatsTODO {G2.10} *Software should ensure that extraction or filtering of
#' single columns from tabular inputs should not presume any particular default
#' behaviour, and should ensure all column-extraction operations behave
#' consistently regardless of the class of tabular data used as input.*
#' @srrstatsTODO {G2.11} *Software should ensure that `data.frame`-like tabular
#' objects which have columns which do not themselves have standard class
#' attributes (typically, `vector`) are appropriately processed, and do not
#' error without reason. This behaviour should be tested. Again, columns created
#' by the [`units` package](https://github.com/r-quantities/units/) provide a
#' good test case.*
#' @srrstatsTODO {G2.12} *Software should ensure that `data.frame`-like tabular
#' objects which have list columns should ensure that those columns are
#' appropriately pre-processed either through being removed, converted to
#' equivalent vector columns where appropriate, or some other appropriate
#' treatment such as an informative error. This behaviour should be tested.*

#' @srrstatsTODO {G5.4} **Correctness tests** *to test that statistical
#' algorithms produce expected results to some fixed test data sets (potentially
#' through comparisons using binding frameworks such as
#' [RStata](https://github.com/lbraglia/RStata)).*
#' @srrstatsTODO {G5.4a} *For new methods, it can be difficult to separate out
#' correctness of the method from the correctness of the implementation, as
#' there may not be reference for comparison. In this case, testing may be
#' implemented against simple, trivial cases or against multiple implementations
#' such as an initial R implementation compared with results from a C/C++
#' implementation.*
#' @srrstatsTODO {G5.4b} *For new implementations of existing methods,
#' correctness tests should include tests against previous implementations. Such
#' testing may explicitly call those implementations in testing, preferably from
#' fixed-versions of other software, or use stored outputs from those where that
#' is not possible.*
#' @srrstatsTODO {G5.5} *Correctness tests should be run with a fixed random
#' seed*

#' @srrstatsTODO {G5.8} **Edge condition tests** *to test that these conditions
#' produce expected behaviour such as clear warnings or errors when confronted
#' with data with extreme properties including but not limited to:*
#' @srrstatsTODO {G5.8a} *Zero-length data*
#' @srrstatsTODO {G5.8b} *Data of unsupported types (e.g., character or complex
#' numbers in for functions designed only for numeric data)*
#' @srrstatsTODO {G5.8c} *Data with all-`NA` fields or columns or all identical
#' fields or columns*
#' @srrstatsTODO {G5.8d} *Data outside the scope of the algorithm (for example,
#' data with more fields (columns) than observations (rows) for some regression
#' algorithms)*

#' @srrstatsTODO {G5.9} **Noise susceptibility tests** *Packages should test for
#' expected stochastic behaviour, such as through the following conditions:*
#' @srrstatsTODO {G5.9a} *Adding trivial noise (for example, at the scale of
#' `.Machine$double.eps`) to data does not meaningfully change results*
#' @srrstatsTODO {G5.9b} *Running under different random seeds or initial
#' conditions does not meaningfully change results*
#' @noRd
NULL

#' NA_standards
#'
#' @noRd
#' @srrstatsΝΑ {G2.4e, G2.5} No inputs are expected to be of `factor` type.
#' @srrstatsNA {G3.1, G3.1a} The package does not include ambiguous covariance
#' calculations. All variance-covariance matrices are appropriately defined in
#' the relevant literature.
#' @srrstatsNA {G5.0} Tests rely on simulated data so that the true parameter
#' values are known. Using real data sets is applicable but not practicable.
#' @srrstatsNA {G5.4c} Where applicable, stored values may be drawn from
#' published paper outputs when applicable and where code from original
#' implementations is not available
#' @srrstatsNA {G5.11, G5.11a} No data download or relevant procedure is
#' included in the package tests.
