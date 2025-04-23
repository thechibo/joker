# Testing

The `joker` package is extensively tested in order to make sure everything works
as intended and meet the rOpenSci statistical software development standards.
More than 1000 tests are implemented in the package.

# Extended Tests

Tests concerning the consistency and asymptotic normality of the estimators,
i.e. tests that verify that the estimators and their asymptotic variance are
coded correctly in the package can be switched on/off by using the
`JOKER_EXTENDED_TESTS="true"` environmental variable. These tests perform Monte
Carlo simulations and therefore require some not-insignificant amount of time.
A rough estimate would be that the (approximately 1000) basic tests take about
30 seconds, while the (approximately 50) extended tests require 5 extra minutes.
