require(asreml)
require(foreach)
require(testthat)
require(agridat)
data(SorghumYield)
data(SorghumCvGroup)

set.seed(123)

test_that("cv_groups works", {
  data("ortiz.tomato.yield")
  # Subset to obtain a character vector of unique environments
  unique_envs <- as.character(levels(ortiz.tomato.yield$env))

  # Partition the environments into 6 groups/folds for k-folds cross validation
  cv_df <- cv_groups(unique_envs, folds=6)
  expect_equal( as.integer(cv_df$cv_group),
                c(3, 3, 2, 4, 1, 3, 6, 5, 4,
                  1, 5, 2, 5, 1, 4, 2, 6, 6))
})






test_that("cv_groups works", {
  # Run baseline model
  baseline_asr <- asreml::asreml( Yld ~ Genotype + density + Genotype:density,
                                  random =~ at(Trial):Rep  + at(Trial):MainPlot +
                                    at(Trial,c('Breeza 1', 'Breeza 2', 'Emerald', 'Moree')):SubPlot +
                                    at(Trial,'Breeza 1'):Column +
                                    Trial + Env +
                                    spl(density, k=6) + spl(density, k=6):Genotype +
                                    str(~Trial:Genotype + Trial:Genotype:density,
                                        ~corh(2):id(48)) +
                                    str(~Env:Genotype + Env:Genotype:density,
                                        ~corh(2):id(136)),
                                  residual=~ dsum(~units|ResidualTOS),
                                  data = SorghumYield,
                                  na.action=na.method(x='include'),
                                  maxit=30, workspace="1Gb")
  ec_rmse_summary <- ec_cv_full(.fm=baseline_asr, .ec=rlang::expr(PrePAW), .G =rlang::expr(Genotype), .E = rlang::expr(Env),
                                .M = rlang::expr(density), .trial=rlang::expr(Trial), .env_cv_df=SorghumCvGroup)

  expect_equal( ec_rmse_summary$Rmse, 3.08660939)
})




