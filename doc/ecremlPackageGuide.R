## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=300,
  dev="jpeg",dev.args=list(type="cairo")
)
# opts_chunk$set(comment=NA, warning=FALSE,message=FALSE,fig.align='center',background='white',
#               dev="jpeg", dev.args=list(type="cairo"), dpi=300)




## ----include = TRUE, echo=FALSE, message=FALSE, warning=FALSE, results="markup", out.width='100%'----


knitr::include_graphics(path="ecFlowChart.png")


## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----load-packages, include = TRUE, message=FALSE, warning=FALSE, results="hide"----
# Load required R packages
library(asreml)
library(ecreml)
library(foreach)
library(rlang)
data(SorghumYield)


## ----load-packageBs, echo=FALSE, include = TRUE, message=FALSE, warning=FALSE, results="hide"----
# Load required R packages
library(bookdown)
library(bslib)
library(tidyverse)
library(viridis)
library(knitr)
library(quarto)


## ----baseline-model, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
baseline_asr <- asreml(
  fixed = Yld ~ Genotype + TargetPop + Genotype:TargetPop,
  random =~ at(Trial):Rep  + at(Trial):MainPlot +
   at(Trial,c('Breeza 1', 'Breeza 2', 'Emerald', 'Moree')):SubPlot +
   Trial + Trial:Genotype + Trial:Genotype:TargetPop +
   Env + Env:Genotype + Env:Genotype:TargetPop,
  residual=~ dsum(~units|Env),
  data = SorghumYield)


## ----include = TRUE, echo=FALSE, message=FALSE, warning=FALSE, results="markup", out.width='100%'----

knitr::include_graphics(path="ec_single.png")


## ----ec-single, eval=FALSE, include=TRUE , message=FALSE, warning=FALSE, results="hide"----
# data(SorghumCvGroup)
# prePAW_rmsep_asr <- ec_single(.fm = baseline_asr,
#                               .ec = PrePAW,
#                               .G = Genotype,
#                               .E = Env,
#                               .M = TargetPop,
#                               .trial= Trial,
#                         .env_cv_df=SorghumCvGroup)


## ----prepaw-full, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
prePAW_full_asr <-  ec_full_model_constructor(.fm = baseline_asr,
   .ec= PrePAW,
   .G = Genotype,
   .M = TargetPop)


## ----prepaw-full-call, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
prePAW_full_asr$call


## ----SorghumCvGroups, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
SorghumCvGroup


## ----prepaw-rmsep, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# prePAW_rmsep_asr$Rmsep


## ----include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, results="markup", out.width='100%'----

knitr::include_graphics(path="ec_finder.png")


## ----ec-finder, eval=FALSE, include=TRUE , message=FALSE, warning=FALSE, results="hide"----
# ec_search <- ec_finder(fm = baseline_asr,
#                        ECs = c(PrePAW:PostPAW, PreFlwEvap),
#                        G = Genotype,
#                        E = Env,
#                        M = TargetPop,
#                        trial = Trial,
#                        env_cv_df = SorghumCvGroup)


## ----ec-finder-results, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# ec_search$summary_ecs
# ec_search$ec_selected


## ----include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, results="markup", out.width='100%'----

knitr::include_graphics(path="simplify_ec_model.png")


## ----postpaw-full, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# postPAW_full_asr <-  ec_full_model_constructor(.fm = baseline_asr,
#                                                .ec = PostPAW,
#                                                .G = Genotype,
#                                                .M = TargetPop)
#


## ----postpaw-full-call, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# postPAW_full_asr$call


## ----postpaw-random, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# random_simplify_asr <- ec_random_model( .fm = postPAW_full_asr,
#                                         .ecs_in_model = PostPAW,
#                                         .G = Genotype,
#                                         .M = TargetPop)


## ----postpaw-random-call, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# random_simplify_asr$call


## ----postpaw-fixed, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# fixed_simplify_asr <- ec_fixed_model(.fm = random_simplify_asr,
#                                      .ecs_in_model = PostPAW,
#                                      .G = Genotype,
#                                      .M = TargetPop)


## ----postpaw-fixed-call, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# fixed_simplify_asr$call


## ----postpaw-simplify, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# simplify_asr <- simplify_ec_model(.fm = postPAW_full_asr,
#                                   .ecs_in_model = PostPAW,
#                                   .G = Genotype,
#                                   .M = TargetPop)


## ----postpaw-simplify-call, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# simplify_asr$call


## ----eval=FALSE, include=TRUE , message=FALSE, warning=FALSE, results="hide", fig.width=7, fig.height=7, out.width='100%'----
#
# # Load required graphics packages
# library(tidyverse)
# library(viridis)
#
# # Create a data frame consisting of all the unique combinations
# #   that we want to obtain predictions for
# PostPAWgeno_df <- unique(SorghumYield[ ,
#                     c('Env', 'Trial', 'Genotype', 'PostPAW')])
#
# # Obtain predictions for the data frame created previously
# PostPAWgeno_pred <- predict.asreml( simplify_asr, classify='Genotype:PostPAW',
#     levels=list('PostPAW'=PostPAWgeno_df$PostPAW,
#     'Genotype'=PostPAWgeno_df$Genotype), parallel=T)
#
# PostPAWgeno_pred$pvals
#
# # Plot the predicted values and 95% credible intervals
# EC_1_geno <-  ggplot( PostPAWgeno_pred$pvals , aes( x=PostPAW,
#   y = predicted.value, group = Genotype )) +
#   geom_ribbon( aes(ymin=predicted.value-2*std.error, #95% prediction interval
#       ymax=predicted.value+2*std.error, group=Genotype),
#       alpha=0.05, linetype=0) +
#   geom_line(colour="Blue", size=0.5) +
#   geom_point(data = SorghumYield , aes(x=PostPAW, y=Yld,colour=Env) ,
#       stat = "identity" , size=0.8) +
#   facet_wrap(~Genotype, ncol=2) +
#   xlab("Post-flowering plant avaliable water (mm)") +
#   ylab("Grain yield (t/ha)") +
#   theme_bw() +
#   theme(  panel.grid.minor = element_blank(),
#           axis.line.x = element_line(colour = "black"))
#
# # Print the ggplot object
# EC_1_geno
#


## ----postpaw-fixed-kr, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# fixed_simplify_asr <- ec_fixed_model( .fm = random_simplify_asr,
#                                       .ecs_in_model = PostPAW,
#                                       .G = Genotype, .M = TargetPop,
#                                       denDF = "numeric")


## ----include=TRUE, echo=FALSE, message=FALSE, warning=FALSE, results="markup", out.width='100%'----

knitr::include_graphics(path="ec_iteration.png")


## ----ec-iteration, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# ec1_model_asr <- ec_iteration(fm = baseline_asr,
#                               ECs = c(PrePAW:PostPAW, PreFlwEvap),
#                               G = Genotype,
#                               E = Env,
#                               M = TargetPop,
#                               trial = Trial,
#                               env_cv_df = SorghumCvGroup)


## ----ec-iteration-selected-ec, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# ec1_model_asr$selected_ec


## ----ec-iteration-call, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# ec1_model_asr$fm$call


## ----ec-iteration-rmsep, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# ec1_model_asr$rmsep


## ----ec-all, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="hide"----
# ec_model_asr <- ec_all(fm = baseline_asr,
#                        ECs = c(PrePAW:PostPAW, PreFlwEvap),
#                        G = Genotype,
#                        E = Env,
#                        M = TargetPop,
#                        trial = Trial,
#                        env_cv_df = SorghumCvGroup)


## ----ec-all-call, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# ec_model_asr$fm$call
# ec_model_asr$rmsep


## ----ec-all-varcomp, eval=FALSE, include = TRUE , message=FALSE, warning=FALSE, results="markup"----
# summary(ec_model_asr$fm)$varcomp

