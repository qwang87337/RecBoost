# RecBoost
This is an example code of RecBoost, a component-wise gradient boost method for recurrent event analysis. It integrates the component-wise gradient boosting method with the EM algorithm, 
and thus enables both variable selection and enhanced model interpretability. For methodological details, please refer to Wang et al., 202x.


## How to use
Please download all .R files in the R folder and run the __Main__ file. Make sure to install the required packages and import all the auxiliary functions. A simulated sample data is also provided for the purpose of demonstration.

### required packages
```{r}
# survival package
install.packages("survival")

# R package for component-wise gradient boost
install.packages("mboost")

# other helpful packages
install.packages("dplyr")
install.packages("Rcpp")
install.packages("MASS")
```

### description of each .R file
`__main__` the main code to be run for demonstration.

`em_arguments.R` contains functions that are used for .

`em_aux.R` contains some auxiliary functions needed for performing the method.

`emboost_fit.R` contains the function that does the EM algorithm and boosting iteration.

`emboost_fit_nl.R` similar to `emboost_fit.R`, but designed explicitly for nonlinear settings.

`fast_Estep.R` Under the assumption of Gamma-distributed frailty, a closed form is available for conditional complete-data log likelihood in the E-step. This file contains the code to calculate conditional complete-data log likelihood in the E-step. 

`plloss` is the loss function for recurrent event data.

`rCoxPH` contains the fundemental functions for performing component-wise gradient boost for recurrent event data, including the loss function and negative gradient for shared frailty model for recurrent event data.
