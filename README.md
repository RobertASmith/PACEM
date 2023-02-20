[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![GitHub last commit](https://img.shields.io/github/last-commit/RobertASmith/PACEM?color=red&style=plastic)
![GitHub top language](https://img.shields.io/github/last-commit/RobertASmith/PACEM?style=plastic)
![GitHub repo size](https://img.shields.io/github/last-commit/RobertASmith/PACEM?style=plastic)
![GitHub forks](https://img.shields.io/github/forks/RobertASmith/PACEM?style=social&label=Fork&maxAge=2592000)

# PACEM

Physical Activity in Children Economic Model repository used to store model code and input data.

Please cite this work as:

Smith, Robert (2023) *Improving the validity and usability of decision models: case studies with a focus on physical activity*. PhD thesis, University of Sheffield.

### Installation

1. Install rStudio
2. Clone this repository `git clone https://github.com/RobertASmith/PACEM.git`
3. Open `PACEM.Rproj` in rStudio
4. Run `install.packages("renv")` in the `R` terminal
5. Run `renv::restore()` in the `R` terminal

### Running the code

The code can only be run with access to the underlying model data. This can be obtained upon request.

Once you have the data, save it in 'Data' and then run the following scripts:

1. `mirosim_MAIN.R`
 Purpose: To run the model (as long as you have access to data).
 Outputs: a set of 9 .rds files saved in 'Results' folder containing model and sensitivity analysis results.
2. `scripts\create_visuals.R`
 Purpose: Create visualisations for reporting   
 Outputs: a set of figures and tables saved in 'Figures' folder.
   


