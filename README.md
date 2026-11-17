Repository that accompanies:

EVALUATION AND ADVICE ON ALTERNATIVE HATCHERY MANAGEMENT SCENARIOS TO INFORM RECOVERY FEASIBILITY ASSESSMENT OF SAKINAW SOCKEYE

CSAS Science Response

Data and code to recreate the analyses and explore further simulations are provided in this repository. We will give an overview of each component and instructions.

data - 2 datasets\
sakinaw_abundances.csv\
-column headers should be self-explanatory, 'surv' means survival and are expressed as percentages (not proportions)\
-note that marine survivals are indexed by the out-migrating smolt cohort year, referring to adult returns from this cohort +2 years\

smolt_broodtable.csv\
-same data as in 'sakinaw_abundances', just with natural-origin smolts indexed by brood year (ie. NO smolts counted +2 years from the returning brood cohort year)\

R - 2 scripts\
funcs.R\
-functions used in the main code, including functions for simulating the Sakinaw captive brood program ('cbp_sim') and the ex situ Big Qualicum River broodtake scenarios ('bq_sim')\
sakinaw_pva.R\
-code that processes the empirical data and runs the simulations for the 8 scenarios presented in the CSAS document, as well as all figures\


Customized simulations for captive brood program or ex situ broodtake program scenario planning\
-the 'cbp_sim' function has 2 main input parameters: 'mean.fry' (the mean number of fry released from the program) and 'mean.smolt' (the mean number of smolts released from the program), with options to also specify a specific standard deviation for each parameter - if SDs are unspecified it defaults to the empirical standard deviations for both parameters. Additional options for 'n.sims' (number of simulations, default 1000) and 'sim.length' (number of years to forecast, default = 24)\
-the 'bq_sim' function has 1 main input parameter: 'bq.release.target', the target number of smolts released to the ex situ site (Big Qualicum River)\
-outputs from either function can be quickly plotted, using the same figures as the CSAS document (2 and 3), using the function 'plot_sims', which only requires the output from either 'cbp_sim' or 'bq_sim'\
