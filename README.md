# Eye movement analysis

## Introduction
This project investigates the relationship between memory for familiar scenes and eye movement behavior, with a specific focus on the persistence of the image familiarity effect over longer time scales. The image familiarity effect refers to the phenomenon where individuals exhibit longer fixation duration and shorter saccade amplitude when viewing familiar images. The main objective of this study is to determine if this effect persists over longer time scales, such as weeks or months, and to clarify if the effect emerges due to visual long-term memory or other mechanisms. The project includes an eye-tracking experiment where participants observe images while their eye movements are recorded, and the data collected is analyzed using statistical methods. The findings of this study have potential implications for fields such as visual perception and cognitive psychology.

## Tools

The project was implemented using Matlab R and the following libraries:

* `Microsaccade Toolbox`
* `spatstat`
* `dplyr` for data cleaning and manipulation
* `tidyr` for data reshaping
* `ggplot2` for data visualization
* `lme4` for mixed-effects modeling
* `car` for model diagnostic and post-hoc analysis

## Data Collection

The eye-tracking data was collected using an Eyelink 1000 desktop mounted eye tracker. Participants observed new and familiar images in two experimental sessions, one week apart from each other. The experiment was programmed in MATLAB and consisted of two blocks: memorization and recall.

## Data Preprocessing

In order to prepare the data for further analysis, we first conducted saccade detection using the velocity-based algorithm proposed by Engbert & Kliegl (2003), which utilizes the saccades’ specific property - linear relationship between peak velocity and amplitude. In order to conduct the saccade detection, we used the Microsaccade Toolbox for R developed by Engbert et al.(2015). 
Another step of the data preprocessing was elimination of outliers for dependent variables. As the result, we set limit for minimal and maximal values for dependent variables.

## Analysis

The data collected was analyzed using a combination of statistical methods, including pair correlation function (PCF) analysis and linear mixed-effects modeling (LMM).

### Pair Correlation Function

In order to investigate the influence of image familiarity effect on the eye movement behavior, the pair correlation function analysis was applied to the eye movement data. The pair correlation function (PCF) is a method from spatial statistics and can be used to study and compare spatial correlations between fixations under different viewing conditions. The PCF analysis for fixation locations was implemented in R using the `spatstat` package.

First, we simulated two control point processes (homogeneous and inhomogeneous) to ensure that correlations in the PCF arise from the empirical data. To do so, we estimated fixation density for each image using the data from all participants. Once the fixation density was estimated, individual scanpaths for the inhomogeneous and homogeneous point processes were simulated using the computed density.

IMAGE HERE
    
## Acknowledgements

I would like to extend my sincere appreciation to my supervisors, Prof. Dr. Ralf Engbert, Lisa Schwetlick, and Daniel Backhaus for their guidance, mentorship, and support throughout this project. Their valuable insights, expertise, and constructive feedback have been instrumental in shaping the direction and outcome of this study. Their unwavering encouragement and support have been a source of inspiration and motivation. I am deeply grateful for the opportunity to work under their guidance and for the invaluable experience gained from this project. I would like to thank them for their dedication, time, and patience throughout the project.
