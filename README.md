# tstat-wildfire

Using smart thermostats to reduce indoor exposure to wildfire fine particulate matter (PM2.5)  
Cite: [https://doi.org/10.1016/j.indenv.2025.100088](https://doi.org/10.1016/j.indenv.2025.100088)  
---  

## 🧠 Summary

Exposure to fine particulate matter (PM2.5) is responsible for millions of premature deaths globally each year. Wildfires are a major source of PM2.5, creating dangerously high levels of air pollution across extensive regions.

Current public health recommendations for wildfire-related PM2.5 exposure include staying indoors and using portable air cleaners or central air systems with adequate filtration. However, central systems are often underutilized. In this study, we used smart thermostat data from ~5000 California homes during the 2020 wildfire peak to show that central systems are not effectively operated to improve indoor air quality.

We simulated the potential health benefits of optimizing central air system usage based on smart thermostat control and outdoor PM2.5 data. Key findings include:

- Automated HVAC control reduced indoor PM2.5 exposure by **up to 54 ± 5 %**, and **up to 61 ± 5 % during peak wildfire days**.
- The increased energy cost (~$5/month per home) is offset by an estimated **53 ± 5 % reduction in premature mortality**, corresponding to **$29M in monetized health benefits**.
- Combining **MERV 13 filters** and **low leakage** designs further improved protection.
- One central air system with proper filtering can be as effective as four portable air cleaners on average.
- The largest health impact gains were found in **lower-income areas**.

This work demonstrates how **existing, often overlooked infrastructure and technologies** can be better leveraged to protect public health during wildfire events.

---

## 📁 Repository Contents

This repository contains the R scripts required to:

- Preprocess housing, weather, and air quality data
- Run dynamic infiltration and indoor air PM2.5 simulations
- Model HVAC control scenarios and compare their impact
- Summarize results at hourly intervals across California census tracts

> See [`1_simulation.R`](./1_simulation.R) for the main simulation.
> See [`2_health.R`](./2_health.R) for the health impact assessment.
> See [`3_sensitivity_analysis.R`](./3_sensitivity_analysis.R) for the sensitivity analysis.

---

## 🛠 Requirements

This project uses the R ecosystem. To run the simulation, install the following packages:

```r
pacman::p_load(tidyverse, here, lubridate, zoo, pbapply, furrr, future.apply)
```

It is recommended to use an RStudio project with the [`here`](https://github.com/jennybc/here_here) package to manage paths.

---

## 📊 Data Access

The simulation requires input data including:

- Smart thermostat indoor temperature and runtime (Ecobee)
- Outdoor PM2.5 concentration (EPA)
- Weather data (ISD)
- Housing characteristics by census tract

> 🔒 These datasets are not included in this repository due to size and privacy constraints.  
> 📥 **To request access**, please contact: [federico.dallo@cnr.it](mailto:federico.dallo@cnr.it)

Data will be made available soon via a dedicated repository.

---

## 📜 Citation

If you use this code or its outputs, please cite the following publication:

> Dallo, F. et al. (2025). *Using smart thermostats to reduce indoor exposure to wildfire PM2.5*. *Indoor Environment*, [https://doi.org/10.1016/j.indenv.2025.100088](https://doi.org/10.1016/j.indenv.2025.100088)

---

## 👤 Contact

**Federico Dallo**  
Researcher, CNR  
📧 [federico.dallo@cnr.it](mailto:federico.dallo@cnr.it)

