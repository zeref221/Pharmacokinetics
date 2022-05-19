# Pharmacokinetics

*Objectives* 
1. Pharmacokinetic simulations of a simple 2-compartment model with oral or I.V. administration 
2. Manually specify the parameters estimates (theta's, omega's and sigma's)
3. Generate a figure of the simulation 

Results:

1) _PK model_

![Output ](https://user-images.githubusercontent.com/70428805/168469143-539a283c-f30d-4492-bed7-3884fae3d3bd.jpeg)

>Fig: Time profile of a drug over several hours vs Concentratation of blood plasma 

2) _For monte carlo integration generate the equation of the above curve using trendline feature in excel_

![Rplot02](https://user-images.githubusercontent.com/70428805/169224552-8c301211-7986-4a42-af2e-f95240b17139.png)

> PLotting the curve in excel

![WhatsApp Image 2022-05-19 at 12 07 05 PM](https://user-images.githubusercontent.com/70428805/169226451-2d0e86e3-a482-4164-aa3a-5cfd9702e418.jpeg)

>  Using trendline to generate the polynomial equation for our curve

3) _Monte carlo intergration using R (sample size= 2000)_

![Rplot01](https://user-images.githubusercontent.com/70428805/169224794-d4343a30-725e-4f49-99c0-ef6031a8b5a6.png)

> Area of the shaded region is the value of our intergral

4) _Checking accuracy of our estimate by varying the sample size (True value: 93.47426)_

![Rplot](https://user-images.githubusercontent.com/70428805/169225075-eabde576-a02a-43df-a0f7-1f78aeafe11b.png)

> For larger value of sample size the value of our simulation tends to be close to the true value (93.47426)

Readings:
https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4394611/
 
http://webpopix.org/shiny/ShinyExamples.html

Further work
1. Make an app via shiny 
2. Have Real Time simulation of the model



