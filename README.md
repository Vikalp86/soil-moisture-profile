### Soil Moisture Profile 

Virtical soil moisture profile development using statistical models:
  - The Exponential Model and 
  - The Principle of Maximum Entropy (POME) model 

The Exponentail model is based on the exponential decay function, mimicing the the time-lag relationship between the soil mositure at the surface with the bottom layers. The model depends upon a pre-defined relation (T) between the two layers and performs well under moist conditions (low evapotranspiration effects) and shallow depths. The model captures a temporal variations in moisture conidtions. The only input needed for the model is the surface moisture content value and a long-term (at least one year) record to compute 'T'. For details refers to [Wagner et al., (1999)](https://www.sciencedirect.com/science/article/abs/pii/S003442579900036X) and [Albergel et al., (2008)](https://d-nb.info/114976970X/34)

The POME model on the other hand develops a vertical soil moisture profile is based on maximising the entropy. The maximization of entropy characterizes the diffusion of moisture through the soil column as a function of time. Unlike exponential model, the POME does not depend upon the any a-priori information about the nature or kind of relation between the layers. The method guarantees the minimum variance unbiased profile subject to bopundary conditions. The model needs three inputs: the surface and bottom moisture content (boundary conditions) and the soil column mean soil moisture value. For more details please refer to Al-Hamdan and Cruise ([2010])(https://ascelibrary.org/doi/abs/10.1061/(ASCE)HE.1943-5584.0000196); Mishra et al., ([2015](https://www.mdpi.com/1099-4300/17/6/4454/htm);[2020](https://www.tandfonline.com/doi/full/10.1080/02626667.2020.1730846))

