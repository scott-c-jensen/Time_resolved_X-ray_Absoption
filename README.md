# Time_resolved_X-ray_Absoption

## Experimental Setup in Brief
I was looking at the changes in Photosystem II in spinach to understand the Kok cycle. **Basically what do plants do with the light and store it's energy during photosynthesis.** We examined this using x-ray absorption of the Manganese atoms which are central to this process.

PSII has 4 semi-stable states (S0-S3) which it goes through when illuminated by light. Different changes occur during each light-driven change between states. One of the most interesting transitions is the S3-S0. We looked at this at the microsecond timescale and examined which of the S3-S0 models would best represent the data.

## Software Description
This code base takes raw dectector data, recorded at 1MHz, and separates the data based on where laser flashed the sample which induced changes. The changes are extracted, background subtracted and combined. 

The data are then fit in a way that may seem complicated because, frankly, the system is difficult to study. Basically 4 states exist. Changing states changes the signal (x-ray absorption). But this is convoluted further with imperfect advancement (see below for an large aside). However, the advancement is unknown so it is globally fit along with the kinetics (time and amplitude of changes) as they are the same in all transitions.

## Dependencies
Python 2.7.11
Numpy 1.11
LMFIT 0.98

### More than you want to know on the state transitons
We are looking at transitions between states, and while we know that we have a nearly pure starting state, there is incomplete advancement and uncertain kinetics to deal with. So each laser flash (of 5 total) creates a more complicated mixture of initial and final states. A fairly good approximation is to assume a single advancement percentage.

There are different oxidation states (S-states) in photosystem II. They cycle from S0-S4. The S4-S0 transition doesn't require light to transition, returning to the most reduced state. Not only are there multiple changes in the system at each state, but only the S1 (the dark adapted state) is a pure state, others are mixed. This program therefore, has global fitting of all s-state transitions for the advancement and individual state transitions with the percent of each state.  

So for example:  
If we have a 80% advancement then after each flash (F) we have  
F0 = 100% S1  
F1 = 20% S1, 80% S2  

But what we are looking at is the change in absorption... so we only care about the changes. So if we make a matrix for the S-state transitions for each flash, such as the first transition S1 to S2 (or S1-S2). This can be written for each flash as an array [S0-S1,S1-S2,S2-S3,S3-S0]:  

1F = [0.00, 0.80, 0.00, 0.00,]  
2F = [0.00, 0.16, 0.64, 0.00,]  
3F = [0.00, 0.03, 0.26, 0.51,]  
4F = [0.41, 0.01, 0.08, 0.31,]  
5F = [0.33, 0.33, 0.02, 0.12,]  


## Experimental Setup in Brief
I was looking at the changes in Photosystem II in spinach to understand the Kok cycle. **Basically what do plants do with the light and store it's energy during photosynthesis.** We examined this using x-ray absorption of the Manganese atoms which are central to this process.

PSII has 4 semi-stable states (S0-S3) which it goes through when illuminated by light. Different changes occur during each light-driven change between states. One of the most interesting transitions is the S3-S0. We looked at this at the microsecond timescale and examined which of the S3-S0 models would best represent the data.

Through a complicated time sequence, the data was captured

