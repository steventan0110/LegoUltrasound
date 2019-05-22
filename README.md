# LegoUltrasound by Weiting Tan 
**from 2019/05/21
Summer Research at Hopkins under Professor Emad Boctor**

The proeject starts with Matlab code, completing simulation based on different phantoms (using K-wave software).
With data retrieved from the simulation, use python to process and analyze with neural network
Try to decrease the # of elements used in the transducer to increase the flexibility of current probes => "LEGO" like ultrasound probe

Process:
1. write bone generator for lumbar puncture based the position of the needle
2. process the phantom with K-wave simulation, get the RF data and write beamforming code to process the data and transform it into the data that reflects the intensity at different intersection of the grid
3. process the intensity value with supervised learning to train the neural network and the test data would be beamformed data from real ultrasound devices in the lab
