# Focused-Ultrasound-FUS-k-Wave

This folder contains the codes for Focused Ultrasound (FUS) in k-Wave MATLAB toolbox. A single element bowel-shaped focused transducer (f-number = 0.9, diameter = 50 mm, and focal length = 40 mm) working at 0.5 MHz and 0.9 MHz was used. Four simulations were done to study the pressure distribution in single layer brain medium (SLB) and four layer brain medium (FLB). Moreover, the focal pressure was kept at 0.5 MPa in all our studies and the computational grid size was taken as 124 x 124 x 164 with ∆x, ∆y, and ∆z all equal to 0.5 mm and width of the perfectly matched layer (PML) was twenty grid points. Inside the computational domain as shown in Figure 1, FLB medium is modelled as concentric semicircles and the background medium is water. The distance between the transducer and the brain is 28 mm while the centre of the FLB medium is 13 mm away from the midpoint of the computational domain. The ultrasound focus was positioned at 10 mm from the CSF layer and the background temperature was about 22 °C. Furthermore, the exposure time for the ultrasound was 2379 steps with a temporal resolution (∆t) of 18.5 ns and for the stability the Courant–Friedrichs–Lewy (CFL) number (va x ∆t/∆x) was 0.1 where va is the speed of US in the skull medium [1].


<img width="699" alt="Screenshot 2022-04-20 at 7 18 07 PM" src="https://user-images.githubusercontent.com/71398563/174534862-9ffe32ce-cc09-4e91-ad8b-dfc729d0353d.png">






[1] http://urn.fi/URN:NBN:fi:oulu-202206152859
