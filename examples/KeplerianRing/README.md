Keplerian Ring Test
===================

This test uses the '2D Keplerian Ring' to assess the accuracy of SPH
calculations. For more information on the test, please see Cullen & Dehnen 2009
(http://arxiv.org/abs/1006.1524) Section 4.3.

It is key that for this test in particular that the initial conditions are
perfectly arranged (here we use the same methodology as described in the paper
referenced above).

The test uses:

+ GM = 1000 for a central point mass
+ r = 10 with sd = 2.5 gaussian distribution for particles
+ 9745 particles
+ c_s = 0.01 << v_\phi = 10.

Please note that the initial condition generator **requires python3 rather than
python2**.