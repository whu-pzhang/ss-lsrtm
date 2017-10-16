#!/usr/bin/env python2
# coding:utf-8

"""
Generate a spike file
according backgroud density/velocity model
"""

import numpy as np
import rsf.api as rsf

par = rsf.Par()

inp = rsf.Input()  # backgroud density/velocity model
spk = rsf.Input("spk")  # coordinates of spikes
out = rsf.Output() # Output spikes model

assert inp.type == 'float', "Need float input!"

n1 = inp.int("n1")
n2 = inp.size(1)

delta = par.float("delta") # perturbation of backgroud
assert delta, "Need delta="

den = np.zeros((n2, n1), dtype=np.float32)
inp.read(den)   # read in backgroud density file

spikes = np.zeros_like(den)

spk.read(spikes)

index = np.argwhere(spikes > 0)

# # 散射点处的背景场大小
# bak_den_in_spk = den[index[:,0], index[:,1]]

delta = 1 # 扰动百分比 100%

den[index[:,0], index[:,1]] *= (1+delta)

# plt.imshow(den.transpose())
# plt.show()

out.write(den)

