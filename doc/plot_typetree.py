import os

from matplotlib import pyplot as plt

from gbtask.utils import TypeTree

fig, ax = plt.subplots()
TypeTree.plot(TypeTree.DefaultTypeTree, ax)
scriptdir = os.path.dirname(__file__)
fig.savefig(os.path.join(scriptdir, "typetree.png"))
