from matplotlib import pyplot as plt
import seaborn as sns

from typing import List

sns.set_theme()

def display_windows(window_value:List[float], title:str=None, ylabel:str=None, ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    plt.plot(window_value)
    plt.xlabel("Window position")
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    plt.show(block=False)
    plt.waitforbuttonpress(timeout=2)
    return 0