from matplotlib import pyplot as plt
import seaborn as sns

sns.set_theme()

from typing import List


def display_appearance_freq(dico_windows : dict, ax : plt.Axes = None):
    plt.plot(list(range(0, len(dico_windows))), list(dico_windows.values()))
    plt.xlabel("Key index")
    plt.ylabel("kmer count in the genome")
    return 0

def display_freq(freq_km_avg : List[float], ax : plt.Axes = None):
    plt.plot(list(range(0, len(freq_km_avg))), freq_km_avg)
    plt.xlabel("Window position")
    plt.ylabel("Average frequency")
    plt.title("Average kmer frequency in the window over the genome length")
    plt.show(block=False)
    plt.waitforbuttonpress(timeout=2)
    return 0