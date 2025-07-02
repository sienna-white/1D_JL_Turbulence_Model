import xarray as xr
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

class ModelResult:
    def __init__(self, fn, start_date, end_date, isave):
        self.ds = xr.open_dataset(fn, decode_timedelta=False)
        self.nt = len(self.ds.time.values)
        self.start_day = pd.to_datetime(start_date)
        self.end_day = pd.to_datetime(end_date)
        self.dates = pd.date_range(start_date, end_date, freq="10s")
        self.dates = self.dates[::isave]
        self.dates = self.dates[0:self.nt]
        print("There are %d saved profiles in this model..." % self.nt)

colors = ["#a4c639", "#73c2fb", "#d1001c"]

# Define some plotting functions
def format_date_ax(ax, days=2):
    """Format the x-axis of a plot with date labels."""
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=days))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b/%d %H:%M"))
    ax.grid(alpha=0.25, which="both")
    # plt.setp(ax.get_xticklabels(), rotation=10, ha='right')   

def add_colorbar(fig, c):
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    cb = fig.colorbar(c, cax=cbar_ax, shrink=0.4)
    return cb 

def filter_to_daytime(df, hr1, hr2): 
    df = df[df.index.hour >= hr1]
    df = df[df.index.hour <= hr2]
    return df 