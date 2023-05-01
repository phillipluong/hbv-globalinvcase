import atomica as at
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def plot_output(df, age, run_id, output):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Plots a time series of a specific output within an age group of a specific run_id among the stochastic outputs

    Inputs:
    - df:       Stochastic Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    - run_id:   Run_id of result to plot (int)
    - output:   Output column of result to plot - either 'cohort_size', 'cases', 'dalys', 'deaths' (str)
    '''
    plt.plot(np.arange(2000,2101), df[(df.age == age) & (df.run_id == run_id)][output])

def plot_cen_output(df, age, output):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Plots a time series of a specific output within an age group of among the central outputs

    Inputs:
    - df:       Central Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    - output:   Output column of result to plot - either 'cohort_size', 'cases', 'dalys', 'deaths' (str)
    '''
    plt.plot(np.arange(2000,2101), df[(df.age == age)][output],color = 'black',linewidth=2)

def plot_all(cen_df, sto_df, age):
    '''
    Created by: Phillip Luong (https://github.com/phillipluong)
    Last Updated: 31/01/23
    Plots time series of all outputs from the central run and all stochastic runs. The central run is specified by a thicker black line, while all stochastic outputs are characterised by thinner time series. Ideally up to 100 stochastic outputs will generate the best idea of the distribution from the uncertainty analysis.

    Inputs:
    - cen_df:   Central Result Dataframe (DataFrame)
    - sto_df:   Stochastic Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    '''
    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'cohort_size')
        plot_cen_output(cen_df, age, 'cohort_size')
        plt.title(f'Cohort Size, age {age}')

    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'cases')
        plot_cen_output(cen_df, age, 'cases')
        plt.title(f'Cases, age {age}')

    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'dalys')
        plot_cen_output(cen_df, age, 'dalys')
        plt.title(f'DALYs, age {age}')


    plt.figure()
    for i in range(1,31):
        plot_output(sto_df, age, i, 'deaths')
        plot_cen_output(cen_df, age, 'deaths')
        plt.title(f'Total Deaths, age {age}')

def plot_bounds(cen_df, sto_df, age, output):
    '''
    At the moment, it's strictly the 20th and 80th percentile, but we can change this to cater to what we need

    Inputs:
    - cen_df:   Central Result Dataframe (DataFrame)
    - sto_df:   Stochastic Result Dataframe (DataFrame)
    - age:      Age of result to plot (int)
    - output:   Output column of result to plot - either 'cohort_size', 'cases', 'dalys', 'deaths' (str)
    '''
    lwr = []
    hgr = []
    for year in range(2000,2101):
        lwr.append(sto_df[(sto_df.age == age) & (sto_df.year == year)][output].quantile(0.2))
        hgr.append(sto_df[(sto_df.age == age) & (sto_df.year == year)][output].quantile(0.8))

    fig, ax = plt.subplots()
    plt.plot(list(range(2000,2101)), lwr, label = '20 percentile')
    plt.plot(list(range(2000,2101)), hgr, label = '80 percentile')
    ax.fill_between(list(range(2000,2101)), lwr, hgr, alpha = 0.3)
    plt.plot(np.arange(2000,2101), cen_df[(cen_df.age == age)][output],color = 'black',linewidth=1.5)
