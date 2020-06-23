import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats
import matplotlib.pylab as pylab
import ipywidgets as widgets
import sys
from os import path
from sklearn.decomposition import PCA
from matplotlib.offsetbox import AnchoredText
from IPython.display import clear_output
from IPython.display import display


class DataLoader:
    def __init__(self, trtmnt, normalized):
        """
        input:
        @trtmnt: a string, must be one of the four: high, low, veh, pt
        @normalized: bool, True or False
        """
        self.__norm = normalized
        self.__kall = None
        self.__rsem = {}
        self.__dir = {'high': ['/home/lima/Project/simulation/Comparision/high',
                               '/home/lima/Project/simulation/Comparision/high/Normalized'],
                      'low': ['/home/lima/Project/simulation/Comparision/low',
                             '/home/lima/Project/simulation/Comparision/low/Normalized'],
                       'veh': ['/home/lima/Project/simulation/Comparision/veh',
                              '/home/lima/Project/simulation/Comparision/veh/Normalized'],
                       'pt': ['/home/lima/Project/simulation/Comparision/pt',
                             '/home/lima/Project/simulation/Comparision/pt/Normalized']
                      }
        if normalized == True:
            if trtmnt == 'pt':
                print("Tumor has only one sample, R cannot perform normalization")
                sys.exit()
            self.__path = self.__dir[trtmnt][1]
            self.__kallPrefix = 'NormalizedKallistoReadCounts'
            self.__rsemPrefix = 'NormalizedCountsCutOffValue'
        else:
            self.__path = self.__dir[trtmnt][0]
            self.__kallPrefix = 'KallistoReadCounts'
            self.__rsemPrefix = 'CountsCutOffValue'

    def __loadKall(self):
        kallDF = pd.read_csv("{path}/{kallPrefix}.csv".format(path = self.__path, kallPrefix = self.__kallPrefix))
        #kallDF = pd.read_csv('/home/lima/Project/simulation/Comparision/high/KallistoReadCounts.csv')
        kallDF = kallDF.rename(columns = {"Unnamed: 0": "GeneID"})
        kallDF.set_index("GeneID", inplace = True)
        # print(kallDF)
        if self.__norm == False:
            kallDF = np.log2(kallDF + 1)
        self.__kall = kallDF

    def __loadRSEM(self):
        cutOff = [0, 2.5, 5, 7.5, 10]
        # key of the tempDict is the cut off threshould, value is the cnt matrix accordingly
        tempDict = {}
        dictList = []
        for cut in cutOff:
            fullPath = "{path}/{rsemPrefix}{cut}.csv".format(cut = cut, path = self.__path, rsemPrefix = self.__rsemPrefix)
            # fullPath = "/home/lima/Project/simulation/Comparision/high/CountsCutOffValue{cut}.csv".format(cut = cut)
            dataFrame = pd.read_csv(fullPath)
            dataFrame = dataFrame.rename(columns = {"Unnamed: 0": "GeneID"})
            dataFrame.set_index("GeneID", inplace = True)
            if self.__norm == False:
                dataFrame = np.log2(dataFrame + 1)
            tempDict[cut] = dataFrame
        self.__rsem = tempDict

    def getKall(self):
        self.__loadKall()
        return self.__kall

    def getRSEM(self):
        self.__loadRSEM()
        return self.__rsem

class ConcateData:
    def __init__(self, rsem, kallisto):
        # input:
        # @groupList: a list contains the names of samples. there are four groups of samples, high, low, veh, pt
        # @rsem: rsem counts
        # @kallisto: kallisto counts
        # both counts are in gene level
        # self.__list = groupList
        self.__rsem = rsem
        self.__kallisto = kallisto
        self.__result = None

    def __concat(self):

        #for sample in self.__list:

        rsemLen = len(self.__rsem)
        kallLen = len(self.__kallisto)
        name = list(self.__rsem.columns.values)
        sample = name[0]

        print("{num} genes have been detected in {sample} by RSEM".format(num = rsemLen, sample = sample))
        print("{num} genes have been detected in {sample} by kallisto".format(num = kallLen, sample = sample))

        # convert the index, gene ids, to the first column for merge two dataframes
        self.__rsem.reset_index(inplace = True)
        self.__kallisto.reset_index(inplace = True)
        # print(targetKall)
        # print(targetRSEM)
        df = pd.merge(self.__rsem, self.__kallisto, how = 'outer', on = 'GeneID', suffixes = ('_RSEM', '_Kallisto'))
        print("{num} genes in both dataframes".format(num = len(df)))
        # take out the GeneIDs in each data frame and put them into two sets
        rsemIdx = set(self.__rsem['GeneID'])
        kallIdx = set(self.__kallisto['GeneID'])

        inBoth = rsemIdx & kallIdx
        rsemOnly = rsemIdx - inBoth
        kallOnly = kallIdx - inBoth

        inBoth = list(inBoth)
        rsemOnly = list(rsemOnly)
        kallOnly = list(kallOnly)

        bothDF = pd.DataFrame(list(zip(inBoth, ["DetectedByBoth"] * len(inBoth))), columns = ["GeneID", "Labels"])
        rsemOnlyDF = pd.DataFrame(list(zip(rsemOnly, ["OnlyDetectedByRSEM"] * len(rsemOnly))), columns = ["GeneID", "Labels"])
        kallOnlyDF = pd.DataFrame(list(zip(kallOnly, ["OnlyDetectedByKallisto"] * len(kallOnly))), columns = ["GeneID", "Labels"])
        labelDF = pd.concat([bothDF, rsemOnlyDF, kallOnlyDF])
        print("{num} labels have been attached".format(num = len(labelDF)))
        # fill nans in df with 0
        df.fillna(0, inplace = True)
        df = pd.merge(df, labelDF, how = 'inner', on = 'GeneID')
        # print(df)
        self.__result = df

    def setData(self):
        self.__concat()

    def getData(self):
        return self.__result

class Plotter:
    def __init__(self, cutOff):
        self.__cutOff = cutOff

    def plotPlotBlock(self, dfs, sample):

        params = {'legend.fontsize': 8,
                    'font.family': 'serif',
                    'font.weight': 'medium',
                    'font.variant': 'normal',
                    'figure.figsize': (10, 10),
                    'figure.edgecolor': '#04253a',
                    'figure.titlesize': 8,
                    'axes.labelsize': 6,
                    'axes.titlesize': 6,
                    'xtick.labelsize': 6,
                    'ytick.labelsize': 6}
        pylab.rcParams.update(params)
        figureOne, axes = plt.subplots(2, 2, constrained_layout = True)
        # self.buildPlotBlock(df, axes, cutOff, names)
        self.buildPlotBlock(dfs, axes)
        plt.suptitle('Statistical plots for comparing counts for {sample} based on RSEM & Kallisto with reads shorter than {cut} being removed'.format(sample = sample, cut = self.__cutOff))
        # title = 'Real data and synthetic data of {realData} with reads shorter than {cut} being removed.png'.format(realData = names[0], cut = cutOff)
        # plt.savefig(title, dpi = 300, bbox_inches = 'tight')
        plt.plot()

    def buildPlotBlock(self, dfs, axes):
        """
        buildPlotBlock is a helper function helping putting the desired plots, such as the scatter plot, box plot,
        violin plot and the distribution plot into a grided pyplot

        Input:
        @df: the dataframe which provides the data
        @axes: a 2x2 subplot matrix
        """
        # df = np.log2(df[names] + 1)
        rsem = dfs[0].copy()
        # print(rsem)
        kall = dfs[1].copy()
        # print(kall)
        rsem['Labels'] = ['RSEM'] * len(rsem)

        kall['Labels'] = ['Kallisto'] * len(kall)

        df = pd.concat([rsem, kall])

        # df = pd.DataFrame([df.iloc[ : , 0], df.iloc[ : , 1]] , columns = ['Counts', 'Labels'])
        colNames = list(df.columns.values)
        sns.set(style = "whitegrid", palette = "muted", color_codes = True)
        # put the violin plot to axes[0, 0]
        axes[0, 0].clear()
        sns.violinplot(x = colNames[1], y = colNames[0], data = df, linewidth = 0.2, ax = axes[0, 0])

        # put the boxplot to axes[1, 1]
        axes[1, 1].clear()
        sns.boxplot(x = colNames[1], y = colNames[0], data = df, linewidth = 0.2, ax = axes[1, 1])

        # put the distribution plot to axes[0, 1]
        axes[0, 1].clear()
        sns.distplot(rsem.iloc[ : , 0], kde = True, color = "b", label = 'RSEM', ax = axes[0, 1])
        sns.distplot(kall.iloc[ : , 0], kde = True, color = "r", label = 'Kallisto', ax = axes[0, 1])
        axes[0, 1].legend(loc = 'best', prop = {'size': 8})
        axes[0, 1].set_xlabel("Distribution plot for real data and simulated data")

        # put the scatter plot to axes [1, 0]
        axes[1, 0].clear()
        self.scatterPlot(dfs, axes[1, 0])

    def scatterPlot(self, dfs, ax):
        """
        the scatterPlot function plots a pairwise scatter plot for the simulated data and the real data for the same
        sample, with the x-coordinates being the counts number of the real data and the y-coordinates being the counts
        number of the simulated data

        input data:
        @df: a pandas dataframe contains the counts for the real and the simulated data
        @ax: an axes object
        @names: a list contains the sample names
        """

        # notInrealData, notInSimulatedData, inNeither are lists stores the row indecies for data as theirs names
        # notInRealData = df.index[df.iloc[ : , 0] == 0].tolist()
        # notInSimulatedData = df.index[df.iloc[ : , 1] == 0].tolist()
        # inNeither = list(set(notInRealData) & set(notInSimulatedData))

        # inBoth, simOnly, realOnly, neither are integers for
        # genes found in both samples, one of the samples or none of them
        # neither = len(inNeither)
        # simOnly = len(notInRealData) - neither
        # realOnly = len(notInSimulatedData) - neither
        # inBoth = len(df) - simOnly - realOnly + neither

        # category the genes into 4 classes
        # df["ExpressionStatus"] = ["BothSamples"] * len(df)
        # df.loc[notInRealData, "ExpressionStatus"] = "SimulatedOnly"
        # df.loc[notInSimulatedData, "ExpressionStatus"] = "RealOnly"
        # df.loc[inNeither, "ExpressionStatus"] = "NeithSample"
        # print(df)
        # set up the plotting
        # print(dfs[0])
        # print(dfs[1])
        conc = ConcateData(dfs[0], dfs[1])
        conc.setData()
        df = conc.getData()
        # print(df)
        sns.set(style = "whitegrid", palette = "muted", color_codes = True)
        # sns.set(style = "dark", palette = "inferno", color_codes = True)
        # calculate the correlation between the two sets of data.
        # print(df.iloc[ : , 0])
        names = list(df.columns)
        # print(names)
        # col_one = df.iloc[ : , 1]
        # col_two = df.iloc[ : , 2]
        # corr = scipy.stats.pearsonr(col_one, col_two)[0]

        detectedByBoth = len(df[df['Labels'] == 'DetectedByBoth'])
        detectedByRSEMOnly = len(df[df['Labels'] == 'OnlyDetectedByRSEM'])
        detectedByKallistoOnly = len(df[df['Labels'] == 'OnlyDetectedByKallisto'])

        sns.scatterplot(x = names[1], y = names[2], hue = 'Labels', marker = '+', style = 'Labels', data = df, ax = ax, palette = 'inferno', s = 2.0)
        at = AnchoredText("{inBoth} were found by both RSEM & Kallisto\n{RSEM} were detected by RSEM only\n{Kallisto} were detected by Kallisto only\nCut off threshould: {cut}".format(
            inBoth = detectedByBoth, RSEM = detectedByRSEMOnly, Kallisto = detectedByKallistoOnly, cut = self.__cutOff),
                          prop = dict(size = 8), frameon = True, loc = 'upper left')
        at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")

        ax.add_artist(at)
        ax.set_xlabel("RSEM")
        ax.set_ylabel("Kallisto")
        ax.legend(loc = 'lower right', fancybox = True, shadow = False, scatterpoints = 200, prop = {'size': 8})

class Pen:
    __trtmnt = {'high': pd.read_csv('/home/lima/Project/simulation/Comparision/high.txt', skiprows = 0).values.tolist(),
                     'low': pd.read_csv('/home/lima/Project/simulation/Comparision/low.txt', skiprows = 0).values.tolist(),
                     'veh': pd.read_csv('/home/lima/Project/simulation/Comparision/veh.txt', skiprows = 0).values.tolist(),
                     'pt': pd.read_csv('/home/lima/Project/simulation/Comparision/pt.txt', skiprows = 0).values.tolist()}
    def __init__(self, trtmnt, normalized):
        self.__loader = DataLoader(trtmnt, normalized)
        self.__kall = self.__loader.getKall()
        self.__rsem = self.__loader.getRSEM()
        self.__cut = [0, 2.5, 5, 7.5, 10]
        self.__samples = self.__trtmnt[trtmnt]

    def draw(self):
        for sample in self.__samples:
            for cut in self.__cut:
                rsemCut = rsem[cut]
                plotter = Plotter(cut)
                plotter.plotPlotBlock([rsemCut[sample], self.__kall[sample]], sample)
