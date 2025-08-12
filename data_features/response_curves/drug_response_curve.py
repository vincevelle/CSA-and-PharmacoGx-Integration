#! /usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager

import sys
import argparse
import numpy as np
import pandas as pd
import re
import _pickle as cp
import os
import copy
import itertools

# from tqdm import tqdm
# from itertools import islice
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit
from scipy import stats



def z_prime_score(pos, neg):
    z_prime = 1 - (3 * np.std(pos) + 3 * np.std(neg)) / np.abs(np.mean(pos) - np.mean(neg))
    return z_prime



def z_prime_robust_score(pos, neg):
    z_prime = 1 - (3 * np.median(np.abs(pos - np.median(pos))) + 3 * np.median(np.abs(neg - np.median(neg)))) / np.abs(np.median(pos) - np.median(neg))
    return z_prime
# z_robust_factor <- 1 - (3 * (mad_stauro + mad_dmso) / abs(median_stauro - median_dmso))


def IQR_outlier_detection(x, lower_quantile=0.25, upper_quantile=0.75, factor=1.5):

    # x should be a 1-D data array whose outliers need to be detected and removed.
    Q1 = np.quantile(x, lower_quantile)
    Q3 = np.quantile(x, upper_quantile)
    IQR = Q3 - Q1
    lower = Q1 - factor * IQR
    upper = Q3 + factor * IQR
    outlier_id = np.sort(np.union1d(np.where(x < lower)[0], np.where(x > upper)[0]))

    return outlier_id, lower, upper



class drug_screening_dataframe:

    def __init__(self, res, control_types, model_col, plate_col, row_col, column_col, comp_col, dose_col, solvent_col,
                 fea_col, across_plate_col, study_id, row_edge_index=None, col_edge_index=None, edge_col=None,
                 exp_col=None, norm_fea_col=None, control_col=None):
        self.res = res                          # Drug screening dataframe
        self.control_types = control_types      # A list of all negative and positive controls in the compound name column, not including control compounds
        self.model_col = model_col              # A string, the name of the column giving the models
        self.plate_col = plate_col              # A string, the name of the column giving the plates
        self.row_col = row_col                  # A string, the name of the column giving the well row indices
        self.column_col = column_col            # A string, the name of the column giving the well column indices
        self.comp_col = comp_col                # A string, the name of the column giving the compound names
        self.dose_col = dose_col                # A string, the name of the column giving the dose of treatment. Concentration must be in log10(M)
        self.solvent_col = solvent_col          # A string, the name of the column providing solvent to be used for normalization
        self.fea_col = fea_col                  # A string, the name of the column of feature values based on which calculation should be performed.
        self.across_plate_col = across_plate_col# A string, the name of the flag column indicating whether a well should be considered within or across plates.
        self.study_id = study_id                # A string, indicating the name of the drug screening study
        self.row_edge_index = row_edge_index    # A dictionary giving the row edge index in every plate.
        self.col_edge_index = col_edge_index    # A dictionary giving the column edge index in every plate.
        self.edge_col = edge_col                # A string, the name of the column indicating whether a well is on the edge of a plate
        self.exp_col = exp_col                  # A string, the name of the column providing unique experiment IDs
        self.norm_fea_col = norm_fea_col        # A string, the name of the column providing normalized feature values
        self.control_col = control_col          # A string, the name of the column indicating whether a well is control or not



    def add_edge_flag_column(self, row_edge_index, col_edge_index, edge_col='edge_flag'):
        self.row_edge_index = row_edge_index
        self.col_edge_index = col_edge_index
        edge = []
        for i in range(self.res.shape[0]):
            if self.res[self.row_col].iloc[i] in self.row_edge_index[self.res[self.plate_col].iloc[i]] \
                    or self.res[self.column_col].iloc[i] in self.col_edge_index[self.res[self.plate_col].iloc[i]]:
                edge.append(True)
            else:
                edge.append(False)
        self.res[edge_col] = edge
        self.edge_col = edge_col



    def add_experiment_id_column(self, exp_start_id, exp_col='exp.id'):
        exp_id = pd.Series(['' for i in range(self.res.shape[0])])
        for i in range(self.res.shape[0]):
            if self.res.iloc[i, :][self.comp_col] in self.control_types:
                continue
            if exp_id[i] == '':
                plate = self.res.iloc[i, :][self.plate_col]
                compound = self.res.iloc[i, :][self.comp_col]
                cell = self.res.iloc[i, :][self.model_col]
                if self.res.iloc[i, :][self.across_plate_col]:
                    # For compounds, all wells of which across plates are taken as one experiment
                    exp_index = np.intersect1d(np.intersect1d(np.where(self.res[self.comp_col] == compound)[0],
                                               np.where(self.res[self.model_col] == cell)[0]),
                                               np.where(self.res[self.across_plate_col])[0])
                else:
                    # For compounds, all wells of which on a single plate form an experiment.
                    exp_index = np.intersect1d(np.intersect1d(np.where(self.res[self.plate_col] == plate)[0],
                                                              np.where(self.res[self.comp_col] == compound)[0]),
                                               np.intersect1d(np.where(self.res[self.model_col] == cell)[0],
                                                              np.where(np.invert(self.res[self.across_plate_col]))[0]))
                exp_start_id = exp_start_id + 1
                exp_id[exp_index] = str(exp_start_id)
            else:
                continue
        self.res[exp_col] = exp_id.values
        self.exp_col = exp_col



    def edge_outlier_detection(self, threshold=0.3):
        '''
        :param threshold: threshold used for outlier definition
        :return:
            The row indices of edge outliers
        '''

        flag = [False for i in range(self.res.shape[0])]
        for i in range(self.res.shape[0]):
            if self.res.iloc[i, :][self.edge_col]:
                id_i = np.intersect1d(np.intersect1d(np.where(self.res[self.comp_col] == self.res.iloc[i, :][self.comp_col])[0],
                                                     np.where(self.res[self.dose_col] == self.res.iloc[i, :][self.dose_col])[0]),
                                      np.intersect1d(np.where(self.res[self.plate_col] == self.res.iloc[i, :][self.plate_col])[0],
                                                     np.where(self.res[self.model_col] == self.res.iloc[i, :][self.model_col])[0]))
                id_i = id_i[np.where(np.invert(self.res.iloc[id_i, :][self.edge_col]))[0]]
                if len(id_i) >= 1:
                    non_edge_mean = np.mean(self.res.iloc[id_i, :][self.fea_col].values)
                    if np.abs((self.res.iloc[i, :][self.fea_col] - non_edge_mean) / non_edge_mean) > threshold:
                        flag[i] = True

        return np.where(flag)[0]



    def IQR_outlier_detection_within_plate(self, comp, dose=None, factor=1.5):
        '''
        :param comp: a string or a list of strings of compounds for which outliers should be detected.
                    These compounds should be included in the self.comp_col column of the dataframe
        :param dose: a number or a list of numbers giving the doses within which outliers should be detected.
                    The length of dose should be the same as the length of comp, so that each compound is evaluated
                    at a specific dose.
        :param factor: the factor used in IQR outlier detection
        :return:
            outlier_id: a list of row indices of outliers
            outlier_per_plate: a nested dictionary of outliers in each plate. In the first tier, keys are plate IDs.
                In the second tier, keys are compound names. In the third tier, keys are 'outlier_id', 'lower', 'upper',
                and 'values'. 'values' includes all values of the feature of the compound in the plate. 'lower' is
                the lower bound. 'upper' is the upper bound. 'outlier_id' is a list of row indices of outliers for
                the compound in the plate.
        '''

        if isinstance(comp, str):
            comp = [comp]

        idna = np.where(pd.isna(self.res[self.fea_col]))[0]

        outlier_per_plate = {}

        plates = np.unique(self.res[self.plate_col].iloc[np.where(np.invert(pd.isna(self.res[self.plate_col])))[0]])
        for p in plates:
            idp = np.where(self.res[self.plate_col] == p)[0]

            for i_c in range(len(comp)):
                c = comp[i_c]
                idc = np.where(self.res[self.comp_col] == c)[0]
                if dose is not None:
                    d = dose[i_c]
                    idc = np.intersect1d(idc, np.where(self.res[self.dose_col] == d)[0])
                    c = c + '___' + str(d)
                idpc = np.setdiff1d(np.intersect1d(idp, idc), idna)

                if len(idpc) > 0:
                    if p not in outlier_per_plate.keys():
                        outlier_per_plate[p] = {}
                    outlier_per_plate[p][c] = {}
                    outlier_per_plate[p][c]['values'] = self.res.iloc[idpc, :][self.fea_col].values
                    outlier_idpc, lower_idpc, upper_idpc = IQR_outlier_detection(x=outlier_per_plate[p][c]['values'],
                                                                                 factor=factor)
                    outlier_per_plate[p][c]['outlier_id'] = idpc[outlier_idpc]
                    outlier_per_plate[p][c]['lower'] = lower_idpc
                    outlier_per_plate[p][c]['upper'] = upper_idpc

        outlier_id = []
        for p in outlier_per_plate.keys():
            for c in outlier_per_plate[p].keys():
                outlier_id = np.concatenate((outlier_id, outlier_per_plate[p][c]['outlier_id']))
        outlier_id = np.sort([int(i) for i in outlier_id])

        return outlier_id, outlier_per_plate



    def calculate_z_prime_per_plate(self, control, dose=None, ttest_flag=False, plot_file=None, robust_plot_file=None,
                                    plot_size=(10, 6)):
        '''
        :param control: a list of two strings giving the names of positive and negative controls. These control names
                        should be included in the self.comp_col column of the dataframe
        :param dose: a list of two numbers giving specific doses of the two controls.
        :param plot_file: the file path to save the plot of z' values. If plot_file is None, do not save the plot.
        :param robust_plot_file: the file path to save the plot of z' robust values. If robust_plot_file is None,
                                do not save the plot.
        :param plot_size: size of the plot
        :return:
            z_prime_result: a dataframe of two columns. The first column is plates. The second column is z' values.
        '''

        if len(control) != 2:
            sys.exit('z_prime can only be calculated for two controls')

        idna = np.where(pd.isna(self.res[self.fea_col]))[0]
        p_list = []
        z_prime = []
        z_prime_robust = []
        if ttest_flag:
            ttest_list = []
            pvalue_list = []
        plates = np.unique(self.res[self.plate_col].iloc[np.where(np.invert(pd.isna(self.res[self.plate_col])))[0]])

        for p in plates:
            idp = np.where(self.res[self.plate_col] == p)[0]
            con_v = {}
            flag_v = True
            c_list = []

            for i_c in range(len(control)):
                c = control[i_c]
                idc = np.where(self.res[self.comp_col] == c)[0]
                if dose is not None:
                    d = dose[i_c]
                    idc = np.intersect1d(idc, np.where(self.res[self.dose_col] == d)[0])
                    c = c + '___' + str(d)
                idpc = np.setdiff1d(np.intersect1d(idp, idc), idna)

                if len(idpc) <= 1:
                    print('z_prime can not be calculated for plate ' + p)
                    flag_v = False
                    break

                con_v[c] = self.res.iloc[idpc, :][self.fea_col].values
                c_list.append(c)

            if flag_v:
                p_list.append(p)
                z_prime.append(z_prime_score(con_v[c_list[0]], con_v[c_list[1]]))
                z_prime_robust.append(z_prime_robust_score(con_v[c_list[0]], con_v[c_list[1]]))
                if ttest_flag:
                    t_stat, p_val = stats.ttest_ind(con_v[c_list[0]], con_v[c_list[1]], equal_var=False)
                    ttest_list.append(t_stat)
                    pvalue_list.append(p_val)

        # Create the plot
        plt.figure(figsize=plot_size)
        plt.plot(p_list, z_prime_robust, marker='o', linestyle='none', markersize=6)
        plt.xlabel('Plate')
        if dose is None:
            plt.ylabel('Z\' robust of ' + control[0] + ' and ' + control[1])
        else:
            plt.ylabel('Z\' robust of ' + control[0] + ' ' + str(dose[0]) + ' and ' + control[1] + ' ' + str(dose[1]))
        plt.title('Z\' robust by Plate')
        plt.xticks(rotation=90)
        plt.ylim(np.min([np.min(z_prime_robust) - 0.1, 0]), 1)
        plt.grid(True)
        plt.tight_layout()
#        plt.show()
        if robust_plot_file is not None:
            plt.savefig(robust_plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        plt.figure(figsize=plot_size)
        plt.plot(p_list, z_prime, marker='o', linestyle='none', markersize=6)
        plt.xlabel('Plate')
        if dose is None:
            plt.ylabel('Z\' of ' + control[0] + ' and ' + control[1])
        else:
            plt.ylabel('Z\' of ' + control[0] + ' ' + str(dose[0]) + ' and ' + control[1] + ' ' + str(dose[1]))
        plt.title('Z\' by Plate')
        plt.xticks(rotation=90)
        plt.ylim(np.min([np.min(z_prime) - 0.1, 0]), 1)
        plt.grid(True)
        plt.tight_layout()
        #        plt.show()
        if plot_file is not None:
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        if ttest_flag:
            z_prime_result = pd.DataFrame({'Plate': p_list, 'Z\'': z_prime, 'Z\'_robust': z_prime_robust, 't-statistic': ttest_list, 'pvalue': pvalue_list})
        else:
            z_prime_result = pd.DataFrame({'Plate': p_list, 'Z\'': z_prime, 'Z\'_robust': z_prime_robust})
        return z_prime_result


    def calculate_CV_per_plate(self, control, dose=None, plot_file=None, plot_size=(10, 6), font_size=10, markersize=6):
        '''
        :param control: a list of compounds for which CV should be calculated. The compound names
                        should be included in the self.comp_col column of the dataframe
        :param dose: a list of numbers giving specific doses of the compounds.
        :param plot_file: the file path to save the plot of CV values. If plot_file is None, do not save the plot.
        :param plot_size: size of the plot
        :return:
            cv_result: a dataframe of CV values. The first column is plates. The second column is compounds. The third
                        column gives doses and is optional. The last column is cv values.
        '''

        idna = np.where(pd.isna(self.res[self.fea_col]))[0]
        cv_list = []
        comp_list = []
        plate_list = []
        plates = np.unique(self.res[self.plate_col].iloc[np.where(np.invert(pd.isna(self.res[self.plate_col])))[0]])
        for p in plates:
            idp = np.where(self.res[self.plate_col] == p)[0]

            for i_c in range(len(control)):
                c = control[i_c]
                idc = np.where(self.res[self.comp_col] == c)[0]
                if dose is not None:
                    d = dose[i_c]
                    idc = np.intersect1d(idc, np.where(self.res[self.dose_col] == d)[0])
                    c = c + '___' + str(d)

                idpc = np.setdiff1d(np.intersect1d(idp, idc), idna)
                if len(idpc) <= 1:
                    print('CV can not be calculated for ' + c + ' on plate ' + p)
                    continue

                tmp_data = self.res.iloc[idpc, :][self.fea_col].values
                if np.std(tmp_data) == 0:
                    cv_list.append(0)
                else:
                    cv_list.append(np.std(tmp_data) / np.mean(tmp_data))
                comp_list.append(c)
                plate_list.append(p)
        cv_list = np.array(cv_list)
        comp_list = np.array(comp_list)
        plate_list = np.array(plate_list)

        colors = get_colors()
        font = {'size': font_size}
        matplotlib.rc('font', **font)
        plt.figure(figsize=plot_size)
        plt.ylim(0, np.max([np.max(cv_list) + 0.1, 1]))
        plt.xlabel('Plate')
        plt.ylabel('Coefficient of Variation')
        plt.title('Coefficient of variation plot by plate')
        rank = 0
        for con in np.unique(comp_list):
            id_con = np.where(comp_list == con)[0]
            color = colors[rank]
            rank = int(np.mod(rank + 1, 20))
            plt.plot(plate_list[id_con], cv_list[id_con], marker='.', color=color, label=con, markersize=markersize, linestyle='none')
        plt.grid(True)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.legend()
#        plt.show()
        if plot_file is not None:
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        if dose is None:
            cv_result = pd.DataFrame({'Plate': plate_list, 'Compound': comp_list, 'CV': cv_list})
        else:
            cv_result = pd.DataFrame({'Plate': plate_list, 'Compound': [i.split('___')[0] for i in comp_list],
                                      'Dose': [float(i.split('___')[1]) for i in comp_list], 'CV': cv_list})
        return cv_result



    def plot_wells_across_plates(self, compound, dose=None, plot_file=None, plot_size=(10, 6), font_size=10, markersize=6):
        '''
        :param compound: a list of compounds that should be plotted. The compound names
                        should be included in the self.comp_col column of the dataframe
        :param dose: a list of numbers giving specific doses of the compounds.
        :param plot_file: the file path to save the plot of CV values. If plot_file is None, do not save the plot.
        :param plot_size: size of the plot
        :return:
            data_result: a dataframe of extracted data values. The first column is plates. The second column is compounds. The third
                        column gives doses and is optional. The last column is data values.
        '''

        idna = np.where(pd.isna(self.res[self.fea_col]))[0]
        data_list = []
        comp_list = []
        plate_list = []
        plates = np.sort(np.unique(self.res[self.plate_col].iloc[np.where(np.invert(pd.isna(self.res[self.plate_col])))[0]]))
        for p in plates:
            idp = np.where(self.res[self.plate_col] == p)[0]

            for i_c in range(len(compound)):
                c = compound[i_c]
                idc = np.where(self.res[self.comp_col] == c)[0]
                if dose is not None:
                    d = dose[i_c]
                    idc = np.intersect1d(idc, np.where(self.res[self.dose_col] == d)[0])
                    c = c + '___' + str(d)
                idpc = np.setdiff1d(np.intersect1d(idp, idc), idna)

                if len(idpc) == 0:
                    print('No data available for ' + c + ' on plate ' + p)

                tmp_data = self.res.iloc[idpc, :][self.fea_col].values
                data_list = data_list + list(tmp_data)
                comp_list = comp_list + [c for i in range(len(tmp_data))]
                plate_list = plate_list + [p for i in range(len(tmp_data))]
        data_list = np.array(data_list)
        comp_list = np.array(comp_list)
        plate_list = np.array(plate_list)

        colors = get_colors()
        font = {'size': font_size}
        matplotlib.rc('font', **font)
        plt.figure(figsize=plot_size)
        plt.ylim(np.min(data_list), np.max(data_list))
        plt.xlabel('Plate')
        plt.ylabel(self.fea_col)
        plt.title(self.fea_col + ' plot by plate')
        rank = 0
        for con in np.unique(comp_list):
            id_con = np.where(comp_list == con)[0]
            color = colors[rank]
            rank = int(np.mod(rank + 1, 20))
            plt.plot(plate_list[id_con], data_list[id_con], marker='.', color=color, label=con, markersize=markersize, linestyle='none')
        plt.grid(True)
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.legend()
#        plt.show()
        if plot_file is not None:
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
        plt.close()

        if dose is None:
            data_result = pd.DataFrame({'Plate': plate_list, 'Compound': comp_list, 'Data': data_list})
        else:
            data_result = pd.DataFrame({'Plate': plate_list, 'Compound': [i.split('___')[0] for i in comp_list],
                                        'Dose': [float(i.split('___')[1]) for i in comp_list], 'Data': data_list})
        return data_result



    def normalization_within_plate(self, norm_method='division'):
        '''
            This function adds a column to the res dataframe, which is the normalized viability data.
            The column name is 'Normalized_' + self.res.fea_col
            :param norm_method: can be 'division', 'add1subtract'
        '''

        x = np.array([np.nan for i in range(self.res.shape[0])])

        idna = np.where(pd.isna(self.res[self.fea_col]))[0]

        norm_v = {}
        for i in range(self.res.shape[0]):
            if pd.isna(self.res.iloc[i, :][self.fea_col]):
                continue

            s = self.res.iloc[i, :][self.solvent_col]
            p = self.res.iloc[i, :][self.plate_col]

            if p not in norm_v.keys():
                norm_v[p] = {}

            if s not in norm_v[p].keys():
                idps = np.intersect1d(np.where(self.res[self.plate_col] == p)[0], np.where(self.res[self.comp_col] == s)[0])
                idps = np.setdiff1d(idps, idna)
                if len(idps) == 0:
                    print('No value available for plate ' + p + ' and solvent ' + s)
                    continue
                norm_v[p][s] = {'norm_value': np.mean(self.res.iloc[idps, :][self.fea_col]), 'num_norm': len(idps)}

            if norm_method == 'division':
                x[i] = self.res.iloc[i, :][self.fea_col] / norm_v[p][s]['norm_value']
            elif norm_method == 'add1subtract':
                x[i] = self.res.iloc[i, :][self.fea_col] + 1 - norm_v[p][s]['norm_value']

        self.res['Normalized_' + self.fea_col] = x
        self.norm_fea_col = 'Normalized_' + self.fea_col



    def identify_extreme_values(self, abnormal_pattern='high', method='mean', factor=0.05, plot_folder=None,
                                plot_file_str='', plot_size=(10, 6), font_size=10, markersize=6):
        '''
        This function should be run based on normalized data.
        :param abnormal_pattern: 'high' or 'low', indicating whether we shall identify extreme values higher or lower
                                than solvent controls
        :param method: 'mean' or 'std', indicating what method is used to detect extreme values
        :param factor: a factor number used to detect extreme values
        :param plot_folder: a string giving the folder to save the plots of extreme values
        :param plot_file_str: a string to form the file names of extreme value plots
        :param plot_size:
        :param font_size:
        :param markersize:
        :return:
            result: a dataframe of detected extreme values. The first column is models. The second column is compounds.
                The third column is plates. The fourth column is the number of wells with extreme values. The fifth column is
                the proportion of wells with extreme values. The sixth column is the average portion beyond the average
                solvent value.
            extreme_id: indices of all the extreme values that are detected.
        '''

        # idna = np.where(pd.isna(self.res[self.fea_col]))[0]
        # all_drug = np.setdiff1d(list(pd.unique(self.res[self.comp_col])), ['', 'nan'] + self.control_types)
        # all_cell = np.setdiff1d(list(pd.unique(self.res[self.model_col])), ['', 'nan'])
        # cell_l = []
        # drug_l = []
        # plate_l = []
        # number_exceed_l = []    # Number of wells exceeding the average solvent value
        # proportion_exceed_l = []
        # average_exceed_l = []   # Average portion (decimal) exceeding the average solvent value
        # extreme_id = []
        #
        # for cell in all_cell:
        #     id_c = np.setdiff1d(np.where(self.res[self.model_col] == cell)[0], idna)
        #
        #     for plate in pd.unique(self.res.iloc[id_c, :][self.plate_col]):
        #         id_c_p = np.intersect1d(id_c, np.where(self.res[self.plate_col] == plate)[0])
        #
        #         for drug in all_drug:
        #             if drug in self.res.iloc[id_c_p, :][self.comp_col].values:
        #                 id_c_p_d = np.intersect1d(id_c_p, np.where(self.res[self.comp_col] == drug)[0])
        #                 drug_value = self.res.iloc[id_c_p_d, :][self.fea_col].values
        #
        #                 if abnormal_pattern == 'high':
        #                     id = np.where(drug_value > 1 + factor)[0]
        #                 if abnormal_pattern == 'low':
        #                     id = np.where(drug_value < 1 - factor)[0]
        #                 if len(id) > 0:
        #                     number_exceed_l.append(int(len(id)))
        #                     proportion_exceed_l.append(round(int(len(id)) / len(drug_value), 3))
        #                     average_exceed_l.append(round(np.mean(drug_value[id] - 1), 3))
        #                     cell_l.append(cell)
        #                     drug_l.append(drug)
        #                     plate_l.append(plate)
        #                     extreme_id = extreme_id + list(id_c_p_d[id])
        #
        # result = pd.DataFrame({'Cell': cell_l, 'Drug': drug_l, 'Plate': plate_l, 'Number_Exceed': number_exceed_l,
        #                        'Proportion_Exceed': proportion_exceed_l, 'Average_Exceed_Portion': average_exceed_l})
        # extreme_id = np.sort(extreme_id)
        # return result, extreme_id



        idna = np.where(pd.isna(self.res[self.norm_fea_col]))[0]
        all_drug = np.setdiff1d(list(pd.unique(self.res[self.comp_col])), ['', 'nan'] + self.control_types)
        all_cell = np.setdiff1d(list(pd.unique(self.res[self.model_col])), ['', 'nan'])
        cell_l = []
        drug_l = []
        plate_l = []
        number_exceed_l = []    # Number of wells exceeding the average solvent value
        proportion_exceed_l = []
        average_exceed_l = []   # Average portion (decimal) exceeding the average solvent value
        extreme_id = []
        solvent_v = {}

        for cell in all_cell:
            id_c = np.setdiff1d(np.where(self.res[self.model_col] == cell)[0], idna)
            if cell not in solvent_v.keys():
                solvent_v[cell] = {}

            for plate in pd.unique(self.res.iloc[id_c, :][self.plate_col]):
                id_c_p = np.intersect1d(id_c, np.where(self.res[self.plate_col] == plate)[0])
                if plate not in solvent_v[cell].keys():
                    solvent_v[cell][plate] = {}

                for drug in all_drug:
                    if drug in self.res.iloc[id_c_p, :][self.comp_col].values:
                        id_c_p_d = np.intersect1d(id_c_p, np.where(self.res[self.comp_col] == drug)[0])
                        solvent = self.res.iloc[id_c_p_d, :][self.solvent_col].values[0]
                        if solvent not in solvent_v[cell][plate].keys():
                            id_c_p_s = np.intersect1d(id_c_p, np.where(self.res[self.comp_col] == solvent)[0])
                            solvent_v[cell][plate][solvent] = {'mean': np.mean(self.res.iloc[id_c_p_s, :][self.norm_fea_col]),
                                                               'std': np.std(self.res.iloc[id_c_p_s, :][self.norm_fea_col])}
                        drug_value = self.res.iloc[id_c_p_d, :][self.norm_fea_col].values

                        if abnormal_pattern == 'high':
                            if method == 'mean':
                                id = np.where(drug_value > (1 + factor) * solvent_v[cell][plate][solvent]['mean'])[0]
                            if method == 'std':
                                id = np.where(drug_value > solvent_v[cell][plate][solvent]['mean'] +
                                              factor * solvent_v[cell][plate][solvent]['std'])[0]

                        if abnormal_pattern == 'low':
                            if method == 'mean':
                                id = np.where(drug_value < (1 - factor) * solvent_v[cell][plate][solvent]['mean'])[0]
                            if method == 'std':
                                id = np.where(drug_value < solvent_v[cell][plate][solvent]['mean'] -
                                              factor * solvent_v[cell][plate][solvent]['std'])[0]

                        if len(id) > 0:
                            number_exceed_l.append(int(len(id)))
                            proportion_exceed_l.append(round(int(len(id)) / len(drug_value), 3))
                            if method == 'mean':
                                average_exceed_l.append(round(np.mean((drug_value[id] - solvent_v[cell][plate][solvent]['mean']) / solvent_v[cell][plate][solvent]['mean']), 3))
                            if method == 'std':
                                average_exceed_l.append(round(np.mean((drug_value[id] - solvent_v[cell][plate][solvent]['mean']) / solvent_v[cell][plate][solvent]['std']), 3))
                            cell_l.append(cell)
                            drug_l.append(drug)
                            plate_l.append(plate)
                            extreme_id = extreme_id + list(id_c_p_d[id])

        result = pd.DataFrame({'Cell': cell_l, 'Drug': drug_l, 'Plate': plate_l, 'Number_Exceed': number_exceed_l,
                               'Portion_Exceed': proportion_exceed_l, 'Average_Exceed_Proportion': average_exceed_l})
        extreme_id = np.sort(extreme_id)

        all_plates = np.sort(pd.unique(self.res.iloc[np.setdiff1d(list(range(self.res.shape[0])), idna), :][self.plate_col]))
        num_exceed = []
        tot_exceed_proportion = []
        for p in all_plates:
            id_p = np.where(result['Plate'] == p)[0]
            num_exceed.append(np.sum(result.iloc[id_p, :]['Number_Exceed']))
            if len(id_p) == 0:
                tot_exceed_proportion.append(0)
            else:
                tot_exceed_proportion.append(np.sum(result.iloc[id_p, :]['Number_Exceed'] * result.iloc[id_p, :]['Average_Exceed_Proportion']))

        if plot_folder is not None:
            # Plot the numbers of extreme values
            font = {'size': font_size}
            matplotlib.rc('font', **font)
            plt.figure(figsize=plot_size)
            # plt.ylim(np.min(num_exceed) - 2, np.max(num_exceed) + 2)
            plt.xlabel('Plate')
            plt.ylabel('Number of extreme values')
            plt.title('Number of extreme values by plates')
            plt.plot(all_plates, num_exceed, marker='.', color='blue', markersize=markersize, linestyle='none')
            plt.grid(True)
            plt.xticks(rotation=90)
            plt.tight_layout()
            # plt.legend()
            plt.savefig(plot_folder + '/' + plot_file_str + '_Number_Of_Extreme_Values.png', dpi=300,
                        bbox_inches='tight')
            plt.close()

            # Plot the exceeding proportions of extreme values
            font = {'size': font_size}
            matplotlib.rc('font', **font)
            plt.figure(figsize=plot_size)
            # plt.ylim(np.min(data_list), np.max(data_list))
            plt.xlabel('Plate')
            plt.ylabel('Total exceeding proportion')
            plt.title('Total exceeding proportion by plates')
            plt.plot(all_plates, tot_exceed_proportion, marker='.', color='blue', markersize=markersize, linestyle='none')
            plt.grid(True)
            plt.xticks(rotation=90)
            plt.tight_layout()
            # plt.legend()
            plt.savefig(plot_folder + '/' + plot_file_str + '_Total_Exceeding_Proportion.png', dpi=300,
                        bbox_inches='tight')
            plt.close()

        return result, extreme_id



        # res = copy.deepcopy(self.res[np.invert(pd.isna(self.res[self.fea_col]))])
        # all_drug = np.setdiff1d(list(pd.unique(res[self.comp_col])), ['', 'nan'] + self.control_types)
        # all_cell = np.setdiff1d(list(pd.unique(res[self.model_col])), ['', 'nan'])
        # cell_l = []
        # drug_l = []
        # plate_l = []
        # number_exceed_l = []    # Number of wells exceeding the average solvent value
        # proportion_exceed_l = []
        # average_exceed_l = []   # Average portion (decimal) exceeding the average solvent value
        # solvent_v = {}
        #
        # for cell in all_cell:
        #     res_c = res[res[self.model_col] == cell].copy()
        #     if cell not in solvent_v.keys():
        #         solvent_v[cell] = {}
        #
        #     for plate in pd.unique(res_c[self.plate_col]):
        #         res_c_p = res_c[res_c[self.plate_col] == plate].copy()
        #         if plate not in solvent_v[cell].keys():
        #             solvent_v[cell][plate] = {}
        #
        #         for drug in all_drug:
        #             if drug in res_c_p[self.comp_col].values:
        #                 res_c_p_d = res_c_p[res_c_p[self.comp_col] == drug].copy()
        #                 solvent = res_c_p_d[self.solvent_col].values[0]
        #                 if solvent not in solvent_v[cell][plate].keys():
        #                     res_c_p_s = res_c_p[res_c_p[self.comp_col] == solvent].copy()
        #                     solvent_v[cell][plate][solvent] = np.mean(res_c_p_s[self.fea_col])
        #
        #                 drug_value = res_c_p_d[self.fea_col].values
        #                 num_exceed = 0
        #                 if abnormal_pattern == 'high':
        #                     id = np.where(drug_value > (1 + factor) * solvent_v[cell][plate][solvent])[0]
        #                     if len(id) > 0:
        #                         number_exceed_l.append(int(len(id)))
        #                         proportion_exceed_l.append(round(int(len(id)) / len(drug_value), 3))
        #                         average_exceed_l.append(round(np.mean((drug_value[id] - solvent_v[cell][plate][solvent]) / solvent_v[cell][plate][solvent]), 3))
        #                 if abnormal_pattern == 'low':
        #                     id = np.where(drug_value < (1 - factor) * solvent_v[cell][plate][solvent])[0]
        #                     if len(id) > 0:
        #                         number_exceed_l.append(int(len(id)))
        #                         proportion_exceed_l.append(round(int(len(id)) / len(drug_value), 3))
        #                         average_exceed_l.append(round(np.mean((solvent_v[cell][plate][solvent] - drug_value[id]) / solvent_v[cell][plate][solvent]), 3))
        #                 if num_exceed > 0:
        #                     cell_l.append(cell)
        #                     drug_l.append(drug)
        #                     plate_l.append(plate)
        #
        # result = pd.DataFrame({'Cell': cell_l, 'Drug': drug_l, 'Plate': plate_l, 'Number_Exceed': number_exceed_l,
        #                        'Proportion_Exceed': proportion_exceed_l, 'Average_Exceed_Portion': average_exceed_l})
        # return result, extreme_id



    def identify_reverse_pattern(self, abnormal_pattern='high', factor=0.05, plot_save_folder=None, plot_size=(10, 6),
                                 font_size=10, markersize=3):
        '''
        This function should be run on normalized data.
        :param abnormal_pattern: 'high' or 'low', indicating whether an increasing or decreasing pattern is abnormal.
        :param factor: a factor number used to detect reverse patterns
        :return:
            result: a dataframe of detected reverse patterns. The first column is models. The second column is compounds.
                The third column is experiment IDs.
        '''

        idna = np.where(pd.isna(self.res[self.norm_fea_col]))[0]
        cells = np.setdiff1d(list(pd.unique(self.res[self.model_col])), ['', 'nan'])
        cell_l = []
        drug_l = []
        exp_l = []

        for cell in cells:
            id_c = np.setdiff1d(np.where(self.res[self.model_col] == cell)[0], idna)
            drugs = np.setdiff1d(list(pd.unique(self.res.iloc[id_c, :][self.comp_col])), ['', 'nan'] + self.control_types)

            for drug in drugs:
                id_c_d = np.intersect1d(id_c, np.where(self.res[self.comp_col] == drug)[0])
                exps = np.setdiff1d(list(pd.unique(self.res.iloc[id_c_d, :][self.exp_col])), ['', 'nan'])

                for exp in exps:
                    id_c_d_e = np.intersect1d(id_c_d, np.where(self.res[self.exp_col] == exp)[0])
                    data = self.res.iloc[id_c_d_e, :][self.norm_fea_col].values
                    dose = self.res.iloc[id_c_d_e, :][self.dose_col].values

                    unique_dose = np.sort(np.unique(dose))
                    mean_measure = []
                    for i in range(len(unique_dose)):
                        dose_id_i = np.where(dose == unique_dose[i])[0]
                        mean_measure.append(np.mean(data[dose_id_i]))
                    abnormal_flag = []
                    for i in range(len(mean_measure) - 1):
                        if abnormal_pattern == 'high':
                            abnormal_flag.append(mean_measure[i + 1] - mean_measure[i] >= factor)
                        if abnormal_pattern == 'low':
                            abnormal_flag.append(mean_measure[i] - mean_measure[i + 1] >= factor)
                    if 'TrueTrue' in ''.join([str(i) for i in abnormal_flag]):
                        flag = True
                    else:
                        flag = False
                    if flag:
                        cell_l.append(cell)
                        drug_l.append(drug)
                        exp_l.append(exp)

                        if plot_save_folder is not None:
                            if not os.path.exists(plot_save_folder):
                                os.makedirs(plot_save_folder)

                            font = {'size': font_size}
                            matplotlib.rc('font', **font)
                            plt.figure(figsize=plot_size)
                            plt.xlim(np.min(dose) - 0.5, np.max(dose) + 0.5)
                            plt.ylim(np.min(data) - 0.2, np.max(data) + 0.2)
                            plt.xlabel('Dose (log10(M))')
                            plt.ylabel(self.norm_fea_col)
                            plt.title(cell + '--' + drug + '--' + exp)
                            plt.plot(dose, data, marker='.', markersize=markersize, linestyle='none')
                            plt.grid(True)
                            plt.tight_layout()
                            file_name = plot_save_folder + '/' + cell + '--' + drug + '--' + exp + '.png'
                            plt.savefig(file_name, dpi=360)
                            plt.close()

        result = pd.DataFrame({'Cell': cell_l, 'Drug': drug_l, 'Experiment_ID': exp_l})
        return result



    def fit_dose_response_curve(self, param_bound, plot_save_folder):
        """
        :param param_bound:
        :param plot_dir: directory to store dose response curve plots
        :return:
        """

        if not os.path.exists(plot_save_folder):
            os.makedirs(plot_save_folder)

        idna = np.where(pd.isna(self.res[self.norm_fea_col]))[0]
        plate_id = np.setdiff1d(np.where(np.invert(pd.isna(self.res[self.plate_col])))[0], idna)
        plates = np.unique(self.res[self.plate_col].iloc[plate_id])
        model_id = np.setdiff1d(np.where(np.invert(pd.isna(self.res[self.model_col])))[0], idna)
        models = np.unique(self.res[self.model_col].iloc[model_id])

        control = {}
        for model in models:
            m_id = np.where(self.res[self.model_col] == model)[0]
            for plate in plates:
                p_id = np.intersect1d(m_id, np.where(self.res[self.plate_col] == plate)[0])
                for comp in self.control_types:
                    c_id = np.intersect1d(p_id, np.where(self.res[self.comp_col] == comp)[0])
                    c_id = np.setdiff1d(c_id, idna)
                    if len(c_id) > 0:
                        if model not in control.keys():
                            control[model] = {}
                        if plate not in control[model].keys():
                            control[model][plate] = {}
                        control[model][plate][comp] = self.res[self.fea_col].iloc[c_id].values

        data = {}
        all_exp = np.setdiff1d(list(pd.unique(self.res[self.exp_col])), ['', 'nan'])
        for exp in all_exp:

            res_e = self.res[self.res[self.exp_col] == exp].copy()

            if len(pd.unique(res_e[self.model_col])) > 1 or len(pd.unique(res_e[self.comp_col])) > 1:
                print('More than one specimens or drugs in experiment ' + exp)
                print(res_e)
                continue

            fit = dose_response_curve(compound=res_e[self.comp_col].values[0],
                                      specimen=res_e[self.model_col].values[0],
                                      study_id=self.study_id, exp_id=exp, model_type='4_parameter_hillslop',
                                      control_exp_flag=res_e[self.control_col].values[0], feature=self.norm_fea_col)

            d_v = res_e[self.dose_col].values
            m_v = res_e[self.norm_fea_col].values
            id_d_m = np.intersect1d(np.intersect1d(np.where(np.invert(np.isnan(d_v)))[0],
                                                   np.where(np.invert(np.isinf(d_v)))[0]),
                                    np.intersect1d(np.where(np.invert(np.isnan(m_v)))[0],
                                                   np.where(np.invert(np.isinf(m_v)))[0]))
            d_v = d_v[id_d_m]
            m_v = m_v[id_d_m]
            metrics = fit.compute_fit_metrics(dose=d_v, measure=m_v, param_bound=param_bound, mode='trapz')

            metrics['R2_orig'] = metrics['R2fit'] 
            metrics['flipped'] = False

            # Notice that param_bound = ([low_Einf, low_EC50, low_HS, low_E0], [up_Einf, up_EC50, up_HS, up_E0])
            # The logic here is:
            # The scale of doses is log10(M)
            # If we use a parameter bound of param_bound = ([0, -19, -5, 1], [1, 5, 5, 1 + 0.000001]) for fitting and get
            # a negative HS, the negative HS means that E_inf is used to fit low doses and E_0 is used to fit high doses.
            # Because E_0 is actually bounded at 1 and E_inf is bounded within [0, 1]. The low doses must have viability values
            # smaller than high doses, which is not normal. In this case, we relax the bound to param_bound = ([0, -19, 0, 0], [inf, 5, 5, inf]),
            # and refit the data. After refitting, if R2 is improved, we will use the refitting result instead of the original result.
            # The positive HS will guarantee E_inf is used to fit high doses and E_0 is used to fit low doses

            if metrics['HS'] < 0:
                fit_neg_HS = copy.deepcopy(fit)
                param_bound_neg_HS = copy.deepcopy(param_bound)

                # Note: Unbounded E0/Einf leads to some cases assuming a very high E0/Einf, resulting in extreme AUC values

                param_bound_neg_HS[0][3] = 0
                param_bound_neg_HS[1][0] = np.inf
            
                param_bound_neg_HS[0][2] = 0
                param_bound_neg_HS[1][3] = np.inf
        
                metrics_neg_HS = fit_neg_HS.compute_fit_metrics(dose=d_v, measure=m_v, param_bound=param_bound_neg_HS, mode='trapz')

                # Create boolean flag to keep track of which curves were inverted 
                
                metrics_neg_HS['R2_orig'] = metrics['R2fit']
                metrics_neg_HS['R2_flipped'] = metrics_neg_HS['R2fit']
                metrics_neg_HS['flipped'] = True

                if metrics_neg_HS['R2fit'] > metrics['R2fit']:
                    fit = copy.deepcopy(fit_neg_HS)
                    metrics = metrics_neg_HS


                # if metrics_neg_HS['R2fit'] > metrics['R2fit']:
                #     fit = copy.deepcopy(fit_neg_HS)
                #     metrics = copy.deepcopy(metrics_neg_HS)

            # Plots are saved by experiment ID, and some have problematic characters for file naming so this helper function cleans the names before saving plots

            def safe_filename(s):
                return re.sub(r'[^A-Za-z0-9_.-]+', '_', s)



            fit.plot_curve(directory=plot_save_folder, save_file_name=safe_filename(res_e[self.model_col].values[0] + '--' +
                res_e[self.comp_col].values[0] + '--' + exp + '--' + self.norm_fea_col))

            if res_e[self.control_col].values[0]:
                cell = res_e[self.model_col].values[0]
                plate = res_e[self.plate_col].values[0]
                if res_e[self.comp_col].values[0] not in control[cell][plate].keys():
                    control[cell][plate][res_e[self.comp_col].values[0]] = {}
                control[cell][plate][res_e[self.comp_col].values[0]][self.norm_fea_col] = fit
            else:
                data[res_e[self.model_col].values[0] + '--' + res_e[self.comp_col].values[0] + '--' +
                     exp + '--' + self.norm_fea_col] = fit

        return data, control



    def copy(self):
        return copy.deepcopy(self)






def update_dict_recursively(base_dict, update_dict):
    """
    Recursively updates a dictionary with values from another dictionary.

    Args:
        base_dict (dict): The dictionary to be updated.
        update_dict (dict): The dictionary containing the updates.
    """

    for key, value in update_dict.items():
        if key in base_dict and isinstance(base_dict[key], dict) and isinstance(value, dict):
            # If both are dictionaries, recursively update
            update_dict_recursively(base_dict[key], value)
        else:
            # Otherwise, update the value directly
            base_dict[key] = value



def read_feature_category_file(file):
    file_handle = open(file, 'r')
    content = file_handle.read()
    file_handle.close()
    content = content.split('\n')
    params = {}
    for i in range(len(content)):
        params[content[i].split(' = ')[0]] = content[i].split(' = ')[1]

    return params



def dose_combination_fitting(ori_drc, plot_dir):

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    cell_l = []
    drug_l = []
    exp_l = []
    f_l = []
    dose_l = []
    num_dose_l = []
    ori_num_dose_l = []
    auc_l = []
    ori_auc_l = []

    drc_dict = {}
    dose_number_evaluation_result = None
    unique_dose = np.sort(pd.unique(ori_drc.dose))
    if len(unique_dose) > 4:
        for num_dose in range(4, len(unique_dose)):
            combinations = itertools.combinations(unique_dose, num_dose)
            for dose_com in combinations:
                dose_com = np.sort(dose_com)
                drc = copy.deepcopy(ori_drc)
                dose_id = np.where(np.isin(drc.dose, dose_com))[0]

                metrics = drc.compute_fit_metrics(dose=drc.dose[dose_id],
                                                  measure=drc.measure[dose_id],
                                                  param_bound=drc.param_bound, mode='trapz')

                cell_l.append(drc.specimen)
                drug_l.append(drc.compound)
                exp_l.append(drc.exp_id)
                f_l.append(drc.feature)
                dose_com_str = ','.join([str(d) for d in dose_com])
                dose_l.append(dose_com_str)
                num_dose_l.append(num_dose)
                ori_num_dose_l.append(len(unique_dose))
                auc_l.append(drc.fit_metrics.AUC)
                ori_auc_l.append(ori_drc.fit_metrics.AUC)
                fitting_label = drc.specimen + '--' + drc.compound + '--' + drc.exp_id + '--' + drc.feature + '||' + \
                                dose_com_str
                drc_dict[fitting_label] = drc

                drc.plot_curve(directory=plot_dir, save_file_name=fitting_label)

        dose_number_evaluation_result = pd.DataFrame({'Cell': cell_l, 'Drug': drug_l, 'Experiment_ID': exp_l,
                                                      'Feature': f_l, 'Dose': dose_l, 'Number_Dose': num_dose_l,
                                                      'Original_Number_Dose': ori_num_dose_l, 'AUC': auc_l,
                                                      'Original_AUC': ori_auc_l})
    return dose_number_evaluation_result, drc_dict



# Fit dose response curves to experiments
def fit_dose_response_curve(study_id, res, info_col, fc, ec50_bound, hs_bound, all_exp, ori_control, plot_dir):
    """
    study_id: dataset/study name
    res: response dataframe
    info_col: names of columns that are information not data
    fc: feature category information
    ec50_bound: bounds on ec50 parameter
    hs_bound: bounds on hill slope parameter
    all_exp: all experimental IDs
    ori_control: all original readouts of control wells
    plot_dir: directory to store dose response curve plots
    """



    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    control = copy.deepcopy(ori_control)

    data = {}
    for f in fc.keys():
        if f not in res.columns:
            continue

        c = fc[f]

        # The feature category label has three parts, separated by '___'.
        # The first part is either 'NoNorm' or 'Norm'. 'Norm' indicates normalization to solvent is needed,
        # while 'NoNorm' indicates not.
        # The second and third parts are for the parameters of efficacies at 0 and infinity concentrations, respectively.
        # Sub-elements for specifying the bounds of these two parameters are separated by '_'.
        # The first sub-element is either 'E0' or 'Einf' indicating which parameter it is for.
        # If there are two sub-elements, the second sub-element is either 'Free' or a number.
        # 'Free' means there is no bound on the parameter. A number indicates the parameter value is fixed at this number.
        # If there are three sub-elements, the second and third sub-elements are the lower and upper bounds of the parameter.

        norm_flag = c.split('___')[0]
        if norm_flag == 'NoNorm':
            norm_flag = False
        elif norm_flag == 'Norm':
            norm_flag = True
        else:
            print('Wrong normalization flag for ' + f)

        min_E0, max_E0 = generate_parameter_bound(c.split('___')[1])
        min_Einf, max_Einf = generate_parameter_bound(c.split('___')[2])
        min_ec50 = ec50_bound[0]
        max_ec50 = ec50_bound[1]
        min_hs = hs_bound[0]
        max_hs = hs_bound[1]
        param_bound = ([min_Einf, min_ec50, min_hs, min_E0], [max_Einf, max_ec50, max_hs, max_E0])
        res_f = res.loc[:, np.concatenate((info_col, [f]))].copy()

        for exp in all_exp:
            res_f_e = res_f[res_f['exp.id'] == exp].copy()

            if len(pd.unique(res_f_e['plate.layout.info.Model'])) > 1 or len(pd.unique(res_f_e['compound.name'])) > 1:
                print('More than one cell lines or drugs in the experiment')
                print(pd.unique(res_f_e['plate.layout.info.Model']))
                print(pd.unique(res_f_e['compound.name']))
                continue

            control_exp_flag = (res_f_e['plate.layout.info.Treatment'].values[0] == 'control')
            if norm_flag:
                if 'Standard deviation' in f:
                    norm_f = f.split('Standard deviation')[0] + 'Mean' + f.split('Standard deviation')[1]
                else:
                    norm_f = f
                norm_f_flag = True
                for plate in pd.unique(res_f_e['plate.name']):
                    cell = res_f_e['plate.layout.info.Model'].values[0]
                    solvent = res_f_e['Solvent'].values[0]
                    if norm_f not in control[cell][plate][solvent].keys():
                        measure_before_norm = None
                        norm_f_flag = False
                        break
                if norm_f_flag:
                    measure_before_norm = res_f_e[f].copy().values
                    for plate in pd.unique(res_f_e['plate.name']):
                        cell = res_f_e['plate.layout.info.Model'].values[0]
                        solvent = res_f_e['Solvent'].values[0]
                        well_id = np.where(res_f_e['plate.name'] == plate)[0]
                        if f == 'Fraction_living_cells':
                            res_f_e.iloc[well_id, -1] = res_f_e.iloc[well_id, -1] + (1 - np.mean(control[cell][plate][solvent][norm_f]))
                        else:
                            res_f_e.iloc[well_id, -1] = res_f_e.iloc[well_id, -1] / np.mean(control[cell][plate][solvent][norm_f])
            else:
                measure_before_norm = None

            fit = dose_response_curve(compound=res_f_e['compound.name'].values[0],
                                      specimen=res_f_e['plate.layout.info.Model'].values[0],
                                      study_id=study_id, exp_id=exp, model_type='4_parameter_hillslop',
                                      control_exp_flag=control_exp_flag, feature=f, measure_before_norm=measure_before_norm)
            metrics = fit.compute_fit_metrics(dose=res_f_e['plate.layout.info.Treatment dose LOG(M)'].values,
                                              measure=res_f_e[f].values,
                                              param_bound=param_bound, mode='trapz')

            fit.plot_curve(directory=plot_dir, save_file_name=res_f_e['plate.layout.info.Model'].values[0] + '--' +
                                                              res_f_e['compound.name'].values[0] + '--' + exp + '--' + f)

            if control_exp_flag:
                cell = res_f_e['plate.layout.info.Model'].values[0]
                plate = res_f_e['plate.name'].values[0]
                if res_f_e['compound.name'].values[0] not in control[cell][plate].keys():
                    control[cell][plate][res_f_e['compound.name'].values[0]] = {}
                control[cell][plate][res_f_e['compound.name'].values[0]][f] = fit
            else:
                data[res_f_e['plate.layout.info.Model'].values[0] + '--' + res_f_e['compound.name'].values[0] + '--' +
                     exp + '--' + f] = fit

    return data, control


    # output = open(data_file_path, 'wb')
    # cp.dump(data, output)
    # cp.dump(control, output)
    # cp.dump(res, output)
    # cp.dump(fc, output)
    # cp.dump(all_exp, output)
    # output.close()



def response_curve_4p(x, einf, ec50, hs, e0):
    """ transformed the original function with ec50 in log10(M) instead of M
    """
    return einf + (e0 - einf) / (1 + 10 ** ((x - ec50) * hs))

def response_integral_4p(x, einf, ec50, hs, e0):
    return (einf - e0) * np.log10(1 + 10 ** ((x - ec50) * hs)) / hs + e0 * x

def dose_at_response_4p(y, einf, ec50, hs, e0):
    return ec50 + np.log10((e0 - y) / (y - einf)) / hs

def compute_area_4p(x1, x2, einf, ec50, hs, e0, mode='trapz'):
    popt = (einf, ec50, hs, e0)
    if mode == 'trapz':
        # trapezoidal numerical integrationcollapse
        xx = np.linspace(x1, x2, 100)
        yy = response_curve_4p(xx, *popt)
        area = np.trapz(yy, xx, dx=0.01)
    else:
        # the integral function can be expressed analytically
        # but sometimes less accurate due to float precision issues
        area = response_integral_4p(x2, *popt) - response_integral_4p(x1, *popt)
    return area



def curve_fit_4p(dose, measure, param_bound, method='trf'):
    """
    :param dose: should be in log10(M)
    :param measure:
    :return:
    """
    param = None
    param_cov = None
    if dose is not None and measure is not None and len(dose) == len(measure):
        nfev = 100 * len(dose)
        while param is None and nfev < 10000:
            try:
                # print(f"DOSE: {dose}")
                # print(f"MEASURE: {measure}")
                # print(f"PARAM_BOUND {param_bound}")
                param, param_cov = curve_fit(response_curve_4p, dose, measure,
                                             bounds=param_bound, max_nfev=nfev, method=method)
                # print(f"PARAM: {param}")
                # print(f"PARAM_COV: {param_cov}")
            except RuntimeError:
                pass
            nfev *= 2
    else:
        print('Input dose and measure are wrong or mismatched.')
    return param, param_cov


def generate_parameter_bound(param_flag, min_diff=0.000001):

    param_flag_split = param_flag.split('_')
    len_param_flag = len(param_flag_split)
    if len_param_flag == 2:
        flag = param_flag_split[1]
        if flag == 'Free':
            lower = -np.inf
            upper = np.inf
        else:
            lower = float(flag)
            upper = float(flag) + min_diff
    if len_param_flag == 3:
        if param_flag_split[1] == '-inf':
            lower = -np.inf
        else:
            lower = float(param_flag_split[1])
        if param_flag_split[2] == 'inf':
            upper = np.inf
        else:
            upper = float(param_flag_split[2])

    return lower, upper



class dose_response_curve:
    def __init__(self, compound, specimen, study_id, exp_id, model_type, control_exp_flag, feature,
                 measure_before_norm=None, num_sec=5, R2t=0.5):
        self.compound = compound    # Compound name
        self.specimen = specimen    # Specimen name
        self.study_id = study_id    # Study name
        self.exp_id = exp_id        # A generated ID to indicate different experiments. One dose response curve is fitted for each experiment
        self.dose = None            # Log10 dose
        self.measure = None         # Measurements
        self.model_type = model_type    # The model used for dose response curve fitting
        self.param = None           # Estimated parameters of einf, ec50, hs, e0
        self.param_cov = None       # Covariance matrix of parameter estimates for einf, ec50, hs, e0
        self.param_bound = None  # Bounds of parameters used in model fitting
        self.fit_metrics = None     # Response metrics obtained from model fitting
        self.control_exp_flag = control_exp_flag    # True or False, indicating whether it is a control experiment of compound
        self.feature = feature      # The feature name of readout
        self.measure_before_norm = measure_before_norm  # Feature/readout values before normalization, if the feature needs to be normalized; otherwise, None.
        self.num_sec = num_sec      # If R2 of fit < 0.5, divide the range of EC50 into num_sec sections, and refit the model to test whether a better model can be obtained.
        self.R2t = R2t


    def compute_area(self, d1=-10, d2=-4, mode='trapz'):
        area = None
        if self.param is not None:
            if mode == 'trapz':
                xx = np.linspace(d1, d2, 100)
                yy = response_curve_4p(xx, *self.param)
                area = np.trapz(yy, xx, dx=0.01)
            elif mode == 'integral':
                # the integral function can be expressed analytically
                # but sometimes less accurate due to float precision issues
                area = response_integral_4p(d2, *self.param) - \
                       response_integral_4p(d1, *self.param)
        return area


    def compute_fit_metrics(self, dose, measure, param_bound, d1=-10, d2=-4, mode='trapz'):
        self.dose = copy.deepcopy(dose)
        self.measure = copy.deepcopy(measure)
        self.param_bound = copy.deepcopy(param_bound)
        self.ori_param_bound = copy.deepcopy(param_bound)

        self.param, self.param_cov = curve_fit_4p(dose=self.dose, measure=self.measure, param_bound=self.param_bound)

        if self.param is not None:
            ypred = response_curve_4p(self.dose, *self.param)
            r2 = r2_score(self.measure, ypred)
        else:
            r2 = None

        if self.param is None or r2 < self.R2t:
            param_bound_sec = copy.deepcopy(self.ori_param_bound)
            sorted_unique_dose = np.sort(np.unique(self.dose))
            for sec_ec50 in range(self.num_sec):
                # reset the bounds on ec50
                min_ec50_sec = self.ori_param_bound[0][1] + (self.ori_param_bound[1][1] - self.ori_param_bound[0][1]) * sec_ec50 / self.num_sec
                max_ec50_sec = self.ori_param_bound[0][1] + (self.ori_param_bound[1][1] - self.ori_param_bound[0][1]) * (sec_ec50 + 1) / self.num_sec
                param_bound_sec[0][1] = min_ec50_sec
                param_bound_sec[1][1] = max_ec50_sec

                for sec_e in range(self.num_sec):
                    # reset bounds on einf and e0
                    low_id = np.where(np.isin(self.dose, sorted_unique_dose[:2]))[0]
                    high_id = np.where(np.isin(self.dose, sorted_unique_dose[-2:]))[0]
                    min_einf_sec = np.mean(self.measure[high_id]) + (np.min(self.measure[high_id]) - np.mean(self.measure[high_id])) * (sec_e + 1)
                    max_einf_sec = np.mean(self.measure[high_id]) + (np.max(self.measure[high_id]) - np.mean(self.measure[high_id])) * (sec_e + 1)
                    min_e0_sec = np.mean(self.measure[low_id]) + (np.min(self.measure[low_id]) - np.mean(self.measure[low_id])) * (sec_e + 1)
                    max_e0_sec = np.mean(self.measure[low_id]) + (np.max(self.measure[low_id]) - np.mean(self.measure[low_id])) * (sec_e + 1)
                    # bounds should not exceed the original bounds, but if there is no overlap with the original bound, then set the bounds to the orignal ones.
                    param_bound_sec[0][0] = np.max((min_einf_sec, self.ori_param_bound[0][0]))
                    param_bound_sec[1][0] = np.min((max_einf_sec, self.ori_param_bound[1][0]))
                    if param_bound_sec[0][0] >= param_bound_sec[1][0]:
                        param_bound_sec[0][0] = self.ori_param_bound[0][0]
                        param_bound_sec[1][0] = self.ori_param_bound[1][0]
                    param_bound_sec[0][3] = np.max((min_e0_sec, self.ori_param_bound[0][3]))
                    param_bound_sec[1][3] = np.min((max_e0_sec, self.ori_param_bound[1][3]))
                    if param_bound_sec[0][3] >= param_bound_sec[1][3]:
                        param_bound_sec[0][3] = self.ori_param_bound[0][3]
                        param_bound_sec[1][3] = self.ori_param_bound[1][3]

                    param_sec, param_cov_sec = curve_fit_4p(dose=self.dose, measure=self.measure, param_bound=param_bound_sec)
                    if param_sec is not None:
                        ypred_sec = response_curve_4p(self.dose, *param_sec)
                        r2_sec = r2_score(self.measure, ypred_sec)
                        if self.param is None or r2_sec > r2:
                            self.param_bound = copy.deepcopy(param_bound_sec)
                            self.param = copy.deepcopy(param_sec)
                            self.param_cov = copy.deepcopy(param_cov_sec)
                            r2 = r2_sec

        if self.param is None:
            cols = 'Einf Einf_se EC50 EC50_se HS HS_se E0 E0_se R2fit AUC IC50 AAC1 AUC1'.split(' ')
            self.fit_metrics = pd.Series([np.nan] * len(cols), index=cols)
            print('An experiment is not fitted')
            print(self.compound)
            print(self.specimen)
            print(self.exp_id)
            print('****************************')
            return self.fit_metrics

        # Check whether we need to swtich einf and e0 and take a negative transformation on hs
        sorted_unique_dose = np.sort(np.unique(self.dose))
        if len(sorted_unique_dose) >= 4:
            low_id = np.where(np.isin(self.dose, sorted_unique_dose[:2]))[0]
            high_id = np.where(np.isin(self.dose, sorted_unique_dose[-2:]))[0]
            # Check for reversed pattern of einf and e0
            if (np.mean(self.measure[low_id]) > np.mean(self.measure[high_id]) and self.param[3] < self.param[0]) or \
                    (np.mean(self.measure[low_id]) < np.mean(self.measure[high_id]) and self.param[3] > self.param[0]):
                # Check whether current einf and e0 fall into the original bounds of e0 and einf, respectively.
                if (self.param[0] >= self.ori_param_bound[0][3] and self.param[0] <= self.ori_param_bound[1][3]) and \
                        (self.param[3] >= self.ori_param_bound[0][0] and self.param[3] <= self.ori_param_bound[1][0]):
                    self.param = self.param[[3, 1, 2, 0]]
                    self.param[2] = -self.param[2]
                    self.param_cov = self.param_cov[[3, 1, 2, 0], :]
                    self.param_cov = self.param_cov[:, [3, 1, 2, 0]]
                    self.param_cov[2, :] = -self.param_cov[2, :]
                    self.param_cov[:, 2] = -self.param_cov[:, 2]

        einf = self.param[0]
        ec50 = self.param[1]
        hs = self.param[2]
        e0 = self.param[3]

        perr = np.sqrt(np.diag(self.param_cov))
        einf_se = perr[0]
        ec50_se = perr[1]
        hs_se = perr[2]
        e0_se = perr[3]

        xmin = np.min(self.dose)
        xmax = np.max(self.dose)

        auc1 = self.compute_area(d1=xmin, d2=xmax, mode=mode) / (xmax - xmin)
        aac1 = 1 - auc1

        ic50 = ec50 + np.log10((e0 - 0.5) / (0.5 - einf)) / hs if einf < 0.5 else np.nan

        auc = self.compute_area(d1=d1, d2=d2, mode=mode) / (d2 - d1)

        self.fit_metrics = pd.Series({'Einf':einf, 'Einf_se':einf_se, 'EC50':ec50, 'EC50_se':ec50_se, 'HS':hs,
                                      'HS_se':hs_se, 'E0':e0, 'E0_se':e0_se, 'R2fit':r2, 'AUC':auc, 'IC50':ic50,
                                       'AAC1':aac1, 'AUC1':auc1}).round(4)

        return self.fit_metrics


    def plot_curve(self, show_flag=False, directory='./', save_file_name=None):
        # Call this plot curve function after calling compute_fit_metrics

        dmax = np.max(self.dose)
        dmin = np.min(self.dose)
        dmin = dmin - (dmax - dmin) / 50
        dmax = dmax + (dmax - dmin) / 50
        mmax = np.max(self.measure)
        mmin = np.min(self.measure)
        if self.param is not None:
            xx = np.linspace(dmin, dmax, 100)
            yy = response_curve_4p(xx, *self.param)
            mmax = np.max((mmax, np.max(yy)))
            mmin = np.min((mmin, np.min(yy)))
        mmin = mmin - (mmax - mmin) / 50
        mmax = mmax + (mmax - mmin) / 50

        title = self.specimen + '--' + self.compound + '--' + self.exp_id
        if save_file_name is not None:
            font = {'size': 14}
            matplotlib.rc('font', **font)
            plt.figure(figsize=(12.0, 6.0))
        plt.xlim(dmin, dmax)
        plt.ylim(mmin, mmax)

        if self.param is not None:
            if self.fit_metrics is None:
                plt.plot(xx, yy, 'r-', label='Einf=%.3f, EC50=%.3f, HS=%.3f, E0=%.3f' % tuple(self.param))
            else:
                plt.plot(xx, yy, 'r-', label='Einf=%.3f, EC50=%.3f, HS=%.3f, E0=%.3f, R2=%.3f, AUC=%.3f' % tuple(np.concatenate((self.param, self.fit_metrics[['R2fit', 'AUC']].values))))
        plt.plot(self.dose, self.measure, 'b*')
        plt.xlabel('Dose (log10(M))')
        plt.ylabel(self.feature)
        plt.title(title)
        plt.tight_layout()
        plt.legend()

        if show_flag:
            plt.show()

        if save_file_name is not None:
            # file_name = directory + save_file_name + '.png'
            file_name = os.path.join(directory, save_file_name + '.png') 
            plt.savefig(file_name, dpi=360)
            plt.close()

        return


    def copy(self):
        return copy.deepcopy(self)



def get_colors():

    # if num_color < 1:
    #     return None
    # num_l = int(np.ceil((num_color + 2) ** (1/3)))
    # can = list(range(num_l))
    # pre_colors = [[i] for i in can]
    # for s in range(2):
    #     colors = []
    #     for i in range(len(pre_colors)):
    #         for j in range(len(can)):
    #             colors.append(pre_colors[i] + [can[j]])
    #     pre_colors = colors
    # colors.pop(0)
    # colors.pop(-1)
    # for i in range(len(colors)):
    #     r, g, b = colors[i]
    #     colors[i] = (r / float(num_l - 1), g / float(num_l - 1), b / float(num_l - 1))

    colors = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120), (148, 103, 189), (44, 160, 44),
              (214, 39, 40), (255, 152, 150), (152, 223, 138), (197, 176, 213), (140, 86, 75), (196, 156, 148),
              (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199), (188, 189, 34), (219, 219, 141),
              (23, 190, 207), (158, 218, 229)]
    for i in range(len(colors)):
        r, g, b = colors[i]
        colors[i] = (r / 255., g / 255., b / 255.)

    return colors



def plot_curves(data, show_flag=False, directory='./', save_file_name=None):

    features = []
    dmin = []
    dmax = []
    mmin = []
    mmax = []
    for k in data.keys():
        features.append(k.split('--')[-1])
        dmin.append(np.min(data[k].dose))
        dmax.append(np.max(data[k].dose))
        mmin.append(np.min(data[k].measure))
        mmax.append(np.max(data[k].measure))
    dmin = np.min(dmin)
    dmax = np.max(dmax)
    mmin = np.min(mmin)
    mmax = np.max(mmax)
    dmin = dmin - (dmax - dmin) / 50
    dmax = dmax + (dmax - dmin) / 50
    xx = np.linspace(dmin, dmax, 100)
    for k in data.keys():
        yy = response_curve_4p(xx, *data[k].param)
        mmin = np.min((mmin, np.min(yy)))
        mmax = np.max((mmax, np.max(yy)))
    mmin = mmin - (mmax - mmin) / 50
    mmax = mmax + (mmax - mmin) / 50

    if len(np.unique(features)) > 1:
        print('There are more than 1 features in data.')
        return
    feature = np.unique(features)[0]

    font = {'size': 14}
    matplotlib.rc('font', **font)
    plt.figure(figsize=(12, 6))

    colors = get_colors()
    plt.xlim(dmin, dmax)
    plt.ylim(mmin, mmax)
    plt.xlabel('Dose (log10(M))')
    plt.ylabel(feature)
    rank = 0
    for k in data.keys():
        yy = response_curve_4p(xx, *data[k].param)
        color = colors[rank]
        # rank = rank + 1
        rank = int(np.mod(rank + 1, 20))
        plt.plot(xx, yy, '-', color=color, label=data[k].specimen + '--' + data[k].compound + '--' + data[k].exp_id +
                                                 '--AUC:' + str(round(data[k].fit_metrics['AUC'], 3)))
        plt.plot(data[k].dose, data[k].measure, '.', color=color, label='')
    plt.tight_layout()
    plt.legend()

    if show_flag:
        plt.show()

    if save_file_name is not None:
        file_name = directory + save_file_name + '.png'
        plt.savefig(file_name, dpi=360)
        plt.close()



def split(x):
    return x.split('--')



def generate_dataframe(data, content='fit_metrics'):
    if content == 'fit_metrics':
        dataframe = pd.DataFrame({})
        for k in data.keys():
            dataframe[data[k].study_id + '--' + data[k].specimen + '--' + data[k].compound + '--' + data[k].exp_id +
                      '--' + data[k].feature] = data[k].fit_metrics
        dataframe = dataframe.transpose()
        dataframe['study'] = [i[0] for i in map(split, dataframe.index)]
        dataframe['cell'] = [i[1] for i in map(split, dataframe.index)]
        dataframe['drug'] = [i[2] for i in map(split, dataframe.index)]
        dataframe['exp_id'] = [i[3] for i in map(split, dataframe.index)]
        dataframe['feature'] = [i[4] for i in map(split, dataframe.index)]
        dataframe = dataframe.iloc[:, [-5, -4, -3, -2, -1] + list(range(dataframe.shape[1] - 5))]

    return dataframe




# HS_BOUNDS_ORIG = ([0, 10**-12, 0], [1, 1, 4])
#
# def hs_response_curve_original(x, einf, ec50, hs):
#     """ from PharmacoDB supp. https://doi.org/10.1093/nar/gkx911
#         bounds:
#           einf: [0, 1]       # fraction of cells not susceptible to drug
#           ec50: [10^-12, 1]  # concentration to have half target receptors bound: [1pM, 1M]
#           hs:   [0, 4]       # hill slope binding cooperativity
#     """
#     return einf + (1 - einf) / (1 + np.power(x/ec50, hs))
#
#
# HS_BOUNDS = ([0, 0, 0], [1, 14, 4])
#
# def response_curve(x, einf, ec50, hs):
#     """ transformed the original function with ec50 in -log10(M) instead of M
#     """
#     return einf + (1 - einf) / (1 + 10 ** ((ec50 - x) * hs))
#
#
# def response_integral(x, einf, ec50, hs):
#     return (1 - einf) * np.log10(1 + 10 ** ((ec50 - x) * hs)) / hs + x
#
#
# def compute_area(x1, x2, einf, ec50, hs, mode='trapz'):
#     popt = (einf, ec50, hs)
#     if mode == 'trapz':
#         # trapezoidal numerical integrationcollapse
#         xx = np.linspace(x1, x2, 100)
#         yy = response_curve(xx, *popt)
#         area = np.trapz(yy, xx, dx=0.01)
#     else:
#         # the integral function can be expressed analytically
#         # but sometimes less accurate due to float precision issues
#         area = response_integral(x2, *popt) - response_integral(x1, *popt)
#     return area
#
#
#
# HS_BOUNDS_2 = ([0, -14, 0], [1, 0, 4])
#
# def response_curve_2(x, einf, ec50, hs):
#     """ transformed the original function with ec50 in log10(M) instead of M
#     """
#     return einf + (1 - einf) / (1 + 10 ** ((x - ec50) * hs))
#
# def response_integral_2(x, einf, ec50, hs):
#     return (einf - 1) * np.log10(1 + 10 ** ((x - ec50) * hs)) / hs + x
#
# def compute_area_2(x1, x2, einf, ec50, hs, mode='trapz'):
#     popt = (einf, ec50, hs)
#     if mode == 'trapz':
#         # trapezoidal numerical integrationcollapse
#         xx = np.linspace(x1, x2, 100)
#         yy = response_curve_2(xx, *popt)
#         area = np.trapz(yy, xx, dx=0.01)
#     else:
#         # the integral function can be expressed analytically
#         # but sometimes less accurate due to float precision issues
#         area = response_integral_2(x2, *popt) - response_integral_2(x1, *popt)
#     return area
#
#
#
# HS_4P_BOUNDS_ORIG = ([0, 10**-12, 0, 0.999999], [1, 1, 4, 1.000001])
#
# def hs_curve_4p_original(x, einf, ec50, hs, e0):
#     """ from PharmacoDB supp. https://doi.org/10.1093/nar/gkx911
#         bounds:
#           einf: [0, 1]       # fraction of cells not susceptible to drug
#           ec50: [10^-12, 1]  # concentration to have half target receptors bound: [1pM, 1M]
#           hs:   [0, 4]       # hill slope binding cooperativity
#           e0: []             # Need to define and adjust
#     """
#     return einf + (e0 - einf) / (1 + np.power(x/ec50, hs))
#
# HS_4P_BOUNDS = ([0, -14, -4, 1], [1, 0, 4, 1.0001])
#
#
# def compute_fit_metrics(xdata, ydata, popt, pcov, d1=4, d2=10):
#     if popt is None:
#         cols = 'AUC IC50 EC50 EC50se R2fit Einf HS AAC1 AUC1 DSS1'.split(' ')
#         return pd.Series([np.nan] * len(cols), index=cols)
#
#     einf, ec50, hs = popt
#     perr = np.sqrt(np.diag(pcov))
#     ec50se = perr[1]
#
#     xmin = xdata.min()
#     xmax = xdata.max()
#
#     ypred = response_curve(xdata, *popt)
#     r2 = r2_score(ydata, ypred)
#
#     auc1 = compute_area(xmin, xmax, *popt) / (xmax - xmin)
#     aac1 = 1 - auc1
#
#     ic50 = ec50 - np.log10(0.5/(0.5-einf)) / hs if einf < 0.5 else np.nan
#     ic90 = ec50 - np.log10(0.9/(0.1-einf)) / hs if einf < 0.1 else np.nan
#     ic10 = ec50 - np.log10(0.1/(0.9-einf)) / hs if einf < 0.9 else np.nan
#
#     ic10x = min(ic10, xmax)
#     int10x = compute_area(xmin, ic10x, *popt)
#     dss1 = (0.9 * (ic10x - xmin) - int10x) / (0.9 * (xmax - xmin)) if xmin < ic10x else 0
#     auc = (response_integral(d2, *popt) - response_integral(d1, *popt)) / (d2 - d1)
#
#     metrics = pd.Series({'AUC':auc, 'IC50':ic50, 'EC50':ec50,
#                          'EC50se':ec50se, 'R2fit':r2, 'Einf':einf, 'HS':hs,
#                          'AAC1':aac1, 'AUC1':auc1, 'DSS1':dss1}).round(4)
#
#     return metrics


# def response_curve_fit(xdata, ydata, bounds=HS_BOUNDS):
#     ydata = ydata.clip(lower=0, upper=1.0)
#     popt, pcov = None, None
#     nfev = 100 * 3
#     while popt is None and nfev < 10000:
#         # print(nfev)
#         try:
#             popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev)
#             # popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev, method='dogbox')
#         except RuntimeError:
#             pass
#         nfev *= 2
#     return popt, pcov


# def response_curve_fit(xdata, ydata, fit_function, flag_yclip=True, bounds=HS_BOUNDS):
#     if flag_yclip:
#         ydata = ydata.clip(lower=0, upper=1.0)
#     popt, pcov = None, None
#     nfev = 100 * 3
#     while popt is None and nfev < 10000:
#         # print(nfev)
#         try:
#             popt, pcov = curve_fit(fit_function, xdata, ydata, bounds=bounds, max_nfev=nfev)
#             # popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev, method='dogbox')
#         except RuntimeError:
#             pass
#         nfev *= 2
#     # print(popt)
#     # print(pcov)
#     return popt, pcov
#
#
#
#
# def fit_exp(df_exp, title=None, dmin=None, dmax=None, save=False):
#     if save:
#         font = {'family' : 'normal',
#                 # 'weight' : 'bold',
#                 'size'   : 14}
#         matplotlib.rc('font', **font)
#         plt.figure(figsize=(12, 6))
#
#     print(df_exp)
#     xdata = df_exp.DOSE.astype(np.float)
#     ydata = df_exp.GROWTH.astype(np.float)
#     # ydata = df_exp.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
#
#     # print(xdata)
#     # print(ydata)
#
#     popt, pcov = response_curve_fit(xdata, ydata)
#     metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
#
#     if popt is None:
#         return metrics
#
#     dmin = dmin or xdata.min()
#     dmax = dmax or xdata.max()
#     xx = np.linspace(dmin, dmax, 100)
#     yy = response_curve(xx, *popt)
#
#     plt.xlim(dmax, dmin)
#     plt.ylim(0, np.max([105, np.max(yy)]))
#     plt.plot(xx, yy*100, 'r-', label='fit: Einf=%.3f, EC50=%.3f, HS=%.3f' % tuple(popt))
#     plt.plot(xdata, ydata.clip(lower=0, upper=1.0)*100, 'b*', label='')
#     plt.xlabel('Dose (-log10(M))')
#     plt.ylabel('Growth%')
#     plt.title(title)
#     plt.tight_layout()
#     plt.legend()
#     if save:
#         plt.savefig('exp.png', dpi=360)
#         plt.close()
#     else:
#         plt.show()
#
#     return metrics.to_frame(name='metrics').T
#
#
# def fit_response(df_all, cell, drug, source, study=None, save=False):
#     cell_ids = ud.cell_name_to_ids(cell) or [cell]
#     drug_ids = ud.drug_name_to_ids(drug) or [drug]
#
#     df_exp = df_all[df_all.CELL.isin(cell_ids) & df_all.DRUG.isin(drug_ids)].copy()
#     df_exp.GROWTH = (df_exp.GROWTH/2 + 0.5)
#     df_exp = df_exp[df_exp.SOURCE == source]
#
#     title = f'{cell} treated with {drug} in {source}'
#
#     studies = df_exp.STUDY.unique()
#     if len(studies) > 1:
#         study = studies[study] if type(study) == int else study or studies[0]
#         title += f' study {study}'
#         df_exp = df_exp[df_exp.STUDY == study]
#
#     return fit_exp(df_exp, title, save=save)
#
#
# def show_dose_distribution(df_all):
#     sources = df_all.SOURCE.unique()
#     qs = [0, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 1]
#     series = []
#     for src in sources:
#         s = df_all[df_all.SOURCE == src].DOSE.quantile(qs)
#         s.name = src
#         series.append(s)
#     df_dose = pd.concat(series, axis=1)
#     return df_dose
#
#
# def process_df(df, fname, sep='\t', ngroups=None):
#     # df = df1.copy()
#     i = 0
#     header = None
#     cols = ['SOURCE', 'CELL', 'DRUG', 'STUDY']
#     groups = df.groupby(cols)
#     f = open(fname, 'w')
#     for name, group in tqdm(groups):
#         # print(name)
#         xdata = group.DOSE.astype(np.float)
#         ydata = group.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
#         popt, pcov = response_curve_fit(xdata, ydata)
#         metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
#         if header is None:
#             header = cols + metrics.index.tolist()
#             print(sep.join(header), file=f)
#         print(sep.join(name), end=sep, file=f)
#         print(sep.join([f'{x:.4g}' for x in metrics]), file=f)
#         i += 1
#         if ngroups and i >= ngroups:
#             break
#     f.close()
#
#
# def process_df_part(df, fname, sep='\t', start=0, count=None):
#     header = None
#     cols = ['SOURCE', 'CELL', 'DRUG', 'STUDY']
#     groups = df.groupby(cols)
#     # count = count or (len(groups) - start)
#     count = count or (4484081 - start)
#     groups = islice(groups, start, start+count)
#     f = open(f'{fname}.{start}', 'w')
#     for name, group in tqdm(groups):
#         # print(name)
#         xdata = group.DOSE.astype(np.float)
#         ydata = group.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
#         popt, pcov = response_curve_fit(xdata, ydata)
#         metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
#         if start == 0 and header is None:
#             header = cols + metrics.index.tolist()
#             print(sep.join(header), file=f)
#         print(sep.join(name), end=sep, file=f)
#         print(sep.join([f'{x:.4g}' for x in metrics]), file=f)
#     f.close()
#
#
# def test():
#     df0 = ud.load_single_dose_response(fraction=True)
#
#     cell_ids = ud.cell_name_to_ids('LOXIMVI')
#     drug_ids = ud.drug_name_to_ids('paclitaxel')
#
#     df1 = df0[df0.CELL.isin(cell_ids) & df0.DRUG.isin(drug_ids)].copy()
#     df1.GROWTH = df1.GROWTH/2 + 0.5
#     df2 = df1[df1.SOURCE == 'NCI60']
#
#     fit_exp(df2)
#
#
# def process_chem_partner_data():
#     df_cp = pd.read_csv('curve/ChemPartner_dose_response', sep='\t')
#     df_cp = df_cp[df_cp.DRUG2.isnull() & df_cp.DOSE2.isnull()].drop(['DRUG2', 'DOSE2'], axis=1)
#     df_cp = df_cp.rename(columns={'DRUG1':'DRUG', 'DOSE1':'DOSE'})
#     df_cp.DOSE = -df_cp.DOSE
#     # df_cp.GROWTH = df_cp.GROWTH/100
#     df_cp.GROWTH = df_cp.GROWTH/200 + 0.5
#
#     # process_df(df_cp, 'curve/ChemPartner_single_response_agg', ngroups=10)
#     process_df(df_cp, 'curve/ChemPartner_single_response_agg.new')
#
#
# def fix_auc_gt_one():
#     dfx = pd.read_table('curve/combined_single_response_agg.0', engine='c', low_memory=False)
#
#
# def notebook():
#     d = pd.read_csv('./curve/combined_single_response_agg', engine='c', sep='\t', low_memory=False)
#     m = 1e-3
#     print(d[(d.AUC < m) & (d.R2fit < m) & (d.EC50se < m)].shape)
#     print(d[(d.AUC < m) & (d.R2fit < m) & (d.EC50se < m)].head())
#
#
# def fit_exp(df_exp, title=None, dmin=None, dmax=None, save=False):
#     if save:
#         font = {'family' : 'normal',
#                 # 'weight' : 'bold',
#                 'size'   : 14}
#         matplotlib.rc('font', **font)
#         plt.figure(figsize=(12, 6))
#
#     print(df_exp)
#     xdata = df_exp.DOSE.astype(np.float)
#     ydata = df_exp.GROWTH.astype(np.float)
#     # ydata = df_exp.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
#
#     # print(xdata)
#     # print(ydata)
#
#     popt, pcov = response_curve_fit(xdata, ydata)
#     metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
#
#     if popt is None:
#         return metrics
#
#     dmin = dmin or xdata.min()
#     dmax = dmax or xdata.max()
#     xx = np.linspace(dmin, dmax, 100)
#     yy = response_curve(xx, *popt)
#
#     plt.xlim(dmax, dmin)
#     plt.ylim(0, np.max([105, np.max(yy)]))
#     plt.plot(xx, yy*100, 'r-', label='fit: Einf=%.3f, EC50=%.3f, HS=%.3f' % tuple(popt))
#     plt.plot(xdata, ydata.clip(lower=0, upper=1.0)*100, 'b*', label='')
#     plt.xlabel('Dose (-log10(M))')
#     plt.ylabel('Growth%')
#     plt.title(title)
#     plt.tight_layout()
#     plt.legend()
#     if save:
#         plt.savefig('exp.png', dpi=360)
#         plt.close()
#     else:
#         plt.show()
#
#     return metrics.to_frame(name='metrics').T
#
#
# def get_tableau20_colors():
#     # tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
#     #          (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),
#     #          (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),
#     #          (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
#     #          (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
#     tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),
#                  (148, 103, 189), (44, 160, 44), (214, 39, 40), (255, 152, 150),
#                  (152, 223, 138), (197, 176, 213), (140, 86, 75), (196, 156, 148),
#                  (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),
#                  (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]
#     # Scale the RGB values to the [0, 1] range, which is the format matplotlib accepts.
#     for i in range(len(tableau20)):
#         r, g, b = tableau20[i]
#         tableau20[i] = (r / 255., g / 255., b / 255.)
#     return tableau20
#
#
# def plot_curves(df_all, cell='LOXIMVI', drug='paclitaxel', study=None, max_reps=2, dmin=4, dmax=10, out=None):
#     cell_ids = ud.cell_name_to_ids(cell)
#     drug_ids = ud.drug_name_to_ids(drug)
#
#     df_exps = df_all[df_all.CELL.isin(cell_ids) & df_all.DRUG.isin(drug_ids)].copy()
#     df_exps.GROWTH = (df_exps.GROWTH/2 + 0.5)
#
#     title = f'{cell} treated with {drug}'
#     out = out or f'{cell}-{drug}'
#
#     # font = {'family': 'normal', 'size': 14}
#     font = {'size': 14}
#     matplotlib.rc('font', **font)
#     plt.figure(figsize=(12, 6))
#     colors = get_tableau20_colors()
#
#     dmin = dmin or df_exps.DOSE.min()
#     dmax = dmax or df_exps.DOSE.max()
#     xx = np.linspace(dmin-0.1, dmax+0.1, 100)
#
#     plt.xlim(dmax+0.1, dmin-0.1)
#     plt.ylim(0, 105)
#     plt.xlabel('Dose (-log10(M))')
#     plt.ylabel('Growth%')
#     plt.title(title)
#
#     df_metrics = None
#     rank = 0
#     order = ['NCI60', 'CTRP', 'GDSC', 'CCLE', 'gCSI']
#     sources = df_exps.SOURCE.unique().tolist() if study is None else study
#     sources = sorted(sources, key=lambda x:order.index(x))
#
#     for source in sources:
#         studies = df_exps[df_exps.SOURCE == source].STUDY.unique()
#         for i, study in enumerate(studies[:max_reps]):
#             df_exp = df_exps[(df_exps.SOURCE == source) & (df_exps.STUDY == study)]
#             xdata = df_exp.DOSE.astype(np.float)
#             ydata = df_exp.GROWTH.astype(np.float)
#             # ydata = df_exp.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
#             popt, pcov = response_curve_fit(xdata, ydata)
#             metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
#             if popt is None:
#                 continue
#             color = colors[rank]
#             rank = (rank + 1) % 20
#             yy = response_curve(xx, *popt)
#             label = source
#             if len(studies) > 1:
#                 label += f' rep {i+1}'
#             plt.plot(xx, yy*100, '-', color=color, label=label)
#             plt.plot(xdata, ydata.clip(lower=0, upper=1.0)*100, '.', color=color, label='')
#             if df_metrics is None:
#                 df_metrics = metrics.to_frame(name=label).T
#             else:
#                 df_metrics = pd.concat([df_metrics, metrics.to_frame(name=label).T])
#
#     plt.tight_layout()
#     plt.legend()
#     plt.savefig(f'{out}.png', dpi=360)
#     plt.close()
#
#     df_metrics.index.name = 'Source'
#     df_metrics.to_csv(f'{out}.csv', float_format='%.5g')
#     print(f'Saved {out}.png and {out}.csv.')
#
#     return df_metrics


# def main():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('-c', '--cell', help='cell line name')
#     parser.add_argument('-d', '--drug', help='drug name')
#     parser.add_argument('-s', '--study', nargs='+', help='list of data sources/studies')
#     parser.add_argument('--reps', type=int, default=2, help='maximum of replicates to include in plots')
#     parser.add_argument('--start', type=int, default=0, help='start index (0-based)')
#     parser.add_argument('--count', type=int, help='number of experiments')
#     parser.add_argument('--out', help='prefix of output file')
#
#     args = parser.parse_args()
#
#     df_all = ud.load_single_dose_response(fraction=True)
#
#     if args.cell and args.drug:
#         plot_curves(df_all, cell=args.cell, drug=args.drug, study=args.study, max_reps=args.reps)
#     else:
#         fname = out or 'combined_single_response_agg'
#         process_df_part(df_all, fname, start=args.start, count=args.count)
#         # 4484081
#         # process_df_part(df_all, 'curve/combined_single_response_agg', start=int(sys.argv[1]), count=int(sys.argv[2]))
#         # process_df(df_all, 'curve/combined_single_response_agg')
#         # process_df_part(df_all, 'curve/combined_single_response_agg.debug', start=0, count=270)
#
#
# if __name__ == '__main__':
#     main()
