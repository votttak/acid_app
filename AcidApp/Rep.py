from numpy.core.function_base import _linspace_dispatcher
import pandas as pd
import sqlite3
import csv
import numpy as np
import yattag
from tkinter import *
from pandas.core.frame import DataFrame

'''
SUBSTRATE_SEARCH_1 = None
SUBSTRATE_SEARCH_2 = None
SUBSTRATE_SEARCH_3 = None
SUBSTRATE_SEARCH_4 = None
SUBSTRATE_SEARCH_5 = None
'''



def clicked():
    title=""
    code, protease_1, protease_2, freq_table1, freq_table2, peptide_table1, peptide_table2, peptide_combination1, peptide_combination2, common_peptides, unique_peptides_1, unique_peptides_2 = execute(txtp1main.get(), txtp1pos.get(), txtp2pos.get(), txtp3pos.get(), txtp4pos.get(), txtp1pos1.get(),txtp2pos2.get(), txtp3pos3.get(), txtp4pos4.get(), txtp1ph.get(), txt2p1main.get(), txt2p1pos.get(), txt2p2pos.get(), txt2p3pos.get(), txt2p4pos.get(), txt2p1pos1.get(), txt2p2pos2.get(), txt2p3pos3.get(), txt2p4pos4.get(), txt2p1ph.get(), txt3p1main.get(), txt3p1pos.get(), txt3p2pos.get(), txt3p3pos.get(), txt3p4pos.get(), txt3p1pos1.get(), txt3p2pos2.get(), txt3p3pos3.get(), txt3p4pos4.get(), txt3p1ph.get(), txt4p1main.get(), txt4p1pos.get(), txt4p2pos.get(), txt4p3pos.get(), txt4p4pos.get(), txt4p1pos1.get(), txt4p2pos2.get(), txt4p3pos3.get(), txt4p4pos4.get(), txt4p1ph.get(), txt5p1main.get(), txt5p1pos.get(), txt5p2pos.get(), txt5p3pos.get(), txt5p4pos.get(), txt5p1pos1.get(), txt5p2pos2.get(), txt5p3pos3.get(), txt5p4pos4.get(), txt5p1ph.get())
    #html_str = "<div style=\"float: right;\"><br><br><br> PEPTIDE_TABLE " + protease_2 + peptide_table2.to_html().replace('nan', '0') + peptide_table2.to_html().replace('nan', '0')+freq_table2.to_html().replace('NaN', '0') +freq_table1.to_html().replace('NaN', '0')+peptide_combination2.to_html()+peptide_combination1.to_html()+common_peptides.to_html()+unique_peptides_2.to_html()+unique_peptides_1.to_html()
    html_str1 = "<div style=\"float: right;\"><br><br><br> PEPTIDE_TABLE " + protease_2 + peptide_table2.to_html().replace('nan', '0') + "</div>" + "\n\n<br><br><br>PEPTIDE_TABLE " + protease_1 + peptide_table1.to_html().replace('nan', '0') + "<div style=\"float: right;\"><br><br><br> FREQUENCY TABLE " + protease_2 + freq_table2.to_html().replace('NaN', '0') + "</div>" + "\n\n<br><br><br> FREQUENCY TABLE " + protease_1 + freq_table1.to_html().replace('NaN', '0') + "<div style=\"float: right;\"><br><br><br> PEPTIDE COMBINATIONS " + protease_2 + " <br>(made out of two best peptides based on frequency for each position)\n" + peptide_combination2.to_html() + "</div>" + "\n\n<br><br><br> PEPTIDE COMBINATIONS " + protease_1 + " <br>(made out of two best peptides based on frequency for each position)\n" + peptide_combination1.to_html() + "<center><br><br><br> COMMON PEPTIDES for " + title + "</center>" + "<center>" + common_peptides.to_html() + "</center>" + "<div style=\"float: right;\"><br><br><br> UNIQUE PEPTIDES for " + protease_2  + unique_peptides_2.to_html() + "</div>" + "\n\n<br><br><br> UNIQUE PEPTIDES for " + protease_1  + unique_peptides_1.to_html()
    Html_file = open('result.html', "w")
    Html_file.write(html_str1)
    Html_file.close()

def create_file_csv(mycursor, table_name, protease_number):
    header = [row[0] for row in mycursor.description]
    rows = mycursor.fetchall()
    f = open(table_name + '.csv', 'w')

    # Write header
    f.write(','.join(header) + '\n')

    for row in rows:
        f.write(','.join(str(r) for r in row) + '\n')
    f.close()

    file = pd.read_csv(table_name + '.csv')

    if protease_number == 1:
        SUBSTRATE_SEARCH_1 = file
    if protease_number == 2:
        SUBSTRATE_SEARCH_2 = file
    if protease_number == 3:
        SUBSTRATE_SEARCH_3 = file
    if protease_number == 4:
        SUBSTRATE_SEARCH_4 = file
    if protease_number == 5:
        SUBSTRATE_SEARCH_5 = file

    return file


def fetch_table_data(protease, cleavage_type_to_exclude, protease_number):
    # conn=sqlite3.connect("sqlDB.db")
    sqlite_connection = sqlite3.connect('sqlDB.db')
    cursor = sqlite_connection.cursor()
    # print("Подключен к SQLite")

    # sqlite_select_query = """SELECT * from sqlitedb_developers"""
    # cursor.execute(sqlite_select_query)
    # records = cursor.fetchall()
    columns = 'Site_P4, Site_P3, Site_P2, Site_P1, Site_P1prime, Site_P2prime, Site_P3prime, Site_P4prime, Protease, cleavage_type '
    table_name = 'Tab1'
    # limit = '20'
    cursor.execute(
        "select " + columns + "from " + table_name + " where Protease = " + "'" + protease + "'" + " AND cleavage_type NOT IN('synthetic') ")

    # create csv file
    return create_file_csv(cursor, protease + "_data", protease_number)


def filter_array(arr, el_to_exclude):
    # Create an empty list
    filter_arr = []

    # go through each element in arr
    for element in arr:
        if element[0:3] != el_to_exclude:
            filter_arr.append(True)
        else:
            filter_arr.append(False)

    newarr = arr[filter_arr]
    return newarr


### the most frequent amino acids for P4-P4' positions ###
def extract_peptide(substrate_search, amino_acid_to_exclude_at_particular_positions=[], frequency_table_name='freq.csv',
                    protease_number=1):
    df = substrate_search

    df_cut = df[
        ['Site_P4', 'Site_P3', 'Site_P2', 'Site_P1', 'Site_P1prime', 'Site_P2prime', 'Site_P3prime', 'Site_P4prime']]
    freq_table = df_cut.apply(pd.value_counts)
    freq_table.to_csv(frequency_table_name, index=True)
    cols = freq_table.columns

    df_peptide = pd.DataFrame()
    number_to_div_by = freq_table.sum(axis=0, skipna=True)[0]

    for i, col in enumerate(cols):
        x = freq_table.sort_values(by=[col], ascending=False)
        x = x[col]

        y = (x.index).to_numpy()  # ['-', 'Ala', ... ]
        y_old = (x.index).to_numpy()
        if ('-' in x.index):
            x = x.drop(['-'], axis=0)
        percent = round((x / (x.sum(axis=0, skipna=True))) * 100, 1)
        percent = percent.apply(str)
        y = filter_array(y, '-')

        for j in range(len(percent)):
            y[j] = y[j] + " (" + percent[y[j]] + "%) "

        if len(amino_acid_to_exclude_at_particular_positions) != 0:
            tuple = amino_acid_to_exclude_at_particular_positions[i]
            for amino_acid in tuple[1]:
                y = filter_array(y, amino_acid)

        df_peptide[col] = pd.Series(y)

    percentage_table = freq_table / number_to_div_by
    return (freq_table, df_peptide)


# def peptide_combinations(table, res = []):
#	if len(res) == 8:
#		print("HERE: ", res)
#	else:
#		#df.loc[row_indexer,column_indexer]
#		res_new = res.copy()
#		res.append(table.loc[0, table.columns[len(res)]])
#		res_new.append(table.loc[1, table.columns[len(res)]])
#		peptide_combinations(table, res)
#		peptide_combinations(table, res_new)


def parse_str(input_str):
    first_part = input_str.split("%)")[0]  # ALA (12.32
    result = first_part.split('(')[1]
    return float(result)


def parse_str_aa(input_str):
    first_part = input_str.split("%)")[0]  # ALA (12.32
    result = first_part.split('(')[0]
    return result


def asd(elem, f_list, s_list):
    return [elem, ]


'''
d = {
   	         'aa' : ['Arg', 'Asp', 'Cys', 'Glu', 'His', 'Lys', 'Tyr'], 
                'pKa' : [12.48, 3.86, 8.33, 4.25, 6.0, 10.53, 10.07],                                # (COOH) 
            }
PKS_TABLE = (pd.DataFrame(data=d)).set_index(['aa'])
d = {
		'aa' : ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg', 'Ser', 'Thr', 'Val', 'Trp', 'Tyr'],
	    'hyd_val' : [-0.5, -1.0, 3.0, 3.0, -2.5, 0.0, -0.5, -1.8, 3.0, -1.8, -1.3, 0.2, 0.0, 0.2, 3.0, 0.3, -0.4, -1.5, -3.4, -2.3]
	    }
HYDROPHILICITY_VALUE_TABLE = (pd.DataFrame(data=d)).set_index(['aa'])

pH1 = 0.0
pH2 = 0.0

peptide_combs_with_letters_1 = []
peptide_combs_with_letters_2 = []


def create_peptide_combinations(f, s, t, e, ph, protease1):
	result = None
	if len(f) == 1:
		result = [[f[0]], [s[0]], [t[0]], [e[0]]]
		print('876543')
		print(result)
	else:
		full_result = []
		temp_res = create_peptide_combinations(f[1:], s[1:], t[1:], e[1:], ph)
		for li in temp_res:
			if not(f[0] is None):
				peptide_list = [f[0]] + li
				#print('543210')
				#print(peptide_list)
				if len(peptide_list) == 8:
					peptide = convert_aa_to_capitals(peptide_list)
					if protease1:
						peptide_combs_with_letters_1.append(peptide)
					else: 
						peptide_combs_with_letters_2.append(peptide)
					pI = predict_isoelectric_point(peptide)
					netCharge = net_charge([f[0]] + li, PKS_TABLE, ph)
					hydr_resid_DIV_tot_n_of_resid = compute_hydrophilic_residues_DIV_tot_num_of_residues(peptide_list, HYDROPHILICITY_VALUE_TABLE) 
					entry = peptide_list + [pI] + [netCharge] + [hydr_resid_DIV_tot_n_of_resid]
					print('654321')
					print(entry)		
					full_result.append(entry)
				else:
					full_result.append(peptide_list)
			if not(s[0] is None):
				peptide_list = [s[0]] + li
				#print('543210')
				#print(peptide_list)
				if len(peptide_list) == 8:
					peptide = convert_aa_to_capitals(peptide_list)
					if protease1:
						peptide_combs_with_letters_1.append(peptide)
					else: 
						peptide_combs_with_letters_2.append(peptide)
					pI = predict_isoelectric_point(peptide)
					netCharge = net_charge([f[0]] + li, PKS_TABLE, ph)
					hydr_resid_DIV_tot_n_of_resid = compute_hydrophilic_residues_DIV_tot_num_of_residues(peptide_list, HYDROPHILICITY_VALUE_TABLE) 
					entry = peptide_list + [pI] + [netCharge] + [hydr_resid_DIV_tot_n_of_resid]
					print('654321')
					print(entry)		
					full_result.append(entry)
				else:
					full_result.append(peptide_list)
			if not(t[0] is None):
				peptide_list = [t[0]] + li
				#print('543210')
				#print(peptide_list)
				if len(peptide_list) == 8:
					peptide = convert_aa_to_capitals(peptide_list)
					if protease1:
						peptide_combs_with_letters_1.append(peptide)
					else: 
						peptide_combs_with_letters_2.append(peptide)
					pI = predict_isoelectric_point(peptide)
					netCharge = net_charge([f[0]] + li, PKS_TABLE, ph)
					hydr_resid_DIV_tot_n_of_resid = compute_hydrophilic_residues_DIV_tot_num_of_residues(peptide_list, HYDROPHILICITY_VALUE_TABLE) 
					entry = peptide_list + [pI] + [netCharge] + [hydr_resid_DIV_tot_n_of_resid]
					print('654321')
					print(entry)		
					full_result.append(entry)
				else:
					full_result.append(peptide_list)
			if not(e[0] is None):
				peptide_list = [e[0]] + li
				#print('543210')
				#print(peptide_list)
				if len(peptide_list) == 8:
					peptide = convert_aa_to_capitals(peptide_list)
					if protease1:
						peptide_combs_with_letters_1.append(peptide)
					else: 
						peptide_combs_with_letters_2.append(peptide)
					pI = predict_isoelectric_point(peptide)
					netCharge = net_charge([f[0]] + li, PKS_TABLE, ph)
					hydr_resid_DIV_tot_n_of_resid = compute_hydrophilic_residues_DIV_tot_num_of_residues(peptide_list, HYDROPHILICITY_VALUE_TABLE) 
					entry = peptide_list + [pI] + [netCharge] + [hydr_resid_DIV_tot_n_of_resid]
					print('654321')
					print(entry)		
					full_result.append(entry)
				else:
					full_result.append(peptide_list)
		result = full_result
		#print('987654')
		#print(result)
	return result
'''


def create_peptide_combinations(f, s, t, e):
    result = None

    if len(f) == 1:
        result = [[f[0]], [s[0]], [t[0]], [e[0]]]
    else:
        full_result = []
        full_result_rejected_peptides = []
        temp_res = create_peptide_combinations(f[1:], s[1:], t[1:], e[1:])
        for li in temp_res:
            if not (f[0] is None):
                full_result.append([f[0]] + li)
            if not (s[0] is None):
                full_result.append([s[0]] + li)
            if not (t[0] is None):
                full_result.append([t[0]] + li)
            if not (e[0] is None):
                full_result.append([e[0]] + li)
        result = full_result
    return result


scales = {
    "EMBOSS": {'Cterm': 3.6, 'pKAsp': 3.9, 'pKGlu': 4.1, 'pKCys': 8.5, 'pKTyr': 10.1, 'pk_his': 6.5, 'Nterm': 8.6,
               'pKLys': 10.8, 'pKArg': 12.5},
    "DTASelect": {'Cterm': 3.1, 'pKAsp': 4.4, 'pKGlu': 4.4, 'pKCys': 8.5, 'pKTyr': 10.0, 'pk_his': 6.5, 'Nterm': 8.0,
                  'pKLys': 10.0, 'pKArg': 12.0},
    "Solomon": {'Cterm': 2.4, 'pKAsp': 3.9, 'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 9.6,
                'pKLys': 10.5, 'pKArg': 12.5},
    "Sillero": {'Cterm': 3.2, 'pKAsp': 4.0, 'pKGlu': 4.5, 'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 6.4, 'Nterm': 8.2,
                'pKLys': 10.4, 'pKArg': 12.0},
    "Rodwell": {'Cterm': 3.1, 'pKAsp': 3.68, 'pKGlu': 4.25, 'pKCys': 8.33, 'pKTyr': 10.07, 'pk_his': 6.0, 'Nterm': 8.0,
                'pKLys': 11.5, 'pKArg': 11.5},
    "Patrickios": {'Cterm': 4.2, 'pKAsp': 4.2, 'pKGlu': 4.2, 'pKCys': 0.0, 'pKTyr': 0.0, 'pk_his': 0.0, 'Nterm': 11.2,
                   'pKLys': 11.2, 'pKArg': 11.2},
    "Wikipedia": {'Cterm': 3.65, 'pKAsp': 3.9, 'pKGlu': 4.07, 'pKCys': 8.18, 'pKTyr': 10.46, 'pk_his': 6.04,
                  'Nterm': 8.2, 'pKLys': 10.54, 'pKArg': 12.48},
    "Grimsley": {'Cterm': 3.3, 'pKAsp': 3.5, 'pKGlu': 4.2, 'pKCys': 6.8, 'pKTyr': 10.3, 'pk_his': 6.6, 'Nterm': 7.7,
                 'pKLys': 10.5, 'pKArg': 12.04},
    'Lehninger': {'Cterm': 2.34, 'pKAsp': 3.86, 'pKGlu': 4.25, 'pKCys': 8.33, 'pKTyr': 10.0, 'pk_his': 6.0,
                  'Nterm': 9.69, 'pKLys': 10.5, 'pKArg': 12.4},
    'Bjellqvist': {'Cterm': 3.55, 'pKAsp': 4.05, 'pKGlu': 4.45, 'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 5.98,
                   'Nterm': 7.5, 'pKLys': 10.0, 'pKArg': 12.0},
    'IPC_peptide': {'Cterm': 2.383, 'pKAsp': 3.887, 'pKGlu': 4.317, 'pKCys': 8.297, 'pKTyr': 10.071, 'pk_his': 6.018,
                    'Nterm': 9.564, 'pKLys': 10.517, 'pKArg': 12.503},  # IPC peptide
    'IPC_protein': {'Cterm': 2.869, 'pKAsp': 3.872, 'pKGlu': 4.412, 'pKCys': 7.555, 'pKTyr': 10.85, 'pk_his': 5.637,
                    'Nterm': 9.094, 'pKLys': 9.052, 'pKArg': 11.84},  # IPC protein
    'Toseland': {'Cterm': 3.19, 'pKAsp': 3.6, 'pKGlu': 4.29, 'pKCys': 6.87, 'pKTyr': 9.61, 'pk_his': 6.33,
                 'Nterm': 8.71, 'pKLys': 10.45, 'pKArg': 12},
    'Thurlkill': {'Cterm': 3.67, 'pKAsp': 3.67, 'pKGlu': 4.25, 'pKCys': 8.55, 'pKTyr': 9.84, 'pk_his': 6.54,
                  'Nterm': 8.0, 'pKLys': 10.4, 'pKArg': 12.0},
    'Nozaki': {'Cterm': 3.8, 'pKAsp': 4.0, 'pKGlu': 4.4, 'pKCys': 9.5, 'pKTyr': 9.6, 'pk_his': 6.3, 'Nterm': 7.5,
               'pKLys': 10.4, 'pKArg': 12},
    'Dawson': {'Cterm': 3.2, 'pKAsp': 3.9, 'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 8.2,
               'pKLys': 10.5, 'pKArg': 12},
}

aaDict = {'Asp': 'D', 'Glu': 'E', 'Cys': 'C', 'Tyr': 'Y', 'His': 'H',
          'Lys': 'K', 'Arg': 'R', 'Met': 'M', 'Phe': 'F', 'Leu': 'L',
          'Val': 'V', 'Ala': 'A', 'Gly': 'G', 'Gln': 'Q', 'Asn': 'N',
          'Ile': 'I', 'Trp': 'W', 'Ser': 'S', 'Thr': 'T', 'Sec': 'U',
          'Pro': 'P', 'Xaa': 'X', 'Sec': 'U', 'Pyl': 'O', 'Asx': 'B',
          'Xle': 'J', }


def predict_isoelectric_point(seq, scale='IPC_peptide'):
    """calculate isoelectric point using 9 pKa set model"""
    pKCterm = scales[scale]['Cterm']
    pKAsp = scales[scale]['pKAsp']
    pKGlu = scales[scale]['pKGlu']
    pKCys = scales[scale]['pKCys']
    pKTyr = scales[scale]['pKTyr']
    pKHis = scales[scale]['pk_his']
    pKNterm = scales[scale]['Nterm']
    pKLys = scales[scale]['pKLys']
    pKArg = scales[scale]['pKArg']
    pH = 6.51  # starting po pI = 6.5 - theoretically it should be 7, but average protein pI is 6.5 so we increase the probability of finding the solution
    pHprev = 0.0
    pHnext = 14.0
    E = 0.01  # epsilon means precision [pI = pH +- E]
    temp = 0.01
    nterm = seq[0]
    if scale == 'Bjellqvist':
        if nterm in pKnterminal.keys():
            pKNterm = pKnterminal[nterm]

    cterm = seq[-1]
    if scale == 'Bjellqvist':
        if cterm in pKcterminal.keys():
            pKCterm = pKcterminal[cterm]

    while 1:  # the infinite loop
        QN1 = -1.0 / (1.0 + pow(10, (pKCterm - pH)))
        QN2 = -seq.count('D') / (1.0 + pow(10, (pKAsp - pH)))
        QN3 = -seq.count('E') / (1.0 + pow(10, (pKGlu - pH)))
        QN4 = -seq.count('C') / (1.0 + pow(10, (pKCys - pH)))
        QN5 = -seq.count('Y') / (1.0 + pow(10, (pKTyr - pH)))
        QP1 = seq.count('H') / (1.0 + pow(10, (pH - pKHis)))
        QP2 = 1.0 / (1.0 + pow(10, (pH - pKNterm)))
        QP3 = seq.count('K') / (1.0 + pow(10, (pH - pKLys)))
        QP4 = seq.count('R') / (1.0 + pow(10, (pH - pKArg)))
        NQ = QN1 + QN2 + QN3 + QN4 + QN5 + QP1 + QP2 + QP3 + QP4
        # print NQ
        # %%%%%%%%%%%%%%%%%%%%%%%%%   BISECTION   %%%%%%%%%%%%%%%%%%%%%%%%
        if NQ < 0.0:  # we are out of range, thus the new pH value must be smaller
            temp = pH
            pH = pH - ((pH - pHprev) / 2.0)
            pHnext = temp
            # print "pH: ", pH, ", \tpHnext: ",pHnext
        else:
            temp = pH
            pH = pH + ((pHnext - pH) / 2.0)
            pHprev = temp
            # print "pH: ", pH, ",\tpHprev: ", pHprev

        if (pH - pHprev < E) and (pHnext - pH < E):  # terminal condition, finding pI with given precision
            return pH


def convert_aa_to_capitals(peptide_list):
    peptide = ''
    num_aas = len(peptide_list)
    for i in range(num_aas):
        aa = aaDict[(peptide_list[i])[0:3]]
        peptide = peptide + aa
    return peptide


def net_charge(peptide_list, pks_table, pH):
    net_charge_val = 0
    Ni = 0
    Nj = 0
    z_left = 0
    z_right = 0
    for i in range(8):
        aa = (peptide_list[i])[0:3]
        if aa in ['Arg', 'Lys', 'His']:
            Ni += 1
            pka = float(pks_table.loc[aa, ['pKa']])
            z_left += pow(10, pka) / (pow(10, pH) + pow(10, pka))
        if aa in ['Asp', 'Glu', 'Cys', 'Tyr']:
            Nj += 1
            pka = float(pks_table.loc[aa, ['pKa']])
            z_right += pow(10, pH) / (pow(10, pH) + pow(10, pka))

    net_charge_val += z_left - z_right
    return net_charge_val


def compute_hydrophilic_residues_DIV_tot_num_of_residues(peptide_list, hydrophilicity_value_table):
    hydrophilic_residues = 0
    for i in range(8):
        if float(hydrophilicity_value_table.loc[(peptide_list[i])[0:3], ['hyd_val']]) > 0:
            hydrophilic_residues += 1
    return abs(hydrophilic_residues / 8) * 100


COMPUTE_ADDITIONAL_PARAMETERS = False


def define_isoelectric_points(peptide_combinations, pks_table, ph, hydrophilicity_value_table,
                              peptides_capital='combinations_cap'):
    num_rows = peptide_combinations.shape[0]
    pIs_list = []
    net_charge_list = []
    hydrophilic_residues_DIV_tot_num_of_residues = []
    list_of_peptides = []
    for i in range(num_rows):
        peptide_list = (peptide_combinations.iloc[i,]).to_numpy()
        # isoelectric point calcul

        peptide = convert_aa_to_capitals(peptide_list)
        list_of_peptides.append(peptide)

        if (COMPUTE_ADDITIONAL_PARAMETERS):
            pIs_list.append(predict_isoelectric_point(peptide))
            net_charge_list.append(net_charge(peptide_list, pks_table, ph))
            hydrophilic_residues_DIV_tot_num_of_residues.append(
                compute_hydrophilic_residues_DIV_tot_num_of_residues(peptide_list, hydrophilicity_value_table))

    peptide_comb_cap = pd.DataFrame(list_of_peptides, columns=['peptide'])

    if (COMPUTE_ADDITIONAL_PARAMETERS):
        pIs_df = pd.DataFrame(pIs_list, columns=['pI'])
        netCharge_df = pd.DataFrame(net_charge_list, columns=['netCharge'])
        hydrophilic_residues_DIV_tot_num_of_residues_df = pd.DataFrame(hydrophilic_residues_DIV_tot_num_of_residues,
                                                                       columns=[
                                                                           'hydrophilic_residues / num_of_residues'])

        peptide_combinations = pd.concat([peptide_combinations, pIs_df], axis=1)
        peptide_combinations = pd.concat([peptide_combinations, netCharge_df], axis=1)
        peptide_combinations = pd.concat([peptide_combinations, hydrophilic_residues_DIV_tot_num_of_residues_df],
                                         axis=1)

    return peptide_combinations, peptide_comb_cap


def execute_2(substrate_search, excludeP4, excludeP3, excludeP2, excludeP1, excludeP1prime, excludeP2prime,
              excludeP3prime, excludeP4prime, pH):
    amino_acid_to_exclude_at_particular_positions = [(1, excludeP4.replace(' ', '').split(',')),
                                                     (2, excludeP3.replace(' ', '').split(',')),
                                                     (3, excludeP2.replace(' ', '').split(',')),
                                                     (4, excludeP1.replace(' ', '').split(',')),
                                                     (5, excludeP1prime.replace(' ', '').split(',')),
                                                     (6, excludeP2prime.replace(' ', '').split(',')),
                                                     (7, excludeP3prime.replace(' ', '').split(',')),
                                                     (8, excludeP4prime.replace(' ', '').split(','))]

    best_peptide_table_AND_freq_table = extract_peptide(substrate_search, amino_acid_to_exclude_at_particular_positions,
                                                        'frequency table_2.csv',
                                                        protease_number=2)  # extract_peptide(True, A, i) exclude amino acids A at position i (i \in [1,...,8], P4=1, ... P4'=8)
    (best_peptide_table_AND_freq_table[1]).to_csv('best_peptide_table_2.csv', index=True)

    list1 = ((best_peptide_table_AND_freq_table[1]).iloc[0,]).values.tolist()
    list2 = ((best_peptide_table_AND_freq_table[1]).iloc[1,]).values.tolist()
    list3 = ((best_peptide_table_AND_freq_table[1]).iloc[2,]).values.tolist()
    list4 = ((best_peptide_table_AND_freq_table[1]).iloc[3,]).values.tolist()

    print('FIRST 4 PEPTIDES FOR PROTEASE_2')
    print(list1)
    print(list2)
    print(list3)
    print(list4)

    # list1 list2 list3 list4 combinations
    for i in range(8):
        if parse_str(list1[i]) > 90.0:
            list2[i] = None
            list3[i] = None
            list4[i] = None
            continue
        if parse_str(list2[i]) > 90.0:
            list1[i] = None
            list3[i] = None
            list4[i] = None
            continue
        if parse_str(list3[i]) > 90.0:
            list1[i] = None
            list2[i] = None
            list4[i] = None
            continue
        if parse_str(list4[i]) > 90.0:
            list1[i] = None
            list2[i] = None
            list3[i] = None
            continue

    peptide_combination_df = pd.DataFrame(create_peptide_combinations(list1, list2, list3, list4))

    d = {
        'aa': ['Arg', 'Asp', 'Cys', 'Glu', 'His', 'Lys', 'Tyr'],
        'pKa': [12.48, 3.86, 8.33, 4.25, 6.0, 10.53, 10.07],  # (COOH)
    }
    pks_table = (pd.DataFrame(data=d)).set_index(['aa'])

    d = {
        'aa': ['Ala', 'Cys', 'Asp', 'Glu', 'Phe', 'Gly', 'His', 'Ile', 'Lys', 'Leu', 'Met', 'Asn', 'Pro', 'Gln', 'Arg',
               'Ser', 'Thr', 'Val', 'Trp', 'Tyr'],
        'hyd_val': [-0.5, -1.0, 3.0, 3.0, -2.5, 0.0, -0.5, -1.8, 3.0, -1.8, -1.3, 0.2, 0.0, 0.2, 3.0, 0.3, -0.4, -1.5,
                    -3.4, -2.3]
    }
    hydrophilicity_value_table = (pd.DataFrame(data=d)).set_index(['aa'])

    peptide_combination, peptide_comb_cap = define_isoelectric_points(peptide_combination_df, pks_table, pH,
                                                                      hydrophilicity_value_table,
                                                                      'peptides_comb_cap_2.csv')
    columns = ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'"]
    if (COMPUTE_ADDITIONAL_PARAMETERS):
        columns = ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'", 'pI',
                   "netChrg", "Hydrophilicity (in %)"]
    peptide_combination.columns = columns

    print("END OF COMMANDS_execute_2")
    return best_peptide_table_AND_freq_table[0], best_peptide_table_AND_freq_table[
        1], peptide_combination, peptide_comb_cap


def execute(protease_1, excludeP4_1, excludeP3_1, excludeP2_1, excludeP1_1, excludeP1prime_1, excludeP2prime_1,
            excludeP3prime_1, excludeP4prime_1, pH_1,
            protease_2, excludeP4_2, excludeP3_2, excludeP2_2, excludeP1_2, excludeP1prime_2, excludeP2prime_2,
            excludeP3prime_2, excludeP4prime_2, pH_2,
            protease_3, excludeP4_3, excludeP3_3, excludeP2_3, excludeP1_3, excludeP1prime_3, excludeP2prime_3,
            excludeP3prime_3, excludeP4prime_3, pH_3,
            protease_4, excludeP4_4, excludeP3_4, excludeP2_4, excludeP1_4, excludeP1prime_4, excludeP2prime_4,
            excludeP3prime_4, excludeP4prime_4, pH_4,
            protease_5, excludeP4_5, excludeP3_5, excludeP2_5, excludeP1_5, excludeP1prime_5, excludeP2prime_5,
            excludeP3prime_5, excludeP4prime_5, pH_5):
    ### query code ###
    # print(protease_1)
    # print(protease_2)
    # print(protease_3)
    # print(protease_4)
    # print(protease_5)
    #
    SUBSTRATE_SEARCH_1 = fetch_table_data(protease_1, ["here is list of cleavage_types to exlude << LATER"], 1)
    SUBSTRATE_SEARCH_2 = fetch_table_data(protease_2, ["here is list of cleavage_types to exlude << LATER"], 2)

    freq_table_1, best_peptide_table_1, peptide_combination_1, peptide_comb_cap_1 = execute_2(SUBSTRATE_SEARCH_1,
                                                                                              excludeP4_1,
                                                                                              excludeP3_1,
                                                                                              excludeP2_1,
                                                                                              excludeP1_1,
                                                                                              excludeP1prime_1,
                                                                                              excludeP2prime_1,
                                                                                              excludeP3prime_1,
                                                                                              excludeP4prime_1,
                                                                                              pH_1)

    freq_table_2, best_peptide_table_2, peptide_combination_2, peptide_comb_cap_2 = execute_2(SUBSTRATE_SEARCH_2,
                                                                                              excludeP4_2,
                                                                                              excludeP3_2,
                                                                                              excludeP2_2,
                                                                                              excludeP1_2,
                                                                                              excludeP1prime_2,
                                                                                              excludeP2prime_2,
                                                                                              excludeP3prime_2,
                                                                                              excludeP4prime_2,
                                                                                              pH_2)

    freq_table_3 = None
    best_peptide_table_3 = None
    peptide_combination_3 = None
    peptide_comb_cap_3 = None

    freq_table_4 = None
    best_peptide_table_4 = None
    peptide_combination_4 = None
    peptide_comb_cap_4 = None

    freq_table_5 = None
    best_peptide_table_5 = None
    peptide_combination_5 = None
    peptide_comb_cap_5 = None

    if protease_3:
        SUBSTRATE_SEARCH_3 = fetch_table_data(protease_3, ["here is list of cleavage_types to exlude << LATER"], 3)
        freq_table_3, best_peptide_table_3, peptide_combination_3, peptide_comb_cap_3 = execute_2(SUBSTRATE_SEARCH_3,
                                                                                                  excludeP4_3,
                                                                                                  excludeP3_3,
                                                                                                  excludeP2_3,
                                                                                                  excludeP1_3,
                                                                                                  excludeP1prime_3,
                                                                                                  excludeP2prime_3,
                                                                                                  excludeP3prime_3,
                                                                                                  excludeP4prime_3,
                                                                                                  pH_3)
    if protease_4:
        SUBSTRATE_SEARCH_4 = fetch_table_data(protease_4, ["here is list of cleavage_types to exlude << LATER"], 4)
        freq_table_4, best_peptide_table_4, peptide_combination_4, peptide_comb_cap_4 = execute_2(SUBSTRATE_SEARCH_4,
                                                                                                  excludeP4_4,
                                                                                                  excludeP3_4,
                                                                                                  excludeP2_4,
                                                                                                  excludeP1_4,
                                                                                                  excludeP1prime_4,
                                                                                                  excludeP2prime_4,
                                                                                                  excludeP3prime_4,
                                                                                                  excludeP4prime_4,
                                                                                                  pH_4)
    if protease_5:
        SUBSTRATE_SEARCH_5 = fetch_table_data(protease_5, ["here is list of cleavage_types to exlude << LATER"], 5)
        freq_table_5, best_peptide_table_5, peptide_combination_5, peptide_comb_cap_5 = execute_2(SUBSTRATE_SEARCH_5,
                                                                                                  excludeP4_5,
                                                                                                  excludeP3_5,
                                                                                                  excludeP2_5,
                                                                                                  excludeP1_5,
                                                                                                  excludeP1prime_5,
                                                                                                  excludeP2prime_5,
                                                                                                  excludeP3prime_5,
                                                                                                  excludeP4prime_5,
                                                                                                  pH_5)

    comb_1 = prepare_combinations_for_matring(peptide_comb_cap_1, peptide_combination_1)
    comb_2 = prepare_combinations_for_matring(peptide_comb_cap_2, peptide_combination_2)
    common_peptides = common_peptide(comb_1, comb_2)
    common_peptides = take_distinct_peptides(common_peptides, allowed_similar_positions_num=4)
    print('distinct_common')
    print(common_peptides)

    unique_peptides_1 = unique_peptide(comb_1, comb_2)
    unique_peptides_1 = take_distinct_peptides(unique_peptides_1, allowed_similar_positions_num=2)
    print('distinct p1')
    print(unique_peptides_1)

    unique_peptides_2 = unique_peptide(comb_2, comb_1)
    unique_peptides_2 = take_distinct_peptides(unique_peptides_2, allowed_similar_positions_num=2)
    print('distinct p2')
    print(unique_peptides_2)

    if protease_3:
        comb_3 = prepare_combinations_for_matring(peptide_comb_cap_3, peptide_combination_3)
        common_peptides = common_peptide(common_peptides, comb_3)
        print("HERE_HERE")
        print(common_peptides)
    if protease_4:
        comb_4 = prepare_combinations_for_matring(peptide_comb_cap_4, peptide_combination_4)
        common_peptides = common_peptide(common_peptides, comb_4)
        print(common_peptides)
    if protease_5:
        comb_5 = prepare_combinations_for_matring(peptide_comb_cap_5, peptide_combination_5)
        common_peptides = common_peptide(common_peptides, comb_5)
        print(common_peptides)

    # clear data a bit
    peptide_combination_1[
        ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'"]] = \
    peptide_combination_1[
        ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'"]].applymap(
        parse_str_aa)
    peptide_combination_2[
        ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'"]] = \
    peptide_combination_2[
        ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'"]].applymap(
        parse_str_aa)

    print("END OF COMMANDS")

    return (
        (protease_1 + "  &  " + protease_2 + "  &  " + protease_3 + "  &  " + protease_4 + "  &  " + protease_5).upper(),
        protease_1.upper(),
        protease_2.upper(),
        freq_table_1,
        freq_table_2,
        best_peptide_table_1,
        best_peptide_table_2,
        peptide_combination_1.head(100),
        peptide_combination_2.head(100),
        common_peptides.head(20),  # .sort_values(by='peptide', ascending=True),
        unique_peptides_1.head(100),  # .sort_values(by='peptide', ascending=True),
        unique_peptides_2.head(100))  # .sort_values(by='peptide', ascending=True))


def prepare_combinations_for_matring(peptide_comb_cap_1, peptide_combination_1):
    comb = peptide_comb_cap_1
    corresponding_percents_1 = peptide_combination_1[
        ["Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1'", "Site_P2'", "Site_P3'", "Site_P4'"]].applymap(
        parse_str)
    corresponding_percents_1 = pd.DataFrame(corresponding_percents_1.sum(axis=1, skipna=True))
    comb = pd.concat([comb, corresponding_percents_1], axis=1)
    comb.columns = ['peptide', 'sum']
    return comb


def common_peptide(comb1, comb2):
    comb1.peptide = comb1.peptide.astype(str)
    comb2.peptide = comb2.peptide.astype(str)
    common_peptides = pd.merge(comb1, comb2, how='inner', on=['peptide'])
    common_peptides.columns = ['peptide', 'sum1', 'sum2']
    sum = pd.DataFrame(common_peptides.sum1 + common_peptides.sum2)
    sum.columns = ['sum']
    common_peptides = pd.concat([common_peptides, sum], axis=1)
    common_peptides = common_peptides.drop(columns=['sum1', 'sum2'])
    common_peptides = common_peptides.sort_values(by='sum', ascending=False)
    return common_peptides


def unique_peptide(comb1, comb2):
    comb1.peptide = comb1.peptide.astype(str)
    comb2.peptide = comb2.peptide.astype(str)
    unique_peptides = comb1[~(comb1.peptide.isin(comb2.peptide))]
    unique_peptides.columns = ['peptide', 'sum']
    unique_peptides = unique_peptides.sort_values(by='sum', ascending=False)
    return unique_peptides


def take_distinct_peptides(peptide_teable, allowed_similar_positions_num):
    res_peptides = peptide_teable.head(1)
    # print(res_peptides)
    num_rows = peptide_teable.shape[0]
    for i in range(num_rows):
        peptide = (peptide_teable.iloc[i, 0])
        num_rows_new = res_peptides.shape[0]
        no_similar_peptides = True
        for j in range(num_rows_new):
            existing_peptide = (res_peptides.iloc[j, 0])
            distance = hamming_distance(peptide, existing_peptide)

            if distance < (8 - allowed_similar_positions_num):
                no_similar_peptides = False
        if no_similar_peptides:
            row = peptide_teable.iloc[i:i + 1, 0:]
            res_peptides = pd.concat([res_peptides, row])

    return res_peptides


def compare_peptides_with_others(peptide, already_selected_peptides, allowed_similar_positions_num):
    a = 1


def hamming_distance(peptide1, peptide2):
    return sum(c1 != c2 for c1, c2 in zip(peptide1, peptide2))


window = Tk()
window.geometry('1000x1000')
window.title("AcidApp")
lblp1main = Label(window, text="Protease_1:")
lblp1main.grid(column=0, row=0)

lblp1pos = Label(window, text="Exclude from position P4:")
lblp1pos.grid(column=0, row=2)
lblp2pos = Label(window, text="Exclude from position P3:")
lblp2pos.grid(column=0, row=3)
lblp3pos = Label(window, text="Exclude from position P2:")
lblp3pos.grid(column=0, row=4)
lblp4pos = Label(window, text="Exclude from position P1:")
lblp4pos.grid(column=0, row=5)

lblp1pos1 = Label(window, text="Exclude from position P4:")
lblp1pos1.grid(column=0, row=6)
lblp2pos2 = Label(window, text="Exclude from position P3:")
lblp2pos2.grid(column=0, row=7)
lblp3pos3 = Label(window, text="Exclude from position P2:")
lblp3pos3.grid(column=0, row=8)
lblp4pos4 = Label(window, text="Exclude from position P1:")
lblp4pos4.grid(column=0, row=9)

lblp1ph = Label(window, text="pH value for net charge calculation:")
lblp1ph.grid(column=0, row=10)

txtp1main = Entry(window, width=25)
txtp1main.grid(column=1, row=0)

txtp1pos = Entry(window, width=25)
txtp1pos.grid(column=1, row=2)
txtp2pos = Entry(window, width=25)
txtp2pos.grid(column=1, row=3)
txtp3pos = Entry(window, width=25)
txtp3pos.grid(column=1, row=4)
txtp4pos = Entry(window, width=25)
txtp4pos.grid(column=1, row=5)

txtp1pos1 = Entry(window, width=25)
txtp1pos1.grid(column=1, row=6)
txtp2pos2 = Entry(window, width=25)
txtp2pos2.grid(column=1, row=7)
txtp3pos3 = Entry(window, width=25)
txtp3pos3.grid(column=1, row=8)
txtp4pos4 = Entry(window, width=25)
txtp4pos4.grid(column=1, row=9)
txtp1ph = Entry(window, width=25)
txtp1ph.grid(column=1, row=10)

lblp2main = Label(window, text="Protease_2:")
lblp2main.grid(column=0, row=11)

lbl2p1pos = Label(window, text="Exclude from position P4:")
lbl2p1pos.grid(column=0, row=12)
lbl2p2pos = Label(window, text="Exclude from position P3:")
lbl2p2pos.grid(column=0, row=13)
lbl2p3pos = Label(window, text="Exclude from position P2:")
lbl2p3pos.grid(column=0, row=14)
lbl2p4pos = Label(window, text="Exclude from position P1:")
lbl2p4pos.grid(column=0, row=15)

lbl2p1pos1 = Label(window, text="Exclude from position P4:")
lbl2p1pos1.grid(column=0, row=16)
lbl2p2pos2 = Label(window, text="Exclude from position P3:")
lbl2p2pos2.grid(column=0, row=17)
lbl2p3pos3 = Label(window, text="Exclude from position P2:")
lbl2p3pos3.grid(column=0, row=18)
lbl2p4pos4 = Label(window, text="Exclude from position P1:")
lbl2p4pos4.grid(column=0, row=19)

lbl2p1ph = Label(window, text="pH value for net charge calculation:")
lbl2p1ph.grid(column=0, row=20)

txt2p1main = Entry(window, width=25)
txt2p1main.grid(column=1, row=11)

txt2p1pos = Entry(window, width=25)
txt2p1pos.grid(column=1, row=12)
txt2p2pos = Entry(window, width=25)
txt2p2pos.grid(column=1, row=13)
txt2p3pos = Entry(window, width=25)
txt2p3pos.grid(column=1, row=14)
txt2p4pos = Entry(window, width=25)
txt2p4pos.grid(column=1, row=15)

txt2p1pos1 = Entry(window, width=25)
txt2p1pos1.grid(column=1, row=16)
txt2p2pos2 = Entry(window, width=25)
txt2p2pos2.grid(column=1, row=17)
txt2p3pos3 = Entry(window, width=25)
txt2p3pos3.grid(column=1, row=18)
txt2p4pos4 = Entry(window, width=25)
txt2p4pos4.grid(column=1, row=19)
txt2p1ph = Entry(window, width=25)
txt2p1ph.grid(column=1, row=20)


lblp4main = Label(window, text="Protease_3:")
lblp4main.grid(column=0, row=21)

lbl3p1pos = Label(window, text="Exclude from position P4:")
lbl3p1pos.grid(column=0, row=22)
lbl3p2pos = Label(window, text="Exclude from position P3:")
lbl3p2pos.grid(column=0, row=23)
lbl3p3pos = Label(window, text="Exclude from position P2:")
lbl3p3pos.grid(column=0, row=24)
lbl3p4pos = Label(window, text="Exclude from position P1:")
lbl3p4pos.grid(column=0, row=25)

lbl3p1pos1 = Label(window, text="Exclude from position P4:")
lbl3p1pos1.grid(column=0, row=26)
lbl3p2pos2 = Label(window, text="Exclude from position P3:")
lbl3p2pos2.grid(column=0, row=27)
lbl3p3pos3 = Label(window, text="Exclude from position P2:")
lbl3p3pos3.grid(column=0, row=28)
lbl3p4pos4 = Label(window, text="Exclude from position P1:")
lbl3p4pos4.grid(column=0, row=29)

lbl3p1ph = Label(window, text="pH value for net charge calculation:")
lbl3p1ph.grid(column=0, row=30)

txt3p1main = Entry(window, width=25)
txt3p1main.grid(column=1, row=21)

txt3p1pos = Entry(window, width=25)
txt3p1pos.grid(column=1, row=22)
txt3p2pos = Entry(window, width=25)
txt3p2pos.grid(column=1, row=23)
txt3p3pos = Entry(window, width=25)
txt3p3pos.grid(column=1, row=24)
txt3p4pos = Entry(window, width=25)
txt3p4pos.grid(column=1, row=25)

txt3p1pos1 = Entry(window, width=25)
txt3p1pos1.grid(column=1, row=26)
txt3p2pos2 = Entry(window, width=25)
txt3p2pos2.grid(column=1, row=27)
txt3p3pos3 = Entry(window, width=25)
txt3p3pos3.grid(column=1, row=28)
txt3p4pos4 = Entry(window, width=25)
txt3p4pos4.grid(column=1, row=29)
txt3p1ph = Entry(window, width=25)
txt3p1ph.grid(column=1, row=30)


lblp4main = Label(window, text="Protease_4:")
lblp4main.grid(column=0, row=31)

lbl4p1pos = Label(window, text="Exclude from position P4:")
lbl4p1pos.grid(column=0, row=32)
lbl4p2pos = Label(window, text="Exclude from position P3:")
lbl4p2pos.grid(column=0, row=33)
lbl4p3pos = Label(window, text="Exclude from position P2:")
lbl4p3pos.grid(column=0, row=34)
lbl4p4pos = Label(window, text="Exclude from position P1:")
lbl4p4pos.grid(column=0, row=35)

lbl4p1pos1 = Label(window, text="Exclude from position P4:")
lbl4p1pos1.grid(column=0, row=36)
lbl4p2pos2 = Label(window, text="Exclude from position P3:")
lbl4p2pos2.grid(column=0, row=37)
lbl4p3pos3 = Label(window, text="Exclude from position P2:")
lbl4p3pos3.grid(column=0, row=38)
lbl4p4pos4 = Label(window, text="Exclude from position P1:")
lbl4p4pos4.grid(column=0, row=39)

lbl4p1ph = Label(window, text="pH value for net charge calculation:")
lbl4p1ph.grid(column=0, row=40)

txt4p1main = Entry(window, width=25)
txt4p1main.grid(column=1, row=31)

txt4p1pos = Entry(window, width=25)
txt4p1pos.grid(column=1, row=32)
txt4p2pos = Entry(window, width=25)
txt4p2pos.grid(column=1, row=33)
txt4p3pos = Entry(window, width=25)
txt4p3pos.grid(column=1, row=34)
txt4p4pos = Entry(window, width=25)
txt4p4pos.grid(column=1, row=35)

txt4p1pos1 = Entry(window, width=25)
txt4p1pos1.grid(column=1, row=36)
txt4p2pos2 = Entry(window, width=25)
txt4p2pos2.grid(column=1, row=37)
txt4p3pos3 = Entry(window, width=25)
txt4p3pos3.grid(column=1, row=38)
txt4p4pos4 = Entry(window, width=25)
txt4p4pos4.grid(column=1, row=39)
txt4p1ph = Entry(window, width=25)
txt4p1ph.grid(column=1, row=40)

lblp5main = Label(window, text="Protease_5:")
lblp5main.grid(column=2, row=1)

lbl5p1pos = Label(window, text="Exclude from position P4:")
lbl5p1pos.grid(column=2, row=2)
lbl5p2pos = Label(window, text="Exclude from position P3:")
lbl5p2pos.grid(column=2, row=3)
lbl5p3pos = Label(window, text="Exclude from position P2:")
lbl5p3pos.grid(column=2, row=4)
lbl5p4pos = Label(window, text="Exclude from position P1:")
lbl5p4pos.grid(column=2, row=5)

lbl5p1pos1 = Label(window, text="Exclude from position P4:")
lbl5p1pos1.grid(column=2, row=6)
lbl5p2pos2 = Label(window, text="Exclude from position P3:")
lbl5p2pos2.grid(column=2, row=7)
lbl5p3pos3 = Label(window, text="Exclude from position P2:")
lbl5p3pos3.grid(column=2, row=8)
lbl5p4pos4 = Label(window, text="Exclude from position P1:")
lbl5p4pos4.grid(column=2, row=9)

lbl5p1ph = Label(window, text="pH value for net charge calculation:")
lbl5p1ph.grid(column=2, row=10)

txt5p1main = Entry(window, width=25)
txt5p1main.grid(column=3, row=1)

txt5p1pos = Entry(window, width=25)
txt5p1pos.grid(column=3, row=2)
txt5p2pos = Entry(window, width=25)
txt5p2pos.grid(column=3, row=3)
txt5p3pos = Entry(window, width=25)
txt5p3pos.grid(column=3, row=4)
txt5p4pos = Entry(window, width=25)
txt5p4pos.grid(column=3, row=5)

txt5p1pos1 = Entry(window, width=25)
txt5p1pos1.grid(column=3, row=6)
txt5p2pos2 = Entry(window, width=25)
txt5p2pos2.grid(column=3, row=7)
txt5p3pos3 = Entry(window, width=25)
txt5p3pos3.grid(column=3, row=8)
txt5p4pos4 = Entry(window, width=25)
txt5p4pos4.grid(column=3, row=9)
txt5p1ph = Entry(window, width=25)
txt5p1ph.grid(column=3, row=10)

btn = Button(window, text="Submit!", command=clicked)
btn.grid(column=2, row=41)
lblres = Label(window, text="")
lblres.grid(column=3, row=11)

window.mainloop()

'''

positions excluded in the sum
veriety
peptides clastering: 
check peptides based on a particular aa on a position 


fisrt: 
unique 
which positions to 

semimanual: inputs: wich position to select: how many residues to check (e.g. 3 or 4)
lets make P3 ala: combination 

selection: top peptide with 


'''