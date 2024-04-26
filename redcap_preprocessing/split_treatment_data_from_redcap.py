import pandas as pd
import numpy as np
import os
import chardet

from redcap_preprocessing.utils import get_cell_line_code, get_content_matching_type_1, get_content_matching_type_2, get_content_matching_type_3, get_content_matching_type_4, add_content
from redcap_preprocessing.utils import get_delimiter
from redcap_preprocessing.utils import standardize_code
from redcap_preprocessing.utils import format_dates

def get_single_patient_treatment_data(patient_treatment_data,
                                      redcap_CRC_conversion_table,
                                      disease_type,):

    cleaned_patient_treatment_data = pd.DataFrame('',
                                                index = patient_treatment_data.index,
                                                columns = redcap_CRC_conversion_table['orakloncology_name'].unique(),)
    
    if disease_type == 'PDAC':
        slice_df = patient_treatment_data[patient_treatment_data['bras_essai_clinique_met'].isna()]
        patient_treatment_data.loc[slice_df.index, 'bras_essai_clinique_met'] = patient_treatment_data['bras_essai_clinique_met'].fillna('')

        slice_df = patient_treatment_data[patient_treatment_data['spe_autre_met'].isna()]
        patient_treatment_data.loc[slice_df.index, 'spe_autre_met'] = patient_treatment_data['spe_autre_met'].fillna('')

    elif disease_type == 'CRC':
        slice_df = patient_treatment_data[patient_treatment_data['other_iv'].isna()]
        patient_treatment_data.loc[slice_df.index, 'other_iv'] = patient_treatment_data['other_iv'].fillna('')

    # Select single therapy row
    for index, row in patient_treatment_data.iterrows():

        # loop through all the columns in the redcap table
        for column_name in cleaned_patient_treatment_data.columns:

            matching_types = redcap_CRC_conversion_table[redcap_CRC_conversion_table['orakloncology_name'] == column_name]['matching_type'].unique()

            for matching_type in matching_types:

                if matching_type == 1:

                    content = get_content_matching_type_1(row, redcap_CRC_conversion_table, column_name)

                if matching_type == 2:

                    content = get_content_matching_type_2(row, redcap_CRC_conversion_table, column_name)

                if matching_type == 3:

                    content = get_content_matching_type_3(row, redcap_CRC_conversion_table, column_name)

                if matching_type == 4:

                    content = get_content_matching_type_4(row, redcap_CRC_conversion_table, column_name)

                # special case for chemotherapy_type
                if disease_type == 'CRC':
                    if (column_name == 'chemotherapy_type')&(content == 'other'):
                        content = row['other_iv'].lower()
                elif disease_type == 'PDAC':
                    if (column_name == 'chemotherapy_type')&(content == 'other'):
                        content = row['spe_autre_met'].lower()
                    elif (column_name == 'chemotherapy_type')&(content == 'clinical_trial'):
                        content = row['bras_essai_clinique_met'].lower()

                # formating
                elif (column_name == 'targeted_therapy_type')&(content != '')&(content == content):
                    content = content.lower()
                
                content= add_content(content, cleaned_patient_treatment_data.loc[index, column_name])
                cleaned_patient_treatment_data.loc[index, column_name] = content

    # add treatment type            
    cleaned_patient_treatment_data['treatment_type'] = patient_treatment_data['redcap_repeat_instrument'].replace('ligne_mtastatique_de_traitement', 'metastatic').replace('hai_chemotherapy', 'hai').replace('', 'neo_adjuvant')

    # if neo_adjuvant, change stop cause
    cleaned_patient_treatment_data.loc[cleaned_patient_treatment_data['treatment_type'] == 'neo_adjuvant', 'stop_chemotherapy_cause'] = 'surgery'

    return cleaned_patient_treatment_data
    
def split_treatment_data_from_redcap(redcap_path: str,
                       redcap_conversion_table_path: str,
                       output_dir: str,
                       disease_type: str,
                       save_as_single_file: bool = False,
                       ):
    """
    Get the treatment data from the redcap data set.

    Parameters
    ----------
    redcap_path: str
        path to redcap data set
    redcap_conversion_table_path: str
        path to redcap CRC conversion table
    output_dir: str
        output directory
    """

    def detect_encoding(file_path):
        with open(file_path, 'rb') as file:
            result = chardet.detect(file.read())
        return result['encoding']

    file_encoding = detect_encoding(redcap_path)
    print("Detected encoding:", file_encoding)

    # detect the separator type and read the data
    redcap_delimiter = get_delimiter(redcap_path)
    redcap = pd.read_csv(redcap_path, delimiter=redcap_delimiter, encoding="unicode_escape")

    # some of the redcap files are not encoded in unicode_escape
    try:
        redcap = pd.read_csv(redcap_path, delimiter=redcap_delimiter, encoding="unicode_escape")
        # get the unique record ids
        record_ids = redcap.record_id.unique()
    except:
        redcap = pd.read_csv(redcap_path, delimiter=redcap_delimiter, encoding=file_encoding)
        # get the unique record ids
        record_ids = redcap.record_id.unique()

    # Read in the data
    redcap['redcap_repeat_instrument'] = redcap['redcap_repeat_instrument'].fillna('')

    # check that the table isn't empty
    if len(redcap) == 0:
        raise ValueError('The redcap table is empty.')
    
    # verify that the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # column name mapping
    conversion_table_delimiter = get_delimiter(redcap_conversion_table_path)
    redcap_CRC_conversion_table = pd.read_csv(redcap_conversion_table_path, delimiter=conversion_table_delimiter)    
    redcap_CRC_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'treatment']
    redcap_CRC_conversion_table['redcap_name'] = redcap_CRC_conversion_table['redcap_name'].str.strip()

    # recap dataframe
    cleaned_patient_treatment_data = pd.DataFrame()
    # index for patients w/o cell lines
    i = 0

    # loop through all the record ids
    for record_id in record_ids:

        if True:

            # get the cell line code
            cell_line_code, date_cell_line = get_cell_line_code(disease_type, redcap, record_id)

            # select the patient data
            patient_treatment_data = redcap[(redcap.record_id == record_id) &
                                            (redcap['redcap_repeat_instrument'].isin(['','ligne_mtastatique_de_traitement', 'hai_chemotherapy']))]
            
            # test if there is any data, otherwise skip
            if len(patient_treatment_data) > 0:

                # get the single patient treatment data
                cleaned_single_patient_treatment_data = get_single_patient_treatment_data(patient_treatment_data,
                                                                                    redcap_CRC_conversion_table,
                                                                                    disease_type)
                # TO IMPROVE
                if 'CGR' in cell_line_code:
                    cell_line_code = cell_line_code.replace('CGR', 'GR')
                elif 'PGR' in cell_line_code:
                    cell_line_code = cell_line_code.replace('PGR', 'GR')

                extracted_cell_line_codes = cell_line_code.split(';')
                extracted_cell_line_dates = date_cell_line.split(';')



                # if multiple treatment lines are present, save twice under different names
                for cell_line_code, cell_line_date in zip(extracted_cell_line_codes, extracted_cell_line_dates):
                    
                    # we had bugs when single patient had multiple cell lines
                    # needed to correct this by running all the normalizations on a virgin table
                    # for that patient
                    cleaned_single_patient_treatment_data_unique_cell_line = cleaned_single_patient_treatment_data.copy()

                    if len(cell_line_code)>0:
                        cell_line_code = standardize_code(cell_line_code, 'GR')
                    else:
                        cell_line_code = 'XX' + f'{i:04d}'
                        i += 1

                    # add the cell line code and date
                    cleaned_single_patient_treatment_data_unique_cell_line['cell_line_code'] = cell_line_code
                    cleaned_single_patient_treatment_data_unique_cell_line['sister_cell_line_codes'] = ';'.join([element for element in extracted_cell_line_codes if element != cell_line_code])
                    cleaned_single_patient_treatment_data_unique_cell_line['date_cell_line'] = cell_line_date
                    cleaned_single_patient_treatment_data_unique_cell_line['record_id'] = record_id

                    cleaned_single_patient_treatment_data_unique_cell_line = format_dates(cleaned_single_patient_treatment_data_unique_cell_line)

                    if save_as_single_file:
                        cleaned_patient_treatment_data = pd.concat([cleaned_patient_treatment_data, cleaned_single_patient_treatment_data_unique_cell_line])
                    else:
                        # create a file name
                        if disease_type == 'CRC':
                            filename = f'{output_dir}/TTR_C_PID_{cell_line_code}_SID_0001.csv'
                        elif disease_type == 'PDAC':
                            filename = f'{output_dir}/TTR_P_PID_{cell_line_code}_SID_0001.csv'
                        cleaned_single_patient_treatment_data_unique_cell_line.to_csv(filename, sep=';')

        #except:
        #    print(f'Error for {record_id}')

    if save_as_single_file:

        # reorder the columns with record_id, cell_line_code and date_cell_line in front
        all_columns = list(cleaned_patient_treatment_data.columns)

        # Columns to be placed at the beginning
        selected_columns = ['record_id', 'cell_line_code', 'date_cell_line']

        # Filter out the selected columns from all_columns
        other_columns = [col for col in all_columns if col not in selected_columns]

        # Concatenate selected columns with other columns
        new_column_order = selected_columns + other_columns

        # Reorder the DataFrame
        cleaned_patient_treatment_data = cleaned_patient_treatment_data[new_column_order]

        # go and get the date of death from clincial data frame
        if disease_type == 'CRC':
            clinical_data = pd.read_csv(f'{output_dir}/CLI_C_PID_ALL.csv', sep=';')
        elif disease_type == 'PDAC':
            clinical_data = pd.read_csv(f'{output_dir}/CLI_P_PID_ALL.csv', sep=';')

        # do something cleaner for data coming from multiple modalities
        cleaned_patient_treatment_data = add_death_date(cleaned_patient_treatment_data, clinical_data)

        if disease_type == 'CRC':
            filename = f'{output_dir}/TTR_C_PID_ALL.csv'
        elif disease_type == 'PDAC':
            filename = f'{output_dir}/TTR_P_PID_ALL.csv'

        cleaned_patient_treatment_data.to_csv(filename, sep=';')
    return None


def add_death_date(treatment_data, clinical_data):
    
    # get the death date
    death_date = clinical_data[['record_id', 'date_death']]

    # merge the data
    treatment_data = pd.merge(treatment_data, death_date, on='record_id', how='left')



    return treatment_data