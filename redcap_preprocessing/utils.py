
import pandas as pd
import numpy as np
import csv

def standardize_code(code, prefix='GR'):

    # Extract the numeric part of the code
    numeric_part = code[len(prefix):]

    # Standardize the numeric part to have 4 digits, padding with zeros if necessary
    standardized_numeric = numeric_part.zfill(4)

    # Combine back with the prefix
    return prefix + standardized_numeric

def get_cell_line_code(disease_type: str,
                       redcap: pd.DataFrame,
                       record_id: str):
    
    if disease_type == 'CRC':
        cell_line_code, date_cell_line = get_cell_line_code_CRC(redcap, record_id)
        prefix = 'CGR'
    elif disease_type == 'PDAC':
        cell_line_code, date_cell_line = get_cell_line_code_PDAC(redcap, record_id)
        prefix = 'PGR'

    # reformat cell line codes
    cell_line_code = standardize_code(cell_line_code, prefix)

    return cell_line_code, date_cell_line

def get_cell_line_code_CRC(redcap: pd.DataFrame,
                       record_id: str):
    """
    Get the cell line code from the redcap data set.add

    Parameters
    ----------
    redcap: pd.DataFrame
        redcap data set
    record_id: str
    """

    selected_row = redcap[(redcap['record_id'] == record_id) & 
                          (redcap['redcap_repeat_instrument'] == 'organoides')]
        
    if len(selected_row) == 0:
        print(f'There are no PDO cell lines for record_id {record_id}.')
        cell_line_code, date_cell_line = '', ''
    else:
        cell_line_code = ';'.join(selected_row['nom_lign_e'].unique())
        date_cell_line = ';'.join(selected_row['date_sample'].unique())

    return cell_line_code, date_cell_line

def get_cell_line_code_PDAC(redcap: pd.DataFrame,
                       record_id: str):
    """
    Get the cell line code from the redcap data set.add

    Parameters
    ----------
    redcap: pd.DataFrame
        redcap data set
    record_id: str
    """

    selected_row = redcap[(redcap['record_id'] == record_id) & 
                          (redcap['redcap_repeat_instrument'] == 'organodes')]
        
    if len(selected_row) == 0:
        print(f'There are no PDO cell lines for record_id {record_id}.')
        cell_line_code, date_cell_line = '', ''
    else:
        cell_line_code = ';'.join(selected_row['namepdo'].unique())
        date_cell_line = ';'.join(selected_row['date_pdo'].fillna('').unique())

    return cell_line_code, date_cell_line

def get_content_matching_type_1(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    assert column_name in redcap_CRC_conversion_table.orakloncology_name.values, 'column name not found in conversion table'
    assert len(redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 1)]['redcap_name']) == 1, 'multiple redcap names found for column name in a type 1 match'

    redcap_column_name = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 1)]['redcap_name'].values[0]
    content = single_patient_clinical_row[redcap_column_name]

    return content

def get_content_matching_type_2(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    content = []
    redcap_column_rows = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 2)]

    for index, row in redcap_column_rows.iterrows():

        if single_patient_clinical_row[row['redcap_name']] == 1:
            content.append(row['orakloncology_options'])

    return ';'.join(content)

def get_content_matching_type_3(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    content = []
    redcap_column_rows = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 3)]

    for index, row in redcap_column_rows.iterrows():

        if row['redcap_options'] == single_patient_clinical_row[row['redcap_name']]:
            content.append(row['orakloncology_options'])

    return ';'.join(content)

def get_content_matching_type_4(single_patient_clinical_row,
                                redcap_CRC_conversion_table,
                                column_name,
                                ):
    
    content = []
    redcap_column_rows = redcap_CRC_conversion_table[(redcap_CRC_conversion_table.orakloncology_name == column_name)&(redcap_CRC_conversion_table.matching_type == 4)]

    for index, row in redcap_column_rows.iterrows():
        
        # test if not nan, then add to content
        if single_patient_clinical_row[row['redcap_name']] == single_patient_clinical_row[row['redcap_name']]:
            content.append(single_patient_clinical_row[row['redcap_name']])

    return ';'.join(content)

def add_content(content, prior_content):

    # test if empty str
    if prior_content == '':
        return content
    
    # test if nan
    elif prior_content != prior_content:
        return content
    
    # test if empty str
    elif content == '':
        return prior_content
    
    # test if nan
    elif content != content:
        return prior_content
    
    else:
        # test if the content and prior content are of the same type
        assert type(content) == type(prior_content), 'content and prior content are not of the same type'
        
        if type(content) == str:
            return prior_content + ';' + content
        
        else:
            print(content)


def get_delimiter(csv_path):

    with open(csv_path, 'r') as file:
        sample = file.read(1024)

    sniffer = csv.Sniffer()

    return sniffer.sniff(sample).delimiter