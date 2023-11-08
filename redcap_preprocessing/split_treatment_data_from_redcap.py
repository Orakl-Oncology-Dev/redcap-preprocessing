import pandas as pd
import numpy as np

from utils import get_cell_line_code, get_content_matching_type_1, get_content_matching_type_2, get_content_matching_type_3, get_content_matching_type_4, add_content

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

def get_single_patient_treatment_data(patient_treatment_data,
                                      redcap_CRC_conversion_table):

    cleaned_patient_treatment_data = pd.DataFrame('',
                                                index = patient_treatment_data.index,
                                                columns = redcap_CRC_conversion_table['orakloncology_name'].unique(),)

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
                if (column_name == 'chemotherapy_type')&(content == 'other'):
                    content = row['other_iv'].lower()

                # formating
                elif (column_name == 'targeted_therapy_type')&(content != '')&(content == content):
                    content = content.lower()
                
                content= add_content(content, cleaned_patient_treatment_data.loc[index, column_name])
                cleaned_patient_treatment_data.loc[index, column_name] = content

    return cleaned_patient_treatment_data
    
def split_treatment_data_from_redcap(redcap_path: str,
                       redcap_conversion_table_path: str,
                       output_dir: str,
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

    # Read in the data
    redcap = pd.read_csv(redcap_path, sep=';')

    # check that the table isn't empty
    if len(redcap) == 0:
        raise ValueError('The redcap table is empty.')

    # column name mapping
    redcap_CRC_conversion_table = pd.read_csv(redcap_conversion_table_path, sep=';')
    redcap_CRC_conversion_table = redcap_CRC_conversion_table[redcap_CRC_conversion_table.data_type == 'treatment']
    redcap_CRC_conversion_table['redcap_name'] = redcap_CRC_conversion_table['redcap_name'].str.strip()

    # get the unique record ids
    record_ids = redcap.record_id.unique()

    # loop through all the record ids
    for record_id in record_ids:

        try:

            # get the cell line code
            cell_line_code, date_cell_line = get_cell_line_code(redcap, record_id)

            # select the patient data
            patient_treatment_data = redcap[(redcap.record_id == record_id) &
                                            (redcap['redcap_repeat_instrument'].isin(['ligne_mtastatique_de_traitement', 'hai_chemotherapy']))]

            # get the single patient treatment data
            cleaned_patient_treatment_data = get_single_patient_treatment_data(patient_treatment_data,
                                                                                redcap_CRC_conversion_table)

            if len(cleaned_patient_treatment_data) > 0:

                # add the cell line code and date
                cleaned_patient_treatment_data['cell_line_code'] = cell_line_code
                cleaned_patient_treatment_data['date_cell_line'] = date_cell_line

                # create a file name
                filename = f'{output_dir}/TR_C_PID_{cell_line_code}_SID_0001.csv'   

                # save the data
                cleaned_patient_treatment_data.to_csv(filename, sep=';')

        except:
            print(f'Error for {record_id}')

    return None