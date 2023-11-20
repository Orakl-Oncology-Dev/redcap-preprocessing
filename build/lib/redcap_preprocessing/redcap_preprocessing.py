
import os

from redcap_preprocessing import split_clinical_data_from_redcap
from redcap_preprocessing import split_treatment_data_from_redcap
from redcap_preprocessing import split_molecular_data_from_redcap

def preprocess_redcap_data(redcap_filepath: str,
                           disease_type: str,
                           save_as_single_file: bool = False,
                           output_dir = None,
                           ):
    
    # check that the disease type is valid
    assert disease_type in ['CRC', 'PDAC']

    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(redcap_filepath), 'preprocessed_redcap_data')

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if disease_type == 'CRC':
        conversion_table_filepath = os.path.join('conversion_table', 'redcap_CRC_conversion_table.csv')
    elif disease_type == 'PDAC':
        conversion_table_filepath = os.path.join('conversion_table', 'redcap_PDAC_conversion_table.csv')

    utils_preprocess_redcap_data(redcap_filepath,
                           conversion_table_filepath,
                           output_dir,
                           disease_type,
                           save_as_single_file
                           )

    return 


def utils_preprocess_redcap_data(redcap_filepath,
                           conversion_table_filepath,
                           output_dir,
                           disease_type,
                           save_as_single_file
                           ):
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    split_clinical_data_from_redcap.split_clinical_data_from_redcap_directory(redcap_filepath,
                                                        conversion_table_filepath,
                                                        output_dir,
                                                        disease_type,
                                                        save_as_single_file)
    
    split_treatment_data_from_redcap.split_treatment_data_from_redcap(redcap_filepath,
                                                        conversion_table_filepath,
                                                        output_dir,
                                                        disease_type,
                                                        save_as_single_file)
    
    split_molecular_data_from_redcap.split_molecular_data_from_redcap(redcap_filepath,
                                                        conversion_table_filepath,
                                                        output_dir,
                                                        disease_type,
                                                        save_as_single_file)
    
    return None