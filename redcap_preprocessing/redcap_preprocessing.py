
import os

from redcap_preprocessing import split_clinical_data_from_redcap
from redcap_preprocessing import split_treatment_data_from_redcap
from redcap_preprocessing import split_molecular_data_from_redcap


def preprocess_redcap_data(redcap_filepath,
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