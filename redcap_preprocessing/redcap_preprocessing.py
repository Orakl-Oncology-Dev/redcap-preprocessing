
import os

from redcap_preprocessing import split_clinical_data_from_redcap
from redcap_preprocessing import split_treatment_data_from_redcap
from redcap_preprocessing import split_molecular_data_from_redcap


def preprocess_redcap_data(redcap_filepath,
                           conversion_table_filepath,
                           output_dir,
                           ):
    

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    split_clinical_data_from_redcap.split_clinical_data_from_redcap_directory(redcap_filepath,
                                                        conversion_table_filepath,
                                                        output_dir)
    
    split_treatment_data_from_redcap.split_treatment_data_from_redcap(redcap_filepath,
                                                        conversion_table_filepath,
                                                        output_dir)
    
    split_molecular_data_from_redcap.split_molecular_data_from_redcap(redcap_filepath,
                                                        conversion_table_filepath,
                                                        output_dir)
    
    return None