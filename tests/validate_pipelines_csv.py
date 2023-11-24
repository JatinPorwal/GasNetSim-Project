#!/usr/bin/env python
# coding: utf-8

# **********************************************************************************************************************
# This script ensures the pipelines.csv file conforms to the established format guidelines.
# **********************************************************************************************************************

import csv
import os
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def validate_pipeline_csv_format(file_path):
    expected_headers = [
        'pipeline_index',
        'inlet_index',
        'outlet_index',
        'diameter_m',
        'length_m',
        ' is_bothDirection',
        'max_cap_M_m3_per_d',
        'max_pressure_bar',
        'remarks'
    ]
    
    mandatory_headers = [
        'pipeline_index',
        'inlet_index',
        'outlet_index',
        'diameter_m',
        'length_m'
    ]
    if not (os.path.isfile(file_path)):
        logger.info(f"The file is either missing or not readable")
        return False

    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter=';')
        headers = reader.fieldnames

        # Check if headers match the expected format
        if headers != expected_headers:
            logger.info(f"The CSV file does not match the expected format.")
            return False
        
        # Check each row for missing values in required columns
        for row in reader:
            for key in mandatory_headers:
                if row[key] == '' or row[key] is None or row[key].lower() == 'nan':  # Check for None values or 'nan'
                    logger.info(f"Missing value in column '{key}' in one or more rows.")
                    return False

    logger.info(f"The CSV file follows the expected format.")
    return True


# Check if the script is being run directly
if __name__ == "__main__":
    # Input the 'File_path' with the path of the pipeline CSV file
    # File_path = '../examples/Irish13/Irish13_pipelines.csv'
    File_path = input("Enter the path of the pipelines.csv.\nFor example "
                      "'../examples/Irish13/Irish13_pipelines.csv':")
    validate_pipeline_csv_format(File_path)
