#import required libraries
import os
import pandas as pd
from PyPDF2 import PdfReader


#Process pt_src files
def process_ptsrc_pdf(file_path, filename):
    try:
        reader = PdfReader(file_path)
        if len(reader.pages) < 3:
            print(f"PDF file '{filename}' has less than 3 pages. Skipping.")
            return None
        page = reader.pages[2]
        text = page.extract_text()
        lines = text.split('\n')
        for line in lines:
            if 'Pb-210' in line:
                ptsrc_pb210, PtSrc_Pb210error = line.split()[-2:]
                return {
                    'File': filename,
                    'Pb-210': float(ptsrc_pb210),
                    'Pb-210 error': float(PtSrc_Pb210error)
                }
        print(f"Pb-210 not found in '{filename}'. Skipping.")
        return None
    except Exception as e:
        print(f"Error processing PDF file '{filename}': {e}")
        return None

#Process regular files
def process_regular_pdf(file_path, filename):
    pb210 = pb210error = Bi214 = Bi214error = Pb214 = Pb214error = None
    try:
        reader = PdfReader(file_path)
        if len(reader.pages) < 3:
            print(f"PDF file '{filename}' has less than 3 pages. Skipping.")
            return None
        page = reader.pages[2]
        text = page.extract_text()
        lines = text.split('\n')
        for line in lines:
            if 'Pb-210' in line:
                pb210, pb210error = line.split()[-2:]
            elif 'Bi-214' in line:
                Bi214, Bi214error = line.split()[-2:]
            elif 'Pb-214' in line:
                Pb214, Pb214error = line.split()[-2:]
        if pb210 is None or pb210error is None:
            print(f"Pb-210 not found in '{filename}'. Skipping.")
            return None
        if Bi214 is None or Bi214error is None:
            print(f"Bi-214 not found in '{filename}'. Skipping.")
            return None
        if Pb214 is None or Pb214error is None:
            print(f"Pb-214 not found in '{filename}'. Skipping.")
            return None
        return {
            'File': filename,
            'Pb-210': float(pb210),
            'Pb-210 error': float(pb210error),
            'Bi-214': float(Bi214),
            'Bi-214 error': float(Bi214error),
            'Pb-214': float(Pb214),
            'Pb-214 error': float(Pb214error)
        }
    except Exception as e:
        print(f"Error processing PDF file '{filename}': {e}")
        return None

#define PDF processing
def process_pdf_file(file_path, filename):
    if filename.startswith("PtSrc_"):
        return process_ptsrc_pdf(file_path, filename)
    else:
        return process_regular_pdf(file_path, filename)

#Extract PDF data from folder
def extract_pdf_data(folder_path):
    combined_data = []
    for filename in os.listdir(folder_path):
        # Check for PDF extension (case-insensitive)
        if filename.lower().endswith(".pdf"):
            file_path = os.path.join(folder_path, filename)
            data = process_pdf_file(file_path, filename)
            if data is not None:
                combined_data.append(data)
    return combined_data

#Sort and process PDFs
def sort_pdf_data(combined_data, parse_numbers=False):
    combined_df = pd.DataFrame(combined_data)

    def extract_numeric_suffix(file_name):
        try:
            parts = file_name.split('_')[-1].split('.')[0]
            return int(parts)
        except ValueError:
            return float('nan')

    combined_df['File_order'] = combined_df['File'].apply(extract_numeric_suffix)
    combined_df = combined_df.sort_values(by='File_order').drop(columns=['File_order'])

    if parse_numbers:
        for col in ['Pb-210', 'Bi-214', 'Pb-214']:
            if col in combined_df.columns:
                combined_df[col] = pd.to_numeric(combined_df[col], errors='coerce')
    return combined_df

#Export to csv
def extract_pdf_values(folder_path, output_csv_path, parse_numbers=False):
    combined_data = extract_pdf_data(folder_path)
    combined_df = sort_pdf_data(combined_data, parse_numbers)
    combined_df.to_csv(output_csv_path, index=False)


def process_activity(weights_file_path, canberra_file_path, output_file, year_of_core):
    #Load and merge data
    # Load the CSV files
    csv1 = pd.read_csv(weights_file_path)
    csv2 = pd.read_csv(canberra_file_path)

    # Extract 'Center point of interval' from csv2 based on the median of the last digits in 'File'
    csv2['Center point of interval'] = csv2['File'].apply(
        lambda x: np.median([int(num) for num in re.findall(r'\d+', x.split('_')[-1])])
    )

    # Merge CSV files on 'Center point of interval'
    data = pd.merge(csv1, csv2, on='Center point of interval', how='left')
    # Calculate activity and correction factors
    data['Pb-210 activity Uncertainty (Bq-g)'] = data['Pb-210 error']/data['sediment weight (g)']
    data['Pb-210 activity (Bq/g)'] = data['Pb-210'] / data['sediment weight (g)']
    data['Pb-210 correction factor'] = data['ptsrc_pb210'] / 151031.56
    data['Self absorb. Corrected Pb-210 activity (Bq/g)'] = (
        data['Pb-210 activity (Bq/g)'] / data['Pb-210 correction factor']
    )
    data['Bi-214 activity (Bq/g)'] = data['Bi-214'] / data['sediment weight (g)']
    data['Pb-214 activity (Bq/g)'] = data['Pb-214'] / data['sediment weight (g)']
    # Calculate averaged supported activity of Bi-214 and Pb-214 (Bq/g)
    data['Averaged supported activity of Bi-214 and Pb-214 (Bq/g)'] = (
        data['Bi-214 activity (Bq/g)'] + data['Pb-214 activity (Bq/g)']
    ) / 2

    # Calculate background activity uncertainty (Bq/g)
    data['Background activity uncertainty (Bq/g)'] = (
        (data['Bi-214 error'] + data['Pb-214 error']) / 2
    ) / data['sediment weight (g)']
    # Calculate Excess Pb-210 (Bq/g)
    data['Excess Pb-210 (Bq/g)'] = (
        data['Self absorb. Corrected Pb-210 activity (Bq/g)'] -
        data['Averaged supported activity of Bi-214 and Pb-214 (Bq/g)']
    )

    # Determine surface activity (first interval value)
    data['Surface activity'] = data['Excess Pb-210 (Bq/g)'].iloc[0]

    # Calculate Age bp using the natural logarithm
    data['Age bp'] = (1 / 0.03114) * np.log(data['Surface activity'] / data['Excess Pb-210 (Bq/g)'])
    # Calculate 'calendar years pre year of core'
    data['calendar years pre year of core'] = year_of_core - data['Age bp']

    # Save the final DataFrame to a CSV file
    data.to_csv(output_file_name, index=False)

    print(f"Calculations completed, data exported to '{output_file_name}'")
